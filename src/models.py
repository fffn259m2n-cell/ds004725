"""Stage 4: QC aggregation, pre-registered moderation models, sensitivities.

Implements reports/preregistration.md sections 5-7. Usage:
    python -I src/models.py            # writes results/ and reports/figures/
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

N_PERM = 10_000
SEED = 20261108
FRAME = "derivatives/phenotype_audit/analysis_frame.tsv"
QC_DIR = Path("derivatives/oxtfc")

# ----------------------------------------------------------------- QC ------


def load_qc(qc_dir: Path = QC_DIR) -> pd.DataFrame:
    rows = []
    for f in sorted(qc_dir.glob("sub-*/sub-*_qc.json")):
        q = json.loads(f.read_text())
        r = {k: q[k] for k in ("participant_id", "tr", "mean_fd", "max_fd",
                               "n_keep", "n_keep_strict", "norm_r",
                               "whole_brain_tsnr")}
        for k, v in q["roi_coverage"].items():
            r[f"cov_{k}"] = v
        for k, v in q["roi_tsnr"].items():
            r[f"tsnr_{k}"] = v
        for tag in ("fc", "fc_gsr", "fc_strict"):
            for k, v in q.get(tag, {}).items():
                r[f"{tag}_{k}"] = v
        rows.append(r)
    return pd.DataFrame(rows)


def qc_exclusions(qc: pd.DataFrame) -> pd.DataFrame:
    """Pre-registered rules (s5, D5). Adds boolean flag columns + 'qc_pass'."""
    q = qc.copy()
    med = q.norm_r.median()
    mad = (q.norm_r - med).abs().median() * 1.4826
    q["flag_fd"] = q.mean_fd > 0.5
    q["flag_nkeep"] = q.n_keep < 150
    q["flag_cov"] = q[[c for c in q if c.startswith("cov_")]].min(1) < 0.5
    q["flag_norm_review"] = q.norm_r < med - 3 * mad
    q["qc_pass"] = ~(q.flag_fd | q.flag_nkeep | q.flag_cov)
    return q


# ------------------------------------------------------------- models ------

PRIMARY = "FC ~ T + zmod + T:zmod + age_c + T:age_c + sex + fd_c"
REDUCED = "FC ~ T + zmod + age_c + T:age_c + sex + fd_c"


def prepare(frame: pd.DataFrame, qc: pd.DataFrame, fc_col: str, mod: str,
            transform: str = "log") -> pd.DataFrame:
    d = frame.merge(qc, on="participant_id", how="inner", validate="1:1")
    raw = d[mod.replace("log_", "") + "_pg_ml"] if transform == "raw" else d[mod]
    d = d.assign(FC=d[fc_col], T=(d.treatment_group == "OT").astype(float),
                 modv=raw)
    d = d.dropna(subset=["FC", "modv", "age", "sex", "mean_fd"])
    d["zmod"] = (d.modv - d.modv.mean()) / d.modv.std(ddof=1)
    d["age_c"] = (d.age - d.age.mean()) / 10
    d["sex"] = (d.sex == "F").astype(float)
    d["fd_c"] = d.mean_fd - d.mean_fd.mean()
    d["older"] = (d.cohort == "older").astype(float)
    return d


def freedman_lane(d: pd.DataFrame, full: str, reduced: str, term: str,
                  n_perm: int = N_PERM, seed: int = SEED) -> float:
    """Two-sided Freedman-Lane permutation p for `term` (OLS t statistic)."""
    rng = np.random.default_rng(seed)
    t_obs = smf.ols(full, d).fit().tvalues[term]
    red = smf.ols(reduced, d).fit()
    fit, res = red.fittedvalues.to_numpy(), red.resid.to_numpy()
    import patsy
    y, X = patsy.dmatrices(full, d, return_type="dataframe")
    X = X.to_numpy()
    j = list(patsy.dmatrices(full, d, return_type="dataframe")[1].columns).index(term)
    XtX_inv = np.linalg.pinv(X.T @ X)
    n, k = X.shape
    count = 0
    for _ in range(n_perm):
        ys = fit + rng.permutation(res)
        b = XtX_inv @ X.T @ ys
        r = ys - X @ b
        se = np.sqrt((r @ r) / (n - k) * XtX_inv[j, j])
        count += abs(b[j] / se) >= abs(t_obs)
    return (count + 1) / (n_perm + 1)


def fit_interaction(d: pd.DataFrame, full: str = PRIMARY,
                    reduced: str = REDUCED, term: str = "T:zmod",
                    n_perm: int = N_PERM) -> dict:
    m = smf.ols(full, d).fit(cov_type="HC3")
    ci = m.conf_int().loc[term]
    out = {"n": int(m.nobs), "n_OT": int(d["T"].sum()),
           "n_PL": int((1 - d["T"]).sum()), "beta": m.params[term],
           "se_hc3": m.bse[term], "ci_low": ci[0], "ci_high": ci[1],
           "p_hc3": m.pvalues[term],
           "p_perm": freedman_lane(d, full, reduced, term, n_perm) if n_perm else np.nan,
           "beta_T": m.params["T"], "p_T": m.pvalues["T"],
           "sd_FC": d.FC.std(ddof=1),
           "beta_std": m.params[term] / d.FC.std(ddof=1)}
    # conditional OT-PL effect at zmod = -1, 0, +1 (other terms at centre)
    names = list(m.params.index)
    for z in (-1, 0, 1):
        L = np.zeros(len(names))
        L[names.index("T")] = 1
        if "T:zmod" in names:
            L[names.index("T:zmod")] = z
        est = float(L @ m.params.to_numpy())
        se = float(np.sqrt(L @ m.cov_params().to_numpy() @ L))
        out[f"cond_T_at_z{z:+d}"] = est
        out[f"cond_T_at_z{z:+d}_lo"] = est - 1.96 * se
        out[f"cond_T_at_z{z:+d}_hi"] = est + 1.96 * se
    return out


def diagnostics(d: pd.DataFrame, full: str = PRIMARY) -> dict:
    m = smf.ols(full, d).fit()
    infl = m.get_influence()
    cooks = infl.cooks_distance[0]
    q = smf.ols(full + " + I(zmod**2) + T:I(zmod**2)", d).fit(cov_type="HC3")
    from scipy import stats
    return {"n_cooks_gt_4n": int((cooks > 4 / len(d)).sum()),
            "max_cooks": float(cooks.max()),
            "max_leverage": float(infl.hat_matrix_diag.max()),
            "shapiro_p_resid": float(stats.shapiro(m.resid).pvalue),
            "quad_TxZ2_beta": float(q.params["T:I(zmod ** 2)"]),
            "quad_TxZ2_p": float(q.pvalues["T:I(zmod ** 2)"]),
            "zmod_range_OT": [float(d.loc[d["T"] == 1, "zmod"].min()),
                              float(d.loc[d["T"] == 1, "zmod"].max())],
            "zmod_range_PL": [float(d.loc[d["T"] == 0, "zmod"].min()),
                              float(d.loc[d["T"] == 0, "zmod"].max())]}


def holm(p: list[float]) -> list[float]:
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for rank, i in enumerate(o):
        run = max(run, (len(p) - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj.tolist()


def robust_fit(d: pd.DataFrame, full: str = PRIMARY) -> dict:
    import statsmodels.api as sm
    m = smf.rlm(full, d, M=sm.robust.norms.HuberT()).fit()
    ci = m.conf_int().loc["T:zmod"]
    return {"beta": m.params["T:zmod"], "ci_low": ci[0], "ci_high": ci[1],
            "p": m.pvalues["T:zmod"]}


# --------------------------------------------------------------- main ------


def main() -> None:
    out = Path("results")
    figs = Path("reports/figures")
    frame = pd.read_csv(FRAME, sep="\t", na_values="n/a")
    qc = qc_exclusions(load_qc())
    qc.drop(columns=[c for c in qc if c.startswith("fc")]).to_csv(
        out / "imaging_qc.csv", index=False, float_format="%.4f")
    base = frame[(frame.days_blood_to_scan > 0) & frame.avp_pg_ml.notna()]
    flow = [("Stage-1 eligible", len(base))]
    have = base.merge(qc, on="participant_id")
    flow.append(("preprocessed (qc.json present)", len(have)))
    for f in ("flag_fd", "flag_nkeep", "flag_cov"):
        flow.append((f"excluded by {f}", int(have[f].sum())))
    ana = have[have.qc_pass]
    flow.append(("analysed (primary)", len(ana)))
    flow.append(("  of which OT / PL", f"{(ana.treatment_group == 'OT').sum()} / "
                 f"{(ana.treatment_group == 'PL').sum()}"))
    flow.append(("flagged for normalization visual review", int(have.flag_norm_review.sum())))
    pd.DataFrame(flow, columns=["step", "n"]).to_csv(out / "imaging_flow.csv", index=False)
    qc_pass = qc[qc.qc_pass]
    elig = qc_pass[qc_pass.participant_id.isin(base.participant_id)]

    rows = []

    def add(label, family, fc_col, mod, qcd=elig, frame_=frame, transform="log",
            full=PRIMARY, reduced=REDUCED, n_perm=N_PERM, extra=None):
        d = prepare(frame_, qcd, fc_col, mod, transform)
        if extra is not None:
            d = extra(d)
        r = fit_interaction(d, full, reduced, n_perm=n_perm)
        r.update(label=label, family=family, fc=fc_col, moderator=mod)
        rows.append(r)
        return d

    # Primary + secondary family
    d_primary = add("H1 primary: T x AVP, amygdala-mPFC", "primary", "fc_amygdala_mpfc", "log_avp")
    add("S1: T x OT, amygdala-mPFC", "secondary", "fc_amygdala_mpfc", "log_ot")
    add("S2: T x AVP, amygdala-antInsula", "secondary", "fc_amygdala_antinsula", "log_avp")
    add("S3: T x OT, amygdala-antInsula", "secondary", "fc_amygdala_antinsula", "log_ot")
    # Sensitivities (primary estimand only)
    add("Sens1: GSR", "sensitivity", "fc_gsr_amygdala_mpfc", "log_avp")
    add("Sens2: raw AVP", "sensitivity", "fc_amygdala_mpfc", "log_avp", transform="raw")
    add("Sens3: + days & T x days", "sensitivity", "fc_amygdala_mpfc", "log_avp",
        full=PRIMARY + " + days_c + T:days_c", reduced=REDUCED + " + days_c + T:days_c",
        extra=lambda d: d.assign(days_c=(d.days_blood_to_scan - d.days_blood_to_scan.mean()) / 30))
    add("Sens4: cohort factor (T x older)", "sensitivity", "fc_amygdala_mpfc", "log_avp",
        full="FC ~ T + zmod + T:zmod + older + T:older + sex + fd_c",
        reduced="FC ~ T + zmod + older + T:older + sex + fd_c")
    add("Sens5: exclude TR=2.5s", "sensitivity", "fc_amygdala_mpfc", "log_avp",
        qcd=elig[elig.tr == 2.0])
    strict = elig[(elig.mean_fd <= 0.3) & (elig.n_keep_strict >= 150)]
    add("Sens6: strict motion (FD<=0.3; censor 0.3)", "sensitivity",
        "fc_strict_amygdala_mpfc", "log_avp", qcd=strict)
    add("Sens8: joint model (+ T x OT)", "sensitivity", "fc_amygdala_mpfc", "log_avp",
        full=PRIMARY + " + zot + T:zot", reduced=REDUCED + " + zot + T:zot",
        extra=lambda d: d.assign(zot=(d.log_ot - d.log_ot.mean()) / d.log_ot.std(ddof=1)))
    unk = frame[frame.days_blood_to_scan.isna() & frame.avp_pg_ml.notna()]
    add("Sens9: + unknown-interval participant", "sensitivity", "fc_amygdala_mpfc",
        "log_avp", qcd=qc_pass[qc_pass.participant_id.isin(
            pd.concat([base, unk]).participant_id)])
    res = pd.DataFrame(rows)
    sec = res.family == "secondary"
    res.loc[sec, "p_perm_holm"] = holm(res.loc[sec, "p_perm"].tolist())
    rob = robust_fit(d_primary)
    res = pd.concat([res, pd.DataFrame([{
        "label": "Sens7: Huber robust regression", "family": "sensitivity",
        "fc": "fc_amygdala_mpfc", "moderator": "log_avp", "n": len(d_primary),
        "beta": rob["beta"], "ci_low": rob["ci_low"], "ci_high": rob["ci_high"],
        "p_hc3": rob["p"]}])], ignore_index=True)
    cols = ["label", "family", "n", "n_OT", "n_PL", "beta", "ci_low", "ci_high",
            "p_hc3", "p_perm", "p_perm_holm", "beta_std", "sd_FC",
            "cond_T_at_z-1", "cond_T_at_z-1_lo", "cond_T_at_z-1_hi",
            "cond_T_at_z+0", "cond_T_at_z+0_lo", "cond_T_at_z+0_hi",
            "cond_T_at_z+1", "cond_T_at_z+1_lo", "cond_T_at_z+1_hi", "beta_T", "p_T"]
    res[cols].to_csv(out / "moderation_results.csv", index=False, float_format="%.5g")
    diag = diagnostics(d_primary)
    (out / "primary_diagnostics.json").write_text(json.dumps(diag, indent=2))

    # Replication anchor (descriptive)
    anc = []
    for fcol in ("fc_amygdala_mpfc", "fc_amygdala_antinsula"):
        m = smf.ols("FC ~ T + age_c + sex + fd_c", d_primary.assign(FC=d_primary[fcol])).fit(cov_type="HC3")
        m3 = smf.ols("FC ~ T * older * sex + fd_c", d_primary.assign(FC=d_primary[fcol])).fit(cov_type="HC3")
        anc.append({"fc": fcol, "beta_T": m.params["T"], "ci_low": m.conf_int().loc["T", 0],
                    "ci_high": m.conf_int().loc["T", 1], "p": m.pvalues["T"],
                    "TxOlderxSex_beta": m3.params["T:older:sex"],
                    "TxOlderxSex_p": m3.pvalues["T:older:sex"]})
        cell = d_primary.assign(FC=d_primary[fcol]).groupby(
            ["older", "sex", "treatment_group"]).FC.agg(["mean", "std", "count"])
        cell.to_csv(out / f"cell_means_{fcol}.csv")
    pd.DataFrame(anc).to_csv(out / "replication_anchor.csv", index=False, float_format="%.5g")
    motion = d_primary.groupby("treatment_group")[["mean_fd", "n_keep"]].describe().T
    motion.to_csv(out / "motion_by_arm.csv")
    plot_primary(d_primary, figs / "primary_interaction.png")
    plot_forest(res, figs / "forest_interactions.png")
    print(pd.DataFrame(flow).to_string(index=False))
    print(res[["label", "n", "beta", "ci_low", "ci_high", "p_hc3", "p_perm",
               "p_perm_holm"]].to_string(index=False))
    print(json.dumps(diag, indent=1))
    print(pd.DataFrame(anc).to_string(index=False))


def plot_primary(d: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = smf.ols(PRIMARY, d).fit(cov_type="HC3")
    colors = {0: "#2a78d6", 1: "#eb6834"}
    surf, ink, muted = "#fcfcfb", "#0b0b0b", "#52514e"
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=surf)
    zs = np.linspace(d.zmod.min(), d.zmod.max(), 60)
    for t, lab in ((0, "PL"), (1, "OT")):
        g = d[d["T"] == t]
        ax[0].scatter(g.zmod, g.FC, s=24, color=colors[t], edgecolors=surf, lw=1,
                      label=f"{lab} (n={len(g)})", zorder=3,
                      marker="o" if t == 0 else "s")
        nd = pd.DataFrame({"T": t, "zmod": zs, "age_c": 0.0, "sex": d.sex.mean(), "fd_c": 0.0})
        pr = m.get_prediction(nd).summary_frame()
        ax[0].plot(zs, pr["mean"], color=colors[t], lw=2)
        ax[0].fill_between(zs, pr.mean_ci_lower, pr.mean_ci_upper, color=colors[t], alpha=.15, lw=0)
    ax[0].set_xlabel("screening plasma AVP, z(ln pg/mL)", color=muted, fontsize=9)
    ax[0].set_ylabel("amygdala–mPFC FC (Fisher z)", color=muted, fontsize=9)
    ax[0].set_title("Predicted FC by arm (age, sex, FD at mean; 95% CI)", color=ink, fontsize=9.5, loc="left")
    ax[0].legend(frameon=False, fontsize=8)
    names = list(m.params.index)
    est, lo, hi = [], [], []
    for z in zs:
        L = np.zeros(len(names)); L[names.index("T")] = 1; L[names.index("T:zmod")] = z
        e = L @ m.params.to_numpy(); s = np.sqrt(L @ m.cov_params().to_numpy() @ L)
        est.append(e); lo.append(e - 1.96 * s); hi.append(e + 1.96 * s)
    ax[1].plot(zs, est, color=ink, lw=2)
    ax[1].fill_between(zs, lo, hi, color="#86b6ef", alpha=.35, lw=0)
    ax[1].axhline(0, color="#bdbcb6", lw=1)
    ax[1].set_xlabel("screening plasma AVP, z(ln pg/mL)", color=muted, fontsize=9)
    ax[1].set_ylabel("OT − PL difference in FC (Fisher z)", color=muted, fontsize=9)
    ax[1].set_title("Conditional treatment effect (95% CI, HC3)", color=ink, fontsize=9.5, loc="left")
    for a in ax:
        a.set_facecolor(surf)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(colors=muted, labelsize=8)
        a.grid(color="#e6e5e0", lw=.6, zorder=0)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=surf)
    plt.close(fig)


def plot_forest(res: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    surf, ink, muted = "#fcfcfb", "#0b0b0b", "#52514e"
    r = res.iloc[::-1].reset_index(drop=True)
    std = r.sd_FC.fillna(r.sd_FC.dropna().iloc[-1])
    fig, ax = plt.subplots(figsize=(8, 0.38 * len(r) + 1), facecolor=surf)
    ax.set_facecolor(surf)
    col = {"primary": "#eb6834", "secondary": "#2a78d6", "sensitivity": "#52514e"}
    for i, x in r.iterrows():
        ax.errorbar(x.beta, i, xerr=[[x.beta - x.ci_low], [x.ci_high - x.beta]],
                    fmt="o", color=col[x.family], ms=5, lw=1.6, capsize=0)
    ax.axvline(0, color="#bdbcb6", lw=1)
    ax.set_yticks(range(len(r)))
    ax.set_yticklabels([f"{l}  (n={int(n)})" for l, n in zip(r.label, r.n)], fontsize=7.5, color=ink)
    ax.set_xlabel("T × moderator coefficient (Fisher-z FC per SD of moderator), 95% CI",
                  color=muted, fontsize=8.5)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=muted, labelsize=8)
    ax.grid(axis="x", color="#e6e5e0", lw=.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=surf)
    plt.close(fig)


if __name__ == "__main__":
    main()
