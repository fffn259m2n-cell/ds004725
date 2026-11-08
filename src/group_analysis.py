"""Pre-registered group analysis (plan v1.0 sections 5-7; OSF osf.io/8yqun).

Two-phase, blinded:
  python -I src/group_analysis.py --qc-only
      reads ONLY *_qc.json + results/qc_visual.csv, applies exclusions, writes
      results/imaging_qc.csv, results/exclusions_by_cell.csv and freezes them
      (results/qc_frozen.json, SHA-256). No *_fc.json file is opened.
  python -I src/group_analysis.py
      refuses to run unless the frozen QC hash matches; then reads *_fc.json
      for included participants and fits H1, H2, S1-S3 and section-7 sensitivities.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

DER = Path("derivatives/prereg_v1")
FRAME = Path("derivatives/phenotype_audit/analysis_frame.tsv")
RES = Path("results")
FIG = Path("reports/figures")
N_PERM = 10_000
SEED = 20261108

# ---------------------------------------------------------------- QC -------


def qc_table() -> pd.DataFrame:
    rows = []
    for f in sorted(DER.glob("sub-*/sub-*_qc.json")):
        q = json.loads(f.read_text())
        r = {k: q.get(k) for k in ("participant_id", "status", "tr", "mean_fd",
                                   "minutes_retained", "prop_censored", "reg_r",
                                   "epi_vs_template_mask_dice", "signals_finite",
                                   "n_keep_fd03")}
        for k, v in (q.get("coverage") or {}).items():
            r[f"cov_{k}"] = v
        rows.append(r)
    return pd.DataFrame(rows)


def apply_exclusions(q: pd.DataFrame, visual: pd.DataFrame | None) -> pd.DataFrame:
    q = q.copy()
    ok = q.status.eq("ok")
    med = q.loc[ok, "reg_r"].median()
    mad = (q.loc[ok, "reg_r"] - med).abs().median()  # raw MAD (plan: "3 MAD")
    q["ex_failed"] = ~ok | ~q.signals_finite.fillna(False).astype(bool)
    q["ex_fd"] = q.mean_fd > 0.5
    q["ex_minutes"] = q.minutes_retained < 5
    q["ex_reg_metric"] = q.reg_r < med - 3 * mad
    if visual is not None:
        v = visual.set_index("participant_id")["visual_pass"].astype(str).str.lower()
        q["ex_reg_visual"] = q.participant_id.map(lambda s: v.get(s, "missing") != "true")
    else:
        q["ex_reg_visual"] = np.nan
    q["ex_cov_main"] = (q.cov_AMY < 0.5) | (q.cov_MPFC < 0.5)
    q["ex_cov_ains"] = q.cov_AINS < 0.5
    base = q.ex_failed | q.ex_fd | q.ex_minutes | q.ex_reg_metric | q.ex_reg_visual.fillna(False).astype(bool)
    q["pass_main"] = ~(base | q.ex_cov_main)
    q["pass_ains"] = ~(base | q.ex_cov_main | q.ex_cov_ains)
    q.attrs["reg_r_median"], q.attrs["reg_r_mad"] = float(med), float(mad)
    return q


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_qc_only() -> None:
    q = qc_table()
    vis_p = RES / "qc_visual.csv"
    if not vis_p.exists():
        pd.DataFrame({"participant_id": q.participant_id, "visual_pass": "",
                      "note": ""}).to_csv(vis_p, index=False)
        sys.exit(f"Wrote blank {vis_p}: complete the visual review, then rerun --qc-only.")
    vis = pd.read_csv(vis_p, dtype=str, keep_default_na=False)
    if (vis.visual_pass.str.lower().isin(["true", "false"]) == False).any():  # noqa: E712
        sys.exit("qc_visual.csv incomplete: every row needs visual_pass true/false.")
    q = apply_exclusions(q, vis)
    ph = pd.read_csv(FRAME, sep="\t", na_values="n/a")
    q = q.merge(ph[["participant_id", "treatment_group", "cohort", "sex"]], on="participant_id")
    q.to_csv(RES / "imaging_qc.csv", index=False, float_format="%.5f")
    cells = q.groupby(["treatment_group", "cohort", "sex"]).agg(
        n=("participant_id", "size"), pass_main=("pass_main", "sum"),
        pass_ains=("pass_ains", "sum")).reset_index()
    cells.to_csv(RES / "exclusions_by_cell.csv", index=False)
    frozen = {"frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
              "imaging_qc_sha256": sha(RES / "imaging_qc.csv"),
              "qc_visual_sha256": sha(vis_p),
              "reg_r_median": q.attrs.get("reg_r_median"), "reg_r_mad": q.attrs.get("reg_r_mad"),
              "n_processed": int(len(q)), "n_pass_main": int(q.pass_main.sum()),
              "n_pass_ains": int(q.pass_ains.sum()),
              "exclusion_counts": {c: int(q[c].fillna(False).astype(bool).sum())
                                   for c in q if c.startswith("ex_")},
              "fc_files_opened": 0}
    (RES / "qc_frozen.json").write_text(json.dumps(frozen, indent=2))
    print(json.dumps(frozen, indent=2))
    print(cells.to_string(index=False))


# ------------------------------------------------------------ models -------


def design(formula: str, d: pd.DataFrame):
    y, X = patsy.dmatrices(formula, d, return_type="dataframe")
    return y.iloc[:, 0].to_numpy(), X


def hc3_fit(y, X):
    Xa = X.to_numpy()
    XtXi = np.linalg.pinv(Xa.T @ Xa)
    b = XtXi @ Xa.T @ y
    e = y - Xa @ b
    h = np.einsum("ij,jk,ik->i", Xa, XtXi, Xa)
    A = XtXi @ Xa.T
    var = (A ** 2) @ (e ** 2 / (1 - h) ** 2)
    return b, np.sqrt(var), e, A, h


def freedman_lane(y, X: pd.DataFrame, term: str, n_perm=N_PERM, seed=SEED) -> float:
    """Two-sided Freedman-Lane p on the HC3 t statistic for column `term`."""
    j = list(X.columns).index(term)
    b, se, _, A, h = hc3_fit(y, X)
    t_obs = b[j] / se[j]
    Xr = X.drop(columns=term).to_numpy()
    br = np.linalg.lstsq(Xr, y, rcond=None)[0]
    fit_r, res_r = Xr @ br, y - Xr @ br
    rng = np.random.default_rng(seed)
    Xa, aj = X.to_numpy(), A[j]
    hit = 0
    for _ in range(n_perm):
        ys = fit_r + rng.permutation(res_r)
        bs = A @ ys
        e = ys - Xa @ bs
        se_j = np.sqrt(np.sum(aj ** 2 * e ** 2 / (1 - h) ** 2))
        hit += abs(bs[j] / se_j) >= abs(t_obs)
    return (hit + 1) / (n_perm + 1)


def bf01_bic(y, X: pd.DataFrame, term: str) -> float:
    full = sm.OLS(y, X).fit()
    red = sm.OLS(y, X.drop(columns=term)).fit()
    return float(np.exp((full.bic - red.bic) / 2))


def test_term(formula, d, term, perm=True) -> dict:
    y, X = design(formula, d)
    m = sm.OLS(y, X).fit(cov_type="HC3")
    ci = m.conf_int().loc[term]
    out = {"n": int(m.nobs), "term": term, "beta": float(m.params[term]),
           "se_hc3": float(m.bse[term]), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
           "p_hc3": float(m.pvalues[term]), "sd_FC": float(np.std(y, ddof=1)),
           "beta_std": float(m.params[term] / np.std(y, ddof=1))}
    if perm:
        out["p_perm"] = freedman_lane(y, X, term)
        out["bf01"] = bf01_bic(y, X, term)
    for k in ("T", "zM"):
        if k in m.params and k != term:
            out[f"beta_{k}"], out[f"p_{k}"] = float(m.params[k]), float(m.pvalues[k])
    return out


def holm(p):
    p = np.asarray(p, float)
    adj, run = np.empty_like(p), 0.0
    for r, i in enumerate(np.argsort(p)):
        run = max(run, (len(p) - r) * p[i])
        adj[i] = min(1.0, run)
    return adj.tolist()


def load_fc(ids, variant="primary") -> pd.DataFrame:
    rows = []
    for s in ids:
        fc = json.loads((DER / s / f"{s}_fc.json").read_text()).get(variant)
        if fc:
            rows.append({"participant_id": s, **fc})
    return pd.DataFrame(rows)


def coded(ph: pd.DataFrame, qc: pd.DataFrame, fc: pd.DataFrame, mod: str | None) -> pd.DataFrame:
    d = ph.merge(qc[["participant_id", "mean_fd", "prop_censored", "tr"]], on="participant_id") \
          .merge(fc, on="participant_id")
    d["T"] = np.where(d.treatment_group == "OT", 0.5, -0.5)
    d["C"] = np.where(d.cohort == "older", 0.5, -0.5)
    d["S"] = np.where(d.sex == "F", 0.5, -0.5)
    d["age_z"] = (d.age - d.age.mean()) / d.age.std(ddof=1)
    d["age_c"] = d.age - d.age.mean()
    d["sexF"] = (d.sex == "F").astype(float)
    if mod:
        d = d[d[mod].notna()].copy()
        d["zM"] = (d[mod] - d[mod].mean()) / d[mod].std(ddof=1)
    return d


H1_F = "FC ~ T * C * S + mean_fd"
H2_F = "FC ~ T * zM + T * C * S + mean_fd"


def run_models() -> None:
    fz = json.loads((RES / "qc_frozen.json").read_text())
    assert sha(RES / "imaging_qc.csv") == fz["imaging_qc_sha256"], "QC changed after freezing"
    qc = pd.read_csv(RES / "imaging_qc.csv")
    ph = pd.read_csv(FRAME, sep="\t", na_values="n/a")
    pre = ph[(ph.days_blood_to_scan > 0) & ph.avp_pg_ml.notna() & ph.ot_pg_ml.notna()].participant_id
    main_ids = qc.loc[qc.pass_main, "participant_id"]
    ains_ids = qc.loc[qc.pass_ains, "participant_id"]
    res, sens = [], []

    def fc_for(ids, variant="primary"):
        return load_fc(ids, variant)

    # ---------------- H1 (all passing QC)
    d1 = coded(ph, qc, fc_for(main_ids).rename(columns={"AMY_MPFC": "FC"}), None)
    for term in ("T", "T:C:S"):
        r = test_term(H1_F, d1, term)
        r.update(label=f"H1 {term}", family="H1", outcome="AMY_MPFC")
        res.append(r)
    cells = d1.groupby(["treatment_group", "cohort", "sex"]).agg(
        n=("FC", "size"), mean=("FC", "mean"), sd=("FC", "std"),
        mean_L=("AMYL_MPFC", "mean"), mean_R=("AMYR_MPFC", "mean")).reset_index()
    cells.to_csv(RES / "h1_cell_means.csv", index=False, float_format="%.4f")

    # ---------------- H2 + secondary (pre-dose draws)
    def h2_data(outcome, mod, ids, variant="primary", id_filter=pre):
        f = fc_for([i for i in ids if i in set(id_filter)], variant).rename(columns={outcome: "FC"})
        return coded(ph, qc, f, mod)

    fams = [("H2 T:zAVP AMY-MPFC", "H2", "AMY_MPFC", "log_avp", main_ids),
            ("S1 T:zOT AMY-MPFC", "secondary", "AMY_MPFC", "log_ot", main_ids),
            ("S2 T:zAVP AMY-AINS", "secondary", "AMY_AINS", "log_avp", ains_ids),
            ("S3 T:zOT AMY-AINS", "secondary", "AMY_AINS", "log_ot", ains_ids)]
    data = {}
    for lab, fam, outc, mod, ids in fams:
        d = h2_data(outc, mod, ids)
        data[lab] = d
        r = test_term(H2_F, d, "T:zM")
        r.update(label=lab, family=fam, outcome=outc, moderator=mod,
                 n_OT=int((d["T"] > 0).sum()), n_PL=int((d["T"] < 0).sum()))
        res.append(r)
    R = pd.DataFrame(res)
    sec = R.family == "secondary"
    R.loc[sec, "p_perm_holm"] = holm(R.loc[sec, "p_perm"])
    R.to_csv(RES / "registered_tests.csv", index=False, float_format="%.5g")

    # ---------------- Sensitivities (section 7), estimates + HC3 CI only
    d2 = data["H2 T:zAVP AMY-MPFC"]

    def S(label, formula, d, term="T:zM", applies="H2"):
        r = test_term(formula, d, term, perm=False)
        r.update(label=label, applies=applies)
        sens.append(r)

    for var, lab in (("gsr", "GSR"), ("nocensor", "No censoring"), ("fd03", "FD > 0.3 censoring")):
        S(f"{lab} (H1 T)", H1_F, coded(ph, qc, fc_for(main_ids, var).rename(columns={"AMY_MPFC": "FC"}), None), "T", "H1")
        S(f"{lab} (H2)", H2_F, h2_data("AMY_MPFC", "log_avp", main_ids, var))
    S("Draw <= 30 days before dosing", H2_F,
      h2_data("AMY_MPFC", "log_avp", main_ids, id_filter=ph[(ph.days_blood_to_scan > 0) & (ph.days_blood_to_scan <= 30)].participant_id))
    S("Draw interval as covariate", H2_F + " + days_blood_to_scan", d2)
    post = ph[ph.days_blood_to_scan.notna() & ph.avp_pg_ml.notna()].participant_id
    S("Including the two post-dose draws", H2_F, h2_data("AMY_MPFC", "log_avp", main_ids, id_filter=post))
    S("Excluding TR 2.5 s", H2_F, h2_data("AMY_MPFC", "log_avp", main_ids).query("tr == 2.0"))
    S("Continuous age instead of age group", "FC ~ T * zM + T * age_z * S + mean_fd", d2)
    S("Adding AVP x age group", H2_F + " + zM:C", d2)
    S("Adding proportion censored", H2_F + " + prop_censored", d2)
    S("Original handoff model (age, sex main effects)", "FC ~ T * zM + age_c + sexF + mean_fd", d2)
    y, X = design(H2_F, d2)
    hub = sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit()
    ci = hub.conf_int().loc["T:zM"]
    sens.append({"label": "Huber robust regression", "applies": "H2", "n": len(y),
                 "beta": float(hub.params["T:zM"]), "ci_low": float(ci[0]),
                 "ci_high": float(ci[1]), "p_hc3": float(hub.pvalues["T:zM"])})
    d2q = d2.assign(zM2=d2.zM ** 2)
    S("Quadratic moderation (T:zAVP^2)", H2_F + " + zM2 + T:zM2", d2q, term="T:zM2")
    pd.DataFrame(sens).to_csv(RES / "sensitivity_analyses.csv", index=False, float_format="%.5g")

    # ---------------- Diagnostics
    m = sm.OLS(y, X).fit()
    infl = m.get_influence()
    from scipy import stats
    from statsmodels.stats.diagnostic import het_breuschpagan
    diag = {"n": int(m.nobs), "max_leverage": float(infl.hat_matrix_diag.max()),
            "n_cooks_gt_4_over_n": int((infl.cooks_distance[0] > 4 / m.nobs).sum()),
            "max_cooks": float(infl.cooks_distance[0].max()),
            "shapiro_p": float(stats.shapiro(m.resid).pvalue),
            "breusch_pagan_p": float(het_breuschpagan(m.resid, X)[1]),
            "zAVP_range_OT": [float(d2.zM[d2["T"] > 0].min()), float(d2.zM[d2["T"] > 0].max())],
            "zAVP_range_PL": [float(d2.zM[d2["T"] < 0].min()), float(d2.zM[d2["T"] < 0].max())],
            "motion_by_arm": d2.groupby("treatment_group").mean_fd.describe().round(4).to_dict()}
    (RES / "h2_diagnostics.json").write_text(json.dumps(diag, indent=2))
    plot_conditional(d2, FIG / "h2_conditional_effect.png")
    plot_forest(R, pd.DataFrame(sens), FIG / "forest_registered_and_sensitivity.png")
    print(R[["label", "n", "beta", "ci_low", "ci_high", "p_hc3", "p_perm", "p_perm_holm", "bf01"]].to_string(index=False))
    print(pd.DataFrame(sens)[["label", "n", "beta", "ci_low", "ci_high", "p_hc3"]].to_string(index=False))
    print(json.dumps(diag, indent=1))
    print(cells.to_string(index=False))


# ------------------------------------------------------------- plots -------
SURF, INK, MUTED, BLUE, ORANGE = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6", "#eb6834"


def _style(a):
    a.set_facecolor(SURF)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.tick_params(colors=MUTED, labelsize=8)
    a.grid(color="#e6e5e0", lw=.6, zorder=0)


def plot_conditional(d, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    y, X = design(H2_F, d)
    m = sm.OLS(y, X).fit(cov_type="HC3")
    names = list(m.params.index)
    zs = np.linspace(d.zM.min(), d.zM.max(), 80)
    est, lo, hi = [], [], []
    for z in zs:
        L = np.zeros(len(names)); L[names.index("T")] = 1; L[names.index("T:zM")] = z
        e = L @ m.params.to_numpy(); s = np.sqrt(L @ m.cov_params().to_numpy() @ L)
        est.append(e); lo.append(e - 1.96 * s); hi.append(e + 1.96 * s)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURF)
    for t, lab, c, mk in ((-0.5, "Placebo", BLUE, "o"), (0.5, "Oxytocin", ORANGE, "s")):
        g = d[d["T"] == t]
        ax[0].scatter(g.zM, g.FC, s=26, color=c, marker=mk, edgecolors=SURF, lw=1,
                      label=f"{lab} (n={len(g)})", zorder=3)
    ax[0].set_xlabel("screening plasma AVP immunoreactivity, z(ln)", color=MUTED, fontsize=9)
    ax[0].set_ylabel("amygdala–mPFC FC (Fisher z)", color=MUTED, fontsize=9)
    ax[0].set_title("Observed FC by arm", color=INK, fontsize=9.5, loc="left")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].fill_between(zs, lo, hi, color="#86b6ef", alpha=.35, lw=0)
    ax[1].plot(zs, est, color=INK, lw=2)
    ax[1].axhline(0, color="#bdbcb6", lw=1)
    ax[1].set_xlabel("screening plasma AVP immunoreactivity, z(ln)", color=MUTED, fontsize=9)
    ax[1].set_ylabel("OT − placebo difference in FC", color=MUTED, fontsize=9)
    ax[1].set_title("Conditional treatment effect b1 + b3·zAVP (95% CI, HC3)", color=INK, fontsize=9.5, loc="left")
    for a in ax:
        _style(a)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_forest(R, Sx, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [(r.label, r.beta, r.ci_low, r.ci_high, r.n, ORANGE if r.family in ("H2",) else BLUE if r.family == "secondary" else INK)
            for r in R.itertuples() if r.family != "H1"]
    rows += [(f"  sens: {r.label}", r.beta, r.ci_low, r.ci_high, r.n, MUTED)
             for r in Sx.itertuples() if r.applies == "H2"]
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(rows) + 1.2), facecolor=SURF)
    for i, (lab, b, lo, hi, n, c) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=c, lw=1.6)
        ax.plot(b, i, "o", color=c, ms=5)
    ax.axvline(0, color="#bdbcb6", lw=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{l} (n={int(n)})" for l, _, _, _, n, _ in rows], fontsize=7.5, color=INK)
    ax.set_xlabel("interaction coefficient (Fisher-z FC per SD of moderator), 95% CI", color=MUTED, fontsize=8.5)
    _style(ax)
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-only", action="store_true")
    a = ap.parse_args()
    run_qc_only() if a.qc_only else run_models()
