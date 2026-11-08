"""Stage 1 phenotype audit for OpenNeuro ds004725.

Reads every TSV as strings first (BIDS ``n/a`` -> missing), validates the
participant join, parses the free-text / clock-time fields, and writes
de-identified summary outputs. No participant IDs are written to figures, and
row-level outputs stay under ``derivatives/`` (git-ignored for anything bulky).

Usage (from the dataset root):
    python -I src/audit_phenotype.py --bids-root . --out-root .
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

NA_TOKENS = {"n/a", "N/A", "NA", ""}
VALID_ARMS = {"OT", "PL"}
STANDARD_TR = 2.0

# ---------------------------------------------------------------- parsing ---


def read_tsv_as_strings(path: Path) -> pd.DataFrame:
    """Read a BIDS TSV with every cell as stripped str; CRLF tolerated."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     na_filter=False, encoding="utf-8")
    df.columns = [c.strip() for c in df.columns]
    return df.apply(lambda s: s.str.strip())


def to_missing(s: pd.Series) -> pd.Series:
    return s.where(~s.isin(NA_TOKENS), other=pd.NA)


def to_float(s: pd.Series) -> pd.Series:
    return pd.to_numeric(to_missing(s), errors="raise")


_DAYS_RE = re.compile(r"^(\d+)\s+days?\s+(before|after)$", re.IGNORECASE)


def parse_days_before(value) -> float:
    """Days the blood draw PRECEDED dosing/scan.

    '27 days before' -> 27.0 ; '5 days after' -> -5.0 (draw was
    POST-treatment) ; n/a -> NaN ; anything else raises.
    """
    if value is pd.NA or value is None or str(value).strip() in NA_TOKENS:
        return np.nan
    m = _DAYS_RE.match(str(value).strip())
    if not m:
        raise ValueError(f"Unrecognised sampling-interval value: {value!r}")
    d = float(m.group(1))
    return d if m.group(2).lower() == "before" else -d


def parse_clock_hmm(value) -> tuple[float, str]:
    """Parse the dataset's 'H.MM' clock strings into decimal hours.

    The dictionary says "hours and minutes"; values like '9.58' cannot be
    decimal hours of a ~8 AM protocol consistently, and minutes never exceed
    59 in the 2-digit values, so 'H.MM' is assumed. Single-digit fractions
    ('10.1') are read as trailing-zero-dropped minutes (10:10), because '10.01'
    also occurs in the same column. Returns (decimal_hours, flag).
    Flags: 'ok', 'single_digit_minutes', 'implausible', 'missing'.
    """
    if value is pd.NA or value is None or str(value).strip() in NA_TOKENS:
        return np.nan, "missing"
    s = str(value).strip()
    if not re.fullmatch(r"\d{1,2}(\.\d{1,2})?", s):
        raise ValueError(f"Unrecognised clock value: {value!r}")
    hh, _, mm = s.partition(".")
    flag = "ok"
    if len(mm) == 1:
        mm, flag = mm + "0", "single_digit_minutes"
    h, m = int(hh), int(mm or 0)
    if m > 59 or h > 23:
        return np.nan, "implausible"
    dec = h + m / 60.0
    # Protocol: in-person sessions ~08:00; anything before 06:00 is suspect
    if dec < 6.0:
        return np.nan, "implausible"
    return dec, flag


# Allele -> APOE isoform. rs429358 C => e4 ; rs7412 T => e2 ; else e3.
def apoe_from_snps(rs429358: str, rs7412: str) -> str | float:
    if pd.isna(rs429358) or pd.isna(rs7412):
        return np.nan
    a = rs429358.split()[-1].split("/")
    b = rs7412.split()[-1].split("/")
    if len(a) != 2 or len(b) != 2:
        return np.nan
    c4 = a.count("C")
    t2 = b.count("T")
    if c4 + t2 > 2:
        return "ambiguous"
    alleles = ["e4"] * c4 + ["e2"] * t2 + ["e3"] * (2 - c4 - t2)
    return "".join(sorted(alleles))


# --------------------------------------------------------------- auditing ---


def validate_ids(parts: pd.DataFrame, pheno: pd.DataFrame,
                 sub_dirs: set[str]) -> dict:
    p_ids, q_ids = parts["participant_id"], pheno["participant_id"]
    return {
        "n_participants_tsv": int(len(p_ids)),
        "n_phenotype_tsv": int(len(q_ids)),
        "n_subject_dirs": int(len(sub_dirs)),
        "dup_participants_tsv": int(p_ids.duplicated().sum()),
        "dup_phenotype_tsv": int(q_ids.duplicated().sum()),
        "in_participants_not_phenotype": sorted(set(p_ids) - set(q_ids)),
        "in_phenotype_not_participants": sorted(set(q_ids) - set(p_ids)),
        "in_participants_no_dir": sorted(set(p_ids) - sub_dirs),
        "dir_not_in_participants": sorted(sub_dirs - set(p_ids)),
        "bad_id_format": sorted(i for i in p_ids
                                if not re.fullmatch(r"sub-\d{5}", i)),
    }


def annexed_size(path: Path) -> int | None:
    """File size from a git-annex SHA256E key (symlink target), else stat."""
    if path.is_symlink():
        m = re.search(r"-s(\d+)--", str(path.readlink()))
        if m:
            return int(m.group(1))
    return path.stat().st_size if path.exists() else None


def imaging_inventory(bids_root: Path, ids: list[str]) -> pd.DataFrame:
    root_tr = json.loads((bids_root / "task-rest_bold.json").read_text())[
        "RepetitionTime"]
    rows = []
    for sid in ids:
        bold = bids_root / sid / "func" / f"{sid}_task-rest_bold.nii"
        t1 = bids_root / sid / "anat" / f"{sid}_T1w.nii"
        side = bold.with_suffix("").with_suffix(".json")
        side = bids_root / sid / "func" / f"{sid}_task-rest_bold.json"
        tr = (json.loads(side.read_text())["RepetitionTime"]
              if side.exists() else root_tr)
        rows.append({
            "participant_id": sid,
            "bold_present": bold.is_symlink() or bold.exists(),
            "t1_present": t1.is_symlink() or t1.exists(),
            "bold_bytes": annexed_size(bold),
            "t1_bytes": annexed_size(t1),
            "bold_content_local": bold.exists(),  # False => not yet fetched
            "tr_s": float(tr),
            "subject_sidecar": side.exists(),
        })
    return pd.DataFrame(rows)


def build_analysis_frame(parts: pd.DataFrame, pheno: pd.DataFrame,
                         img: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    df = parts.merge(pheno, on="participant_id", how="outer",
                     validate="one_to_one", indicator=True)
    df = df.merge(img, on="participant_id", how="left", validate="one_to_one")
    flags: dict = {}
    df["age"] = to_float(df["age"])
    df["bmi"] = to_float(df["bmi"])
    df["years_of_education"] = to_float(df["years_of_education"])
    df["sex"] = to_missing(df["sex"])
    df["treatment_group"] = to_missing(df["treatment_group"])
    df["avp_pg_ml"] = to_float(df["BLOOD_VASOPRESSIN_LEVELS"])
    df["ot_pg_ml"] = to_float(df["BLOOD_OXYTOCIN_LEVELS"])
    df["days_blood_to_scan"] = df[
        "TIMING_BIOLOGICAL_SAMPLING_RELATIVE_TO_TREATMENT"].map(
            parse_days_before)
    for src, dst in [("TIME_OF_DAY_BLOOD_SAMPLING", "blood_clock_h"),
                     ("TIME_OF_DAY_SALIVARY_SAMPLING", "saliva_clock_h"),
                     ("TIME_OF_DAY_RESTING_SCAN", "rest_clock_h")]:
        parsed = df[src].map(parse_clock_hmm)
        df[dst] = parsed.map(lambda t: t[0])
        df[dst + "_flag"] = parsed.map(lambda t: t[1])
        flags[dst] = df[dst + "_flag"].value_counts().to_dict()
    # Age cohort: the design is bimodal (18-31 vs 63-81)
    df["cohort"] = np.select([df["age"] < 40, df["age"] >= 60],
                             ["younger", "older"], default="unclassified")
    # Sub-ID prefix convention (inferred, verified below, not relied on)
    df["id_digit1"] = df["participant_id"].str[4]
    df["id_digit2"] = df["participant_id"].str[5]
    df["apoe_recomputed"] = [
        apoe_from_snps(a, b) for a, b in zip(
            to_missing(df["SALIVARY_SNP_RS429358"]),
            to_missing(df["SALIVARY_SNP_RS7412"]))]
    df["apoe_reported"] = to_missing(df["SALIVARY_APOE_GENOTYPE"])
    df["log_avp"] = np.log(df["avp_pg_ml"])
    df["log_ot"] = np.log(df["ot_pg_ml"])
    df["med_at_scan"] = to_float(df["MEDICATION_STATUS_AT_SCAN"])
    flags["contraception_codes"] = to_missing(
        df["CONTRACEPTION_USE_AT_SCAN"]).value_counts(dropna=False).to_dict()
    flags["merge_indicator"] = df["_merge"].value_counts().to_dict()
    return df.drop(columns="_merge"), flags


def missingness_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = {"avp_pg_ml": "plasma AVP", "ot_pg_ml": "plasma OT",
            "days_blood_to_scan": "blood-to-scan days",
            "blood_clock_h": "blood draw clock time (missing/implausible)",
            "rest_clock_h": "rest scan clock time (missing/implausible)",
            "apoe_reported": "APOE genotype", "bmi": "BMI",
            "years_of_education": "education"}
    rows = []
    groups = [("all", df)] + [
        (f"{a}|{c}", g) for (a, c), g in df.groupby(
            ["treatment_group", "cohort"])] + [
        (a, g) for a, g in df.groupby("treatment_group")]
    for gname, g in groups:
        for col, label in cols.items():
            rows.append({"group": gname, "variable": label,
                         "n": len(g), "n_missing": int(g[col].isna().sum()),
                         "pct_missing": round(100 * g[col].isna().mean(), 1)})
    return pd.DataFrame(rows)


def eligibility_flow(df: pd.DataFrame) -> pd.DataFrame:
    steps = [
        ("listed in participants.tsv", pd.Series(True, index=df.index)),
        ("valid treatment code (OT/PL)", df["treatment_group"].isin(VALID_ARMS)),
        ("in phenotype table", df["TICS_SCORE"].notna()),
        ("non-missing plasma AVP", df["avp_pg_ml"].notna()),
        ("non-missing plasma OT", df["ot_pg_ml"].notna()),
        ("blood draw known to PRECEDE dosing (days > 0)",
         df["days_blood_to_scan"] > 0),
        ("resting BOLD file in dataset", df["bold_present"].fillna(False)),
        ("T1w file in dataset", df["t1_present"].fillna(False)),
        ("age and sex non-missing", df["age"].notna() & df["sex"].notna()),
    ]
    keep = pd.Series(True, index=df.index)
    rows = []
    for name, cond in steps:
        keep = keep & cond.astype(bool)
        sub = df[keep]
        rows.append({"step": name, "n_remaining": int(keep.sum()),
                     "n_OT": int((sub["treatment_group"] == "OT").sum()),
                     "n_PL": int((sub["treatment_group"] == "PL").sum()),
                     "n_younger": int((sub["cohort"] == "younger").sum()),
                     "n_older": int((sub["cohort"] == "older").sum())})
    # Sensitivity filter (not a primary exclusion): standard TR only
    std = keep & (df["tr_s"] == STANDARD_TR)
    sub = df[std]
    rows.append({"step": f"[sensitivity] TR == {STANDARD_TR}s",
                 "n_remaining": int(std.sum()),
                 "n_OT": int((sub["treatment_group"] == "OT").sum()),
                 "n_PL": int((sub["treatment_group"] == "PL").sum()),
                 "n_younger": int((sub["cohort"] == "younger").sum()),
                 "n_older": int((sub["cohort"] == "older").sum())})
    return pd.DataFrame(rows)


def smd(a: pd.Series, b: pd.Series) -> float:
    a, b = a.dropna(), b.dropna()
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / sp) if sp > 0 else np.nan


def describe_biomarkers(df: pd.DataFrame) -> pd.DataFrame:
    from scipy import stats
    rows = []
    for col in ["avp_pg_ml", "ot_pg_ml", "log_avp", "log_ot",
                "days_blood_to_scan", "blood_clock_h", "rest_clock_h", "age"]:
        for gname, g in [("all", df)] + list(df.groupby("treatment_group")) \
                + list(df.groupby("cohort")):
            x = g[col].dropna()
            if len(x) == 0:
                continue
            rows.append({
                "variable": col, "group": gname, "n": len(x),
                "mean": x.mean(), "sd": x.std(ddof=1), "min": x.min(),
                "q25": x.quantile(.25), "median": x.median(),
                "q75": x.quantile(.75), "max": x.max(),
                "skew": stats.skew(x, bias=False) if len(x) > 2 else np.nan,
                "n_unique": x.nunique(),
                "n_beyond_3sd": int((np.abs((x - x.mean()) / x.std()) > 3).sum()),
            })
    return pd.DataFrame(rows)


def balance_table(df: pd.DataFrame) -> pd.DataFrame:
    """Descriptive randomisation-balance check (SMD OT - PL). No p-values:
    covariate choice must not be driven by significance tests."""
    ot, pl = df[df.treatment_group == "OT"], df[df.treatment_group == "PL"]
    rows = []
    for col in ["age", "log_avp", "log_ot", "days_blood_to_scan",
                "blood_clock_h", "rest_clock_h", "bmi", "years_of_education"]:
        rows.append({"variable": col, "mean_OT": ot[col].mean(),
                     "mean_PL": pl[col].mean(), "smd_OT_minus_PL": smd(
                         ot[col], pl[col])})
    for lvl in ["F"]:
        rows.append({"variable": "prop_female", "mean_OT": (ot.sex == lvl).mean(),
                     "mean_PL": (pl.sex == lvl).mean(), "smd_OT_minus_PL": np.nan})
    rows.append({"variable": "prop_older",
                 "mean_OT": (ot.cohort == "older").mean(),
                 "mean_PL": (pl.cohort == "older").mean(),
                 "smd_OT_minus_PL": np.nan})
    return pd.DataFrame(rows)


def plot_biomarkers(df: pd.DataFrame, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"PL": "#2a78d6", "OT": "#eb6834"}
    ink, muted, surf = "#0b0b0b", "#52514e", "#fcfcfb"
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), facecolor=surf)
    panels = [("log_avp", "log plasma AVP (ln pg/mL)"),
              ("log_ot", "log plasma OT, unextracted (ln pg/mL)"),
              ("days_blood_to_scan", "days from blood draw to dosing/scan"),
              ("age", "age (years)")]
    rng = np.random.default_rng(0)
    for ax, (col, label) in zip(axes.ravel(), panels):
        ax.set_facecolor(surf)
        for i, arm in enumerate(["PL", "OT"]):
            for j, coh in enumerate(["younger", "older"]):
                x = df.loc[(df.treatment_group == arm) & (df.cohort == coh), col].dropna()
                pos = j * 2.6 + i
                ax.scatter(np.full(len(x), pos) + rng.uniform(-.18, .18, len(x)),
                           x, s=22, color=colors[arm], alpha=.85,
                           edgecolors=surf, linewidths=1.0, zorder=3)
                if len(x):
                    ax.hlines(x.median(), pos - .3, pos + .3, color=ink, lw=2, zorder=4)
        ax.set_xticks([0, 1, 2.6, 3.6])
        ax.set_xticklabels(["PL\nyounger", "OT\nyounger", "PL\nolder", "OT\nolder"],
                           color=muted, fontsize=8)
        ax.set_title(label, color=ink, fontsize=9.5, loc="left")
        ax.tick_params(colors=muted, labelsize=8)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        for s in ["left", "bottom"]:
            ax.spines[s].set_color("#d0cfca")
        ax.grid(axis="y", color="#e6e5e0", lw=.6, zorder=0)
    fig.suptitle("ds004725 screening biomarkers by randomized arm and age cohort "
                 "(points = participants, bar = median)", color=ink, fontsize=10.5,
                 x=.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(out_png, dpi=150, facecolor=surf)
    plt.close(fig)


def plot_avp_ot(df: pd.DataFrame, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"PL": "#2a78d6", "OT": "#eb6834"}
    surf, ink, muted = "#fcfcfb", "#0b0b0b", "#52514e"
    fig, ax = plt.subplots(figsize=(5.5, 4.6), facecolor=surf)
    ax.set_facecolor(surf)
    for arm, mk in [("PL", "o"), ("OT", "s")]:
        g = df[df.treatment_group == arm]
        ax.scatter(g.log_ot, g.log_avp, s=30, marker=mk, color=colors[arm],
                   edgecolors=surf, linewidths=1, label=f"{arm} arm", zorder=3)
    ax.set_xlabel("log plasma OT (unextracted)", color=muted, fontsize=9)
    ax.set_ylabel("log plasma AVP", color=muted, fontsize=9)
    ax.set_title("Screening plasma OT vs AVP", color=ink, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=ink)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=muted, labelsize=8)
    ax.grid(color="#e6e5e0", lw=.6, zorder=0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=surf)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bids-root", type=Path, default=Path("."))
    ap.add_argument("--out-root", type=Path, default=Path("."))
    a = ap.parse_args(argv)
    root, out = a.bids_root.resolve(), a.out_root.resolve()
    (out / "results").mkdir(exist_ok=True)
    (out / "reports" / "figures").mkdir(parents=True, exist_ok=True)

    parts = read_tsv_as_strings(root / "participants.tsv")
    pheno_dir = [d for d in root.iterdir()
                 if d.is_dir() and d.name.lower() in {"phenotype", "phenotypes"}]
    if len(pheno_dir) != 1:
        raise SystemExit(f"Expected one phenotype dir, found {pheno_dir}")
    pheno_tsvs = sorted(pheno_dir[0].glob("*.tsv"))
    if len(pheno_tsvs) != 1:
        raise SystemExit(f"Expected one phenotype TSV, found {pheno_tsvs}")
    pheno = read_tsv_as_strings(pheno_tsvs[0])
    sub_dirs = {d.name for d in root.glob("sub-*") if d.is_dir()}

    ids = validate_ids(parts, pheno, sub_dirs)
    img = imaging_inventory(root, sorted(sub_dirs))
    df, flags = build_analysis_frame(parts, pheno, img)

    bad_arm = sorted(set(df["treatment_group"].dropna()) - VALID_ARMS)
    apoe_ok = (df["apoe_recomputed"].dropna()
               == df.loc[df["apoe_recomputed"].notna(), "apoe_reported"]).all()
    prefix = pd.crosstab([df.id_digit1, df.id_digit2], [df.cohort, df.sex])

    missingness_table(df).to_csv(out / "results/phenotype_missingness.csv",
                                 index=False)
    flow = eligibility_flow(df)
    flow.to_csv(out / "results/eligibility_flow.csv", index=False)
    desc = describe_biomarkers(df)
    desc.round(4).to_csv(out / "results/biomarker_descriptives.csv", index=False)
    bal = balance_table(df)
    bal.round(4).to_csv(out / "results/arm_balance.csv", index=False)
    pd.crosstab([df.treatment_group, df.cohort], df.sex, margins=True).to_csv(
        out / "results/design_cells.csv")
    img_summary = {
        "n_bold_symlinks": int(img.bold_present.sum()),
        "n_t1_symlinks": int(img.t1_present.sum()),
        "n_bold_content_local": int(img.bold_content_local.sum()),
        "total_bold_GB": round(img.bold_bytes.sum() / 1e9, 3),
        "total_t1_GB": round(img.t1_bytes.sum() / 1e9, 3),
        "tr_counts": {str(k): int(v) for k, v in img.tr_s.value_counts().items()},
    }
    from scipy import stats
    cc = df[["log_avp", "log_ot"]].dropna()
    rho = stats.spearmanr(cc.log_avp, cc.log_ot)
    plot_biomarkers(df, out / "reports/figures/biomarkers_by_arm_cohort.png")
    plot_avp_ot(df, out / "reports/figures/avp_vs_ot.png")

    audit = {
        "id_validation": ids, "invalid_treatment_codes": bad_arm,
        "treatment_counts": df.treatment_group.value_counts().to_dict(),
        "cohort_counts": df.cohort.value_counts().to_dict(),
        "flags": flags, "apoe_reported_matches_snps": bool(apoe_ok),
        "apoe_reported_levels": df.apoe_reported.value_counts().to_dict(),
        "imaging": img_summary,
        "spearman_logAVP_logOT": {"rho": float(rho.statistic),
                                  "p": float(rho.pvalue), "n": len(cc)},
        "id_prefix_vs_cohort_sex": {
            "|".join(map(str, k)): {"|".join(map(str, c)): int(v)
                                    for c, v in row.items()}
            for k, row in prefix.to_dict(orient="index").items()},
        "missing_biomarker_rows": df.loc[df.avp_pg_ml.isna() | df.ot_pg_ml.isna(),
                                         ["treatment_group", "cohort"]
                                         ].to_dict(orient="records"),
        "implausible_clock_values": {
            "blood": df.loc[df.blood_clock_h_flag == "implausible",
                            "TIME_OF_DAY_BLOOD_SAMPLING"].tolist(),
            "rest": df.loc[df.rest_clock_h_flag == "implausible",
                           "TIME_OF_DAY_RESTING_SCAN"].tolist()},
        "blood_draw_after_treatment": df.loc[df.days_blood_to_scan <= 0,
            ["treatment_group", "cohort", "days_blood_to_scan"]
            ].to_dict(orient="records"),
        "blood_interval_missing": int(df.days_blood_to_scan.isna().sum()),
    }
    (out / "results/phenotype_audit.json").write_text(
        json.dumps(audit, indent=2, default=str))
    # Row-level analysis frame (contains IDs; stays in derivatives/)
    dpath = out / "derivatives" / "phenotype_audit"
    dpath.mkdir(parents=True, exist_ok=True)
    keep = ["participant_id", "age", "cohort", "sex", "treatment_group",
            "avp_pg_ml", "ot_pg_ml", "log_avp", "log_ot", "days_blood_to_scan",
            "blood_clock_h", "blood_clock_h_flag", "rest_clock_h", "bmi",
            "years_of_education", "med_at_scan", "apoe_reported",
            "bold_present", "t1_present", "tr_s"]
    df[keep].to_csv(dpath / "analysis_frame.tsv", sep="\t", index=False,
                    na_rep="n/a")
    print(json.dumps(audit, indent=2, default=str))
    print(flow.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
