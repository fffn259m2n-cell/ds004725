"""Sensitivity/power simulation for the Treatment x z(log AVP) interaction.

Uses the OBSERVED eligible design matrix (arm, z(log AVP), centred age, sex)
from the Stage 1 analysis frame; simulates only the outcome. No FC data used.
Effect sizes are in outcome-SD units per 1 SD of log AVP (difference in
OT-vs-PL slopes). Usage: python -I src/power_sim.py [--moderator log_avp]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def design(df, mod):
    d = df[(df.days_blood_to_scan > 0) & df[mod].notna() & df.bold_present]
    t = (d.treatment_group == "OT").astype(float).to_numpy()
    z = ((d[mod] - d[mod].mean()) / d[mod].std(ddof=1)).to_numpy()
    age = (d.age - d.age.mean()).to_numpy() / 10
    sex = (d.sex == "F").astype(float).to_numpy()
    X = np.column_stack([np.ones(len(d)), t, z, t * z, age, sex])
    return X, len(d)


def power(X, beta_int, n_sim=4000, alpha=0.05, seed=1):
    from scipy import stats
    rng = np.random.default_rng(seed)
    n, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    se_fac = np.sqrt(XtX_inv[3, 3])
    crit = stats.t.ppf(1 - alpha / 2, n - k)
    mu = X[:, 3] * beta_int
    hits = 0
    for _ in range(n_sim):
        y = mu + rng.standard_normal(n)
        b = XtX_inv @ X.T @ y
        resid = y - X @ b
        s2 = resid @ resid / (n - k)
        hits += abs(b[3] / (np.sqrt(s2) * se_fac)) > crit
    return hits / n_sim


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", default="derivatives/phenotype_audit/analysis_frame.tsv")
    ap.add_argument("--out", default="results/power_interaction.csv")
    a = ap.parse_args()
    df = pd.read_csv(a.frame, sep="\t", na_values="n/a")
    rows = []
    for mod in ["log_avp", "log_ot"]:
        X, n = design(df, mod)
        for b in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
            rows.append({"moderator": mod, "n": n, "beta_int_sd": b,
                         "power": power(X, b)})
    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(exist_ok=True)
    out.to_csv(a.out, index=False)
    print(out.to_string(index=False))
