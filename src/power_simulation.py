"""Sensitivity simulation for H2 under the registered model (plan v1.0 s8; OSF).
Observed design (N = 83: T, zAVP, C, S effect-coded, full T*C*S), HC3 two-sided
alpha .05, 2,000 simulated outcomes per effect size. Outcome-only simulation.
Usage: python -I src/power_simulation.py  -> results/power_simulation.csv"""
import numpy as np
import pandas as pd
import patsy

d = pd.read_csv("derivatives/phenotype_audit/analysis_frame.tsv", sep="\t", na_values="n/a")
d = d[(d.days_blood_to_scan > 0) & d.avp_pg_ml.notna() & d.ot_pg_ml.notna()].copy()
d["T"] = np.where(d.treatment_group == "OT", .5, -.5)
d["C"] = np.where(d.cohort == "older", .5, -.5)
d["S"] = np.where(d.sex == "F", .5, -.5)
d["zM"] = (d.log_avp - d.log_avp.mean()) / d.log_avp.std(ddof=1)
X = patsy.dmatrix("T * zM + T * C * S", d, return_type="dataframe")
j = list(X.columns).index("T:zM")
Xa = X.to_numpy(); XtXi = np.linalg.pinv(Xa.T @ Xa); A = XtXi @ Xa.T
h = np.einsum("ij,jk,ik->i", Xa, XtXi, Xa)
from scipy import stats
crit = stats.norm.ppf(.975)
rng = np.random.default_rng(20261108)
rows = []
for b in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.69, 0.7, 0.8]:
    hits = 0
    for _ in range(2000):
        y = b * Xa[:, j] + rng.standard_normal(len(d))
        e = y - Xa @ (A @ y)
        se = np.sqrt(np.sum(A[j] ** 2 * e ** 2 / (1 - h) ** 2))
        hits += abs((A[j] @ y) / se) > crit
    rows.append({"n": len(d), "beta_int_sd": b, "power": hits / 2000})
out = pd.DataFrame(rows); out.to_csv("results/power_simulation.csv", index=False)
print(out.to_string(index=False))
