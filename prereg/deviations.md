# Deviations log: plan v1.0 / OSF osf.io/8yqun

Governing documents:
- OSF registration osf.io/8yqun, registered 2026-11-08 03:30 Beirut time (**01:30 UTC**).
- Attached plan `prereg/analysis_plan_v1.0.md` (sha256 `a09d703f2f6b8598c4dd94656be4639c0561891724d17e1dd8aeba0843e751cc`). It was added to this repository on 2026-11-08 after registration, as a copy of the file supplied by the author.

Every entry records whether connectivity (FC) values had been seen at the time. **For every entry below, no FC value from any pipeline had been opened.**

| # | Date (UTC) | Item | Deviation / clarification | Reason | FC seen? |
|---|---|---|---|---|---|
| D1 | 11-08 | Foreknowledge (plan line 5; OSF "Foreknowledge" field) | The pilot participant's BOLD and T1w were downloaded into this analysis environment at 00:37 UTC, **before** registration at 01:30 UTC. Headers, motion, coverage and registration QC were viewed. A different (unregistered) pipeline then processed **25 participants**, from 00:56 until it was stopped at about 02:15 UTC. Those outputs contained FC values that were never opened, and they were pushed to GitHub (commits `9d7b5bf`, `0f954d8`, `70ff33f`, `f4ee629`). They live in `derivatives/oxtfc/` and are **excluded from all registered analyses**. Every registered result is recomputed from raw data under the registered pipeline | Analysis started before the author's registration existed. Disclosed | No |
| D2 | 11-08 | Dataset version | Uses tag 1.0.2 (commit `8372e9f`) instead of 1.0.1 (`d40af5c`). `git diff 1.0.1 1.0.2` touches only README and dataset_description.json; all imaging and phenotype files are identical | Fork was at 1.0.2 | No |
| D3 | 11-08 | Git tag `prereg-v1.0` | Not present in this repository. The plan was written elsewhere. The OSF registration is the timestamp of record | Provenance | No |
| D4 | 11-08 | "MNI152 WM probability" / "lateral ventricles at 50%" / "excluding grey matter" (plan §3) | Implemented with the TemplateFlow MNI152NLin6Asym HOSPA probabilistic atlas. Stored values are probability/100 and were rescaled ×100 (verified: channel maxima ≈ 1, voxel sums = 1). WM = L+R cerebral WM > 0.9, eroded once. CSF = L+R lateral ventricle > 0.5, minus HOSPA th25 grey-matter labels (cortex and subcortical nuclei) | The plan does not name a source file | No |
| D5 | 11-08 | HOCPA labels | TemplateFlow's HOCPA TSV indices are 0-based; image label = index + 1 (Insular Cortex = 2), verified by centroid and covered by a unit test | Atlas metadata quirk | No |
| D6 | 11-08 | Registration target | "MNI152 2 mm template" = tpl-MNI152NLin6Asym_res-02_T1w (whole head). ANTs `SyN` (rigid + affine + SyN), Mattes MI for linear and SyN stages, random_seed = 1. Mean EPI N4-corrected before registration | Specification detail | No |
| D7 | 11-08 | "EPI mask" | nilearn `compute_epi_mask` on the native mean EPI | Specification detail | No |
| D8 | 11-08 | GSR sensitivity | Global signal = mean within the native EPI mask, plus its derivative | Specification detail | No |
| D9 | 11-08 | Registration QC metric | `reg_r` (warped mean EPI vs T1w template inside the template brain mask) is low by construction (pilot r = 0.06) because of the T2* vs T1 contrast. Applied exactly as registered (median − 3 × raw MAD), with the visual check as the main safeguard. The pilot's native EPI vs inverse-warped template mask Dice = 0.95 is logged as an extra QC field (not used for exclusion) | Disclosed so readers can judge the metric | No |
| D10 | 11-08 | Voxel size | NIfTI headers give 3.15 × 3.15 × 3.5 mm, not "3 mm isotropic" | Data fact | No |
| D11 | 11-08 | Power figure | Re-simulated under the registered model (`results/power_simulation.csv`, 2,000 sims, HC3): 0.3 SD → 21%, 0.5 SD → 48%, 0.69 SD → 76%. 80% is reached at ≈ 0.72 SD (registered: ≈ 0.69) | Simulation error / setup detail | No |
| D12 | 11-08 | Implementation | Code named in the plan (`src/build_rois.py`, `src/preprocess_subject.py`, `src/group_analysis.py`) was written after registration, to the plan. Implementation tag `prereg-v1.0-impl1` is recorded in every qc.json | Provenance | No |
| D13 | 11-08 | "3 MAD" | Raw MAD (unscaled) | The plan does not specify scaling | No |
