# Pre-registered analysis plan: frozen before any group-level FC is computed

> **SUPERSEDED (2026-11-08):** the governing pre-registration is OSF osf.io/8yqun. This file was an internal working plan. See `reports/osf_reconciliation.md` for the differences and the foreknowledge disclosure.

**Frozen:** 2026-11-08, in the git commit that adds this file. All later deviations are logged in §9.
**Already known when this was written:** the Stage 1 phenotype audit (`reports/feasibility.md`) and the literature (`literature/novelty_matrix.md`). The author has seen BOLD/T1 headers for one pilot participant (sub-11002). No functional-connectivity values from any participant had been computed when this file was frozen.
**Status:** gate AMBER, so this is an **exploratory estimation study**.

## 1. Estimands

- **Primary (H1):** the interaction between randomized treatment (OT = 1, PL = 0) and z-scored ln(screening plasma AVP) on Fisher-z resting-state FC between the **bilateral amygdala** and the **medial prefrontal cortex (Frontal Medial Cortex)**.
- **Secondary (Holm-corrected as a family of 3):**
  - S1: T × z(ln OT) on amygdala–mPFC FC.
  - S2: T × z(ln AVP) on amygdala–anterior insula FC.
  - S3: T × z(ln OT) on amygdala–anterior insula FC.
- **Replication anchor (descriptive):**
  - The main effect of T on amygdala–mPFC FC.
  - The T × age-cohort × sex pattern reported by Ebner et al. (2016).
  - Neither is a new claim. Both check the pipeline.

## 2. Sample (decided at Stage 1)

The sample is the 83 participants who have a valid treatment code, plasma OT and AVP, a blood draw confirmed **before** dosing, and BOLD and T1w images (OT 44, PL 39).

Excluded before imaging:
- 1 participant with missing assays.
- 2 participants whose blood was drawn after dosing.
- 1 participant with an unknown draw interval.

## 3. Regions of interest

All ROIs are in MNI152NLin6Asym at 2 mm, taken from TemplateFlow (SHA-256 sums in `derivatives/templates/SHA256SUMS`). They come from the Harvard–Oxford atlas, thresholded at 25%. This atlas is independent of these data.

| ROI | Atlas | Image label(s) | Rule |
|---|---|---|---|
| Bilateral amygdala | HOSPA th25 | 10 (L) + 20 (R) | Labels combined; mean signal over the voxels |
| mPFC | HOCPA th25 | **25** (Frontal Medial Cortex) | Whole label |
| Anterior insula | HOCPA th25 | **2** (Insular Cortex) | Voxels with MNI y ≥ 0 only |

The HOCPA TSV in TemplateFlow lists FSL's 0-based indices, so image label = TSV index + 1. This was verified by centroid: label 25 is at (0, 43, −19) and label 2 is the insula. A unit test checks it. Each ROI is intersected with the participant's EPI brain-coverage mask before averaging.

## 4. Preprocessing (ANTsPy 0.6.3 / nilearn 0.14.1; fMRIPrep unavailable, see §9)

1. **T1w.**
   - N4 bias correction.
   - SyN registration to the MNI152NLin6Asym 2 mm template.
   - Brain mask: the template mask warped back to native space.
   - Tissue segmentation: Atropos 3-class (CSF/GM/WM) within the mask.
   - WM and CSF masks: posterior > 0.9, eroded by one voxel.
2. **BOLD.**
   - Discard the first 4 volumes (non-steady state). Of the remaining 236 volumes, TR = 2 s.
   - No slice-timing correction: the JSON has no SliceTiming and TR is 2 s.
   - Rigid motion correction to the temporal mean (`ants.motion_correction`).
   - Rigid boundary-free coregistration of the mean EPI to the T1w (mutual information).
   - Concatenate with the T1→MNI warp and resample to 2 mm MNI in a single step from the motion-corrected series.
   - No distortion correction (no fieldmaps).
3. **Framewise displacement (FD).** Power FD from the 6 rigid parameters, using a 50 mm head radius.
4. **Denoising (primary, applied with nilearn `signal.clean`).**
   - Nuisance regressors: 6 motion parameters + their derivatives; mean WM and CSF signals + their derivatives (16 regressors).
   - Linear detrending.
   - Band-pass filter 0.01–0.08 Hz.
   - Censoring: drop volumes with FD > 0.5 mm, plus the volume before and the two after each. Censored volumes are interpolated before filtering and removed afterwards.
5. **Global signal regression (GSR)** is run only as a sensitivity analysis: the 16 regressors + global signal and its derivative.
6. **FC.** Pearson r between ROI mean time series over the retained volumes, then Fisher z = atanh(r).

## 5. Quality-control exclusions

These rules are applied automatically and blind to FC values:
- Mean FD > 0.5 mm.
- Fewer than 150 retained volumes (5 min).
- Less than 50% of any ROI's voxels inside the EPI coverage mask (signal dropout).
- Failed normalization: Pearson r between the warped T1w and the template brain < 0.80. Any participant flagged here gets a visual check, with the decision logged.

## 6. Statistical model

```
FCz ~ T + zAVP + T:zAVP + age_c + T:age_c + sex + meanFD
```

- `zAVP` = z-scored ln(AVP) in the analysed sample.
- `age_c` = age centred, in decades.
- `sex` = F = 1.
- `meanFD` = centred.

Inference:
- OLS with HC3 standard errors; report the 95% CI for the T:zAVP coefficient.
- Freedman–Lane permutation p-value with 10,000 permutations; the permutation p is the confirmatory p.
- Two-sided α = .05 for H1. Holm correction across S1–S3.
- Report the conditional OT-minus-PL effect at zAVP = −1, 0 and +1 with 95% CIs, and plot predicted FC by arm across AVP.

Diagnostics:
- Overlap of AVP by arm.
- Leverage and Cook's distance (report any participant with D > 4/n, without excluding them).
- Residual vs fitted plots.
- A quadratic term in zAVP, reported only as a diagnostic.

## 7. Sensitivity analyses (reported in full, never used to replace the primary)

1. GSR.
2. Raw (untransformed) z(AVP).
3. Adding days from blood draw to scan and its interaction with T.
4. Cohort factor instead of continuous age (T × cohort).
5. Excluding the TR = 2.5 s participant(s).
6. Stricter motion threshold (mean FD > 0.3 mm excluded; censoring at FD > 0.3 mm).
7. Huber robust regression.
8. Joint model with both T:zAVP and T:zOT. Interpret with caution, because ρ(OT, AVP) ≈ −0.64.
9. Adding the participant with an unknown draw interval.

## 8. Interpretation rules

- A significant H1 means only that the OT-vs-PL difference in amygdala–mPFC FC varies with a screening-day peripheral immunoassay value, after allowing the effect to vary with age.
- It does **not** show that AVP causes the difference. It says nothing about central AVP or receptor occupancy, and nothing about behaviour, morality or clinical status.
- A null result is reported with its CI and the power curve (`results/power_interaction.csv`). Never report it as "no effect".

## 9. Deviations from `CLAUDE.md` and this plan

| # | Deviation | Reason |
|---|---|---|
| D1 | ANTsPy/nilearn pipeline instead of fMRIPrep | No Docker daemon or Singularity in the session; ~1 h per participant for fMRIPrep would also exceed the compute available. The pipeline is documented, with its versions pinned in `results/` |
| D2 | Model adds T × age_c | Required by the Stage 1 finding r(AVP, age) = .40 and the published T × age effects |
| D3 | Dataset v1.0.2 rather than 1.0.1 | Only the README/DOI changed between versions |
| D4 | No OT/AVP ratio | Dropped because of the strong anticorrelation and the unextracted assay |
| D5 | Normalization QC rule changed from "r < 0.80" to "r < sample median − 3·MAD, plus visual review of every participant's QC montage" | Pilot (sub-11002) was well registered on visual inspection but had r = 0.69, because the metric compares a skull-stripped T1 with the template brain across contrasts. Logged 2026-11-08 after the pilot only, before any group FC was read |
| D6 | Pipeline parameters checked on the pilot only: Power FD from ANTs rigid parameters (mean 0.30–0.33 mm across 4 estimator variants, ~7% of frames > 0.5 mm) | Kept the pre-registered estimator. Recorded so readers know the FD scale is higher than typical FSL MCFLIRT values |
