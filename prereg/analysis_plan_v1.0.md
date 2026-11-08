# Preregistered analysis plan, v1.0 (frozen before imaging data were accessed)

Project: Does peripheral neuropeptide immunoreactivity predict the neural response to intranasal oxytocin? A preregistered reproduction and extension in OpenNeuro ds004725.

Frozen: 2026-11-08 (UTC), git tag `prereg-v1.0`. At the time of freezing, no BOLD or T1w image from ds004725 had been downloaded into, or viewed in, the analysis environment. Only metadata and phenotype files were examined (`reports/feasibility.md`). The plan was written by the analysis team (F. Wakim). This plan should be deposited on OSF; the git tag and commit hash give the time-stamped record until then.

## 1. Data
OpenNeuro ds004725, tag 1.0.1 (commit d40af5c2689a9a71baf97875e9001219033e1c12): 87 healthy adults (44 aged 18 to 31; 43 aged 63 to 81). Each received a single 24 IU dose of intranasal oxytocin (OT, n = 45) or placebo (PL, n = 42), randomised and double-blind (NCT01823146). Resting-state fMRI was acquired 70 to 90 min after dosing (3T Philips Achieva, TR 2 s, 240 volumes, 3 mm isotropic). Plasma OT and AVP were measured at a screening visit with an Enzo EIA on unextracted plasma.

## 2. Questions and hypotheses
- **H1 (reproduction).** Ebner et al. (2016, doi:10.1016/j.psyneuen.2016.03.013) reported that OT increased amygdala–mPFC resting connectivity (main effect), with a Treatment × Age × Sex interaction (largest increase in young women). We test both terms with a different, fully documented pipeline and the full released sample. Prediction: OT > PL (positive T coefficient). The three-way term is tested two-sided.
- **H2 (primary, new).** Screening plasma AVP immunoreactivity moderates the OT vs PL difference in amygdala–mPFC connectivity. Two-sided, because no directional prior is strong enough.
- **Secondary (Holm-corrected as one family of three):**
  - S1: plasma OT immunoreactivity as moderator of amygdala–mPFC connectivity.
  - S2: AVP as moderator of amygdala–anterior insula connectivity (exploratory edge).
  - S3: OT as moderator of amygdala–anterior insula connectivity (exploratory edge).

## 3. Regions (files in `rois/`, built by `src/build_rois.py`, hashes in `reports/provenance.md`)
- **AMY:** bilateral Harvard-Oxford amygdala, max-probability among subcortical labels, threshold 25%.
- **MPFC:** 10 mm sphere at MNI (−1, 49, −5).
- **AINS:** Harvard-Oxford insular cortex (max-probability among cortical labels, threshold 25%), y > 0.
- **Nuisance masks:** WM, from MNI152 WM probability > 0.9, eroded once; CSF, from the lateral ventricles at 50% probability, excluding grey matter.

## 4. Preprocessing (`src/preprocess_subject.py`, v1.0)
1. Discard 4 volumes.
2. ANTs rigid realignment to the mean volume.
3. Power framewise displacement (FD).
4. Mean EPI to MNI152 2 mm template: ANTs rigid + affine + SyN, Mattes MI.
5. Masks carried to native space with the inverse transform and intersected with an EPI mask.
6. Simultaneous nuisance regression and band-pass filtering:
   - band-pass 0.01–0.08 Hz with linear detrending;
   - 24 motion regressors;
   - WM and CSF signals and their derivatives;
   - censoring of volumes with FD > 0.5 mm, plus the following volume.
7. Mean ROI time series; Pearson r converted to Fisher z.

No slice-timing correction and no spatial smoothing. The repetition time is short relative to the band-pass filter, and analysis is at ROI level. The sidecars contain no slice-timing information.

## 5. Exclusions (applied from `qc.json` files before any connectivity value is read; `--qc-only` mode)
- Mean FD > 0.5 mm.
- Fewer than 5 minutes of retained data.
- Registration failure. This is defined as a correlation between the warped mean EPI and the template that falls more than 3 MAD below the sample median, or a failed visual check of `qc_reg.png`. The visual checks are recorded in `results/qc_visual.csv` before connectivity is read.
- AMY or MPFC coverage < 50% of mask voxels inside the EPI mask. For S2/S3 only, the same rule applies to AINS.
- For H2 and the secondary family, participants must also have a blood draw dated before dosing and non-missing OT and AVP (analytic set from `reports/feasibility.md`: 83 before imaging QC).

## 6. Models (OLS, HC3 standard errors; `src/group_analysis.py`)
Treatment (T), age group (C) and sex (S) are effect-coded ±0.5.
- **H1:** `FC ~ T * C * S + mean_FD`. Test T and T:C:S.
- **H2:** `FC ~ T * zAVP + T * C * S + mean_FD`, where zAVP is the standardised log of plasma AVP within the analytic set. Test T:zAVP. S1–S3 use the same form.

Adjusting for T × C × S is required. AVP rises with age (ρ = 0.39), and Ebner et al. reported a Treatment × Age × Sex effect on this edge. Without the adjustment, an age-dependent OT effect would appear as AVP moderation (`results/power_simulation.csv`).

Inference for H1 (T, T:C:S), H2 and S1–S3:
- 95% HC3 confidence interval;
- Freedman–Lane permutation p (10,000 permutations, seed 20261108), which is the primary p-value;
- BIC-approximated BF01 for the interaction term.

Alpha: H1 at .05 per term; H2 at .05; secondary family Holm-adjusted at .05.

## 7. Sensitivity analyses (H2; H1 where marked)
- GSR (H1 and H2).
- No censoring (H1 and H2).
- FD > 0.3 mm censoring (H1 and H2).
- Blood draw ≤ 30 days before dosing.
- Draw interval as covariate.
- Including the two post-dose draws.
- Excluding the participant with TR 2.5 s.
- Continuous age instead of age group.
- Adding AVP × age group.
- Adding proportion of censored volumes.
- The original handoff model (age and sex main effects only).
- Huber robust regression.
- Quadratic moderation (T:zAVP²).
- Diagnostics: leverage, Cook's distance, residual normality, heteroscedasticity, moderator overlap between arms.

## 8. Interpretation rules (fixed now)
- **Detection limit.** With N ≈ 83, the study has 80% power only for interactions of about 0.69 SD of FC per SD of log AVP. A non-significant H2 will be reported as a confidence interval with the effect sizes it excludes. It will not be reported as evidence of no moderation unless BF01 > 3.
- **What the moderator measures.** The moderator is described as unextracted EIA immunoreactivity from a single screening draw (median 16 days before dosing). It is not described as circulating or central AVP.
- **Plasma OT is observational.** Moderation is an observational association within a randomised trial and supports no causal claim about AVP.
- **Scope of claims.** No claim about social behaviour, empathy, morality, clinical traits or receptor occupancy.

## 9. Deviations
Any change after this tag is logged in `prereg/deviations.md` with date, reason and whether connectivity values had been seen.
