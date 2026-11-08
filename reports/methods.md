# Methods

## Data

- **Source.** OpenNeuro ds004725 v1.0.2 (git `8372e9f`, CC0; Horta, Polk & Ebner, 2023, *Data in Brief*, doi:10.1016/j.dib.2023.109669). The dataset contains 87 generally healthy younger (18–31 y) and older (63–81 y) adults who took part in a randomized, double-blind, between-participant trial (NCT01823146).
- **Treatment.** Participants self-administered 24 IU of intranasal oxytocin (OT, n = 45) or placebo (PL, n = 42).
- **Scans.** A T1w MP-RAGE scan and an eyes-open resting-state BOLD scan were acquired 70–90 min after the spray (3T Philips Achieva; TR = 2 s, TE = 30 ms, 80 × 80 × 38 matrix, header voxel size 3.15 × 3.15 × 3.5 mm, 240 volumes). Two participants have TR = 2.5 s.
- **Blood measures.** Plasma OT and arginine vasopressin (AVP) came from a separate screening visit, measured by Enzo enzyme immunoassay (EIA) **without extraction**, in pg/mL (Plasencia et al., 2019).

Provenance and checksums: `reports/data_inventory.md`, `results/input_sha256.txt`. Every MRI file was checked against its git-annex SHA-256 key when fetched (`scripts/fetch_mri.sh`).

## Eligibility (Stage 1; `src/audit_phenotype.py`)

Participants were required to have:
- a valid treatment code;
- non-missing plasma OT and AVP;
- a blood draw documented as **preceding** dosing;
- BOLD and T1w images.

Two participants whose blood was drawn 5 and 8 days *after* dosing were excluded, as was one participant with an unknown interval. One participant had missing assays. These steps left n = 83.

## Preprocessing (`src/preprocess.py`; ANTsPy 0.6.3, nilearn 0.14.1)

fMRIPrep could not be run in the compute environment (deviation D1). The pre-registered pipeline was:

**Anatomical**
- N4 bias correction.
- SyN registration to MNI152NLin6Asym (2 mm, TemplateFlow).
- Brain mask from the inverse-warped template mask.
- 3-class Atropos segmentation; WM and CSF masks set at posterior > 0.9 and eroded by one voxel.

**Functional**
- Removal of 4 initial volumes.
- Rigid motion correction to the temporal mean.
- Rigid mutual-information coregistration of the mean EPI to the skull-stripped T1w.
- A single resampling of the motion-corrected series to 2 mm MNI space.
- No slice-timing correction (no slice timing in the metadata) and no susceptibility-distortion correction (no fieldmaps).

**Framewise displacement (FD)**
- Power FD computed with a 50 mm radius.

**Denoising**
- Regressors: 6 motion parameters + their derivatives, and mean WM and CSF signals + their derivatives.
- Linear detrending and a 0.01–0.08 Hz band-pass filter.
- Censoring of frames with FD > 0.5 mm, plus 1 preceding and 2 following frames. Censored frames are interpolated before filtering and dropped afterwards (nilearn `signal.clean`, `sample_mask`).
- Global signal regression (GSR) only as a sensitivity analysis.

## Regions and connectivity

Harvard–Oxford atlas, thresholded at 25%, in MNI152NLin6Asym:
- **Bilateral amygdala:** HOSPA labels 10 + 20.
- **mPFC:** HOCPA label 25, Frontal Medial Cortex.
- **Anterior insula:** HOCPA label 2, voxels with y ≥ 0.

HOCPA image labels equal TemplateFlow TSV indices + 1; this was verified by centroid and is covered by a unit test. Each ROI was intersected with the participant's EPI coverage mask. FC is the Fisher-z transformed Pearson correlation of the ROI mean signals over the retained frames.

## Quality control (pre-registered, blind to FC)

Exclusion criteria:
- Mean FD > 0.5 mm.
- Fewer than 150 retained frames.
- ROI coverage < 50%.

Normalization was reviewed visually for every participant (contact sheets). The relative flag was normalization r < median − 3·MAD (deviation D5).

## Statistics (`src/models.py`)

Primary model:

`FC ~ T + zAVP + T:zAVP + age_c + T:age_c + sex + meanFD`

- zAVP is the z-scored ln AVP.
- The model was fitted by OLS with HC3 standard errors and 95% CIs.
- The confirmatory p is a Freedman–Lane permutation p (10,000 permutations, seed 20261108).
- Conditional OT−PL effects are reported at zAVP = −1, 0 and +1.

Secondary tests (Holm-corrected, family of 3):
- T × OT on amygdala–mPFC.
- T × AVP on amygdala–anterior insula.
- T × OT on amygdala–anterior insula.

Sensitivity analyses are listed in `reports/preregistration.md` §7. Diagnostics covered Cook's distance, leverage, residual normality, a quadratic moderator term and moderator overlap by arm.

The replication anchor reproduces Ebner et al. (2016): the main effect of T, and T × cohort × sex. Power at the observed design is given in `results/power_interaction.csv`.

## Software

Python 3.13 with pandas, numpy, scipy, statsmodels, nibabel, nilearn 0.14.1, antspyx 0.6.3 and matplotlib, plus DataLad 1.8.0 and git-annex 10.20240129. The full reproduction sequence is in `PROJECT_README.md`.
