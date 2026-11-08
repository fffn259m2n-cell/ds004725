# Reconciliation with the official pre-registration (OSF 8yqun)

**The governing plan is** Wakim, F. (2026, November 8), *Peripheral Neuropeptides and Oxytocin Effects on Amygdala-mPFC Connectivity (ds004725)*, osf.io/8yqun, registered "Nov 8, 2026, 3:30 AM" (timezone to be confirmed by the author).
`reports/preregistration.md` (git `8a8b4fb`) was an internal working plan. **It is superseded.** Where the two differ, OSF governs.

## A. Foreknowledge disclosure (for an OSF "Transparent Changes" update)

The registration states that imaging files would be downloaded only after it was submitted. The analysis assistant's logs show the following times (UTC):

| Time (UTC) | Event |
|---|---|
| 00:32 | Phenotype audit committed. This was disclosed in the registration |
| 00:37–00:43 | Pilot participant (sub-11002) T1w + BOLD downloaded. NIfTI headers, motion, ROI-coverage and registration QC were inspected, and the QC registration image was viewed |
| 00:56 onward | Batch preprocessing started under a **different** pipeline. By 01:30 UTC, 12 participants had been processed (25 by the time the batch fully stopped). Their connectivity values were written to `derivatives/oxtfc/*/_qc.json` and pushed to GitHub (commits `9d7b5bf`, `0f954d8`) |
| 01:30 (if 03:30 is UTC+2) | OSF registration |
| ≈ 02:00 | Batch stopped after 15 participants once the OSF plan was received |

**No connectivity value was opened, printed, summarised or modelled at any time.** No statistical model has been run on imaging data. Those pre-registration outputs come from a non-registered pipeline and are **excluded** from every registered analysis. They are kept only as an audit trail. Every registered analysis is recomputed from raw data with the pipeline in §B.

If the registration timestamp is *earlier* than 00:37 UTC (for example, if it is UTC+3 or later), no imaging access preceded registration, and this section only records the order of events.

## B. Element-by-element comparison and actions

| Element | OSF 8yqun (governs) | Previous implementation | Action |
|---|---|---|---|
| Dataset | v1.0.1 | v1.0.2 | Use 1.0.2. `git diff 1.0.1 1.0.2` changes only README and dataset_description (imaging and phenotype files are identical). Disclose |
| H1 sample | All participants passing QC, **up to 87** | 84 processed | **Process all 87** |
| H2/S sample | 83 with a pre-dose draw | Same | No change |
| Amygdala | HO subcortical, 25% max-prob, bilateral | HOSPA th25 labels 10 + 20 | Same |
| **mPFC** | **10 mm sphere at MNI (−1, 49, −5)**, per Ebner et al. 2016 | HO Frontal Medial Cortex label | **Switch to the sphere.** The coordinate is taken from the registration and has not been independently checked against Ebner 2016 |
| Anterior insula | HO insular cortex, 25%, **y > 0** | y ≥ 0 | Switch to strict y > 0 |
| ROI space | **ROIs mapped into native EPI space**; time series extracted natively, restricted to voxels with EPI signal | Data resampled to MNI | **Rebuild:** inverse-warp ROIs to native EPI and extract there |
| Normalization | **Mean EPI registered nonlinearly to MNI152 (ANTs)** | T1-based SyN + rigid EPI→T1 | **Rebuild:** SyN of the mean EPI to the MNI152 T1 template. The T1w is used only for WM/CSF masks (not specified by OSF; disclosed) |
| Realignment | Drop 4 volumes; rigid to mean | Same | Same |
| Nuisance | **24 motion parameters** + WM + CSF + derivatives; detrend + band-pass 0.01–0.08 Hz in one step | 12 motion + WM/CSF/derivatives | **Switch to 24P** (6 + derivatives + squares of both) |
| Censoring | FD > 0.5 mm **and the following volume** | 1 before + 2 after | **Switch** |
| FD | Power 2012, 50 mm | Same | Same |
| QC: motion | Mean FD > 0.5 mm | Same | Same |
| QC: data retained | **< 5 minutes** | < 150 volumes | Switch to minutes (TR-aware) |
| QC: registration | **Registered mean EPI** vs template r < median − 3·MAD, or a failed visual check **recorded in a file before FC is read** | Warped T1 r | Switch the metric; write `results/visual_qc.tsv` before any FC is read |
| QC: coverage | < 50% of **amygdala or mPFC**; insula only for the insula analyses | Min over all ROIs | Switch to analysis-specific rules |
| Preprocessing failure | Exclude and log the reason | Crash, now NaN | Log in `results/imaging_qc.csv` |
| Coding | T, A, S effect-coded ±0.5 | 0/1 and continuous age | **Switch** |
| H1 model | FC ~ T·A·S (full factorial) + meanFD; test b(T) and b(T×A×S) | Simpler anchor | **Implement as registered** |
| H2 model | FC ~ T + zAVP + T×zAVP + A, S, T×A, T×S, A×S, T×A×S + meanFD | Continuous age, T×age only | **Implement as registered** |
| Permutation | Freedman–Lane, 10,000 perms, seed 20261108, **HC3 t-statistic** | OLS t | **Switch to HC3 t** |
| Bayes factor | BIC-approximated BF01 (full vs reduced) | Not implemented | **Add**. Null read as evidence only if BF01 > 3 |
| Secondary | S1–S3, Holm | Same | Same |
| OT and AVP in one model | **Never** | A "joint model" sensitivity was planned | **Drop** |
| Robustness | Continuous age; days interval; proportion censored; TR; including the post-dose draws | GSR, raw AVP, strict FD, Huber, etc. | Run the OSF list as *registered robustness checks*. Any other analysis is labelled **non-registered exploratory** |
| Descriptives | Cell means (8 cells); left and right amygdala–mPFC FC; exclusions by arm × age × sex | Partial | Add |
| Software | statsmodels 0.15.0 | 0.15.0 installed | Same |

## C. Open questions for the author

1. What timezone is the "3:30 AM" registration timestamp in?
2. The registration says preprocessing is "described in full in the attached plan". Is there a file under osf.io/8yqun/files? It is not readable from this session without credentials. Please paste or upload it so the pipeline can follow it exactly.
3. WM/CSF masks are not specified. The plan here is a T1w Atropos segmentation, posterior > 0.9, eroded by one voxel, mapped to native EPI by the inverse EPI→T1 rigid transform. Does the attached plan say otherwise?
