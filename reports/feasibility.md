# Stage 1 feasibility report: phenotype audit

Generated 2026-11-08 by `scripts/run_stage1.sh` from ds004725 **v1.0.2** (git `8372e9f`).
Every number below comes from files in `results/`. No imaging data has been downloaded or analysed.

## Gate decision: **AMBER**

Coverage is fine: 83 participants have a screening AVP value drawn before dosing, and both arms are represented (OT 44, PL 39). The problem is measurement quality. The AVP values look like an assay artifact, not a physiological signal (see §3), and AVP is partly a proxy for age. The primary T × AVP model can go ahead only as a **pre-registered exploratory or estimation analysis**, with the changes in §6. Results must not be presented as confirmatory evidence about "baseline vasopressin biology".

## 1. Verified sample (eligibility flow, `results/eligibility_flow.csv`)

| Step | N | OT | PL | Younger | Older |
|---|---|---|---|---|---|
| Listed in `participants.tsv` | 87 | 45 | 42 | 44 | 43 |
| Valid treatment code (OT/PL) | 87 | 45 | 42 | 44 | 43 |
| Non-missing plasma AVP (and OT) | 86 | 44 | 42 | 44 | 42 |
| Blood draw known to precede dosing | **83** | **44** | **39** | 41 | 42 |
| Resting BOLD + T1w present in dataset | 83 | 44 | 39 | 41 | 42 |
| [sensitivity] standard TR = 2.0 s | 82 | 44 | 38 | 40 | 42 |

- The 87 IDs match one to one across `participants.tsv`, `phenotype/phenotypic_data.tsv` and the `sub-*` directories. There are no duplicates or orphans, and no bad treatment codes.
- Design cells (`results/design_cells.csv`): OT-younger 24 (11 F), OT-older 21 (11 F), PL-younger 20 (10 F), PL-older 22 (13 F).
- Every participant has a resting BOLD and a T1w annex pointer. The reported README attrition (participants who lost their rest scan) happened **before** release, so the release only contains people who have a scan.
- The ID prefix exactly encodes cohort and sex: 11 = younger M, 12 = younger F, 21 = older M, 22 = older F. This is recorded only as a consistency check. Analyses must use the `age` and `sex` columns.

## 2. Exclusions identified by the audit (made before seeing any outcome)

| Issue | n | Arm/cohort | Handling |
|---|---|---|---|
| Plasma OT & AVP `n/a` | 1 | OT, older | Excluded (complete case) |
| Blood drawn **after** dosing (`"5 days after"`, `"8 days after"`) | 2 | PL, younger | **Excluded.** These are post-treatment measurements and cannot act as baseline moderators. |
| Sampling interval `n/a` | 1 | PL, younger | Excluded, because the draw cannot be confirmed as pre-treatment. Sensitivity analysis: include it. |
| Non-standard TR 2.5 s (sub-level JSON) | 2 | 1 PL; 1 also post-dosing | Keep in the primary analysis (FC is a correlation, so it is TR-tolerant after filtering). Exclude in a sensitivity analysis. |

## 3. Plasma biomarkers: units, distributions, validity concerns

`results/biomarker_descriptives.csv`, `reports/figures/biomarkers_by_arm_cohort.png`, `reports/figures/avp_vs_ot.png`

| | n | Median (IQR), pg/mL | Range | Skew (raw → log) | Distinct values |
|---|---|---|---|---|---|
| AVP | 86 | 138.6 (120.4–149.5) | 97.0–195.7 | 0.24 → −0.14 | 53 |
| OT (unextracted) | 86 | 794.8 (705.9–903.6) | 489.8–1207.7 | −0.02 → −0.54 | 58 |

Both are recorded in pg/mL, which the dictionary confirms. Neither distribution is skewed, and the log transform does not improve symmetry. The pre-specified ln transform is therefore optional. **Recommendation: pre-register raw z-scores as primary and ln as a sensitivity analysis**, or the other way round, but decide before any FC data is touched.

**Validity red flags (must be in the limitations and discussed with the supervisor):**

1. **AVP magnitude.** Values of about 100–200 pg/mL are roughly one to two orders of magnitude above typical extracted-plasma AVP in healthy adults (usually single-digit pg/mL). The dictionary says the OT assay was *unextracted* (Enzo EIA). It says nothing about extraction for AVP. Unextracted immunoassays are known to give inflated, poorly specific neuropeptide readings because of matrix interference. *Source details to be verified and cited in Stage 2.*
2. **OT–AVP correlation of ρ = −0.64** (Spearman, n = 86). It holds within each cohort (younger −0.75, older −0.62). A strong *inverse* coupling between two separately secreted peptides is biologically unexpected. A shared matrix or plate artifact is a plausible alternative explanation. As a result, OT and AVP cannot be treated as independent moderators, and a log(OT/AVP) ratio would be dominated by this artifact. **Do not use the ratio.**
3. **Ties.** 86 values contain only 53 (AVP) and 58 (OT) distinct values, consistent with back-calculation from a standard curve. This is fine for linear models but limits the effective resolution.
4. **AVP is associated with age:** r(log AVP, age) = 0.40, with older adults higher (median 145.9 vs 130.4). Ebner et al. (2016) already reported OT × age (× sex) effects on rsFC in this sample, so an unadjusted T × AVP term would partly re-estimate T × age. OT shows no age association (r = −0.06).

## 4. Timing

- Days from blood draw to dosing and scan, for draws that preceded dosing (n = 83): median 15, range 1–138. Eight people had more than 60 days, 6 of them in the OT arm. Because of this spread, "baseline" is a single screening-day value and should not be described as a stable trait. No within-person reliability data exist in the dataset.
- **The interval is imbalanced by arm:** mean 30.9 d for OT vs 17.6 d for PL, SMD 0.48 (`results/arm_balance.csv`). The interval was not randomised, so it should be a pre-registered covariate or sensitivity term.
- Blood-draw clock times are stored as `H.MM` strings, so `9.58` means 09:58. Seven values have one-digit minutes (`10.1`), which are read as 10:10 because `10.01` also occurs in the same column (an assumption, covered by a test). Two values (`0`, `1.05`) are impossible under the ~08:00 protocol and are set to missing.
- Rest-scan clock times: five values (`11.80`, `11.90`) are impossible as H.MM and are set to missing. That leaves 82 usable values. Do not guess a conversion.

## 5. Other data-quality notes

- **The dictionary does not match the data:**
  - `CONTRACEPTION_USE_AT_SCAN` is documented as 0/1, but the data use 1/2 (8 × "1", 13 × "2"). Which code means "yes" cannot be resolved from the files, so the variable must not be used until the curators clarify it.
  - `SALIVARY_APOE_GENOTYPE` is documented as `ε3/ε3`, but the data use `e3e3`. The genotype matches the two SNP columns in all 82 genotyped participants, so it is internally consistent.
- `MEDICATION_STATUS_AT_SCAN`: 33% of OT vs 59% of PL participants reported a medication. This is post-randomisation context, so it is reported but not used as a covariate.
- The phenotype TSV uses CRLF line endings. The parser strips them, and a test covers this.
- Rest BOLD JSON has no `SliceTiming` field and there are no fieldmaps. The README states 38 interleaved ascending slices (1, 3, 5, …), TR 2 s, and 4 dummy scans at the start and end. fMRIPrep will need either `--ignore slicetiming` or a derived `SliceTiming` documented in `derivatives/`. Distortion correction is limited to SyN-SDC or none. Decide this in Stage 3.

## 6. Power and sensitivity (`results/power_interaction.csv`)

The simulation uses the observed design: N = 83, the actual arm split, z(log biomarker), centred age and sex. Only the outcome was simulated, using OLS with α = .05 two-sided.

| Interaction (FC SD per SD of moderator) | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |
|---|---|---|---|---|---|
| Power, AVP | 0.25 | 0.41 | 0.57 | 0.73 | 0.86 |

About 80% power needs an interaction of roughly 0.65 SD of FC per SD of AVP. That is a large effect for peripheral-marker moderation. Adding T × age, which is required (below), will lower power further. **Frame the primary analysis as estimation: report the CI and do not interpret a null as evidence of no effect.**

### Required changes to the pre-specified model before Stage 3

```
FC_z ~ T * z(AVP) + T * age_c + sex + mean_FD [+ days_blood_to_scan as sensitivity]
```

- Add **T × age_c**, using either continuous centred age or a cohort factor; pre-register which. This stops the AVP moderation from proxying age moderation.
- Run OT as a secondary moderator in a **separate** model, with no OT/AVP ratio. Optionally fit a joint model, interpreted with caution because ρ = −0.64.
- Exclude the 2 post-dosing draws and the 1 unknown-interval draw (§2).

## 7. Next executable commands

```bash
scripts/setup_annex.sh                 # once per clone; attaches OpenNeuro S3 remote
scripts/fetch_one_mri.sh sub-11002     # Stage 3 pilot: one T1w + BOLD, ~0.11 GB, header check
# Full eligible fetch is ~11.9 GB raw (BOLD 10.1 GB + T1w 1.8 GB); fMRIPrep needs ~
# 4-8 GB RAM/subject and ~100+ GB working space; not available in this cloud session.
```

Stage 2 (literature novelty audit) has **not been run**. Before any claim of novelty, it must check whether Liu et al. 2022 (salience network), Ebner et al. 2016 (rsFC × age/sex), or the 2026 *Mol Psychiatry* paper already tested peripheral OT/AVP moderation.
