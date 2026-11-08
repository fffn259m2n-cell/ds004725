# Stage 2: literature novelty matrix

Searched 2026-11-08 in PubMed (via MCP) and on the web. Queries:
- `Ebner NC[Author] AND oxytocin`, which returned 33 records; all titles and abstracts were screened.
- `Ebner N[Author] AND oxytocin AND (plasma OR vasopressin)`.
- `oxytocin AND (endogenous OR plasma OR baseline) AND intranasal AND (connectivity OR fMRI) AND (predict* OR moderat* OR individual differences)`.
- A search for the DOI of the 2026 paper.
- A web search for `"ds004725" oxytocin resting-state`.

Limits: the PMC full-text host was not reachable from this session. Full text was read through the PubMed MCP for Plasencia et al. 2019 only. Google Scholar, OSF and preprint servers were **not** searched. **Novelty is therefore provisional.**

## Studies using this cohort (UF single-dose OT study, NCT01823146; 2013–2014)

| Study | Data / N | Measures | Model / ROIs | Main conclusion | Overlap with our estimand |
|---|---|---|---|---|---|
| Ebner et al. 2016, *Psychoneuroendocrinology* 69:50 ([doi](https://doi.org/10.1016/j.psyneuen.2016.03.013)) | rsfMRI, N = 79 | OT vs PL | Bilateral amygdala ↔ mPFC FC; treatment × age × sex | OT > PL and men > women in amygdala–mPFC FC. Three-way interaction: OT increased FC only in young women | **High for circuit and treatment effect** (same ROI pair, same scans). No plasma moderator |
| Liu, Lin, Feifel & Ebner 2022, *NeuroImage* 253:119045 ([doi](https://doi.org/10.1016/j.neuroimage.2022.119045)) | rsfMRI, 42 young + 43 older | OT vs PL | ICA salience network; peaks L amygdala (−18, 0, −15), R insula (39, 6, −6); age moderation | OT lowered within-SN FC; the amygdala effect was present only in older adults | High for amygdala/insula (our exploratory edge). No plasma moderator |
| Xiao, Ebner, … Fischer 2023, *Psychoneuroendocrinology* ([doi](https://doi.org/10.1016/j.psyneuen.2023.106666)) | rsfMRI, 46 young + 44 older | OT vs PL | Connectivity density and strength of OXTR-rich regions (thalamus, pallidum, caudate, putamen, olfactory bulb) | Age-dependent: OT ↓ thalamus (young), ↑ caudate (older) | Moderate (treatment × age). No plasma moderator |
| Klugah-Brown et al. 2026, *Mol Psychiatry* 31:4310 ([doi](https://doi.org/10.1038/s41380-026-03549-9)) | ds004725 rsfMRI | OT vs PL | Co-activation patterns (CAPs); age-stratified; CAP dwell time × DSST | OT altered CAP dynamics involving amygdala, mPFC and SN; age-differential | Moderate (same circuit, dynamic method). The abstract mentions no plasma moderator |
| Fooladi et al. 2025, *Brain Res Bull* 232:111604 ([doi](https://doi.org/10.1016/j.brainresbull.2025.111604)) | rsfMRI, N = 87 (matches ds004725) | OT vs PL | Papez circuit seed-to-voxel, ROI-to-ROI, graph metrics | Group differences in thalamus, PCC, ACC and parahippocampus | Low–moderate. No plasma moderator |
| Plasencia et al. 2019, *Psychoneuroendocrinology* 110:104419 ([doi](https://doi.org/10.1016/j.psyneuen.2019.104419)) | Screening plasma, N = 101 (superset of ds004725) | Plasma OT and AVP (**Enzo EIA, both unextracted**; dilution 1:8 OT, 1:2 AVP; CV < 8%) | ANOVA and regressions of OT/AVP on sex, age, attachment, cognition | Older > young in AVP. **OT and AVP strongly negatively correlated in every sex × age group.** | **Key context.** The authors write that endogenous × exogenous OT interactions "cannot be applied … the current sample is too small to allow testing these interaction effects", and that "baseline OT assessments did not take place on the same day as the OT intervention" |
| Horta et al. 2019, *Neurobiol Aging* ([doi](https://doi.org/10.1016/j.neurobiolaging.2019.01.016)); Lin et al. 2022, *Front Psychol* ([doi](https://doi.org/10.3389/fpsyg.2022.838642)); 2021 *Psychol Aging* ([doi](https://doi.org/10.1037/pag0000545)) | Task fMRI (same visit) | OT vs PL | Task PLS, amygdala | Age × OT task effects | Low (task, not rest). No plasma moderator |
| Horta, Polk & Ebner 2023, *Data in Brief* ([doi](https://doi.org/10.1016/j.dib.2023.109669)) | Data descriptor | – | – | – | Source documentation |

Related work from the same lab but other cohorts: Polk et al. 2022 *Phil Trans B* ([doi](https://doi.org/10.1098/rstb.2021.0048)) and 2025 *Sci Rep* ([doi](https://doi.org/10.1038/s41598-025-07353-3)) used plasma OT/AVP and social cognition in older adults. A 2024 *Neuropharmacology* paper ([doi](https://doi.org/10.1016/j.neuropharm.2024.110130)) studied chronic OT and rsFC in older adults. None tests plasma moderation of an acute OT effect on FC.

## Assessment

- **No record found** of a treatment × screening plasma AVP (or OT) interaction on rsFC in this dataset. The original investigators explicitly declined to run it because the sample is too small. **Novelty is provisional**, pending searches of Google Scholar, OSF and bioRxiv/PsyArXiv and a check of Klugah-Brown 2026's supplementary material.
- The **circuit is not novel**. Ebner 2016 already reports OT → amygdala–mPFC FC, with age × sex moderation, on these scans. Our analysis must be framed as **heterogeneity of a known effect**, and it must adjust for treatment × age. Because AVP rises with age (Plasencia 2019; r = .40 here), T × AVP without a T × age term would largely re-describe Ebner 2016 and Liu 2022.
- The **OT–AVP anticorrelation is published** (Plasencia 2019), so we must not present it as a new finding. Its cause (biology vs shared unextracted-EIA interference) remains contested. Plasencia et al. defend unextracted EIA, while the specificity literature (e.g. Szeto et al. 2011 and later kit comparisons) warns of inflation and non-specificity. *Cite primary sources after full-text verification.*
- Ebner 2016 analysed 79 participants and Liu 2022 analysed 85. Our 83 eligible participants are a different subset: excluding post-dosing draws removes 2 young PL participants. Report how our sample overlaps with theirs.

## Consequences for the pre-registered estimand

The estimand is not already published, so there is no stop condition. Proceed under the AMBER framing:
- Primary analysis: T × z(AVP) on amygdala–mPFC FC, **adjusting for T × age**. This is an estimation study.
- Add a **replication anchor**: the main OT effect on amygdala–mPFC FC, as a check of the pipeline against Ebner 2016. It is not a new claim.
