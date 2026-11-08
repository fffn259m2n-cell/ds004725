# Claude Code: Oxytocin × baseline neuropeptides × resting-state connectivity

## Mission
Build a reproducible, publication-oriented **secondary analysis** of OpenNeuro `ds004725`, version `1.0.1`, to test whether **pre-treatment peripheral oxytocin (OT) or vasopressin (AVP)** moderates the randomized effect of **intranasal OT vs placebo** on prespecified resting-state fMRI connectivity. Start **immediately** with download/metadata auditing. Do not interpret plasma assays as central neurotransmission; do not claim clinical ASPD/psychopathy effects or moral enhancement. No intervention was administered by this code.

## Source-of-truth links
- Dataset: https://openneuro.org/datasets/ds004725/versions/1.0.1
- Datalad/Git-annex mirror: https://github.com/OpenNeuroDatasets/ds004725
- Data descriptor: https://pmc.ncbi.nlm.nih.gov/articles/PMC10630614/
- 2026 secondary connectivity/brain-state paper: https://doi.org/10.1038/s41380-026-03549-9
- Earlier salience-network findings: find and document in `literature/novelty_matrix.md` before claims of novelty.

## Verified background vs required verification
- Reported: 87 younger/older healthy participants; randomized double-blind between-person intranasal OT vs placebo; T1 and resting-state BOLD; screening-visit plasma OT/AVP; cognitive testing; APOE-related measures. Confirm actual files.
- `participants.tsv` contains `participant_id`, `age`, `sex`, `handedness`, `treatment_group`, `years_of_education`, `self_report_race`, `self_report_ethnicity`, `bmi` (check actual headers).
- Data descriptor lists phenotype fields `BLOOD_OXYTOCIN_LEVELS`, `BLOOD_VASOPRESSIN_LEVELS`, `TIMING_BIOLOGICAL_SAMPLING_RELATIVE_TO_TREATMENT`, `TIME_OF_DAY_BLOOD_SAMPLING`, `TIME_OF_DAY_RESTING_SCAN`, `TICS_SCORE`, `DIGIT_TOTAL`, `DIGIT_CORRECT`, `RAVLT_TOTAL`, `RAVLT_CORRECT`; confirm TSV values, units, missingness and dictionary.
- The repository uses a folder named `phenotype/`, but the data descriptor calls it `phenotypes/`. Discover, do not hardcode both as fact.
- Clinical antisocial/psychopathy scores and moral-behavior endpoints **are not verified** in this dataset.
- Scan existence, imaging N after QC, plasma complete-case N, time separation and batch effects **are unknown** until audited.

## Non-negotiable methodology
1. **Audit first.** Never make up a results table, exact count, ROI availability or downloaded file. Log commands, errors, filenames, hashes, disk/compute needs, number of participants, randomized arm counts, plasma missingness by arm and age cohort, assay units and distribution, time gap from screening to dosing/scanning.
2. **Avoid false novelty.** Search published literature, Google Scholar/PubMed/preprints/registered plans, especially analyses of ds004725. Write a novelty matrix with what was already measured, models, ROIs and main conclusions. Mark novelty as provisional; never 'absolutely new'.
3. **One prespecified primary test:** Treatment × baseline **plasma AVP** interaction on **bilateral amygdala-to-medial prefrontal cortex** functional connectivity, with baseline AVP log-transformed if appropriate and standardized, plus main treatment and main moderator. Choose ROIs from a cited independent atlas and pin exact atlas labels/coordinates before seeing group results. Motivate AVP as potential moderator of OT effects, not an AVP intervention.
4. **Secondary:** Treatment × baseline OT on same FC; exploratory amygdala-to-anterior insula FC and separate OT/AVP moderators; if considering log(OT/AVP), confirm positive, comparable assay units and avoid ratio as an unquestioned biological measure. Predefine family-wise/FDR approach for secondary tests.
5. **Confounders:** Age (continuous, centered), sex if valid, mean FD and proportion scrubbed frames; consider acquisition variations, time of scan, blood-to-scan interval, and time of blood draw. Be parsimonious. No post-treatment variables as predictors; clearly classify screening measurements. Treatment randomized; biomarkers observational.
6. **Preprocess reproducibly:** inspect BIDS metadata, field maps, TR, slice timing, acquisition quirks, scan duration, and any existing derivatives. If no approved derivatives, use pinned fMRIPrep/container version, document FreeSurfer/TemplateFlow and resources, and QC. Resting-state nuisance regression (motion, WM/CSF and drift); describe global-signal-regression sensitivity, band-pass filtering, censoring, and minimum usable minutes. Preserve original data and write to `derivatives/`.
7. **FC definition:** extract ROI time series with an independent atlas; Fisher-z Pearson correlations; do not use model-fitted connectivity or same data for choosing ROI and hypothesis. Plot motion by group and assess outliers.
8. **Model:** `FC_z ~ T * z(log(AVP)) + age_centered + sex + mean_FD` (exact final covariates pre-specified after QC without reading outcome-group differences). Test coefficient and 95% CI for interaction; plot predicted marginal difference by biomarker; use permutation / robust sensitivity appropriate to randomization and nuisance structure, and bootstrap/HC3 intervals if warranted. Check overlap, leverage, nonlinearity and small sample biases.
9. **Power and interpretation:** run a sensitivity/power simulation for plausible interaction effects at observed N and covariate structure. If sample too small or biomarker missingness large, report uncertainty; do not hype null results or post hoc significance.
10. **Ethics:** Use only public deidentified dataset; never attempt reidentification. Keep limitations explicit: healthy volunteers, no ASPD, no moral outcomes, peripheral assays, between-subject design and possible age confounding.
11. **No clinical/ethical claims:** A connectivity difference is not proof of changed empathy, morality, antisocial traits, receptor occupancy, or therapeutic benefit.

## Execution stages and deliverables
### Stage 0 — repository bootstrap and immediate download
- Create a git repository with `README.md`, `environment.yml` or `pyproject.toml`, `src/`, `tests/`, `reports/`, `results/`, `literature/`, `data/README.md`, `scripts/` and `.gitignore` excluding participant MRI and private/sensitive paths.
- Check installed tools and disk. Prefer `datalad install -s https://github.com/OpenNeuroDatasets/ds004725.git data/ds004725`; `datalad get` only `participants.tsv`, `participants.json`, and the discovered `phenotype` table/dictionary first. If annex access fails, document errors and try official OpenNeuro supported download mechanisms instead of inventing mirrors.
- Verify the OpenNeuro version/commit, source URLs and file hashes when possible.
- Write `reports/data_inventory.md` with paths, sizes, dataset version and retrievability.

### Stage 1 — phenotype audit (no fMRI download required)
- Write a robust Python audit (`src/audit_phenotype.py`) that reads TSV as strings before type conversion; handles BIDS `n/a`; validates one-to-one participant IDs; checks treatment codes; summarizes age cohorts, sex, missingness by arm, plasma units/distributions and sampling-to-scan gaps.
- Outputs: `results/phenotype_missingness.csv`, `results/eligibility_flow.csv`, `reports/feasibility.md`, plots without identifiable participant IDs.
- Go/no-go after available N and biomarker spread are established. If no usable AVP, propose *explicitly labeled* fallback OT or imaging-only analysis rather than silently switching endpoint.

### Stage 2 — literature novelty audit
- Search `ds004725`, full dataset title, trial `NCT01823146`, PI and authors with combinations of `vasopressin`, `plasma`, `baseline`, `moderation`, `amygdala`, `mPFC`, `oxytocin`, `salience`, `connectivity`, `2026`. Create `literature/novelty_matrix.md` with source URLs/DOIs and overlap. Cite the 2022 salience study and 2026 dynamic connectivity study.
- If exact interaction already published, stop and propose a materially different preregistered estimand. Never claim novelty as certain.

### Stage 3 — image fetch and preprocessing
- After phenotype go/no-go, selectively fetch one MRI scan and validate its NIfTI header, modality paths and dependencies. Estimate disk/CPU before bulk retrieval. Then obtain eligible scans, retaining checksums and logs.
- Generate BIDS validation and preprocessing QC reports. Apply predefined exclusions without reference to treatment×biomarker outcomes.

### Stage 4 — ROI & statistics
- Write versioned ROI extraction scripts and deterministic model code. Add synthetic-data unit tests: subject join, medication labels, biomarker transforms, FC Fisher-z, interaction coefficient, missingness logic and anti-leakage.
- Produce all code, a `methods.md`, visualization/report script, tidy statistical output tables, reproducible command manifest, model diagnostics and limitations.

### Stage 5 — release
- Summarize verified sample, primary point estimate with CI/p, corrected secondary tests, negatives, sensitivity analyses and interpretation in `reports/research_report.md`. Add reproducibility instructions and exact provenance. State unrun steps as unrun.

## First-turn instructions for Claude Code
Start Stage 0 and Stage 1 *now*. Execute download commands if network tools available, inspect real filenames/columns and write the audit. Do **not** download all imaging data until the phenotype audit and storage check are complete. Do not merely outline steps: produce files/code and run tests, then report concrete results, blockers and next executable command.

## Local adaptation (added after Stage 0/1 audit, 2026-11-08)
- This repository IS the BIDS dataset (fork of OpenNeuroDatasets/ds004725 at
  commit `8372e9f`, tag **1.0.2**, not 1.0.1). Project code lives beside the BIDS
  files: `src/ tests/ scripts/ results/ reports/ literature/`, all listed in
  `.bidsignore`. Outputs go to `derivatives/`.
- The fork has no `git-annex` branch. Run `scripts/setup_annex.sh` once per
  clone to attach the OpenNeuro public S3 remote; then `datalad get <file>` works.
- Stage 1 findings that change the plan are in `reports/feasibility.md`
  (gate = AMBER). Read it before Stage 3. In particular: exclude the 2
  post-dosing blood draws, add `T x age` to every moderation model, and treat
  the AVP assay as of uncertain validity.
