# Data inventory and provenance (Stage 0)

| Item | Value |
|---|---|
| Dataset | OpenNeuro ds004725, "Single Dose Intranasal Oxytocin Administration: Data from Healthy Younger and Older Adults" |
| Version found | **1.0.2** (2026-04-20, Data in Brief DOI added). The handoff specified 1.0.1; the only 1.0.1 → 1.0.2 change in `CHANGES` is the README/DOI link |
| Upstream | https://github.com/OpenNeuroDatasets/ds004725 (tags 1.0.0, 1.0.1, 1.0.2) |
| Commit | `8372e9fad9e947131a9e4c76a53cc3284edeb4f4`. The fork `fffn259m2n-cell/ds004725` and upstream HEAD are identical |
| Dataset DOI | 10.18112/openneuro.ds004725.v1.0.2 |
| Descriptor | https://doi.org/10.1016/j.dib.2023.109669 |
| Licence | CC0 |
| Retrieval | `datalad clone` of upstream succeeded (DataLad 1.8.0, git-annex 10.20240129). The fork has no `git-annex` branch, so `scripts/setup_annex.sh` attaches the upstream annex and the `s3-PUBLIC` special remote. HTTP HEAD to the S3 object for `sub-12000` BOLD returned 200 |

## Files (metadata, stored directly in git, not annexed)

| Path | Bytes | SHA-256 |
|---|---|---|
| participants.tsv | 4961 | `187fca13…8887` |
| participants.json | 1383 | `2e141b24…329d` |
| phenotype/phenotypic_data.tsv | 11809 | `81db53e2…1f58` |
| phenotype/phenotypic_data.json | 4328 | `06ae5e30…7f258` |

The full hashes are in `results/input_sha256.txt`. The folder is named `phenotype/` (singular), whereas the descriptor says `phenotypes/`. The audit discovers it at runtime.

## Imaging (annexed, NOT downloaded)

| | Count | Total size (from annex keys) |
|---|---|---|
| `sub-*/func/*_task-rest_bold.nii` (uncompressed) | 87 | 10.11 GB |
| `sub-*/anat/*_T1w.nii` (uncompressed) | 87 | 1.76 GB |
| Subject-level BOLD JSON (TR = 2.5 s) | 2 (sub-11001, sub-12000) | – |

- Other files: root `task-rest_bold.json` (TR = 2 s, Philips Achieva 3T); no SliceTiming, PhaseEncodingDirection or fieldmaps. There are no derivatives in the dataset.
- Disk in this session: 27 GB free. That is enough for raw data, but not for fMRIPrep working directories on all 87 participants.

## Tools

Python 3.13, pandas, numpy, scipy, statsmodels, matplotlib, nibabel, pytest, DataLad 1.8.0, git-annex 10.20240129.
