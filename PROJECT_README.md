# ds004725 secondary analysis: does screening plasma AVP/OT moderate intranasal-OT effects on resting-state FC?

The BIDS `README` belongs to the original dataset authors. This file describes the analysis project that sits on top of it. Read `CLAUDE.md` for the full plan.

| Stage | Status | Main outputs |
|---|---|---|
| 0 Bootstrap/download | done | `reports/data_inventory.md`, `scripts/setup_annex.sh` |
| 1 Phenotype audit | done, gate **AMBER** | `reports/feasibility.md`, `results/*.csv`, `reports/figures/` |
| 2 Literature novelty | **not run** | `literature/novelty_matrix.md` (to create) |
| 3-5 Imaging, stats, report | **not run** | – |

To reproduce Stage 1:

```bash
pip install -e .[test]
scripts/run_stage1.sh
```

This is a secondary analysis of public, de-identified data. It makes no clinical, moral or antisocial-trait claims, and plasma assays are not central neurotransmission.
