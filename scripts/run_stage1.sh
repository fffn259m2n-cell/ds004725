#!/usr/bin/env bash
# Reproduce every Stage 1 output from the dataset root.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
python3 -m pytest -q tests
python3 -I src/audit_phenotype.py --bids-root . --out-root . > results/audit_stdout.txt
python3 -I src/power_sim.py
sha256sum participants.tsv participants.json phenotype/phenotypic_data.tsv phenotype/phenotypic_data.json > results/input_sha256.txt
