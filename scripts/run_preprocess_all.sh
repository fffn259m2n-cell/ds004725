#!/usr/bin/env bash
# Fetch -> verify -> preprocess -> drop raw, for every Stage-1-eligible
# participant (plus the unknown-interval participant for sensitivity 9).
# Usage: scripts/run_preprocess_all.sh [n_parallel] [threads_each]
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
NP=${1:-2}; TH=${2:-2}
mkdir -p derivatives/oxtfc/logs
python3 -I - > derivatives/oxtfc/subjects.txt <<'PY'
import pandas as pd
d = pd.read_csv("derivatives/phenotype_audit/analysis_frame.tsv", sep="\t", na_values="n/a")
ok = d.avp_pg_ml.notna() & d.bold_present & ((d.days_blood_to_scan > 0) | d.days_blood_to_scan.isna())
print("\n".join(sorted(d.loc[ok, "participant_id"])))
PY
one() {
  s=$1; TH=$2
  [ -f "derivatives/oxtfc/$s/${s}_qc.json" ] && { echo "skip $s"; return; }
  scripts/fetch_mri.sh "$s/anat/${s}_T1w.nii" "$s/func/${s}_task-rest_bold.nii" >> "derivatives/oxtfc/logs/$s.log" 2>&1 \
   && python3 src/preprocess.py "$s" --threads "$TH" >> "derivatives/oxtfc/logs/$s.log" 2>&1 \
   && echo "done $s" || echo "FAIL $s"
  git annex drop --force "$s/anat/${s}_T1w.nii" "$s/func/${s}_task-rest_bold.nii" >/dev/null 2>&1
}
export -f one
wc -l < derivatives/oxtfc/subjects.txt
xargs -a derivatives/oxtfc/subjects.txt -P "$NP" -I{} bash -c 'one {} '"$TH"
echo ALL_FINISHED
