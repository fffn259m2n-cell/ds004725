#!/usr/bin/env bash
# Plan v1.0: preprocess ALL 87 participants (H1 uses everyone passing QC).
# Fetch BOLD (HTTPS + SHA-256 verify) -> preprocess_subject.py -> drop raw.
# Usage: scripts/run_prereg_all.sh [n_parallel] [threads_each]
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
NP=${1:-2}; TH=${2:-2}
mkdir -p derivatives/prereg_v1/logs
tail -n +2 participants.tsv | cut -f1 | tr -d '\r' > derivatives/prereg_v1/subjects.txt
one() {
  s=$1; TH=$2
  [ -f "derivatives/prereg_v1/$s/${s}_qc.json" ] && { echo "skip $s"; return; }
  f="$s/func/${s}_task-rest_bold.nii"
  scripts/fetch_mri.sh "$f" >> "derivatives/prereg_v1/logs/$s.log" 2>&1 \
    && python3 src/preprocess_subject.py "$s" --threads "$TH" >> "derivatives/prereg_v1/logs/$s.log" 2>&1 \
    && echo "done $s" || echo "FAIL $s"
  git annex drop --force "$f" >/dev/null 2>&1
}
export -f one
xargs -a derivatives/prereg_v1/subjects.txt -P "$NP" -I{} bash -c 'one {} '"$TH"
echo ALL_FINISHED
