#!/usr/bin/env bash
# Stage 3 step 1: fetch ONE participant's T1w + BOLD and print NIfTI headers.
# Usage: scripts/fetch_one_mri.sh sub-11002
set -euo pipefail
sub="${1:?participant id, e.g. sub-11002}"
cd "$(git rev-parse --show-toplevel)"
datalad get "$sub/anat/${sub}_T1w.nii" "$sub/func/${sub}_task-rest_bold.nii"
git annex fsck --fast "$sub" 
python3 -I - "$sub" <<'PY'
import sys, nibabel as nib
s = sys.argv[1]
for f in (f"{s}/anat/{s}_T1w.nii", f"{s}/func/{s}_task-rest_bold.nii"):
    im = nib.load(f); h = im.header
    print(f, im.shape, h.get_zooms(), h.get_xyzt_units(), nib.aff2axcodes(im.affine))
PY
