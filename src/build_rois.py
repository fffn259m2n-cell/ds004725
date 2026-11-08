"""Build the pre-registered ROI and nuisance masks (plan v1.0 section 3).

All masks are in MNI152NLin6Asym, 2 mm (TemplateFlow). Writes rois/*.nii.gz and
rois/rois.json (voxel counts, centroids, SHA-256). Usage: python -I src/build_rois.py
"""
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion

T = Path("derivatives/templates")
OUT = Path("rois")
HOSPA = T / "tpl-MNI152NLin6Asym_res-02_atlas-HOSPA_desc-th25_dseg.nii.gz"
HOCPA = T / "tpl-MNI152NLin6Asym_res-02_atlas-HOCPA_desc-th25_dseg.nii.gz"
PROB = T / "tpl-MNI152NLin6Asym_res-02_atlas-HOSPA_probseg.nii.gz"
MPFC_XYZ, MPFC_R = (-1.0, 49.0, -5.0), 10.0
# HOSPA image labels (TemplateFlow TSV indices are 1-based for HOSPA)
AMY_L, AMY_R = 10, 20
HOSPA_GM = [2, 4, 5, 6, 7, 9, 10, 11, 13, 15, 16, 17, 18, 19, 20, 21]
# HOCPA image label = TSV index + 1 (verified by centroid): Insular Cortex = 2
INSULA = 2
# HOSPA probseg channel = label - 1; stored values are probability / 100
CH_WM, CH_LV = (0, 11), (2, 13)


def mni_coords(img):
    ijk = np.indices(img.shape[:3]).reshape(3, -1).T
    return nib.affines.apply_affine(img.affine, ijk).reshape(*img.shape[:3], 3)


def build():
    ref = nib.load(HOSPA)
    sub = np.asarray(ref.dataobj)
    cort = np.asarray(nib.load(HOCPA).dataobj)
    xyz = mni_coords(ref)
    prob = np.asarray(nib.load(PROB).dataobj).astype(float) * 100.0
    masks = {
        "AMY": np.isin(sub, [AMY_L, AMY_R]),
        "AMY_L": sub == AMY_L,
        "AMY_R": sub == AMY_R,
        "MPFC": np.linalg.norm(xyz - np.array(MPFC_XYZ), axis=-1) <= MPFC_R,
        "AINS": (cort == INSULA) & (xyz[..., 1] > 0),
        "WM": binary_erosion(prob[..., list(CH_WM)].sum(-1) > 0.9),
        "CSF": (prob[..., list(CH_LV)].sum(-1) > 0.5) & ~np.isin(sub, HOSPA_GM),
    }
    OUT.mkdir(exist_ok=True)
    info = {}
    for k, m in masks.items():
        p = OUT / f"space-MNI152NLin6Asym_res-02_desc-{k}_mask.nii.gz"
        nib.save(nib.Nifti1Image(m.astype(np.uint8), ref.affine), p)
        info[k] = {"file": p.name, "n_vox": int(m.sum()),
                   "centroid_mni": np.round(xyz[m].mean(0), 1).tolist(),
                   "sha256_data": hashlib.sha256(m.astype(np.uint8).tobytes()).hexdigest()}
    (OUT / "rois.json").write_text(json.dumps(info, indent=2))
    return info


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
