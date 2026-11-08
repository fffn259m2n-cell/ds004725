"""Resting-state preprocessing and ROI FC extraction for one participant.

Implements reports/preregistration.md section 4 with ANTsPy + nilearn
(fMRIPrep unavailable in this environment; deviation D1).

Usage:
    python -I src/preprocess.py sub-11002 [--root .] [--threads 2]
Outputs (derivatives/oxtfc/<sub>/):
    <sub>_timeseries.tsv   ROI mean signals in MNI space (raw, pre-cleaning)
    <sub>_confounds.tsv    motion params, FD, WM/CSF/global signals
    <sub>_qc.json          QC metrics + FC (primary pipeline and GSR)
    <sub>_qc.png           registration overlay
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

N_DUMMY = 4
FD_THRESH = 0.5
FD_THRESH_STRICT = 0.3
HEAD_RADIUS = 50.0
BAND = (0.01, 0.08)
TPL = "derivatives/templates/tpl-MNI152NLin6Asym_res-02_{}.nii.gz"

# Harvard-Oxford image labels (HOCPA label = TSV index + 1; see prereg s3)
ROI_DEF = {
    "amygdala": ("HOSPA", [10, 20], None),
    "mpfc": ("HOCPA", [25], None),
    "antinsula": ("HOCPA", [2], "y>=0"),
}


def roi_masks(root: Path) -> dict[str, np.ndarray]:
    import nibabel as nib
    out = {}
    for name, (atlas, labels, rule) in ROI_DEF.items():
        im = nib.load(root / TPL.format(f"atlas-{atlas}_desc-th25_dseg"))
        d = np.asarray(im.dataobj)
        m = np.isin(d, labels)
        if rule == "y>=0":
            ijk = np.indices(d.shape).reshape(3, -1).T
            y = nib.affines.apply_affine(im.affine, ijk)[:, 1].reshape(d.shape)
            m &= y >= 0
        out[name] = m
    return out


def rigid_params(tx_files: list[str]) -> np.ndarray:
    """6 rigid parameters per volume: 3 rotations (rad), 3 translations (mm)."""
    import ants
    from scipy.spatial.transform import Rotation as R
    rows = []
    for f in tx_files:
        t = ants.read_transform(f[0] if isinstance(f, (list, tuple)) else f)
        mat = np.asarray(t.parameters[:9]).reshape(3, 3) \
            if len(t.parameters) == 12 else None
        p = np.asarray(t.parameters)
        if mat is not None:  # AffineTransform / MatrixOffset form
            rot = R.from_matrix(mat).as_euler("xyz")
            trans = p[9:12]
        elif t.type.startswith("Euler3D"):
            rot, trans = p[:3], p[3:6]
        else:  # VersorRigid3D: versor (3) + translation (3)
            v = p[:3]
            w = np.sqrt(max(0.0, 1 - (v ** 2).sum()))
            rot = R.from_quat([*v, w]).as_euler("xyz")
            trans = p[3:6]
        rows.append([*rot, *trans])
    return np.asarray(rows)


def power_fd(params: np.ndarray, radius: float = HEAD_RADIUS) -> np.ndarray:
    d = np.abs(np.diff(params, axis=0))
    fd = d[:, :3].sum(1) * radius + d[:, 3:].sum(1)
    return np.concatenate([[0.0], fd])


def censor_mask(fd: np.ndarray, thresh: float) -> np.ndarray:
    """True = keep. Drop FD>thresh volumes plus 1 before and 2 after."""
    bad = fd > thresh
    out = bad.copy()
    idx = np.where(bad)[0]
    for i in idx:
        out[max(0, i - 1):min(len(fd), i + 3)] = True
    return ~out


def deriv(x: np.ndarray) -> np.ndarray:
    return np.vstack([np.zeros((1, x.shape[1])), np.diff(x, axis=0)])


def fisher_fc(ts: pd.DataFrame, confounds: np.ndarray, keep: np.ndarray,
              tr: float) -> dict:
    from nilearn import signal
    clean = signal.clean(ts.to_numpy(), detrend=True, standardize=False,
                         confounds=confounds, low_pass=BAND[1],
                         high_pass=BAND[0], t_r=tr,
                         sample_mask=np.where(keep)[0])
    c = pd.DataFrame(clean, columns=ts.columns)
    r = c.corr().to_numpy()
    cols = list(ts.columns)
    a = cols.index("amygdala")
    return {f"amygdala_{k}": float(np.arctanh(r[a, cols.index(k)]))
            for k in ("mpfc", "antinsula")}


def run(sub: str, root: Path, threads: int) -> dict:
    import ants
    import nibabel as nib
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(threads)
    t0 = time.time()
    out = root / "derivatives" / "oxtfc" / sub
    out.mkdir(parents=True, exist_ok=True)
    side = root / sub / "func" / f"{sub}_task-rest_bold.json"
    tr = float(json.loads((side if side.exists() else root /
                           "task-rest_bold.json").read_text())["RepetitionTime"])

    tpl = ants.image_read(str(root / TPL.format("T1w")))
    tpl_brain = ants.image_read(str(root / TPL.format("desc-brain_T1w")))
    tpl_mask = ants.image_read(str(root / TPL.format("desc-brain_mask")))

    # ---- anatomical
    t1 = ants.image_read(str(root / sub / "anat" / f"{sub}_T1w.nii"))
    t1 = ants.n4_bias_field_correction(t1)
    reg = ants.registration(fixed=tpl, moving=t1, type_of_transform="SyN",
                            random_seed=1)
    mask = ants.apply_transforms(fixed=t1, moving=tpl_mask,
                                 transformlist=reg["invtransforms"],
                                 interpolator="nearestNeighbor")
    t1b = t1 * mask
    seg = ants.atropos(a=t1, x=mask, i="kmeans[3]", m="[0.2,1x1x1]",
                       c="[5,0]")
    csf_p, _, wm_p = seg["probabilityimages"]
    wm = ants.iMath(ants.threshold_image(wm_p, 0.9, 1.1), "ME", 1)
    csf = ants.iMath(ants.threshold_image(csf_p, 0.9, 1.1), "ME", 1)
    wm_mni = ants.apply_transforms(tpl, wm, reg["fwdtransforms"],
                                   interpolator="nearestNeighbor").numpy() > 0
    csf_mni = ants.apply_transforms(tpl, csf, reg["fwdtransforms"],
                                    interpolator="nearestNeighbor").numpy() > 0
    t1_mni = ants.apply_transforms(tpl, t1b, reg["fwdtransforms"])
    m = tpl_mask.numpy() > 0
    norm_r = float(np.corrcoef(t1_mni.numpy()[m], tpl_brain.numpy()[m])[0, 1])

    # ---- functional
    with tempfile.TemporaryDirectory() as td:
        im = nib.load(str(root / sub / "func" / f"{sub}_task-rest_bold.nii"))
        n_raw = im.shape[3]
        trim = nib.Nifti1Image(np.asarray(im.dataobj)[..., N_DUMMY:]
                               .astype(np.float32), im.affine, im.header)
        trim.header.set_data_dtype(np.float32)
        p = os.path.join(td, "bold.nii.gz")
        nib.save(trim, p)
        bold = ants.image_read(p)
    mean0 = ants.get_average_of_timeseries(bold)
    mc = ants.motion_correction(bold, fixed=mean0, type_of_transform="Rigid")
    params = rigid_params(mc["motion_parameters"])
    fd = power_fd(params)
    mean_mc = ants.get_average_of_timeseries(mc["motion_corrected"])
    co = ants.registration(fixed=t1b, moving=ants.n4_bias_field_correction(mean_mc),
                           type_of_transform="Rigid", aff_metric="mattes",
                           random_seed=1)
    txl = reg["fwdtransforms"] + co["fwdtransforms"]
    bold_mni = ants.apply_transforms(fixed=tpl, moving=mc["motion_corrected"],
                                     transformlist=txl, imagetype=3)
    mean_mni = ants.apply_transforms(fixed=tpl, moving=mean_mc,
                                     transformlist=txl)
    data = bold_mni.numpy()  # x,y,z,t
    del bold_mni

    from nilearn.masking import compute_epi_mask
    mean_nii = nib.Nifti1Image(mean_mni.numpy(), nib.load(
        str(root / TPL.format("T1w"))).affine)
    cover = (np.asarray(compute_epi_mask(mean_nii).dataobj) > 0) & m
    rois = roi_masks(root)
    ts, cov, tsnr = {}, {}, {}
    for k, rm in rois.items():
        use = rm & cover
        cov[k] = float(use.sum() / rm.sum())
        x = data[use]  # V x T
        ts[k] = x.mean(0)
        tsnr[k] = float(np.median(x.mean(1) / (x.std(1) + 1e-6)))
    wm_u, csf_u = wm_mni & cover, csf_mni & cover
    sig = pd.DataFrame(ts)
    conf_df = pd.DataFrame(params, columns=["rot_x", "rot_y", "rot_z",
                                            "trans_x", "trans_y", "trans_z"])
    conf_df["fd"] = fd
    conf_df["wm"] = data[wm_u].mean(0)
    conf_df["csf"] = data[csf_u].mean(0)
    conf_df["global"] = data[cover].mean(0)
    whole_tsnr = float(np.median(data[cover].mean(1) / (data[cover].std(1) + 1e-6)))
    del data

    base = np.column_stack([params, deriv(params),
                            conf_df[["wm", "csf"]].to_numpy(),
                            deriv(conf_df[["wm", "csf"]].to_numpy())])
    gsr = np.column_stack([base, conf_df[["global"]].to_numpy(),
                           deriv(conf_df[["global"]].to_numpy())])
    keep = censor_mask(fd, FD_THRESH)
    keep_s = censor_mask(fd, FD_THRESH_STRICT)
    fc = fisher_fc(sig, base, keep, tr)
    fc_gsr = fisher_fc(sig, gsr, keep, tr)
    fc_strict = fisher_fc(sig, base, keep_s, tr) if keep_s.sum() > 20 else {}

    sig.to_csv(out / f"{sub}_timeseries.tsv", sep="\t", index=False,
               float_format="%.4f")
    conf_df.to_csv(out / f"{sub}_confounds.tsv", sep="\t", index=False,
                   float_format="%.6f")
    qc = {
        "participant_id": sub, "tr": tr, "n_vols_raw": int(n_raw),
        "n_vols_after_dummy": int(len(fd)), "mean_fd": float(fd.mean()),
        "max_fd": float(fd.max()), "n_keep": int(keep.sum()),
        "n_keep_strict": int(keep_s.sum()),
        "mean_fd_strict_flag": bool(fd.mean() > FD_THRESH_STRICT),
        "roi_coverage": cov, "roi_tsnr": tsnr, "whole_brain_tsnr": whole_tsnr,
        "norm_r": norm_r, "n_wm_vox": int(wm_u.sum()), "n_csf_vox": int(csf_u.sum()),
        "fc": fc, "fc_gsr": fc_gsr, "fc_strict": fc_strict,
        "versions": {"antspyx": ants.__version__,
                     "nilearn": __import__("nilearn").__version__},
        "runtime_s": round(time.time() - t0, 1),
    }
    (out / f"{sub}_qc.json").write_text(json.dumps(qc, indent=2))
    _qc_png(t1_mni, mean_mni, tpl_brain, out / f"{sub}_qc.png", sub)
    return qc


def _qc_png(t1_mni, mean_mni, tpl_brain, path, sub):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t1, ep, tp = t1_mni.numpy(), mean_mni.numpy(), tpl_brain.numpy()
    fig, ax = plt.subplots(2, 3, figsize=(9, 6), facecolor="black")
    cuts = [("x", 45), ("y", 63), ("z", 36)]
    for j, (axn, c) in enumerate(cuts):
        for i, img in enumerate([t1, ep]):
            sl = {"x": img[c], "y": img[:, c], "z": img[:, :, c]}[axn]
            tsl = {"x": tp[c], "y": tp[:, c], "z": tp[:, :, c]}[axn]
            a = ax[i, j]
            a.imshow(np.rot90(sl), cmap="gray")
            a.contour(np.rot90(tsl), levels=[np.percentile(tp[tp > 0], 40)],
                      colors="#eb6834", linewidths=0.6)
            a.axis("off")
    ax[0, 0].set_title(f"{sub}: T1w (top), mean EPI (bottom) in MNI; "
                       "template contour", color="white", fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=70, facecolor="black")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sub")
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    qc = run(a.sub, a.root.resolve(), a.threads)
    print(json.dumps({k: qc[k] for k in ("participant_id", "mean_fd", "n_keep",
                                          "roi_coverage", "norm_r", "runtime_s")}))
