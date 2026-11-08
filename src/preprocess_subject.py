"""Pre-registered preprocessing (plan v1.0 section 4; OSF osf.io/8yqun).

Per participant:
 1. discard 4 volumes;  2. ANTs rigid realignment to the mean volume;
 3. Power FD (50 mm);   4. mean EPI -> MNI152 2 mm template, ANTs rigid +
    affine + SyN, Mattes MI;  5. masks (rois/) to native space with the
    inverse transform, intersected with a native EPI mask;
 6. simultaneous nuisance regression + 0.01-0.08 Hz band-pass + linear
    detrend: 24 motion regressors, WM and CSF signals + derivatives;
    censoring FD > 0.5 mm plus the following volume;
 7. mean ROI time series in native space; Fisher-z Pearson FC.
No slice-timing correction, no smoothing.

Outputs in derivatives/prereg_v1/<sub>/:
    <sub>_qc.json      QC metrics ONLY (read by group_analysis --qc-only)
    <sub>_fc.json      connectivity values (read only after QC is frozen)
    <sub>_timeseries.tsv, <sub>_confounds.tsv, <sub>_qc_reg.png
Usage: python src/preprocess_subject.py sub-XXXXX [--threads 2]
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

N_DUMMY = 4
HEAD_RADIUS = 50.0
BAND = (0.01, 0.08)
OUTDIR = Path("derivatives/prereg_v1")
TPL = Path("derivatives/templates/tpl-MNI152NLin6Asym_res-02_T1w.nii.gz")
TPL_MASK = Path("derivatives/templates/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz")
ROI_DIR = Path("rois")
ROIS = ["AMY", "AMY_L", "AMY_R", "MPFC", "AINS"]
NUIS = ["WM", "CSF"]
VERSION = "prereg-v1.0-impl1"


# ------------------------------------------------------------ helpers ------

def rigid_params(tx_files) -> np.ndarray:
    """Per-volume [rot_x, rot_y, rot_z (rad), trans_x, trans_y, trans_z (mm)]."""
    import ants
    from scipy.spatial.transform import Rotation as R
    rows = []
    for f in tx_files:
        t = ants.read_transform(f[0] if isinstance(f, (list, tuple)) else f)
        p = np.asarray(t.parameters)
        if len(p) == 12:
            rot = R.from_matrix(p[:9].reshape(3, 3)).as_euler("xyz")
            trans = p[9:12]
        elif t.type.startswith("Euler3D"):
            rot, trans = p[:3], p[3:6]
        else:
            v = p[:3]
            rot = R.from_quat([*v, np.sqrt(max(0.0, 1 - (v ** 2).sum()))]).as_euler("xyz")
            trans = p[3:6]
        rows.append([*rot, *trans])
    return np.asarray(rows)


def power_fd(params: np.ndarray, radius: float = HEAD_RADIUS) -> np.ndarray:
    d = np.abs(np.diff(params, axis=0))
    return np.concatenate([[0.0], d[:, :3].sum(1) * radius + d[:, 3:].sum(1)])


def censor_keep(fd: np.ndarray, thresh: float | None) -> np.ndarray:
    """True = keep. Censor FD > thresh and the FOLLOWING volume (plan s4.6)."""
    if thresh is None:
        return np.ones(len(fd), bool)
    bad = fd > thresh
    bad[1:] |= bad[:-1].copy()
    return ~bad


def deriv(x: np.ndarray) -> np.ndarray:
    x = np.atleast_2d(x.T).T
    return np.vstack([np.zeros((1, x.shape[1])), np.diff(x, axis=0)])


def motion24(params: np.ndarray) -> np.ndarray:
    d = deriv(params)
    return np.column_stack([params, d, params ** 2, d ** 2])


def clean_fc(ts: pd.DataFrame, conf: np.ndarray, keep: np.ndarray, tr: float) -> dict:
    from nilearn import signal
    c = signal.clean(ts.to_numpy(), detrend=True, standardize=None,
                     confounds=conf, low_pass=BAND[1], high_pass=BAND[0],
                     t_r=tr, sample_mask=np.where(keep)[0])
    r = np.corrcoef(c.T)
    i = {k: n for n, k in enumerate(ts.columns)}
    z = lambda a, b: float(np.arctanh(r[i[a], i[b]]))  # noqa: E731
    return {"AMY_MPFC": z("AMY", "MPFC"), "AMY_AINS": z("AMY", "AINS"),
            "AMYL_MPFC": z("AMY_L", "MPFC"), "AMYR_MPFC": z("AMY_R", "MPFC")}


# --------------------------------------------------------------- main ------

def run(sub: str, threads: int = 2) -> dict:
    import ants
    import nibabel as nib
    from nilearn.masking import compute_epi_mask
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(threads)
    t0 = time.time()
    out = OUTDIR / sub
    out.mkdir(parents=True, exist_ok=True)
    side = Path(sub) / "func" / f"{sub}_task-rest_bold.json"
    tr = float(json.loads((side if side.exists() else Path("task-rest_bold.json"))
                          .read_text())["RepetitionTime"])
    qc = {"participant_id": sub, "version": VERSION, "tr": tr, "status": "ok"}
    try:
        with tempfile.TemporaryDirectory() as td:
            im = nib.load(str(Path(sub) / "func" / f"{sub}_task-rest_bold.nii"))
            qc["n_vols_raw"] = int(im.shape[3])
            qc["voxel_size"] = [float(v) for v in im.header.get_zooms()[:3]]
            trim = nib.Nifti1Image(np.asarray(im.dataobj)[..., N_DUMMY:].astype(np.float32),
                                   im.affine, im.header)
            trim.header.set_data_dtype(np.float32)
            nib.save(trim, os.path.join(td, "bold.nii.gz"))
            bold = ants.image_read(os.path.join(td, "bold.nii.gz"))
        # 2-3 realignment + FD
        mean0 = ants.get_average_of_timeseries(bold)
        mc = ants.motion_correction(bold, fixed=mean0, type_of_transform="Rigid")
        params = rigid_params(mc["motion_parameters"])
        fd = power_fd(params)
        mean_mc = ants.get_average_of_timeseries(mc["motion_corrected"])
        data = mc["motion_corrected"].numpy()  # x,y,z,t native
        # 4 registration mean EPI -> template
        tpl = ants.image_read(str(TPL))
        mean_n4 = ants.n4_bias_field_correction(mean_mc)
        reg = ants.registration(fixed=tpl, moving=mean_n4, type_of_transform="SyN",
                                aff_metric="mattes", syn_metric="mattes", random_seed=1)
        warped = reg["warpedmovout"]
        tm = ants.image_read(str(TPL_MASK)).numpy() > 0
        qc["reg_r"] = float(np.corrcoef(warped.numpy()[tm], tpl.numpy()[tm])[0, 1])
        # 5 masks to native, intersect with EPI mask
        mean_nii = nib.Nifti1Image(mean_mc.numpy(), _ants_affine(mean_mc))
        epi = np.asarray(compute_epi_mask(mean_nii).dataobj) > 0
        tm_native = ants.apply_transforms(fixed=mean_mc, moving=ants.image_read(str(TPL_MASK)),
                                          transformlist=reg["invtransforms"],
                                          interpolator="nearestNeighbor").numpy() > 0
        qc["epi_vs_template_mask_dice"] = float(2 * (epi & tm_native).sum() /
                                                (epi.sum() + tm_native.sum()))
        native, cov = {}, {}
        for k in ROIS + NUIS:
            m = ants.apply_transforms(
                fixed=mean_mc,
                moving=ants.image_read(str(ROI_DIR / f"space-MNI152NLin6Asym_res-02_desc-{k}_mask.nii.gz")),
                transformlist=reg["invtransforms"], interpolator="nearestNeighbor").numpy() > 0
            native[k] = m & epi
            cov[k] = float(native[k].sum() / m.sum()) if m.sum() else 0.0
            qc[f"n_vox_{k}"] = int(native[k].sum())
        qc["coverage"] = cov
        # 6-7 signals
        ts = pd.DataFrame({k: data[native[k]].mean(0) if native[k].any() else np.full(data.shape[3], np.nan)
                           for k in ROIS})
        conf = pd.DataFrame(params, columns=["rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z"])
        conf["fd"] = fd
        for k in NUIS:
            conf[k] = data[native[k]].mean(0) if native[k].any() else np.nan
        conf["GS"] = data[epi].mean(0)
        del data
        keep = censor_keep(fd, 0.5)
        qc.update(mean_fd=float(fd.mean()), max_fd=float(fd.max()),
                  n_vols=int(len(fd)), n_keep=int(keep.sum()),
                  minutes_retained=float(keep.sum() * tr / 60),
                  prop_censored=float(1 - keep.mean()),
                  n_keep_fd03=int(censor_keep(fd, 0.3).sum()))
        finite = bool(np.isfinite(ts.to_numpy()).all() and np.isfinite(conf[NUIS].to_numpy()).all())
        qc["signals_finite"] = finite
        base = np.column_stack([motion24(params), conf[NUIS].to_numpy(), deriv(conf[NUIS].to_numpy())])
        gsr = np.column_stack([base, conf[["GS"]].to_numpy(), deriv(conf[["GS"]].to_numpy())])
        fc = {"participant_id": sub}
        if finite:
            fc["primary"] = clean_fc(ts, base, keep, tr)
            fc["gsr"] = clean_fc(ts, gsr, keep, tr)
            fc["nocensor"] = clean_fc(ts, base, censor_keep(fd, None), tr)
            k3 = censor_keep(fd, 0.3)
            fc["fd03"] = clean_fc(ts, base, k3, tr) if k3.sum() > 30 else None
        ts.to_csv(out / f"{sub}_timeseries.tsv", sep="\t", index=False, float_format="%.4f")
        conf.to_csv(out / f"{sub}_confounds.tsv", sep="\t", index=False, float_format="%.6f")
        (out / f"{sub}_fc.json").write_text(json.dumps(fc, indent=1))
        _qc_png(warped, tpl, out / f"{sub}_qc_reg.png", sub)
    except Exception as e:  # logged as preprocessing failure (plan s5)
        qc["status"] = f"failed: {type(e).__name__}: {e}"
    import ants as _a
    import nilearn as _n
    qc["software"] = {"antspyx": _a.__version__, "nilearn": _n.__version__}
    qc["runtime_s"] = round(time.time() - t0, 1)
    (out / f"{sub}_qc.json").write_text(json.dumps(qc, indent=1))
    return qc


def _ants_affine(img) -> np.ndarray:
    """NIfTI (RAS) affine for an ANTs (LPS) image."""
    D = np.asarray(img.direction) * np.asarray(img.spacing)
    A = np.eye(4)
    A[:3, :3] = D
    A[:3, 3] = img.origin
    return np.diag([-1, -1, 1, 1]) @ A


def _qc_png(warped, tpl, path, sub):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    w, t = warped.numpy(), tpl.numpy()
    lev = np.percentile(t[t > 0], 60)
    fig, ax = plt.subplots(1, 3, figsize=(9, 3.2), facecolor="black")
    for a, sl, tl in zip(ax, (w[45], w[:, 63], w[:, :, 32]), (t[45], t[:, 63], t[:, :, 32])):
        a.imshow(np.rot90(sl), cmap="gray")
        a.contour(np.rot90(tl), levels=[lev], colors="#eb6834", linewidths=0.6)
        a.axis("off")
    ax[0].set_title(f"{sub}: mean EPI in MNI, template contour", color="white", fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=70, facecolor="black")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sub")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    q = run(a.sub, a.threads)
    print(json.dumps({k: q.get(k) for k in ("participant_id", "status", "mean_fd", "minutes_retained",
                                             "reg_r", "epi_vs_template_mask_dice", "coverage", "runtime_s")}))
