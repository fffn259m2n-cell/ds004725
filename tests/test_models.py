"""Synthetic-data tests for preprocessing helpers and moderation models."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import models as M  # noqa: E402
import preprocess as P  # noqa: E402


def synth(n=84, beta_int=0.0, seed=0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({
        "T": rng.permutation(np.r_[np.ones(n // 2), np.zeros(n - n // 2)]),
        "zmod": rng.standard_normal(n), "age_c": rng.standard_normal(n),
        "sex": rng.integers(0, 2, n).astype(float), "fd_c": rng.normal(0, .1, n)})
    d["FC"] = 0.2 * d["T"] + 0.1 * d.zmod + beta_int * d["T"] * d.zmod \
        + 0.1 * d.age_c + rng.standard_normal(n) * 0.3
    return d


def test_interaction_recovered():
    d = synth(n=400, beta_int=0.3, seed=1)
    r = M.fit_interaction(d, n_perm=200)
    assert r["ci_low"] < 0.3 < r["ci_high"]
    assert r["p_perm"] < 0.01
    # conditional effect at z=+1 = beta_T + beta_int
    assert r["cond_T_at_z+1"] == pytest.approx(r["beta_T"] + r["beta"])


def test_freedman_lane_null_calibrated():
    ps = [M.freedman_lane(synth(seed=s), M.PRIMARY, M.REDUCED, "T:zmod",
                          n_perm=199, seed=s) for s in range(40)]
    assert np.mean(np.array(ps) < 0.05) < 0.2
    assert 0.35 < np.mean(ps) < 0.65


def test_holm():
    assert M.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


def test_fisher_z_and_censor():
    fd = np.zeros(20); fd[10] = 0.9
    keep = P.censor_mask(fd, 0.5)
    assert keep.sum() == 16 and not keep[9:13].any() and keep[13]
    assert np.arctanh(0.5) == pytest.approx(0.5493, abs=1e-4)


def test_power_fd():
    params = np.zeros((3, 6)); params[1, 3] = 1.0; params[2, 0] = 0.01
    fd = P.power_fd(params)
    assert fd.tolist() == pytest.approx([0, 1.0, 1.0 + 0.5])


@pytest.mark.skipif(not (ROOT / "derivatives/templates/tpl-MNI152NLin6Asym_res-02_atlas-HOCPA_desc-th25_dseg.nii.gz").exists(),
                    reason="templates not downloaded")
def test_roi_centroids_anatomically_correct():
    import nibabel as nib
    rois = P.roi_masks(ROOT)
    aff = nib.load(ROOT / P.TPL.format("atlas-HOSPA_desc-th25_dseg")).affine
    c = {k: nib.affines.apply_affine(aff, np.argwhere(v)).mean(0) for k, v in rois.items()}
    assert abs(c["amygdala"][1] + 4) < 5 and abs(c["amygdala"][2] + 18) < 5   # y~-4, z~-18
    assert c["mpfc"][1] > 30 and c["mpfc"][2] < 0 and abs(c["mpfc"][0]) < 5   # vmPFC
    assert c["antinsula"][1] >= 0 and (np.argwhere(rois["antinsula"]).shape[0] > 200)


def test_qc_rules():
    q = pd.DataFrame({"participant_id": list("abcd"), "mean_fd": [.2, .6, .2, .2],
                      "n_keep": [200, 200, 100, 200], "norm_r": [.7, .7, .7, .7],
                      "cov_amygdala": [1, 1, 1, .4], "cov_mpfc": [1, 1, 1, 1]})
    r = M.qc_exclusions(q)
    assert r.qc_pass.tolist() == [True, False, False, False]
