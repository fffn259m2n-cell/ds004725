"""Tests for the registered pipeline (preprocess_subject, group_analysis, build_rois)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import group_analysis as G  # noqa: E402
import preprocess_subject as P  # noqa: E402


def test_censor_following_volume_only():
    fd = np.zeros(10); fd[3] = 0.6; fd[9] = 0.7
    keep = P.censor_keep(fd, 0.5)
    assert np.where(~keep)[0].tolist() == [3, 4, 9]
    assert P.censor_keep(fd, None).all()


def test_motion24_shape_and_content():
    p = np.arange(18, dtype=float).reshape(3, 6)
    m = P.motion24(p)
    assert m.shape == (3, 24)
    assert np.allclose(m[:, 12:18], p ** 2) and np.allclose(m[0, 6:12], 0)


def _synth(n=83, b3=0.0, seed=0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({"T": rng.choice([-.5, .5], n), "C": rng.choice([-.5, .5], n),
                      "S": rng.choice([-.5, .5], n), "zM": rng.standard_normal(n),
                      "mean_fd": rng.uniform(.1, .4, n)})
    d["FC"] = .1 * d["T"] + .2 * d.C + b3 * d["T"] * d.zM + rng.standard_normal(n) * .3
    return d


def test_hc3_matches_statsmodels():
    import statsmodels.api as sm
    y, X = G.design(G.H2_F, _synth(seed=3))
    b, se, *_ = G.hc3_fit(y, X)
    m = sm.OLS(y, X).fit(cov_type="HC3")
    assert np.allclose(b, m.params) and np.allclose(se, m.bse)


def test_freedman_lane_power_and_null():
    y, X = G.design(G.H2_F, _synth(n=400, b3=0.4, seed=1))
    assert G.freedman_lane(y, X, "T:zM", n_perm=300) < 0.01
    ps = [G.freedman_lane(*G.design(G.H2_F, _synth(seed=s)), "T:zM", n_perm=199, seed=s)
          for s in range(30)]
    assert 0.3 < np.mean(ps) < 0.7


def test_bf01_direction():
    y, X = G.design(G.H2_F, _synth(n=300, b3=0.0, seed=2))
    assert G.bf01_bic(y, X, "T:zM") > 1
    y, X = G.design(G.H2_F, _synth(n=300, b3=0.5, seed=2))
    assert G.bf01_bic(y, X, "T:zM") < 1 / 3


def test_exclusion_rules():
    q = pd.DataFrame({"participant_id": list("abcdef"), "status": ["ok"] * 5 + ["failed: x"],
                      "mean_fd": [.2, .6, .2, .2, .2, .2], "minutes_retained": [7, 7, 4, 7, 7, 7],
                      "reg_r": [.5, .5, .5, .5, .5, .5], "signals_finite": [True] * 6,
                      "cov_AMY": [1, 1, 1, .4, 1, 1], "cov_MPFC": [1] * 6, "cov_AINS": [1, 1, 1, 1, .3, 1]})
    vis = pd.DataFrame({"participant_id": list("abcdef"), "visual_pass": ["true"] * 6})
    r = G.apply_exclusions(q, vis)
    assert r.pass_main.tolist() == [True, False, False, False, True, False]
    assert r.pass_ains.tolist() == [True, False, False, False, False, False]


def test_holm():
    assert G.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


@pytest.mark.skipif(not (ROOT / "rois/rois.json").exists(), reason="ROIs not built")
def test_rois_match_plan():
    import json
    r = json.loads((ROOT / "rois/rois.json").read_text())
    assert r["MPFC"]["centroid_mni"] == [-1.0, 49.0, -5.0] and 500 < r["MPFC"]["n_vox"] < 560
    assert r["AINS"]["centroid_mni"][1] > 0
    assert abs(r["AMY_L"]["centroid_mni"][0] + 23) < 3 and abs(r["AMY_R"]["centroid_mni"][0] - 23) < 3
