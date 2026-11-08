"""Unit tests for the phenotype audit parsers and joins (synthetic data)."""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import audit_phenotype as ap  # noqa: E402


@pytest.mark.parametrize("raw,expected", [
    ("27 days before", 27.0), ("1 day before", 1.0),
    ("5 days after", -5.0), ("n/a", np.nan)])
def test_parse_days(raw, expected):
    got = ap.parse_days_before(raw)
    assert (math.isnan(got) and math.isnan(expected)) or got == expected


def test_parse_days_rejects_garbage():
    with pytest.raises(ValueError):
        ap.parse_days_before("about a month")


@pytest.mark.parametrize("raw,hours,flag", [
    ("9.58", 9 + 58 / 60, "ok"), ("10.01", 10 + 1 / 60, "ok"),
    ("10.1", 10 + 10 / 60, "single_digit_minutes"), ("10", 10.0, "ok"),
    ("11.90", np.nan, "implausible"), ("0", np.nan, "implausible"),
    ("1.05", np.nan, "implausible"), ("n/a", np.nan, "missing")])
def test_parse_clock(raw, hours, flag):
    h, f = ap.parse_clock_hmm(raw)
    assert f == flag
    assert (math.isnan(h) and math.isnan(hours)) or h == pytest.approx(hours)


@pytest.mark.parametrize("a,b,apoe", [
    ("Homozygous T/T", "Homozygous C/C", "e3e3"),
    ("Homozygous C/C", "Homozygous C/C", "e4e4"),
    ("Heterozygous C/T", "Homozygous C/C", "e3e4"),
    ("Homozygous T/T", "Heterozygous C/T", "e2e3"),
    ("Heterozygous C/T", "Heterozygous C/T", "e2e4")])
def test_apoe(a, b, apoe):
    assert ap.apoe_from_snps(a, b) == apoe


def test_na_and_crlf_roundtrip(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_bytes(b"participant_id\tv\r\nsub-00001\tn/a\r\nsub-00002\t3.5\r\n")
    df = ap.read_tsv_as_strings(p)
    assert list(df.columns) == ["participant_id", "v"]
    v = ap.to_float(df["v"])
    assert v.isna().tolist() == [True, False] and v.iloc[1] == 3.5


def test_validate_ids_detects_mismatch():
    parts = pd.DataFrame({"participant_id": ["sub-00001", "sub-00002", "sub-00002"]})
    pheno = pd.DataFrame({"participant_id": ["sub-00001", "sub-00003"]})
    r = ap.validate_ids(parts, pheno, {"sub-00001"})
    assert r["dup_participants_tsv"] == 1
    assert r["in_participants_not_phenotype"] == ["sub-00002"]
    assert r["in_phenotype_not_participants"] == ["sub-00003"]
    assert r["in_participants_no_dir"] == ["sub-00002"]


def test_post_treatment_draws_excluded_from_flow():
    df = pd.DataFrame({
        "treatment_group": ["OT", "PL", "PL"], "TICS_SCORE": ["n/a"] * 3,
        "avp_pg_ml": [1.0, 1.0, 1.0], "ot_pg_ml": [1.0, 1.0, 1.0],
        "days_blood_to_scan": [10.0, -5.0, np.nan],
        "bold_present": [True] * 3, "t1_present": [True] * 3,
        "age": [20.0] * 3, "sex": ["F"] * 3, "cohort": ["younger"] * 3,
        "tr_s": [2.0] * 3})
    flow = ap.eligibility_flow(df).set_index("step")
    assert flow.loc["blood draw known to PRECEDE dosing (days > 0)",
                    "n_remaining"] == 1
