"""One-micromolar binding verdicts of an ASD row."""

import pytest

from structcooker.instructions.readers.asd_binding import (
    DELTA_G_MAX,
    binding_verdict,
    pair_binds,
)


@pytest.mark.parametrize(
    ("kind", "value", "expected"),
    [
        ("bool", "1", True),
        ("bool", "1.0", True),
        ("bool", "True", True),
        ("bool", "0", False),
        ("bool", "0.0", False),
        ("bool", "False", False),
        ("bool", "nan", None),
        ("alphaseq", "2.07", True),
        ("alphaseq", "3.0", True),
        ("alphaseq", "3.01", False),
        ("alphaseq", "nan", None),
        ("-log Kd (M)", "6.0", True),
        ("-log Kd (M)", "5.99", False),
        ("Kd (nM)", "1000", True),
        ("Kd (nM)", "1830", False),
        ("Kd (M)", "6.35e-08", True),
        ("Kd (M)", "1.39e-05", False),
        ("Kd (M)", "I.C.", None),
        ("delta_g", "-11.0", True),
        ("delta_g", "-8.0", False),
        ("-log Kd (nM)", "8.2", None),
        ("EC50 (nM)", "5", True),
        ("EC50 (nM)", "303", True),
        ("EC50 (nM)", "5000", None),
        ("ddg", "-1", None),
        ("fuzzy", "h", True),
        ("fuzzy", "m", False),
        ("fuzzy", "l", False),
        ("fuzzy", "x", None),
    ],
)
def test_verdict(kind, value, expected):
    assert binding_verdict(kind, value) is expected


def test_delta_g_cut_is_one_micromolar_at_25_celsius():
    assert pytest.approx(-8.185, abs=1e-3) == DELTA_G_MAX


def test_pair_binds_needs_a_yes_and_no_no():
    assert pair_binds([True, None])
    assert not pair_binds([None])
    assert not pair_binds([True, False])
    assert not pair_binds([])


def test_unit_uncertain_dataset_gets_no_verdict():
    assert binding_verdict("-log Kd (M)", "7.0") is True
    assert binding_verdict("-log Kd (M)", "7.0", "flab_warszawski2019") is None
    assert binding_verdict("-log Kd (M)", "2.0", "flab_warszawski2019") is None
