"""The AbRank and POISK verdicts: Kd <= 1 uM, a measured IC50 <= 10 ug/mL, an audit-dropped row, a patent-asserted pair."""

import pytest

from structcooker.instructions.readers.abag import columns_swapped
from structcooker.instructions.readers.abag_binding import abrank_verdict, poisk_verdict


@pytest.mark.parametrize(
    ("kd", "ic50", "op", "expected"),
    [
        ("3.4", "", "=", True),
        ("5000", "", "=", False),
        ("1000", "", "=", True),
        ("", "0.3", "=", True),
        ("", "30", "=", None),
        ("", "0.3", ">", None),
        ("", "", "=", None),
        ("nan", "", "=", None),
        ("", "nan", "=", None),
        ("3.4", "30", "=", True),
    ],
)
def test_abrank(kd, ic50, op, expected):
    assert abrank_verdict(kd, ic50, op) is expected


@pytest.mark.parametrize(
    ("text", "audit", "expected"),
    [
        ("", "CONFIRMED", (True, True)),  # a patent reports the pair, no number
        ("", "PLAUSIBLE", (True, True)),
        ("", "CONTRADICTED", (False, None)),
        ("Kd 3 nM", "NOT_FOUND", (False, None)),
        ("Kd = 1.5 x 10^-9 M", "CONFIRMED", (True, True)),
        ("KD: 5 nM", "CONFIRMED", (True, True)),
        ("Kd = 2 uM", "CONFIRMED", (True, False)),
        ("Kd > 2 uM", "CONFIRMED", (True, False)),
        (
            "Kd > 5 nM",
            "CONFIRMED",
            (True, None),
        ),  # a lower bound below 1 uM decides nothing
        ("< 1 nM (Kd)", "CONFIRMED", (True, True)),
        ("IC50: 0.06 ug/ml", "CONFIRMED", (True, True)),
        ("IC50 = 50 ug/ml", "CONFIRMED", (True, None)),
        ("EC50 = 0.3 nM", "CONFIRMED", (True, True)),
        ("5 nM", "CONFIRMED", (True, None)),  # unlabelled
        ("Ki 2 nM", "CONFIRMED", (True, None)),
        ("8.77x10^9 M-1", "CONFIRMED", (True, None)),  # an association constant
        ("10^-7 to 10^-8 M", "CONFIRMED", (True, None)),
        ("Titer 1:64", "CONFIRMED", (True, None)),
        ("Not Detectable", "CONFIRMED", (True, False)),
    ],
)
def test_poisk(text, audit, expected):
    assert poisk_verdict(text, audit) == expected


def test_columns_swapped_needs_both_fr4s():
    heavy_end, light_end = "XXXXXXXXWGQGTLVTVSS", "XXXXXXXXFGQGTKVEIK"
    assert columns_swapped("EVQL" + light_end, "QVQL" + heavy_end)
    assert not columns_swapped("EVQL" + heavy_end, "QVQL" + light_end)
    assert not columns_swapped(
        "EVQL" + light_end, "QVQL" + light_end
    )  # one chain with a light FR4 is not a swap
