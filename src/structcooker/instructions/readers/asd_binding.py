"""Does an ASD row say the antibody binds? One verdict per row, one rule per ``affinity_type``.

The release mixes scales that cannot be compared (``bool``, a log10 assay score, -log Kd, Kd in several units, a free
energy, ratios), so each type is read on its own and compared with one affinity: Kd <= 1 uM binds. Types that are not
an affinity (EC50, ka, ddG, enrichment, a ratio, IC50, ``fuzzy``) and types whose unit is not known give no verdict.
"""
from __future__ import annotations

import math

KD_MAX_M = 1e-6  # the cut: Kd of 1 uM

# delta G = RT ln(Kd) at 25 C (298.15 K), kcal/mol: the unit and the temperature are assumed; the release gives neither.
GAS_CONSTANT_KCAL = 1.987204e-3
TEMPERATURE_K = 298.15
DELTA_G_MAX = GAS_CONSTANT_KCAL * TEMPERATURE_K * math.log(KD_MAX_M)  # -8.185 kcal/mol

ALPHASEQ_MAX = math.log10(KD_MAX_M * 1e9)  # the score is log10(Kd in nM): 3.0

_BOOL = {"1": True, "1.0": True, "true": True, "0": False, "0.0": False, "false": False}
# ``fuzzy`` is the Buzz trastuzumab-HER2 library (Chinery et al.): FACS sorted into high / medium / low affinity bins, the high
# bin the binders, medium and low the non-binders ("we assign 'high' ... positive binders", "group 'medium' and 'low' ... negative
# binders"). ASD writes the bins h / m / l; that reading is inferred (172,149 h of 524,346 = 32.8%, the paper's share of high-affinity
# binders), the paper never uses the letters. No Kd behind it: "binder" here is the FACS high bin, not Kd <= 1 uM.
_FUZZY = {"h": True, "m": False, "l": False}


def _number(affinity: object) -> float | None:
    try:
        value = float(affinity)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


# Datasets whose affinity_type label does not fit their values: no verdict whatever the value. flab_warszawski2019 is labelled
# ``-log Kd (M)`` but its values (-1.78 .. 4.43) are not pKd in molar (that would be millimolar non-binders) and the unit that
# would make them sensible (nM) is a guess; its 2,048 rows are left out rather than judged on a guess (docs/abag_distillation/affinity-filter.md).
UNIT_UNCERTAIN_DATASETS = frozenset({"flab_warszawski2019"})


def binding_verdict(affinity_type: str | None, affinity: object, dataset: str | None = None) -> bool | None:
    """``True`` binds (Kd <= 1 uM), ``False`` does not, ``None`` the row cannot say.

    ``bool``: the label. ``alphaseq``: log10(Kd nM) <= 3 (a missing score is unknown). ``-log Kd (M)``: >= 6.
    ``fuzzy`` (Buzz): ``h`` binds, ``m`` / ``l`` do not. ``Kd (nM)`` / ``Kd (M)``: <= 1000 nM / 1e-6 M. ``EC50 (nM)``: <= 1000 nM binds, above gives no verdict. ``delta_g``: <= -8.185 kcal/mol (25 C). ``-log Kd (nM)`` is
    unknown: the values (4.5-9.5) read as pKd in molar, not in nanomolar, and the release does not say which.
    """
    if dataset in UNIT_UNCERTAIN_DATASETS:
        return None
    kind = (affinity_type or "").strip()
    if kind == "bool":
        return _BOOL.get(str(affinity).strip().lower())
    if kind == "fuzzy":
        return _FUZZY.get(str(affinity).strip().lower())
    value = _number(affinity)
    if value is None:
        return None
    if kind == "alphaseq":
        return value <= ALPHASEQ_MAX
    if kind == "-log Kd (M)":
        return value >= -math.log10(KD_MAX_M)
    if kind == "Kd (nM)":
        return value <= KD_MAX_M * 1e9
    if kind in ("Kd (M)", "KD (M)"):
        return value <= KD_MAX_M
    if kind == "EC50 (nM)":
        # a half-effect concentration is about a Kd or above it, so <= 1 uM says Kd <= ~1 uM; above says nothing
        return True if value <= KD_MAX_M * 1e9 else None
    if kind == "delta_g":
        return value <= DELTA_G_MAX
    return None


def pair_binds(verdicts: list[bool | None]) -> bool:
    """A pair binds when some row says so and none says it does not."""
    return True in verdicts and False not in verdicts
