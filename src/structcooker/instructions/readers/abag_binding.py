"""Does a row of AbRank or POISK say the antibody binds? One verdict per row, the same rule as ASD (``asd_binding``): Kd <= 1 uM binds.

AbRank carries its unit in the column header: ``Affinity_Kd [nM]`` is a Kd, ``IC50 [ug/mL]`` a neutralisation concentration, ``escape`` /
``log(Kd_ratio)`` the loss of binding of an RBD mutant (no verdict on their own). POISK carries free text in ``affinity_value``.
"""

from __future__ import annotations

import math
import re

from .asd_binding import KD_MAX_M

IC50_MAX_UG_ML = 10.0  # neutralising needs binding and an IC50 is about the Kd or above it; 1 uM is >= 15 ug/mL for any format, so 10 is safe for all


def _number(value: object) -> float | None:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def abrank_verdict(
    kd_nm: object, ic50_ug_ml: object, operator: str = "="
) -> bool | None:
    """Return the verdict of an AbRank row.

    ``Affinity_Kd [nM]`` <= 1000 binds, above does not; without a Kd, a measured ``IC50 [ug/mL]`` (operator ``=``) <= 10 binds (above: no verdict);
    an IC50 beyond the tested range (``>`` / ``<``), ``escape`` and ``log(Kd_ratio)`` alone give no verdict.
    """
    kd = _number(kd_nm)
    if kd is not None:
        return kd <= KD_MAX_M * 1e9
    ic50 = _number(ic50_ug_ml)
    if ic50 is not None and operator in ("=", ""):
        return True if ic50 <= IC50_MAX_UG_ML else None
    return None


_UNITS = {
    "pm": 1e-12,
    "nm": 1e-9,
    "um": 1e-6,
    "µm": 1e-6,
    "μm": 1e-6,
    "mm": 1e-3,
    "m": 1.0,
}
_NUM = r"(\d+(?:[.,]\d+)?(?:\s*[x\u00d7*]\s*10\s*\^?\s*[-\u2212\u2013]?\s*\d+)?)"
_MOLAR = re.compile(
    r"([<>≤≥~=]*)\s*" + _NUM + r"\s*(pM|nM|[uµμ]M|mM|M)(?![\w^-])", re.IGNORECASE
)
_MASS = re.compile(r"([<>≤≥~=]*)\s*" + _NUM + r"\s*[uµμ]g\s*/\s*m[lL]", re.IGNORECASE)
_NEGATIVE = re.compile(
    r"not\s+detect|no\s+(?:detectable\s+)?bind|non[- ]?bind|\bnegative\b|no\s+affinity",
    re.IGNORECASE,
)
_KD = re.compile(r"\bK\s*_?\s*D\b|\bKd\b|dissociation", re.IGNORECASE)
_IC50 = re.compile(r"\bIC\s*50\b", re.IGNORECASE)
_EC50 = re.compile(r"\bEC\s*50\b", re.IGNORECASE)


def _value(text: str) -> float | None:
    number = (
        text.replace(",", ".")
        .replace("\u2212", "-")
        .replace("\u2013", "-")
        .replace(" ", "")
    )
    power = re.match(r"([\d.]+)[x\u00d7*]10\^?(-?\d+)", number)
    try:
        return (
            float(power.group(1)) * 10 ** int(power.group(2))
            if power
            else float(number)
        )
    except ValueError:
        return None


def poisk_verdict(text: str, audit_verdict: str) -> tuple[bool, bool | None]:
    """``(keep_row, verdict)`` of a POISK row.

    A row POISK's own audit calls ``CONTRADICTED`` (the patent text says otherwise) or ``NOT_FOUND`` is dropped. A row **without any affinity text** is a
    pair the patent reports as binding: a binder, like the ASD ``patents`` rows. With text: a negative phrase (``Not Detectable``) does not bind; a value
    **labelled Kd** in a molar unit follows the Kd rule (``<`` and ``>`` count only when they decide it); a measured ``IC50`` / ``EC50`` (molar: <= 1 uM,
    ug/mL: <= 10) binds, above gives no verdict; anything else -- an unlabelled concentration, ``Ki``, an association constant ``M-1``, a range, a ratio to
    another antibody, a titer -- gives no verdict.
    """
    if audit_verdict in ("CONTRADICTED", "NOT_FOUND"):
        return False, None
    text = (text or "").strip()
    if not text:
        return True, True
    if _NEGATIVE.search(text):
        return True, False
    molar, mass = _MOLAR.findall(text), _MASS.findall(text)
    if len(molar) + len(mass) != 1:
        return True, None
    kd, ic50, ec50 = (
        bool(_KD.search(text)),
        bool(_IC50.search(text)),
        bool(_EC50.search(text)),
    )
    if mass:
        operator, number = mass[0]
        value = _value(number)
        if value is None or not (ic50 or ec50) or kd or (ic50 and ec50):
            return True, None
        return True, (
            True
            if value <= IC50_MAX_UG_ML and ">" not in operator and "≥" not in operator
            else None
        )
    operator, number, unit = molar[0]
    value = _value(number)
    if value is None:
        return True, None
    molar_value = value * _UNITS[unit.lower()]
    upper, lower = operator in ("<", "<=", "≤"), operator in (">", ">=", "≥")
    if kd and not (ic50 or ec50):
        if upper:
            return True, (True if molar_value <= KD_MAX_M else None)
        if lower:
            return True, (False if molar_value >= KD_MAX_M else None)
        return True, molar_value <= KD_MAX_M
    if (ic50 or ec50) and not kd and not (ic50 and ec50):
        return True, (True if molar_value <= KD_MAX_M and not lower else None)
    return True, None
