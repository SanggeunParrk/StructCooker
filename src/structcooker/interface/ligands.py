"""Ligand exclusion: which CCD components never count as a binding partner."""

from __future__ import annotations

from pathlib import Path

import yaml

DEFAULT_PATH = Path(__file__).resolve().parents[3] / "configs" / "interface" / "ligand_exclusion.yaml"


def load_excluded(
    categories: list[str],
    path: Path | str = DEFAULT_PATH,
) -> frozenset[str]:
    """Union of the CCD ids listed under ``categories`` (upper-cased)."""
    table = yaml.safe_load(Path(path).read_text())
    missing = [c for c in categories if c not in table]
    if missing:
        msg = f"unknown ligand-exclusion categories {missing}; have {sorted(table)}"
        raise KeyError(msg)
    return frozenset(str(c).upper() for cat in categories for c in table[cat])


def is_binding_ligand(ccd_id: str, n_heavy: int, excluded: frozenset[str], min_heavy: int) -> bool:
    """True if the component may act as a ligand party."""
    return ccd_id.upper() not in excluded and n_heavy >= min_heavy
