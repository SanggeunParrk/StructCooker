"""Unit tests for the CCD per-entry validators (synthetic input, no I/O)."""

from __future__ import annotations

import numpy as np
import pytest
from biomol.core.container import FeatureContainer
from biomol.core.feature import EdgeFeature, NodeFeature

from structcooker.instructions.transforms.ccd import (
    validate_atom_container,
    validate_chem_comp,
)


def make_atom(ids, elements, *, xyz=None, bonds=None) -> FeatureContainer:
    """Build a minimal atom FeatureContainer for testing."""
    feats = {
        "id": NodeFeature(value=np.array(ids)),
        "element": NodeFeature(value=np.array(elements)),
    }
    if xyz is not None:
        feats["model_xyz"] = NodeFeature(value=np.array(xyz))
    if bonds is not None:
        src = np.array([b[0] for b in bonds], dtype=int)
        dst = np.array([b[1] for b in bonds], dtype=int)
        feats["bond_type"] = EdgeFeature(
            value=np.array(["SING"] * len(bonds)),
            src_indices=src,
            dst_indices=dst,
        )
    return FeatureContainer(features=feats)


def codes(atom: FeatureContainer) -> list[str]:
    return [issue["code"] for issue in validate_atom_container(atom)]


def test_clean_entry_has_no_issues() -> None:
    atom = make_atom(
        ["C1", "C2", "O1"], ["C", "C", "O"],
        xyz=[["0", "0", "0"], ["1.5", "0", "0"], ["2.7", "0", "0"]],
        bonds=[(0, 1), (1, 2)],
    )
    assert validate_atom_container(atom) == []


def test_self_bond_is_flagged() -> None:
    atom = make_atom(["C1", "C2"], ["C", "C"],
                     xyz=[["0", "0", "0"], ["1.5", "0", "0"]], bonds=[(0, 0)])
    assert "self_bond" in codes(atom)


def test_duplicate_bond_is_flagged() -> None:
    atom = make_atom(["C1", "C2"], ["C", "C"],
                     xyz=[["0", "0", "0"], ["1.5", "0", "0"]], bonds=[(0, 1), (1, 0)])
    assert "duplicate_bond" in codes(atom)


def test_bond_too_long_is_flagged() -> None:
    atom = make_atom(["C1", "C2"], ["C", "C"],
                     xyz=[["0", "0", "0"], ["5.0", "0", "0"]], bonds=[(0, 1)])
    assert "bond_too_long" in codes(atom)


def test_atom_overlap_is_flagged() -> None:
    atom = make_atom(["C1", "C2"], ["C", "C"],
                     xyz=[["0", "0", "0"], ["0.1", "0", "0"]], bonds=[(0, 1)])
    assert "atom_overlap" in codes(atom)


def test_missing_xyz_is_flagged() -> None:
    atom = make_atom(["C1", "C2"], ["C", "C"],
                     xyz=[["0", "0", "0"], ["?", "?", "?"]], bonds=[(0, 1)])
    assert "missing_xyz" in codes(atom)


def test_unparsed_entry_reports_issue() -> None:
    instruction = validate_chem_comp()
    issues = instruction(None)
    assert [i["code"] for i in issues] == ["unparsed"]


@pytest.mark.parametrize("element", ["XX", "Uuo"])
def test_unknown_element_skips_geometry(element: str) -> None:
    # Unknown elements have no covalent radius, so geometry can't be judged:
    # a far bond must NOT raise a geometry issue (only known structural ones).
    atom = make_atom(["A1", "A2"], [element, element],
                     xyz=[["0", "0", "0"], ["9", "0", "0"]], bonds=[(0, 1)])
    assert "bond_too_long" not in codes(atom)
