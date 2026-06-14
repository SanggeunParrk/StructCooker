"""Unit test: the CCD ingest recipe parses a real component correctly."""

from __future__ import annotations

from pathlib import Path

from datacooker import describe, parse_file

from structcooker.instructions.readers.cif import get_cif_data
from structcooker.workflows.ingest.ccd import RECIPE, TARGETS


def test_recipe_shape() -> None:
    # The recipe needs no external inputs and ends at chem_comp_dict.
    assert TARGETS == ["chem_comp_dict"]
    assert RECIPE.required_inputs(TARGETS) == set()
    assert "chem_comp_dict <- parse_chem_comp(" in describe(
        RECIPE, targets=TARGETS, detail="compact",
    )


def test_parses_alanine(fixture_dir: Path) -> None:
    entry = parse_file(
        RECIPE, fixture_dir / "ala.cif", get_cif_data, targets=TARGETS,
    )["chem_comp_dict"]

    atom = entry["atom"]
    # ALA, hydrogens removed: N, CA, C, O, CB, OXT.
    assert list(atom["element"].value) == ["N", "C", "C", "O", "C", "O"]
    # 5 heavy-atom bonds, none dangling.
    bond = atom["bond_type"]
    assert len(bond.src_indices) == 5
    n = len(atom["id"].value)
    assert all(0 <= i < n for i in bond.src_indices)
    assert all(0 <= i < n for i in bond.dst_indices)
