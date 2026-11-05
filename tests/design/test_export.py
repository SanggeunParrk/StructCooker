import numpy as np
import pytest
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.instructions.design.export import (
    CHAIN_LETTERS,
    backbone_pdb_text,
    n_protein_chains,
    write_backbone_pdb,
)


def _mol():
    lig = chain(
        "L", [residue("XYZ", "X", [("C1", "C", (0, 0, 9))])], etype="non-polymer"
    )
    return make_mol([protein_chain("A", 3), lig, protein_chain("B", 2, z=4.0)])


def test_only_polypeptide_chains_with_complete_backbone_are_counted():
    assert n_protein_chains(_mol()) == 2


def test_residue_with_missing_oxygen_is_left_out():
    a = protein_chain("A", 3)
    a["residues"][1]["atoms"] = [x for x in a["residues"][1]["atoms"] if x[0] != "O"]
    text = backbone_pdb_text(make_mol([a]))
    assert len(text.splitlines()) == 2 * 4 + 1
    assert "   2    " not in text  # residue number 2 (cif_idx 1 + 1) is absent


def test_pdb_text_has_n_ca_c_o_per_residue_chain_letters_and_end():
    lines = backbone_pdb_text(_mol()).splitlines()
    assert lines[-1] == "END"
    atoms = lines[:-1]
    assert len(atoms) == 5 * 4
    assert [ln[12:16].strip() for ln in atoms[:4]] == ["N", "CA", "C", "O"]
    assert {ln[21] for ln in atoms[:12]} == {"A"}
    assert {ln[21] for ln in atoms[12:]} == {"B"}
    first = atoms[0]
    assert first[17:20] == "ALA"
    assert float(first[30:38]) == pytest.approx(-1.3)
    assert float(first[38:46]) == pytest.approx(0.5)
    assert int(first[22:26]) == 1  # cif_idx + 1
    assert int(atoms[-1][6:11]) == 20  # serial numbers run on across chains
    assert first[76:78] == " N"


def test_more_than_62_chains_cannot_be_exported():
    chains = [
        protein_chain(f"c{i}", 1, z=10.0 * i) for i in range(len(CHAIN_LETTERS) + 1)
    ]
    mol = make_mol(chains)
    assert n_protein_chains(mol) == 63
    with pytest.raises(ValueError, match="63 protein chains"):
        backbone_pdb_text(mol)


def test_write_backbone_pdb_returns_the_exported_residue_indices(tmp_path):
    out = tmp_path / "bb.pdb"
    rows = write_backbone_pdb(_mol(), out)
    assert out.read_text() == backbone_pdb_text(_mol())
    assert isinstance(rows, np.ndarray)
    # residues: A0 A1 A2, ligand residue (index 3), B0 B1 (4, 5)
    assert rows.tolist() == [0, 1, 2, 4, 5]
