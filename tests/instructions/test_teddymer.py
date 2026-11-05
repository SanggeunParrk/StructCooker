"""Teddymer PDB reader and bond derivation."""

import numpy as np
import pytest

from structcooker.instructions.readers import teddymer as reader
from structcooker.instructions.transforms import teddymer as bonds


def _atom(
    serial,
    name,
    res,
    chain,
    resid,
    xyz=(1.0, 2.0, 3.0),
    bf=90.0,
    el="",
    rec="ATOM  ",
    alt=" ",
    ins=" ",
):
    return (
        f"{rec}{serial:5d} {name:<4s}{alt}{res:>3s} {chain}{resid:4d}{ins}   "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{bf:6.2f}          {el:>2s}\n"
    )


PDB = [
    _atom(1, "N", "GLY", "A", 1, el="N", bf=91.5),
    _atom(2, "CA", "GLY", "A", 1, el="C"),
    _atom(3, "C", "GLY", "A", 1, el="C"),
    _atom(4, "N", "GLY", "A", 2, el="N"),
    _atom(5, "CA", "GLY", "A", 2, el="C"),
    _atom(6, "N", "GLY", "A", 5, el="N"),  # residue-number jump (domain boundary)
    _atom(7, "N", "ALA", "B", 1, el="N"),
    "TER\n",
    "END\n",
]


def test_entry_key_is_filename_stem():
    assert reader.teddymer_entry_key("/x/12DI_A0A3N5GBR8_TED02_TED03.pdb") == (
        "12DI_A0A3N5GBR8_TED02_TED03"
    )


def test_pdb_atom_table_parses_columns():
    t = reader.pdb_atom_table(PDB, name="x.pdb")
    assert t["atom_name"].tolist() == ["N", "CA", "C", "N", "CA", "N", "N"]
    assert t["res_name"][-1] == "ALA"
    assert t["chain_id"].tolist() == list("AAAAAAB")
    assert t["res_id"].tolist() == [1, 1, 1, 2, 2, 5, 1]
    assert t["b_factor"][0] == pytest.approx(91.5)
    assert t["coord"][0].tolist() == [1.0, 2.0, 3.0]
    assert t["molecule_type_id"].tolist() == [0] * 7
    assert not t["hetero"].any()


def test_pdb_atom_table_skips_non_first_altloc_and_optionally_hydrogens():
    lines = [
        _atom(1, "N", "GLY", "A", 1, el="N"),
        _atom(2, "N", "GLY", "A", 1, el="N", alt="B"),
        _atom(3, "H", "GLY", "A", 1, el="H"),
        _atom(4, "O", "HOH", "A", 9, el="O", rec="HETATM"),
    ]
    assert len(reader.pdb_atom_table(lines, name="x")["atom_name"]) == 3
    t = reader.pdb_atom_table(lines, name="x", drop_hydrogens=True)
    assert t["element"].tolist() == ["N", "O"]
    assert t["hetero"].tolist() == [False, True]


def test_pdb_atom_table_element_falls_back_to_atom_name():
    line = _atom(1, "CA", "GLY", "A", 1, el="")
    assert reader.pdb_atom_table([line], name="x")["element"].tolist() == ["C"]


def test_pdb_atom_table_without_atoms_raises():
    with pytest.raises(ValueError, match="empty.pdb has no ATOM records"):
        reader.pdb_atom_table(["TER\n", "END\n"], name="empty.pdb")


def test_entity_ids_group_identical_sequences_across_chains():
    lines = [
        _atom(1, "CA", "GLY", "A", 1, el="C"),
        _atom(2, "CA", "ALA", "A", 2, el="C"),
        _atom(3, "CA", "GLY", "B", 1, el="C"),
        _atom(4, "CA", "ALA", "B", 2, el="C"),
        _atom(5, "CA", "SER", "C", 1, el="C"),
    ]
    t = reader.pdb_atom_table(lines, name="x")
    assert t["entity_id"].tolist() == ["1", "1", "1", "1", "2"]


def test_get_teddymer_structure_data_reads_file(tmp_path):
    p = tmp_path / "7DI_ACC_TED01_TED02.pdb"
    p.write_text("".join(PDB))
    out = reader.get_teddymer_structure_data(p)
    assert out["entry_id"] == "7DI_ACC_TED01_TED02"
    assert len(out["raw_atom_site_dict"]["atom_name"]) == 7


def _gly_ccd():
    return {
        "atom": {
            "nodes": {"id": {"value": np.array(["N", "CA", "C"])}},
            "edges": {
                "bond_type": {
                    "value": np.array(["SING", "SING"]),
                    "src_indices": np.array([0, 1]),
                    "dst_indices": np.array([1, 2]),
                }
            },
        }
    }


def test_attach_bonds_links_intra_residue_and_consecutive_peptide_only():
    raw = reader.pdb_atom_table(PDB, name="x")
    out = bonds.attach_bonds(raw, {"GLY": _gly_ccd(), "ALA": None})
    pairs = sorted(map(tuple, out["bonds"].tolist()))
    # residue 1: N-CA, CA-C; residue 2: N-CA (C absent); peptide C1-N2 (res 1 -> 2).
    # residue 2 -> 5 is a numbering jump, so no peptide bond; chain B is separate.
    assert pairs == [(0, 1), (1, 2), (2, 3), (3, 4)]
    assert "bonds" not in raw  # input table is not mutated


def test_attach_bonds_without_any_bond_gives_empty_2col_array():
    raw = reader.pdb_atom_table([_atom(1, "N", "ALA", "A", 1, el="N")], name="x")
    out = bonds.attach_bonds(raw, {"ALA": None})
    assert out["bonds"].shape == (0, 2)
    assert out["bonds"].dtype == np.int64
