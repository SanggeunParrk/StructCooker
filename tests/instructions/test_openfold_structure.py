"""Openfold structure steps: hierarchy recovery, CCD-derived features, assembly."""

import numpy as np

from structcooker.instructions.mols import CIFMol
from structcooker.instructions.transforms import openfold_structure as ofs


def _gly_entry(alt: list[str] | None = None) -> dict:
    """Build a tiny CCD-like GLY component: N-CA (SING), CA-C (SING), C=O (DOUB)."""
    names = np.array(["N", "CA", "C", "O"])
    nodes = {
        "id": {"value": names},
        "aromatic": {"value": np.array(["N", "N", "N", "N"])},
        "stereo": {"value": np.array(["N", "N", "N", "N"])},
        "charge": {"value": np.array(["0", "0", "0", "-1"])},
        "model_xyz": {"value": np.array([[f"{i}.0", "0.0", "0.0"] for i in range(4)])},
    }
    if alt is not None:
        nodes["alt_atom_id"] = {"value": np.array(alt)}
    src, dst = np.array([0, 1, 2]), np.array([1, 2, 3])
    return {
        "atom": {
            "nodes": nodes,
            "edges": {
                "bond_type": {
                    "value": np.array(["SING", "SING", "DOUB"]),
                    "src_indices": src,
                    "dst_indices": dst,
                },
                "bond_aromatic": {"value": np.array(["N", "N", "N"])},
                "bond_stereo": {"value": np.array(["N", "N", "N"])},
            },
        },
        "residue": {"nodes": {"formula": {"value": np.array(["C2 H5 N O2"])}}},
    }


def _atom_table() -> dict:
    """Chain A: GLY1, GLY2 (protein); chain B: one ligand-ish residue ZZZ (mol type 3)."""
    rows = [
        ("N", "GLY", "A", 1, "1", 0),
        ("CA", "GLY", "A", 1, "1", 0),
        ("C", "GLY", "A", 1, "1", 0),
        ("N", "GLY", "A", 2, "1", 0),
        ("CA", "GLY", "A", 2, "1", 0),
        ("ZZ", "ZZZ", "B", 1, "2", 3),
    ]
    n = len(rows)
    return {
        "atom_name": np.array([r[0] for r in rows]),
        "res_name": np.array([r[1] for r in rows]),
        "chain_id": np.array([r[2] for r in rows]),
        "res_id": np.array([r[3] for r in rows]),
        "ins_code": np.array([""] * n),
        "entity_id": np.array([r[4] for r in rows]),
        "molecule_type_id": np.array([r[5] for r in rows]),
        "hetero": np.array([False] * 5 + [True]),
        "element": np.array(["N", "C", "C", "N", "C", "X"]),
        "coord": np.arange(n * 3, dtype=float).reshape(n, 3),
        "occupancy": np.ones(n),
        "bonds": np.array([[0, 1], [1, 2], [2, 3], [3, 4]]),
    }


def test_group_ids_increment_on_any_column_change():
    a = np.array(["A", "A", "A", "B"])
    b = np.array([1, 1, 2, 2])
    assert ofs._group_ids(a, b).tolist() == [0, 0, 1, 2]


def test_first_per_group_takes_first_member():
    ids = np.array([0, 0, 1, 1, 2])
    assert ofs._first_per_group(ids, np.array(list("abcde"))).tolist() == list("ace")


def test_build_hierarchy_recovers_residues_and_chains():
    (a2r, r2c, n_chain, rn, rid, het, cid, eid, mol) = ofs.build_hierarchy(
        _atom_table()
    )
    assert a2r.tolist() == [0, 0, 0, 1, 1, 2]
    assert r2c.tolist() == [0, 0, 1]
    assert n_chain == 2
    assert rn.tolist() == ["GLY", "GLY", "ZZZ"]
    assert rid.tolist() == [1, 2, 1]
    assert het.tolist() == [0, 0, 1]
    assert cid.tolist() == ["A", "B"]
    assert eid.tolist() == ["1", "2"]
    assert mol.tolist() == [0, 3]


def test_build_hierarchy_splits_same_chain_label_by_entity():
    t = _atom_table()
    t["chain_id"] = np.array(["A"] * 6)
    n_chain = ofs.build_hierarchy(t)[2]
    assert n_chain == 2  # same chain label, different entity ids -> two chains


def test_ccd_atom_lookup_none_is_empty_and_canonical_wins_over_alt():
    assert ofs._ccd_atom_lookup(None) == {}
    # atom N has alt "CA" which collides with the canonical id of atom CA
    lookup = ofs._ccd_atom_lookup(_gly_entry(alt=["CA", "", "C4A", ""]))
    assert lookup["CA"]["model_xyz"][0] == "1.0"  # CA's own features, not N's
    assert lookup["C4A"]["model_xyz"][0] == "2.0"  # alt alias resolves to atom C
    assert lookup["O"]["charge"] == "-1"


def test_read_ccd_entry_returns_none_when_component_missing(tmp_path):
    assert ofs._read_ccd_entry(tmp_path / "no_such_db", "GLY") is None


def test_load_ccd_entries_reads_each_unique_residue_name(monkeypatch, tmp_path):
    calls = []

    def fake(_path, name):
        calls.append(name)

    monkeypatch.setattr(ofs, "_read_ccd_entry", fake)
    out = ofs.load_ccd_entries(_atom_table(), tmp_path)
    assert sorted(calls) == ["GLY", "ZZZ"]
    assert out == {"GLY": None, "ZZZ": None}


def test_derive_atom_features_fills_known_atoms_and_defaults_unknown():
    t = _atom_table()
    feats = ofs.derive_atom_features(t, {"GLY": _gly_entry(), "ZZZ": None})
    assert feats["model_xyz"][1].tolist() == ["1.0", "0.0", "0.0"]
    assert feats["charge"].tolist() == ["0"] * 6  # O absent from table; ZZ default
    assert feats["aromatic"].tolist() == ["N"] * 6
    assert feats["model_xyz"][5].tolist() == ["", "", ""]


def test_derive_bond_edges_types_intra_bonds_and_collapses_residue_bonds():
    t = _atom_table()
    t["bonds"] = np.array([[0, 1], [1, 2], [2, 3], [2, 4]])
    e = ofs.derive_bond_edges(t, {"GLY": _gly_entry(), "ZZZ": None})
    assert e["bond_type"].tolist() == ["SING", "SING", "SING", "SING"]
    # C(res0)-N(res1) and C(res0)-CA(res1) collapse into one residue-level bond
    assert list(zip(e["res_src"].tolist(), e["res_dst"].tolist(), strict=True)) == [
        (0, 1)
    ]
    assert e["res_bond_value"].tolist() == ["SING"]


def test_derive_bond_edges_uses_ccd_bond_type_and_alt_names():
    t = _atom_table()
    t["atom_name"] = np.array(["N", "CA", "Cx", "N", "CA", "ZZ"])
    t["bonds"] = np.array([[1, 2]])
    # atom C has alt name "Cx": the alt endpoint must still resolve the CCD bond
    e = ofs.derive_bond_edges(t, {"GLY": _gly_entry(alt=["", "", "Cx", ""])})
    assert e["bond_type"].tolist() == ["SING"]
    t["atom_name"][:3] = ["C", "O", "x"]
    t["bonds"] = np.array([[0, 1]])
    e = ofs.derive_bond_edges(t, {"GLY": _gly_entry()})
    assert e["bond_type"].tolist() == ["DOUB"]


def test_derive_residue_features_maps_one_letter_formula_and_unknown():
    f = ofs.derive_residue_features(
        np.array(["GLY", "A", "ZZZ"]),
        np.array([1, 2, 3]),
        np.array([0, 0, 1]),
        {"GLY": _gly_entry(), "A": None, "ZZZ": None},
    )
    assert f["one_letter_code"].tolist() == ["G", "A", "X"]
    assert f["formula"].tolist() == ["C2 H5 N O2", "", ""]
    assert f["cif_idx"].tolist() == ["1", "2", "3"]
    assert f["hetero"].tolist() == [0, 0, 1]


def test_derive_chain_features_maps_molecule_type_to_entity_type():
    f = ofs.derive_chain_features(
        np.array(["A", "B", "C", "D", "E"]),
        np.array(["1", "2", "3", "4", "5"]),
        np.array([0, 1, 2, 3, 9]),
    )
    assert f["entity_type"].tolist() == [
        "polypeptide(L)",
        "polyribonucleotide",
        "polydeoxyribonucleotide",
        "non-polymer",
        "non-polymer",
    ]
    assert f["auth_asym_id"].tolist() == list("ABCDE")


def _build(table, cache):
    (a2r, r2c, n_chain, rn, rid, het, cid, eid, mol) = ofs.build_hierarchy(table)
    return ofs.assemble_cifmol(
        table,
        a2r,
        r2c,
        n_chain,
        ofs.derive_atom_features(table, cache),
        ofs.derive_bond_edges(table, cache),
        ofs.derive_residue_features(rn, rid, het, cache),
        ofs.derive_chain_features(cid, eid, mol),
    )


def test_assemble_cifmol_round_trips_and_defaults_b_factor_to_zero():
    mol = CIFMol.from_dict(_build(_atom_table(), {"GLY": _gly_entry(), "ZZZ": None}))
    assert len(mol.atoms.id.value) == 6
    assert len(mol.chains.chain_id.value) == 2
    assert mol.atoms.b_factor.value.tolist() == [0.0] * 6


def test_assemble_cifmol_keeps_b_factor_when_present():
    t = _atom_table()
    t["b_factor"] = np.linspace(10, 60, 6)
    mol = CIFMol.from_dict(_build(t, {"GLY": _gly_entry(), "ZZZ": None}))
    assert mol.atoms.b_factor.value.tolist() == np.linspace(10, 60, 6).tolist()


def test_build_one_cifmol_uses_ccd_db(monkeypatch, tmp_path):
    monkeypatch.setattr(
        ofs, "_read_ccd_entry", lambda _p, n: _gly_entry() if n == "GLY" else None
    )
    out = ofs._build_one_cifmol(_atom_table(), tmp_path)
    mol = CIFMol.from_dict(out)
    assert mol.residues.one_letter_code_can.value.tolist() == ["G", "G", "X"]


def test_build_disordered_template_mols_skips_malformed_and_reduces_good(
    monkeypatch, tmp_path
):
    seen = {}

    def fake_reduce(cifmol, alignment):
        seen["seq"] = alignment
        return "TEMPLATE"

    import structcooker.instructions.transforms.template as template_mod

    monkeypatch.setattr(template_mod, "to_template_mol", fake_reduce)
    monkeypatch.setattr(
        ofs, "_read_ccd_entry", lambda _p, n: _gly_entry() if n == "GLY" else None
    )
    out = ofs.build_disordered_template_mols(
        {"good_A": _atom_table(), "bad_A": {"nothing": np.array([])}}, tmp_path
    )
    assert out == {"good_A": "TEMPLATE"}
    assert seen["seq"] == ("GGX", "GGX")  # identity alignment of its own sequence


def test_wrap_cifmol_layout():
    assembly, meta = ofs.wrap_cifmol({"x": 1}, "MGYP9")
    assert assembly == {"1_1_.": {"x": 1}}
    assert meta == {"id": ["MGYP9"], "deposition_date": None, "resolution": None}
