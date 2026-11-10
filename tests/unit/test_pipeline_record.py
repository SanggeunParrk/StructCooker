import numpy as np
import pytest
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.instructions.unit.config import load_config
from structcooker.instructions.unit.pipeline import extract_units
from structcooker.instructions.unit.record import (
    SCHEMA,
    asm_fields,
    manifest_row,
    unit_to_record,
)


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def _complex():
    lig = chain(
        "L",
        [
            residue(
                "XYZ",
                "X",
                [(f"C{i}", "C", (3.8 * i + 1, -2.5, 1.0)) for i in range(7)],
            )
        ],
        etype="non-polymer",
    )
    water = chain(
        "W",
        [residue("HOH", "X", [("O", "O", (10.0, -3.0, 2.5))])],
        etype="water",
    )
    far_water = chain(
        "V", [residue("HOH", "X", [("O", "O", (10.0, 40.0, 40.0))])], etype="water"
    )
    return make_mol(
        [
            protein_chain("A", 8),
            protein_chain("B", 8, z=4.0, sign=-1.0),
            lig,
            water,
            far_water,
        ]
    )


@pytest.fixture(scope="module")
def extracted(cfg):
    return extract_units(_complex(), cfg)


def test_extract_units_finds_the_protein_protein_unit_with_ligand(extracted):
    asm, units, duplex = extracted
    assert asm.kinds == ["protein", "protein", "ligand"]
    assert duplex == {}
    assert len(units) >= 1
    assert units[0].nodes[:2] == [0, 1]
    assert asm.stats["ab_class"] == [None, None, None]


def test_ab_lookup_is_called_with_the_key_sequence_for_proteins_only(cfg):
    seen = []

    def lookup(seq):
        seen.append(seq)
        return "H"

    asm, _, _ = extract_units(_complex(), cfg, ab_lookup=lookup)
    assert seen == ["AAAAAAAA", "AAAAAAAA"]
    assert asm.stats["ab_class"] == ["H", "H", None]


def test_empty_assembly_returns_no_units(cfg):
    asm, units, duplex = extract_units(make_mol([]), cfg)
    assert len(asm.kinds) == 0
    assert units == []
    assert duplex == {}


def test_na_chain_triggers_duplex_search_without_pairs_for_a_single_strand(cfg):
    dna = chain(
        "D",
        [
            residue(
                "DA", "A", [("C1'", "C", (3.0 * i, 0, 0)), ("N1", "N", (3.0 * i, 1, 0))]
            )
            for i in range(5)
        ],
        etype="polydeoxyribonucleotide",
    )
    off = load_config(overrides={"monomer.enabled": False})
    asm, units, duplex = extract_units(make_mol([dna]), off)
    assert asm.kinds == ["na"]
    assert units == []
    assert duplex == {}
    # with monomers on, the lone strand is a unit of its own
    _, units, duplex = extract_units(make_mol([dna]), cfg)
    assert [u.nodes for u in units] == [[0]]
    assert duplex == {}


def _record(extracted, cfg, idx=0):
    asm, units, duplex = extracted
    meta = {
        "pdb_id": "1abc",
        "assembly_id": "1",
        "model_id": "1",
        "alt_id": ".",
        "resolution": 2.0,
        "release_date": "2020-01-01",
        "deposition_date": "2019-01-01",
    }
    f = asm_fields(_complex(), asm)
    return unit_to_record(
        asm, f, units[idx], duplex, cfg, meta, asm.stats["ab_class"], idx
    )


def test_record_meta_and_node_tables(extracted, cfg):
    rec = _record(extracted, cfg)
    meta = rec["meta"]
    assert meta["unit_id"] == "1abc_1_1_0"
    assert meta["schema"] == SCHEMA
    assert meta["n_nodes"] == len(rec["nodes"]["kind"]) >= 2
    assert rec["nodes"]["kind"][:2] == ["protein", "protein"]
    assert rec["nodes"]["chain_id"][:2] == ["A", "B"]
    assert rec["nodes"]["auth_asym_id"][:2] == ["Aa", "Ba"]
    assert all(len(h) == 16 for h in rec["nodes"]["hash"])
    assert rec["nodes"]["hash"][0] == rec["nodes"]["hash"][1]  # same sequence
    assert rec["nodes"]["ab_class"][0] == ""
    assert len(meta["signature"]) == 16


def test_record_atoms_residues_and_edges_are_consistent(extracted, cfg):
    rec = _record(extracted, cfg)
    at, rs, ed = rec["atoms"], rec["residues"], rec["edges"]
    n = len(at["xyz"])
    assert at["xyz"].dtype == np.float32
    assert at["xyz"].shape == (n, 3)
    assert at["res"].max() == len(rs["chem_comp_id"]) - 1
    assert at["contact"].any()
    assert rs["interface"].any()
    assert len(ed["a"]) >= 1
    assert ed["res_pair_ptr"][0] == 0
    assert ed["res_pair_ptr"][-1] == len(ed["res_pairs"])
    # every interface residue belongs to some residue pair
    assert set(np.flatnonzero(rs["interface"])) == set(ed["res_pairs"].ravel())
    assert (ed["min_dist"] < 5.0).all()
    assert ed["center"].shape == (len(ed["a"]), 3)
    assert rs["node"].dtype == np.int16
    assert len(rec["duplex"]["a"]) == 0


def test_record_keeps_only_waters_near_the_unit(extracted, cfg):
    w = _record(extracted, cfg)["water"]
    assert len(w["xyz"]) == 1
    assert w["chain_id"].tolist() == ["W"]
    assert len(w["bridging"]) == 1


def test_manifest_row_summarises_the_record(extracted, cfg):
    rec = _record(extracted, cfg)
    row = manifest_row(rec)
    kinds = rec["nodes"]["kind"]
    assert row["unit_id"] == "1abc_1_1_0"
    assert row["n_nodes"] == len(kinds)
    assert row["type"].startswith("2protein")
    assert row["n_protein"] == 2
    assert row["n_atoms"] == len(rec["atoms"]["xyz"])
    assert row["n_interface_res"] == int(rec["residues"]["interface"].sum())
    assert row["node_kinds"] == "|".join(kinds)
    assert row["node_len"].split("|")[:2] == ["8", "8"]
    assert row["has_ab"] == 0
    assert row["min_dist"] == pytest.approx(float(rec["edges"]["min_dist"].min()))
    assert row["n_water"] == 1


def test_manifest_row_min_dist_is_blank_without_edges():
    rec = {
        "meta": {"unit_id": "u"},
        "nodes": {
            "kind": [],
            "key": [],
            "hash": [],
            "ab_class": [],
            "auth_asym_id": [],
        },
        "atoms": {"xyz": np.zeros((0, 3)), "contact": np.zeros(0, bool)},
        "residues": {"chem_comp_id": [], "interface": np.zeros(0, bool)},
        "edges": {"n_atom_pairs": np.zeros(0), "min_dist": np.zeros(0)},
        "water": {"xyz": [], "bridging": np.zeros(0, bool)},
        "duplex": {"a": []},
    }
    assert manifest_row(rec)["min_dist"] == ""
