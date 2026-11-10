"""Monomer units: polymer nodes that no contact unit holds become edge-less units of their own."""

from types import SimpleNamespace

import numpy as np
import pytest
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.build import schemas
from structcooker.instructions.unit import units
from structcooker.instructions.unit.config import load_config
from structcooker.instructions.unit.detect import Unit, monomer_units
from structcooker.instructions.unit.pipeline import extract_units
from structcooker.instructions.unit.record import manifest_row
from structcooker.instructions.unit.tables import CLASS_TOTALS
from structcooker.instructions.unit.units import entry_records, entry_units


@pytest.fixture(scope="module")
def cfg():
    return load_config()


@pytest.fixture(autouse=True)
def fake_cifmol(monkeypatch):
    monkeypatch.setattr(
        units, "CIFMol", SimpleNamespace(from_dict=lambda item: item["mol"])
    )
    monkeypatch.setattr(units, "default_classifier", lambda: lambda _seq: None)


def _entry(chains, assemblies=("1_1_.",)):
    return {
        "assembly_dict": {k: {"mol": make_mol(chains)} for k in assemblies},
        "metadata_dict": {
            "id": ["1ABC"],
            "release_date": "2020-01-01",
            "deposition_date": "2019-01-01",
            "resolution": 2.0,
        },
    }


def _pair():
    """Two 8-mers in contact (one contact unit)."""
    return [protein_chain("A", 8), protein_chain("B", 8, z=4.0, sign=-1.0)]


def _lone(cid, seq="G" * 8, z=100.0):
    return protein_chain(cid, len(seq), z=z, seq=seq)


def _water(x, y, z, cid="W"):
    return chain(cid, [residue("HOH", "X", [("O", "O", (x, y, z))])], etype="water")


def test_a_lone_chain_is_a_monomer_unit(cfg):
    recs = entry_records(_entry([protein_chain("A", 8)]), cfg)
    assert list(recs) == ["1abc_1_1_0"]
    rec = recs["1abc_1_1_0"]
    assert rec["meta"]["n_nodes"] == 1
    assert rec["meta"]["n_seed_edges"] == 0
    assert rec["nodes"]["kind"] == ["protein"]
    assert len(rec["edges"]["a"]) == 0
    assert not rec["atoms"]["contact"].any()
    assert not rec["residues"]["interface"].any()
    assert len(rec["atoms"]["xyz"]) == 8 * 5
    assert not schemas.get("J").validate(rec)


def test_monomer_manifest_row_has_no_contact_columns(cfg):
    rec = entry_records(_entry([protein_chain("A", 8)]), cfg)["1abc_1_1_0"]
    row = manifest_row(rec)
    assert row["type"] == "1protein"
    assert row["n_seed_edges"] == 0
    assert row["n_atom_pairs"] == 0
    assert row["min_dist"] == ""
    assert row["n_interface_res"] == 0


def test_monomer_units_follow_the_contact_units_and_do_not_renumber_them(cfg):
    with_lone = entry_records(_entry([*_pair(), _lone("C")]), cfg)
    without = entry_records(_entry(_pair()), cfg)
    assert list(without) == ["1abc_1_1_0"]
    assert list(with_lone) == ["1abc_1_1_0", "1abc_1_1_1"]
    assert (
        with_lone["1abc_1_1_0"]["meta"]["signature"]
        == without["1abc_1_1_0"]["meta"]["signature"]
    )
    mono = with_lone["1abc_1_1_1"]
    assert mono["nodes"]["key"] == ["P:" + "G" * 8]
    assert mono["meta"]["n_seed_edges"] == 0


def test_chains_of_a_contact_unit_are_not_also_monomers(cfg):
    assert list(entry_records(_entry(_pair()), cfg)) == ["1abc_1_1_0"]


def test_symmetry_copies_of_a_monomer_are_kept_once(cfg):
    chains = [_lone("C", z=100.0), _lone("D", z=200.0), _lone("E", z=300.0)]
    assert len(entry_units(_entry(chains), cfg, None)) == 1
    # a repeated assembly adds nothing either
    subs = entry_units(_entry(chains, ("1_1_.", "2_1_.")), cfg, None)
    assert list(subs) == ["1_1_0"]


def test_different_sequences_are_different_monomers(cfg):
    chains = [_lone("C", "G" * 8, z=100.0), _lone("D", "A" * 8, z=200.0)]
    recs = entry_records(_entry(chains), cfg)
    keys = sorted(r["nodes"]["key"][0] for r in recs.values())
    assert keys == ["P:" + "A" * 8, "P:" + "G" * 8]


def test_min_residues_is_counted_over_resolved_residues(cfg):
    assert entry_units(_entry([protein_chain("A", 3)]), cfg, None) == {}
    assert len(entry_units(_entry([protein_chain("A", 4)]), cfg, None)) == 1
    loose = load_config(overrides={"monomer.min_residues": 1})
    assert len(entry_units(_entry([protein_chain("A", 3)]), loose, None)) == 1


def test_monomers_can_be_switched_off(cfg):
    off = load_config(overrides={"monomer.enabled": False})
    assert entry_units(_entry([protein_chain("A", 8)]), off, None) == {}
    asm, got, _ = extract_units(make_mol([*_pair(), _lone("C")]), off)
    assert len(got) == 1
    assert [len(u.nodes) for u in got] == [2]


def test_the_polymer_chain_cap_still_skips_the_assembly():
    small = load_config(overrides={"assembly.max_polymer_chains": 1})
    chains = [_lone("C", z=100.0), _lone("D", "A" * 8, z=200.0)]
    assert entry_units(_entry(chains), small, None) == {}


def test_waters_around_a_monomer_are_context_and_never_bridging(cfg):
    near = _water(3.0, 3.0, 0.0)
    far = _water(50.0, 50.0, 50.0, "V")
    rec = entry_records(_entry([protein_chain("A", 8), near, far]), cfg)["1abc_1_1_0"]
    assert len(rec["water"]["xyz"]) == 1
    assert not rec["water"]["bridging"].any()


def test_ligand_that_only_touches_a_free_chain_makes_a_contact_unit_not_a_monomer(cfg):
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
    recs = entry_records(_entry([protein_chain("A", 8), lig]), cfg)
    assert [r["nodes"]["kind"] for r in recs.values()] == [["protein", "ligand"]]


# --- the selection rule on its own ---------------------------------------------------------
def test_monomer_units_skips_covered_and_non_polymer_nodes():
    kinds = ["protein", "protein", "ligand", "ion", "na"]
    got = monomer_units([Unit([0, 2])], kinds, [10, 10, 5, 1, 10])
    assert [u.nodes for u in got] == [[1], [4]]
    assert all(u.edges == [] and u.seeds == [] and not u.truncated for u in got)


def test_monomer_units_keeps_base_paired_free_strands_together():
    kinds = ["na", "na", "na", "protein"]
    got = monomer_units([], kinds, [12, 12, 12, 12], {(0, 1): 9})
    assert [u.nodes for u in got] == [[0, 1], [2], [3]]


def test_monomer_units_ignores_a_duplex_partner_that_is_in_a_unit():
    kinds = ["na", "na", "protein"]
    got = monomer_units([Unit([1, 2])], kinds, [12, 12, 12], {(0, 1): 9})
    assert [u.nodes for u in got] == [[0]]


def test_monomer_units_min_residues():
    got = monomer_units([], ["protein", "protein"], np.array([3, 4]), None, 4)
    assert [u.nodes for u in got] == [[1]]


def test_manifest_class_totals_separate_monomers_from_protein_protein():
    totals = dict(CLASS_TOTALS)
    mono = {"n_na": "0", "n_ligand": "0", "n_ion": "0", "n_seed_edges": "0"}
    pair = {**mono, "n_seed_edges": "1"}
    mono_pred = totals["monomer units (no contact edge, n_seed_edges = 0)"]
    pp_pred = totals["protein-protein only (no NA / ligand / ion)"]
    assert mono_pred(mono)
    assert not mono_pred(pair)
    assert pp_pred(pair)
    assert not pp_pred(mono)
