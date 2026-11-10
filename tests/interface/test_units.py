from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from datacooker.config import load_config as load_engine_config
from datacooker.lmdb.core import _emit_rebuilt, rebuild_entry
from datacooker.writers import WriterHooks
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.build import schemas
from structcooker.instructions.interface import units
from structcooker.instructions.interface.config import load_config
from structcooker.instructions.interface.units import (
    UnitInput,
    build_unit_record,
    entry_pdb_id,
    entry_records,
    entry_units,
    explode_units,
    select_assembly_keys,
)
from structcooker.instructions.transforms.codecs import from_bytes, to_bytes

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def _mol():
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
    return make_mol(
        [protein_chain("A", 8), protein_chain("B", 8, z=4.0, sign=-1.0), lig]
    )


def _entry(assemblies=("1_1_.",), **md):
    """Build a decoded cif_pdb entry whose assembly items carry the fake mol (CIFMol is patched out)."""
    return {
        "assembly_dict": {k: {"mol": _mol()} for k in assemblies},
        "metadata_dict": {
            "id": ["1ABC"],
            "release_date": "2020-01-01",
            "deposition_date": "2019-01-01",
            "resolution": 2.0,
            **md,
        },
    }


@pytest.fixture(autouse=True)
def fake_cifmol(monkeypatch):
    monkeypatch.setattr(
        units, "CIFMol", SimpleNamespace(from_dict=lambda item: item["mol"])
    )
    monkeypatch.setattr(units, "default_classifier", lambda: lambda _seq: None)


def test_select_assembly_keys_takes_the_first_model_and_prefers_no_altloc():
    keys = ["1_2_.", "1_1_B", "1_1_A", "2_1_.", "2_1_A", "10_1_."]
    assert select_assembly_keys(keys) == ["1_1_A", "10_1_.", "2_1_."]
    assert select_assembly_keys(["1_1_A", "1_1_."]) == ["1_1_."]
    assert select_assembly_keys([]) == []


def test_entry_pdb_id_accepts_a_list_or_a_string():
    assert entry_pdb_id({"id": ["1TSR"]}) == "1tsr"
    assert entry_pdb_id({"id": "1TSR"}) == "1tsr"


def test_entry_units_are_sub_entries_keyed_assembly_model_index(cfg):
    subs = entry_units(_entry(), cfg, None)
    assert list(subs) == ["1_1_0"]
    ui = subs["1_1_0"]["unit_input"]
    assert isinstance(ui, UnitInput)
    assert ui.index == 0
    assert ui.meta == {
        "pdb_id": "1abc",
        "resolution": 2.0,
        "release_date": "2020-01-01",
        "deposition_date": "2019-01-01",
        "assembly_id": "1",
        "model_id": "1",
        "alt_id": ".",
    }


def test_only_the_first_model_of_an_assembly_is_used(cfg):
    subs = entry_units(_entry(("1_2_.", "1_1_.", "1_3_.")), cfg, None)
    assert {u["unit_input"].meta["model_id"] for u in subs.values()} == {"1"}


def test_units_repeated_in_another_assembly_are_kept_once(cfg):
    subs = entry_units(_entry(("1_1_.", "2_1_.", "3_1_.")), cfg, None)
    assert list(subs) == ["1_1_0"]


def test_assemblies_with_too_many_polymer_chains_are_skipped():
    small = load_config(overrides={"assembly.max_polymer_chains": 1})
    assert entry_units(_entry(), small, None) == {}


def test_assembly_without_unit_gives_no_sub_entry(cfg):
    lone = {
        "assembly_dict": {"1_1_.": {"mol": make_mol([protein_chain("A", 8)])}},
        "metadata_dict": _entry()["metadata_dict"],
    }
    assert entry_units(lone, cfg, None) == {}


def test_missing_or_nan_resolution_is_stored_as_none(cfg):
    for res in (None, float("nan")):
        subs = entry_units(_entry(resolution=res), cfg, None)
        assert subs["1_1_0"]["unit_input"].meta["resolution"] is None


def test_antibody_classes_reach_the_record(cfg):
    recs = entry_records(_entry(), cfg, ab_lookup=lambda _seq: "H")
    assert recs["1abc_1_1_0"]["nodes"]["ab_class"][:2] == ["H", "H"]


def test_build_unit_record_is_unit_to_record(cfg):
    ui = entry_units(_entry(), cfg, None)["1_1_0"]["unit_input"]
    rec = build_unit_record(ui, cfg)
    assert rec["meta"]["unit_id"] == "1abc_1_1_0"
    assert rec["meta"]["pdb_id"] == "1abc"
    assert rec["nodes"]["kind"][:2] == ["protein", "protein"]
    assert not schemas.get("J").validate(rec)


def test_entry_records_are_keyed_by_unit_id(cfg):
    recs = entry_records(_entry(("1_1_.", "2_1_.")), cfg)
    assert list(recs) == ["1abc_1_1_0"]
    assert recs["1abc_1_1_0"]["meta"]["unit_id"] == "1abc_1_1_0"


def test_default_classifier_reads_the_cache_named_by_the_environment(
    monkeypatch, tmp_path
):
    monkeypatch.undo()  # the autouse fixture patched it; this test wants the real one
    cache = tmp_path / "ab.tsv"
    cache.write_text("QVQL\tH\nEVQL\t\n")
    monkeypatch.setenv(units.ANTIBODY_CACHE_ENV, str(cache))
    units.default_classifier.cache_clear()
    try:
        ab = units.default_classifier()
        assert ab("QVQL") == "H"
        assert ab("EVQL") is None
    finally:
        units.default_classifier.cache_clear()


def test_the_stage_config_wires_adapter_recipe_and_explode_into_unit_records():
    """The rebuild engine path of db_configs/interface/interface_pdb.yaml on a synthetic entry."""
    engine = load_engine_config(REPO / "db_configs/interface/interface_pdb.yaml")
    assert engine["split_entries"]
    assert engine["explode_entries"]
    assert engine["reader"]["adapter"] is explode_units
    data = engine["reader"]["adapter"](_entry(("1_1_.", "2_1_.")))
    rebuilt = rebuild_entry(
        recipe=REPO / engine["recipe"],
        data=data,
        metadata_dict={},
        parameter_dict={},
        key_transform=None,
        split_entries=True,
        extra_kwargs={},
    )
    out = _emit_rebuilt(
        "1abc",
        rebuilt,
        writer=WriterHooks(serializer=to_bytes),
        explode_entries=True,
        explode_key_sep="_",
        explode_value_key=engine["explode_value_key"],
    )
    assert [k.decode() for k, _, _ in out] == ["1abc_1_1_0"]
    record = from_bytes(out[0][1])
    assert not schemas.get(str(engine.get("schema", "J"))).validate(record)
    assert record["meta"]["unit_id"] == "1abc_1_1_0"


def test_schema_j_names_what_is_missing():
    assert schemas.get("J").validate([]) == ["not a dict: list"]
    assert "missing water" in schemas.get("J").validate(
        {k: {} for k in ("meta", "nodes", "atoms", "residues", "edges", "duplex")}
    )
    rec = {k: {} for k in ("nodes", "residues", "edges", "duplex", "water")}
    rec["meta"] = {"unit_id": "u"}
    rec["atoms"] = {"xyz": np.zeros((2, 3)), "node": np.zeros(3)}
    assert schemas.get("J").validate(rec) == [
        "atoms.xyz and atoms.node differ in length"
    ]
    rec["atoms"] = {"xyz": np.zeros((0, 3)), "node": np.zeros(0)}
    assert schemas.get("J").validate(rec) == ["no atoms"]


def test_the_manifest_stage_turns_a_stored_record_into_a_manifest_row(tmp_path):
    from datacooker.api import parse_dict

    from structcooker.instructions.interface.tables import read_tsv, write_manifest

    engine = load_engine_config(REPO / "db_configs/interface/manifest.yaml")
    assert engine["writer"]["materializer"] is write_manifest
    record = from_bytes(to_bytes(entry_records(_entry())["1abc_1_1_0"]))
    result = parse_dict(
        recipe_path=REPO / engine["extract_recipe_path"],
        datadict={"key": "1abc_1_1_0", "db_data": record},
    )
    out = tmp_path / "manifest.tsv"
    write_manifest({"1abc_1_1_0": result}, out)
    (row,) = read_tsv(out)
    assert row["unit_id"] == "1abc_1_1_0"
    assert row["type"] == "2protein+1ligand"
    assert row["n_atoms"] == str(len(record["atoms"]["xyz"]))
