"""RFD4-Proteina sets: record keys, the worklist, the shard read and the entity tables derived from a thin CIF."""

from __future__ import annotations

import hashlib

import lmdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import zstandard as zstd

from structcooker.instructions.readers import rfd4
from structcooker.instructions.transforms import rfd4_worklist as wl
from structcooker.instructions.transforms.codecs import from_bytes


def thin_site(residues: list[tuple[str, str, str]]) -> dict:
    """Make an ``_atom_site`` of (asym, comp, label_seq_id) residues, two atoms each, one model, no entity ids."""
    cols: dict[str, list[str]] = {
        c: []
        for c in (
            [
                "id",
                "label_asym_id",
                "label_seq_id",
                "label_comp_id",
                "auth_seq_id",
                "label_entity_id",
                "pdbx_PDB_model_num",
            ]
        )
    }
    for i, (asym, comp, seq) in enumerate(residues):
        for _ in range(2):
            cols["id"].append(str(len(cols["id"]) + 1))
            cols["label_asym_id"].append(asym)
            cols["label_seq_id"].append(seq)
            cols["label_comp_id"].append(comp)
            cols["auth_seq_id"].append(str(i + 1))
            cols["label_entity_id"].append(".")
            cols["pdbx_PDB_model_num"].append("1")
    return cols


def test_protein_dna_entities_are_derived():
    d = {
        "_atom_site": thin_site(
            [
                ("A", "MET", "1"),
                ("A", "ALA", "2"),
                ("B", "DA", "1"),
                ("B", "DG", "2"),
                ("C", "MET", "1"),
                ("C", "ALA", "2"),
            ]
        )
    }
    info = rfd4.synthesize_tables(d)
    assert info == {
        "entities_made": 2,
        "polymers_made": 2,
    }  # A and C share a sequence: one entity
    types = dict(zip(d["_entity"]["id"], d["_entity"]["type"], strict=True))
    assert set(types.values()) == {"polymer"}
    kinds = dict(
        zip(d["_entity_poly"]["entity_id"], d["_entity_poly"]["type"], strict=True)
    )
    assert sorted(kinds.values()) == ["polydeoxyribonucleotide", "polypeptide(L)"]
    codes = dict(
        zip(
            d["_entity_poly"]["entity_id"],
            d["_entity_poly"]["pdbx_seq_one_letter_code_can"],
            strict=True,
        )
    )
    assert sorted(codes.values()) == ["AG", "MA"]
    assert d["_entity_poly_seq"]["num"] == ["1", "2", "1", "2"]
    assert d["_chem_comp"]["id"] == ["MET", "ALA", "DA", "DG"]


def test_unknown_residue_letters_and_stated_sequence():
    d = {
        "_atom_site": thin_site(
            [("A", "MET", "1"), ("A", "OTY", "2"), ("A", "ALA", "3")]
        )
    }
    rfd4.synthesize_tables(d)
    assert d["_entity_poly"]["pdbx_seq_one_letter_code_can"] == ["MXA"]
    d = {
        "_atom_site": thin_site(
            [("A", "MET", "1"), ("A", "OTY", "2"), ("A", "ALA", "3")]
        )
    }
    rfd4.synthesize_tables(d, sequence="MYA")  # the dataset's own canonical sequence
    assert d["_entity_poly"]["pdbx_seq_one_letter_code_can"] == ["MYA"]
    d = {"_atom_site": thin_site([("A", "MET", "1"), ("A", "OTY", "2")])}
    rfd4.synthesize_tables(d, sequence="MYAA")  # wrong length: not used
    assert d["_entity_poly"]["pdbx_seq_one_letter_code_can"] == ["MX"]


def test_existing_tables_are_kept_and_missing_columns_filled():
    d = {
        "_atom_site": thin_site([("A", "G", "1"), ("A", "C", "2")]),
        "_entity": {"id": ["1"], "type": ["polymer"]},
        "_entity_poly": {"entity_id": ["1"], "type": ["polyribonucleotide"]},
        "_entity_poly_seq": {
            "entity_id": ["1", "1"],
            "num": ["1", "2"],
            "mon_id": ["G", "C"],
            "hetero": ["n", "n"],
        },
        "_chem_comp": {"id": ["G", "C"], "name": ["g", "c"]},
    }
    d["_atom_site"]["label_entity_id"] = ["1"] * 4
    assert rfd4.synthesize_tables(d) == {"entities_made": 0, "polymers_made": 0}
    rfd4._ensure_poly_codes(d)
    assert d["_entity_poly"]["pdbx_seq_one_letter_code_can"] == ["GC"]
    assert d["_chem_comp"]["name"] == ["g", "c"]


def test_single_residue_chain_is_not_a_polymer_unless_standard():
    d = {
        "_atom_site": thin_site(
            [("A", "ATP", "1"), ("B", "HOH", "1"), ("C", "ALA", "1")]
        )
    }
    rfd4.synthesize_tables(d)
    types = dict(zip(d["_entity"]["id"], d["_entity"]["type"], strict=True))
    assert sorted(types.values()) == ["non-polymer", "polymer", "water"]


def test_atom_site_defaults_and_entry_date():
    site = thin_site([("A", "ALA", "1")])
    rfd4._fill_atom_site(site)
    assert site["occupancy"] == ["1.0", "1.0"]
    assert site["B_iso_or_equiv"] == ["0.0", "0.0"]
    assert (
        rfd4.entry_date({"_entry": {"date": ["2026-10-01"]}}, "2026-01-01")
        == "2026-10-01"
    )
    assert rfd4.entry_date({"_entry": {"date": ["?"]}}, "2026-01-01") == "2026-01-01"
    assert rfd4.entry_date({}, "2026-01-01") == "2026-01-01"


def test_parse_cif_text_nests_by_category():
    text = "data_x\n#\n_entry.id x\n_entry.date 2026-01-02\n#\nloop_\n_atom_site.id\n_atom_site.label_comp_id\n1 ALA\n2 GLY\n"
    d = rfd4.parse_cif_text(text)
    assert d["_entry"]["id"] == ["x"]
    assert d["_atom_site"]["label_comp_id"] == ["ALA", "GLY"]


# ---------------------------------------------------------------------------------------------- keys and worklist
def test_record_keys():
    assert wl.make_keys(pa.array(["a_1", "b_2"]), "id") == ["a_1", "b_2"]
    assert wl.make_keys(pa.array(["tcr_10.cif", "x/tcr_11.cif"]), "stem") == [
        "tcr_10",
        "tcr_11",
    ]
    assert wl.make_keys(pa.array([7, 12]), "idx") == ["000000007", "000000012"]
    assert wl.make_keys(pa.array(["7"]), "idx") == [
        "000000007"
    ]  # OMol's index holds lmdb_idx as text


def fake_component(
    tmp_path, monkeypatch, *, key_col="example_id", kind="id", sequence=False
):
    root = tmp_path / "raw"
    base = wl.component_dir(
        root, "tcr_distill"
    )  # any registered path; COMPONENTS is patched below
    comp = wl.Component(
        "tcr_distill_fake",
        key_col,
        kind,
        "2026-09-09",
        sequence_col="sequence" if sequence else None,
    )
    monkeypatch.setitem(wl.COMPONENTS, "tcr_distill", comp)
    base = wl.component_dir(root, "tcr_distill")
    (base / "blob/data").mkdir(parents=True)
    ids = ["r1", "r2", "r3"]
    shard, frames, rows, blob = "s0.bin", [], [], b""
    for i in ids:
        cif = f"data_{i}\n#\n".encode()
        fr = zstd.ZstdCompressor().compress(cif)
        rows.append(
            {
                key_col: i,
                "shard": shard,
                "offset": len(blob),
                "length": len(fr),
                "sha256": hashlib.sha256(cif).hexdigest(),
            }
        )
        blob += fr
        frames.append(cif)
    (base / "blob/data" / shard).write_bytes(blob)
    pq.write_table(pa.Table.from_pylist(rows), base / "blob/index.parquet")
    ex = [
        {key_col: i, "score": float(n), "sequence": "AC"} for n, i in enumerate(ids)
    ] + [{key_col: "r2", "score": 9.0, "sequence": "AC"}]
    pq.write_table(pa.Table.from_pylist(ex), base / "examples.parquet")
    return root, base, frames


def test_worklist_rows_locate_the_records(tmp_path, monkeypatch):
    root, base, frames = fake_component(tmp_path, monkeypatch, sequence=True)
    (rows,) = list(wl.iter_worklist_rows("tcr_distill", root))
    assert [r["key"] for r in rows] == ["r1", "r2", "r3"]
    assert rows[0]["shard"] == "data/s0.bin"  # resolved relative to blob/
    assert rows[0]["sequence"] == "AC"
    assert [rfd4.read_blob(base / "blob", r).encode() for r in rows] == frames
    bad = dict(rows[1], sha256="0" * 64)
    with pytest.raises(ValueError, match="sha256"):
        rfd4.read_blob(base / "blob", bad)
    short = dict(rows[2], length=rows[2]["length"] + 50)
    with pytest.raises(ValueError, match="shorter"):
        rfd4.read_blob(base / "blob", short)


def test_write_worklist_stores_rows_and_a_sidecar_with_repeated_keys(
    tmp_path, monkeypatch
):
    root, _, _ = fake_component(tmp_path, monkeypatch)
    out = tmp_path / "out"
    msg = wl.write_worklist(
        "tcr_distill",
        str(root),
        str(out / "w.lmdb"),
        str(out / "side.parquet"),
        str(out / "w.txt"),
    )
    assert "3 records" in msg
    env = lmdb.open(str(out / "w.lmdb"), readonly=True, lock=False)
    with env.begin() as txn:
        row = from_bytes(bytes(txn.get(b"r2")))
    assert row["key"] == "r2"
    assert row["date"] == "2026-09-09"
    assert row["shard"] == "data/s0.bin"
    side = pq.read_table(out / "side.parquet")
    assert side.column("key").to_pylist() == [
        "r1",
        "r2",
        "r2",
        "r3",
    ]  # r2 has two examples rows
    assert (out / "w.txt").read_text().strip() == str(out / "w.lmdb")


def test_missing_shard_fails_loudly(tmp_path, monkeypatch):
    root, base, _ = fake_component(tmp_path, monkeypatch)
    (base / "blob/data/s0.bin").unlink()
    with pytest.raises(FileNotFoundError, match="s0.bin"):
        list(wl.iter_worklist_rows("tcr_distill", root))


def test_every_component_has_a_distinct_blob_dir_and_stage():
    subdirs = [c.subdir for c in wl.COMPONENTS.values()]
    assert len(set(subdirs)) == len(subdirs)
    from pathlib import Path

    for name in wl.COMPONENTS:
        assert (Path("db_configs/rfd4") / f"{name}.yaml").is_file()
        assert (Path("db_configs/rfd4") / f"{name}_worklist.yaml").is_file()


def test_a_component_without_a_definition_fails_by_name(monkeypatch):
    monkeypatch.setattr(rfd4, "_ccd_has", lambda c: c in {"ALA"})
    d = {"_atom_site": thin_site([("A", "ALA", "1"), ("B", "Q:1", "1")])}
    with pytest.raises(ValueError, match=r"Q:1"):
        rfd4._require_definitions(d)
    d["_chem_comp_atom"] = {"comp_id": ["Q:1"]}
    rfd4._require_definitions(d)  # the record defines it itself


def test_synthesized_chem_comp_has_the_columns_the_recipe_reads():
    d = {"_atom_site": thin_site([("A", "A:1", "1")])}
    rfd4.synthesize_tables(d)
    assert {"id", "name", "formula"} <= set(d["_chem_comp"])


def test_blob_dir_follows_biomol_root_and_the_rows_set(tmp_path, monkeypatch):
    monkeypatch.delenv(rfd4.BLOB_ENV, raising=False)
    monkeypatch.setenv("BIOMOL_ROOT", str(tmp_path))
    got = rfd4.blob_dir({"set": "tcr_distill"})
    assert (
        got
        == tmp_path
        / "materials/raw/rfd4_proteina"
        / wl.SYNTHETIC
        / "2026_09_09_tcr_af2/blob"
    )
    monkeypatch.setenv(rfd4.BLOB_ENV, str(tmp_path / "elsewhere"))
    assert rfd4.blob_dir({"set": "tcr_distill"}) == tmp_path / "elsewhere"


def test_worklist_rows_name_their_set(tmp_path, monkeypatch):
    root, _, _ = fake_component(tmp_path, monkeypatch)
    (rows,) = list(wl.iter_worklist_rows("tcr_distill", root))
    assert {r["set"] for r in rows} == {"tcr_distill"}
