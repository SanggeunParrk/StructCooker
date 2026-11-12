"""Re-Glyco reader: worklists, shard reads, the entry adapter and the sidecar."""

from __future__ import annotations

import datetime as dt
import hashlib
import json

import lmdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import zstandard as zstd
from glycan_fakes import PDB, SPEC, make_refs, same, value_for

from structcooker.instructions.readers import reglyco
from structcooker.instructions.transforms.delta_codec import (
    DeltaReader,
    canonical_of,
    is_delta,
)
from structcooker.instructions.transforms.reglyco_worklist import (
    build_worklists,
    write_reglyco_worklists,
)
from structcooker.release import sample_e2e

KEYS = [r.key for r in make_refs(SPEC)]


def make_shard(tmp_path, keys):
    """One shard of zstd frames (a fake CIF per record) and the index rows that locate them."""
    shard = tmp_path / "rs_00001.bin"
    rows, blob = [], b""
    for k in keys:
        cif = f"data_{k}\n#\n".encode() * 3
        frame = zstd.ZstdCompressor().compress(cif)
        rows.append(
            {
                "path": f"{k[:2]}/{k[:4]}/{k}.cif",
                "shard": shard.name,
                "offset": len(blob),
                "length": len(frame),
                "uncompressed_length": len(cif),
                "cif_sha256": hashlib.sha256(cif).hexdigest(),
            }
        )
        blob += frame
    shard.write_bytes(blob)
    return rows


def make_tables(rows):
    index = pa.Table.from_pylist(rows)
    ex = pa.table(
        {
            "path": [r["path"] for r in rows],
            "deposition_date": [dt.date(2020, 1, 1 + i % 20) for i in range(len(rows))],
            "resolution": [2.0 + 0.1 * i for i in range(len(rows))],
            "example_id": [f"{r['path']}|1|X" for r in rows],
            "n_atoms": list(range(len(rows))),
        }
    )
    gly = pa.table(
        {
            "path": ex.column("path"),
            "n_atoms": [i + 1000 for i in range(len(rows))],
            "extra": [1] * len(rows),
        }
    )
    return index, ex, gly


def test_keys():
    assert (
        reglyco.record_key("2p/2pmv/2pmv_A293_GS00178_conf012.cif")
        == "2pmv_A293_GS00178_conf012"
    )
    ref = reglyco.split_key("2pmv_A_293_GS00178_conf012")
    assert (ref.group, ref.conf) == ("A_293_GS00178", 12)
    with pytest.raises(ValueError, match="unexpected"):
        reglyco.split_key("not a record")


def test_worklists_sidecar_and_list(tmp_path):
    rows = make_shard(tmp_path, KEYS)
    index, ex, gly = make_tables(rows)
    pq.write_table(index, tmp_path / "index.parquet")
    pq.write_table(ex, tmp_path / "examples.parquet")
    pq.write_table(gly, tmp_path / "examples_glycan.parquet")
    msg = write_reglyco_worklists(
        str(tmp_path / "index.parquet"),
        str(tmp_path / "examples.parquet"),
        str(tmp_path / "examples_glycan.parquet"),
        str(tmp_path / "entries"),
        str(tmp_path / "entries.txt"),
        str(tmp_path / "side.parquet"),
    )
    assert f"{len(KEYS)} records of 1 entries" in msg
    assert (tmp_path / "entries.txt").read_text().split() == [f"{PDB}.json"]
    wl = json.loads((tmp_path / "entries" / f"{PDB}.json").read_text())
    assert set(wl) == {k.removeprefix(PDB + "_") for k in KEYS}
    rec = wl["A293_GS00001_conf000"]
    assert rec["key"] == f"{PDB}_A293_GS00001_conf000"
    assert rec["deposition_date"] == "2020-01-01"
    assert reglyco.read_worklist(tmp_path / "entries" / f"{PDB}.json") == {
        "records": wl
    }
    side = pq.read_table(tmp_path / "side.parquet")
    assert side.num_rows == len(KEYS)
    assert side.column("key").to_pylist() == sorted(KEYS)
    # the glycan-view column that differs is kept under gv_, a new one as it is, the Hub's path stays
    assert {"path", "n_atoms", "gv_n_atoms", "extra"} <= set(side.column_names)


def test_worklist_needs_an_examples_row(tmp_path):
    rows = make_shard(tmp_path, KEYS)
    index, ex, _ = make_tables(rows)
    with pytest.raises(ValueError, match="no examples row"):
        build_worklists(index, ex.slice(1))
    assert list(build_worklists(index, ex, {"zzzz"})) == []


def test_blob_read_checks_length_and_sha(tmp_path):
    rows = make_shard(tmp_path, KEYS[:2])
    rec = {
        "key": KEYS[0],
        "shard": rows[0]["shard"],
        "offset": rows[0]["offset"],
        "length": rows[0]["length"],
        "sha256": rows[0]["cif_sha256"],
    }
    assert reglyco.read_blob(tmp_path, rec).startswith(b"data_")
    with pytest.raises(ValueError, match="sha256"):
        reglyco.read_blob(tmp_path, {**rec, "sha256": "0" * 64})
    with pytest.raises(ValueError, match="shorter"):
        reglyco.read_blob(tmp_path, {**rec, "offset": 10**9})


def test_entry_adapter_writes_the_records_of_the_entry(tmp_path, monkeypatch):
    rows = make_shard(tmp_path, KEYS)
    index, ex, _ = make_tables(rows)
    wl = build_worklists(index, ex)[PDB]
    monkeypatch.setenv(reglyco.BLOB_ENV, str(tmp_path))
    seen = []

    def fake_convert(cif, key, meta, tmp_dir=None):
        seen.append((key, meta["deposition_date"]))
        assert cif.startswith(f"data_{key}".encode())
        return value_for(key, 300)

    monkeypatch.setattr(reglyco, "convert_record", fake_convert)
    out = reglyco.explode_reglyco_entry({"records": wl})
    assert sorted(seen) == sorted(
        (k, wl[k.removeprefix(PDB + "_")]["deposition_date"]) for k in KEYS
    )  # once each
    assert set(out) == {
        k.removeprefix(PDB + "_") for k in KEYS
    }  # the rebuild prefixes <pdb id>_
    values = {f"{PDB}_{sub}": v["reglyco_record"] for sub, v in out.items()}
    assert sum(is_delta(v) for v in values.values()) == len(KEYS) - 1
    env = lmdb.open(str(tmp_path / "db"), map_size=1 << 26)
    with env.begin(write=True) as txn:
        for k, v in values.items():
            txn.put(k.encode(), v)
    reader = DeltaReader(env)
    for k in KEYS:
        assert same(reader.read(k), value_for(k, 300))


def test_adapter_stops_on_a_record_that_does_not_match_the_index(tmp_path, monkeypatch):
    rows = make_shard(tmp_path, KEYS)
    index, ex, _ = make_tables(rows)
    wl = build_worklists(index, ex)[PDB]
    next(iter(wl.values()))["sha256"] = "0" * 64
    monkeypatch.setenv(reglyco.BLOB_ENV, str(tmp_path))
    monkeypatch.setattr(
        reglyco,
        "convert_record",
        lambda *_a, **_k: value_for("x_A1_GS00001_conf000", 50),
    )
    with pytest.raises(ValueError, match="sha256"):
        reglyco.explode_reglyco_entry({"records": wl})


def test_adapt_annotated_cif_derives_the_missing_tables():
    d = {
        "_entity": {"id": ["1", "2"], "type": ["polymer", "branched"]},
        "_chem_comp": {"id": ["ASN"], "name": ["ASN"], "type": ["L-peptide linking"]},
        "_atom_site": {
            "label_asym_id": ["A", "A", "B", "B"],
            "label_entity_id": ["1", "1", "2", "2"],
            "label_seq_id": ["1", "1", ".", "."],
            "label_comp_id": ["ASN", "ASN", "NAG", "NAG"],
            "auth_seq_id": ["293", "293", "1", "1"],
            "pdbx_PDB_ins_code": ["?"] * 4,
            "pdbx_PDB_model_num": ["1"] * 4,
        },
        "_entity_poly_seq": {"entity_id": ["1"], "num": ["1"], "mon_id": ["ASN"]},
        "_struct_conn": {
            "id": ["c1", "c2"],
            "ptnr1_label_asym_id": ["A", "B"],
            "ptnr2_label_asym_id": ["B", "B"],
            "ptnr1_label_comp_id": ["ASN", "NAG"],
            "ptnr2_label_comp_id": ["NAG", "NAG"],
            "ptnr1_label_atom_id": ["ND2", "O4"],
            "ptnr2_label_atom_id": ["C1", "C1"],
            "ptnr1_label_seq_id": ["1", "."],
            "ptnr2_label_seq_id": [".", "."],
            "pdbx_value_order": ["sing", "sing"],
        },
    }
    diag = reglyco.adapt_annotated_cif(d)
    assert diag["chem_comp_added"] == 1
    assert "NAG" in d["_chem_comp"]["id"]
    assert d["_pdbx_poly_seq_scheme"]["mon_id"] == ["ASN"]
    assert d["_pdbx_branch_scheme"]["num"] == ["1"]
    assert d["_pdbx_entity_branch_list"]["comp_id"] == ["NAG"]
    assert d["_entity_poly_seq"]["hetero"] == ["n"]
    # the sugar-sugar link moved to _pdbx_entity_branch_link; the protein-glycan link stays in _struct_conn
    assert diag["branch_links_moved"] == 1
    assert d["_struct_conn"]["id"] == ["c1"]
    assert d["_pdbx_entity_branch_link"]["leaving_atom_id_1"] == [
        "?"
    ]  # unknown, not guessed
    assert d["_struct_conn"]["ptnr1_auth_seq_id"] == ["293"]


def test_sample_e2e_compares_delta_dbs_on_the_decoded_records(tmp_path):
    from structcooker.instructions.transforms.delta_codec import encode_entry

    refs = make_refs(SPEC)
    prod = encode_entry(refs, lambda k: value_for(k, 300))
    # a sample holding only part of the entry stores other bases: the bytes differ, the records do not
    part = [r for r in refs if not r.key.endswith("GS00001_conf000")]
    sample = encode_entry(part, lambda k: value_for(k, 300))
    dbs = {}
    for name, values in (("sample", sample), ("prod", prod)):
        p = tmp_path / f"{name}.lmdb"
        p.mkdir()
        env = lmdb.open(str(p), map_size=1 << 26)
        with env.begin(write=True) as txn:
            for k, v in values.items():
                txn.put(k.encode(), v)
        env.close()
        dbs[name] = p
    report = sample_e2e.compare_lmdb(dbs["sample"], dbs["prod"], "K")
    assert report["records"] == len(part)
    assert report.get("different", 0) == 0
    assert canonical_of(
        value_for(KEYS[0], 300)
    )  # the helper used by the checks stays importable


def test_a_lone_sugar_is_a_non_polymer_like_in_a_pdb_mmcif():
    d = {
        "_entity": {"id": ["1", "2"], "type": ["polymer", "branched"]},
        "_chem_comp": {"id": ["ASN"], "name": ["ASN"], "type": ["L-peptide linking"]},
        "_atom_site": {
            "label_asym_id": ["A", "B"],
            "label_entity_id": ["1", "2"],
            "label_seq_id": ["1", "."],
            "label_comp_id": ["ASN", "NAG"],
            "auth_seq_id": ["293", "1"],
            "pdbx_PDB_ins_code": ["?"] * 2,
            "pdbx_PDB_model_num": ["1"] * 2,
        },
        "_entity_poly_seq": {"entity_id": ["1"], "num": ["1"], "mon_id": ["ASN"]},
    }
    diag = reglyco.adapt_annotated_cif(d)
    assert diag["lone_sugar_entities"] == 1
    assert d["_entity"]["type"] == ["polymer", "non-polymer"]
    assert d["_pdbx_nonpoly_scheme"]["mon_id"] == ["NAG"]
    assert d["_pdbx_entity_nonpoly"] == {"entity_id": ["2"], "comp_id": ["NAG"]}
    assert d["_pdbx_branch_scheme"]["mon_id"] == []
