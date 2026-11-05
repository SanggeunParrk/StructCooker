"""AFDB multimer key builders and CIF override resolution; LMDB cif loader."""

import json
from pathlib import Path

import pytest

from structcooker.instructions.readers import afdb
from structcooker.instructions.readers import io as reader_io


@pytest.mark.parametrize(
    ("name", "key"),
    [
        ("AF-123-model_v1.pdb", "AF-123"),
        ("AF-7-model_v1.cif.gz", "AF-7"),
        ("AF-42-msa_v1.a3m.zst", "AF-42"),
    ],
)
def test_afdb_entity_key_drops_suffix(name, key):
    assert afdb.afdb_entity_key(Path("/x/y") / name) == key


def test_afdb_entity_key_rejects_foreign_names():
    with pytest.raises(ValueError, match="not an AFDB release file name: foo.pdb"):
        afdb.afdb_entity_key(Path("foo.pdb"))


def test_afm_msa_seqid_key_reads_map_once_and_caches(tmp_path, monkeypatch):
    m = tmp_path / "map.tsv"
    m.write_text("AF-1\tP1\nAF-2\tP2\nbad-line\n")
    monkeypatch.setenv("AFM_MSA_SEQID_MAP", str(m))
    monkeypatch.setattr(afdb, "_AFM_MSA_SEQID", None)
    assert afdb.afm_msa_seqid_key(Path("AF-2-msa_v1.a3m.zst")) == "P2"
    m.unlink()  # a second lookup must come from the process cache
    assert afdb.afm_msa_seqid_key(Path("AF-1-msa_v1.a3m.zst")) == "P1"
    with pytest.raises(KeyError):
        afdb.afm_msa_seqid_key(Path("AF-3-msa_v1.a3m.zst"))
    assert "bad-line" not in afdb._afm_msa_seqid_map()


def test_afm_msa_seqid_map_default_path_is_under_biomol_root(tmp_path, monkeypatch):
    monkeypatch.delenv("AFM_MSA_SEQID_MAP", raising=False)
    monkeypatch.setenv("BIOMOL_ROOT", str(tmp_path))
    monkeypatch.setattr(afdb, "_AFM_MSA_SEQID", None)
    with pytest.raises(FileNotFoundError, match="materials/intermediate/worklists"):
        afdb._afm_msa_seqid_map()


def test_get_afm_cif_data_substitutes_repaired_copy(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(afdb, "get_cif_data", lambda p: seen.append(p) or {"p": p})
    manifest = tmp_path / "ov" / "overrides.json"
    manifest.parent.mkdir()
    manifest.write_text(
        json.dumps({"heterodimer/cif/267/AF-9-model_v1.cif.gz": "fixed/AF-9.cif.gz"})
    )
    monkeypatch.setenv("AFM_CIF_OVERRIDES", str(manifest))
    monkeypatch.setattr(afdb, "_CIF_OVERRIDES", None)

    broken = Path("/raw/afdb/heterodimer/cif/267/AF-9-model_v1.cif.gz")
    assert afdb.get_afm_cif_data(broken) == {"p": manifest.parent / "fixed/AF-9.cif.gz"}
    other = Path("/raw/afdb/heterodimer/cif/267/AF-10-model_v1.cif.gz")
    assert afdb.get_afm_cif_data(other) == {"p": other}
    assert seen[0] == manifest.parent / "fixed/AF-9.cif.gz"


def test_get_afm_cif_data_without_manifest_reads_as_is(monkeypatch):
    monkeypatch.setattr(afdb, "get_cif_data", lambda p: {"p": p})
    monkeypatch.delenv("AFM_CIF_OVERRIDES", raising=False)
    monkeypatch.setattr(afdb, "_CIF_OVERRIDES", None)
    assert afdb.get_afm_cif_data(Path("a/AF-1-model_v1.cif.gz")) == {
        "p": Path("a/AF-1-model_v1.cif.gz")
    }
    assert afdb._CIF_OVERRIDES == {}


def test_load_cif_missing_key_raises_keyerror(monkeypatch, tmp_path):
    monkeypatch.setattr(reader_io, "load_raw_data", lambda _k, _p: None)
    with pytest.raises(KeyError, match="Key '1abc' not found"):
        reader_io.load_cif("1abc", tmp_path)


def test_load_cif_splits_assembly_keys_and_merges_metadata(monkeypatch, tmp_path):
    payload = {
        "assembly_dict": {"1_2_a": {"n": 1}, "3_4_.": {"n": 2}},
        "metadata_dict": {"id": ["1abc"]},
    }
    monkeypatch.setattr(reader_io, "load_raw_data", lambda _k, _p: b"raw")
    monkeypatch.setattr(reader_io, "load_bytes", lambda _b: payload)
    monkeypatch.setattr(reader_io.CIFMol, "from_dict", staticmethod(lambda d: dict(d)))
    out = reader_io.load_cif("1abc", tmp_path)
    assert set(out) == {"1_2_a", "3_4_."}
    md = out["1_2_a"]["cifmol"]["metadata"]
    assert (md["assembly_id"], md["model_id"], md["alt_id"]) == ("1", "2", "a")
    assert md["id"] == ["1abc"]
    assert "assembly_id" not in payload["metadata_dict"]  # shared metadata not mutated


def test_load_raw_data_delegates_to_lmdb(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        reader_io, "read_lmdb_raw", lambda p, k: calls.append((p, k)) or b"x"
    )
    assert reader_io.load_raw_data("k", tmp_path) == b"x"
    assert calls == [(tmp_path, "k")]
