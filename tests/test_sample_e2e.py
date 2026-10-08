"""The sampled end-to-end check's selection, tree building and comparison."""
import json
from pathlib import Path

import lmdb
import numpy as np
import yaml
from datacooker.config import load_config

from structcooker import sample_e2e
from structcooker.instructions.transforms.codecs import to_bytes

REPO = Path(__file__).resolve().parents[1]


def test_pick_is_deterministic_and_order_free():
    items = [f"id{i}" for i in range(1000)]
    chosen = sample_e2e.pick(items, 10)
    assert chosen == sample_e2e.pick(reversed(items), 10)
    assert len(chosen) == 10
    assert set(chosen) <= set(items)


def test_fasta_entry_ids():
    assert sample_e2e.fasta_entry_id("MGYP000000008201_1_. | polypeptide(L) | Auth:1") == "MGYP000000008201"
    assert sample_e2e.fasta_entry_id("AAAA02045115.1_275-1799_1_. | polyribonucleotide") == "AAAA02045115.1_275-1799"
    assert sample_e2e.fasta_entry_id("16vp_1_. | polypeptide(L)") == "16vp"


def test_link_tree_links_files_not_directories(tmp_path):
    src = tmp_path / "src"
    (src / "hhblits").mkdir(parents=True)
    (src / "a.a3m").write_text("x")
    (src / "hhblits" / "b.a3m").write_text("y")
    assert sample_e2e.link_tree(src, tmp_path / "dst") == 2
    assert not (tmp_path / "dst" / "hhblits").is_symlink()
    assert (tmp_path / "dst" / "hhblits" / "b.a3m").resolve() == src / "hhblits" / "b.a3m"
    assert sample_e2e.link_tree(tmp_path / "absent", tmp_path / "dst2") == 0


def test_write_fasta_subset(tmp_path):
    src = tmp_path / "in.fasta"
    src.write_text(">A_1_. | p\nMK\nLV\n>B_1_. | p\nGG\n")
    kept = sample_e2e.write_fasta_subset(src, tmp_path / "out.fasta", {"A"})
    assert kept == {"A_1_. | p": "MKLV"}
    assert (tmp_path / "out.fasta").read_text() == ">A_1_. | p\nMKLV\n"


def test_sample_environment_confines_every_config(tmp_path, monkeypatch):
    """With the sample env, no release config reads or writes outside the sample root."""
    from structcooker.cli import _resolve_db

    monkeypatch.setenv("MMCIF_ROOT", "/elsewhere")
    env = sample_e2e.sample_environment(tmp_path)
    assert "MMCIF_ROOT" not in env
    monkeypatch.delenv("MMCIF_ROOT")
    for key in ("DATA_ROOT", "OUTPUT_ROOT", "STRUCTCOOKER_MAX_NODES"):
        monkeypatch.setenv(key, env[key])
    graph = yaml.safe_load((REPO / "db/MANIFEST_all.yaml").read_text())
    configs = {n: load_config(_resolve_db(n)) for n in graph}
    assert sample_e2e.escaping_paths(configs, tmp_path) == []


def test_diff_names_the_first_differing_field():
    a = {"chains": {"seq_id": np.array(["P1", "P2"]), "xyz": np.zeros((2, 3))}, "meta": {"cluster": "c1"}}
    b = {"chains": {"seq_id": np.array(["P1", "P2"]), "xyz": np.zeros((2, 3))}, "meta": {"cluster": "c2"}}
    assert sample_e2e.diff(a, a) is None
    assert sample_e2e.diff(a, b) == "/meta/cluster"
    c = {"chains": {"seq_id": np.array(["P1", "P2"]), "xyz": np.ones((2, 3))}, "meta": {"cluster": "c1"}}
    assert sample_e2e.diff(a, c) == "/chains/xyz"


def test_compare_lmdb_counts_identical_different_and_absent(tmp_path):
    def write(path, records):
        with lmdb.open(str(path), map_size=1 << 20) as env, env.begin(write=True) as txn:
            for k, v in records.items():
                txn.put(k.encode(), to_bytes(v))

    write(tmp_path / "sample", {"a": {"x": 1}, "b": {"x": 2}, "c": {"x": 3}})
    write(tmp_path / "release", {"a": {"x": 1}, "b": {"x": 9}})
    report = sample_e2e.compare_lmdb(tmp_path / "sample", tmp_path / "release", "B")
    assert (report["records"], report["identical"], report["different"], report["absent_in_release"]) == (3, 1, 1, 1)
    assert report["differing_fields"] == {"/x": 1}
    assert json.dumps(report)


def test_compare_lines(tmp_path):
    (tmp_path / "s.tsv").write_text("P1\tMK\nP9\tGG\n")
    (tmp_path / "r.tsv").write_text("P0\tAA\nP1\tMK\n")
    assert sample_e2e.compare_lines(tmp_path / "s.tsv", tmp_path / "r.tsv") == {"lines": 2, "in_release": 1}


def test_renumbered_seq_ids_are_translated(tmp_path):
    (tmp_path / "sample.tsv").write_text("P1\tMK\nP7\tGG\nR2\tACG\n")
    (tmp_path / "release.tsv").write_text("P1\tMK\nP9\tGG\nP7\tAA\nP5\tACG\nR2\tACG\n")
    table = sample_e2e.seq_id_translation(tmp_path / "sample.tsv", tmp_path / "release.tsv")
    assert table == {"P7": "P9"}
    value = {"seq_id": np.array(["P1", "P7"]), "P7": ["P7", 3]}
    out = sample_e2e.translate(value, table)
    assert list(out["seq_id"]) == ["P1", "P9"]
    assert out["P9"] == ["P9", 3]

    def write(path, records):
        with lmdb.open(str(path), map_size=1 << 20) as env, env.begin(write=True) as txn:
            for k, v in records.items():
                txn.put(k.encode(), to_bytes(v))

    write(tmp_path / "s.lmdb", {"P7": {"seq_id": "P7"}})
    write(tmp_path / "r.lmdb", {"P9": {"seq_id": "P9"}})
    assert sample_e2e.compare_lmdb(tmp_path / "s.lmdb", tmp_path / "r.lmdb", "E", table)["identical"] == 1
    (tmp_path / "lines").write_text("P7\tGG\n")
    assert sample_e2e.compare_lines(tmp_path / "lines", tmp_path / "release.tsv", table) == {"lines": 1, "in_release": 1}
