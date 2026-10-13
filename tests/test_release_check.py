"""release_check: a shard collection must hold its shards; counts must agree."""
import json

import lmdb
from datacooker.lmdb.sharded import publish_shards

from structcooker import release_check
from structcooker.instructions.transforms.codecs import to_bytes


def _write(path, n):
    with lmdb.open(str(path), map_size=1 << 22) as env, env.begin(write=True) as txn:
        for i in range(n):
            txn.put(f"k{i}".encode(), to_bytes({"chem_comp_dict": None}))


def _sidecars(target, n):
    (target.parent / (target.name + ".meta.json")).write_text(json.dumps({"entries": n, "schema": "G", "total_bytes": n, "max_bytes": 1}))
    (target.parent / (target.name + ".index.tsv")).write_text("key\tbytes\n" + "".join(f"k{i}\t1\n" for i in range(n)))


def test_single_lmdb_passes_and_counts_are_checked(tmp_path):
    db = tmp_path / "ccd.lmdb"
    _write(db, 3)
    _sidecars(db, 3)
    assert release_check.check_lmdb("ccd/ccd", str(db), "G", None)["problems"] == []
    _sidecars(db, 4)
    problems = release_check.check_lmdb("ccd/ccd", str(db), "G", None)["problems"]
    assert any("meta.json entries" in p for p in problems)


def test_shards_outside_the_collection_are_reported(tmp_path):
    (tmp_path / "build").mkdir()
    _write(tmp_path / "build" / "s0", 2)
    db = tmp_path / "coll.lmdb"
    publish_shards([tmp_path / "build" / "s0"], db)
    _sidecars(db, 2)
    problems = release_check.check_lmdb("x", str(db), "G", None)["problems"]
    assert any("outside its directory" in p for p in problems)


def test_file_outputs(tmp_path):
    (tmp_path / "a.tsv").write_text("x\n")
    (tmp_path / "empty").mkdir()
    r = release_check.check_files("m", [str(tmp_path / "a.tsv"), str(tmp_path / "empty")])
    assert r["problems"] == [f"empty directory {tmp_path / 'empty'}"]
