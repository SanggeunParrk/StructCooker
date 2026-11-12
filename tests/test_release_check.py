"""release_check: a shard collection must hold its shards; counts must agree."""

import json

import lmdb
from datacooker.lmdb.sharded import publish_shards

from structcooker.instructions.transforms.codecs import to_bytes
from structcooker.release import check as release_check


def _write(path, n):
    with lmdb.open(str(path), map_size=1 << 22) as env, env.begin(write=True) as txn:
        for i in range(n):
            txn.put(f"k{i}".encode(), to_bytes({"chem_comp_dict": None}))


def _sidecars(target, n):
    (target.parent / (target.name + ".meta.json")).write_text(
        json.dumps({"entries": n, "schema": "G", "total_bytes": n, "max_bytes": 1})
    )
    (target.parent / (target.name + ".index.tsv")).write_text(
        "key\tbytes\n" + "".join(f"k{i}\t1\n" for i in range(n))
    )


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
    r = release_check.check_files(
        "m", [str(tmp_path / "a.tsv"), str(tmp_path / "empty")]
    )
    assert r["problems"] == [f"empty directory {tmp_path / 'empty'}"]


def test_a_delta_coded_db_must_say_so_in_meta_json(tmp_path):
    import json

    import lmdb

    from structcooker.release.check import check_lmdb

    path = tmp_path / "reglyco.lmdb"
    path.mkdir()
    with lmdb.open(str(path), map_size=1 << 24) as env, env.begin(write=True) as txn:
        txn.put(b"2pmv_A1_GS00001_conf000", b"x")
    (tmp_path / "reglyco.lmdb.index.tsv").write_text(
        "key\tsize\n2pmv_A1_GS00001_conf000\t1\n"
    )
    meta = {
        "schema": "K",
        "codec": "biomol_zstd",
        "entries": 1,
        "total_bytes": 1,
        "max_bytes": 1,
    }
    meta_path = tmp_path / "reglyco.lmdb.meta.json"
    meta_path.write_text(json.dumps(meta))
    report = check_lmdb("glycan/reglyco", str(path), "K", None)
    assert any(
        "codec biomol_zstd != biomol_zstd_delta" in p for p in report["problems"]
    )
    meta["codec"] = "biomol_zstd_delta"
    meta_path.write_text(json.dumps(meta))
    report = check_lmdb("glycan/reglyco", str(path), "K", None)
    assert not any("codec" in p for p in report["problems"])
