import json

import lmdb
import pytest
from datacooker.lmdb.core import (
    count_lmdb_entries,
    extract_lmdb_keys,
    merge_lmdb_shards,
)
from datacooker.lmdb.index import scan_sizes, write_index
from datacooker.lmdb.sharded import open_env, publish_shards, should_merge

from structcooker.build_state import Fingerprints


def make_db(path, items):
    with lmdb.open(str(path), map_size=1024**2) as env, env.begin(write=True) as txn:
        for key, value in items:
            txn.put(key, value)
    return path


@pytest.mark.parametrize(("count", "merged"), [(4_999_999, True), (5_000_000, True), (5_000_001, False)])
def test_five_million_boundary(count, merged):
    assert should_merge(count) is merged


def test_collection_reads_indexes_and_fingerprints(tmp_path):
    a = make_db(tmp_path / "a", [(b"a", b"123"), (b"c", b"45")])
    b = make_db(tmp_path / "b", [(b"b", b"6789")])
    output = tmp_path / "logical"
    assert publish_shards([a, b], output) == 3
    assert not (output / "data.mdb").exists()
    assert count_lmdb_entries(output) == 3
    assert extract_lmdb_keys(output) == ["a", "b", "c"]
    with open_env(output, readonly=True, lock=False) as env, env.begin(buffers=True) as txn:
        assert bytes(txn.get(b"c")) == b"45"
        assert txn.get(b"missing") is None
        assert [bytes(v) for v in txn.cursor().iternext(keys=False)] == [b"123", b"6789", b"45"]
    assert {s.key: s.value_bytes for s in scan_sizes(output)} == {"a": 3, "b": 4, "c": 2}
    assert write_index(output).entries == 3
    from datacooker.lmdb.core import rekey_lmdb

    assert rekey_lmdb(output, tmp_path / "copied", map_size=1024**2).written == 3
    previous = Fingerprints().path(output)
    with lmdb.open(str(a)) as env, env.begin(write=True) as txn:
        txn.put(b"a", b"tampered")
    assert Fingerprints().path(output) != previous
    with pytest.raises(ValueError, match="immutable"):
        open_env(output)


def test_duplicate_publication_preserves_manifest(tmp_path):
    a = make_db(tmp_path / "a", [(b"a", b"123")])
    b = make_db(tmp_path / "b", [(b"a", b"456")])
    output = tmp_path / "logical"
    publish_shards([a], output)
    previous = (output / "shards.json").read_bytes()
    with pytest.raises(ValueError, match="Duplicate"):
        publish_shards([a, b], output, overwrite=True)
    assert (output / "shards.json").read_bytes() == previous
    with pytest.raises(ValueError, match="unique"):
        publish_shards([a, a], tmp_path / "bad")
    (a / "data.mdb").rename(a / "missing")
    with pytest.raises(FileNotFoundError, match="Missing"):
        open_env(output, readonly=True)


def test_merge_policy_uses_actual_outputs(tmp_path, monkeypatch):
    from datacooker.lmdb import sharded

    a = make_db(tmp_path / "a", [(b"a", b"123")])
    b = make_db(tmp_path / "b", [(b"b", b"456")])
    monkeypatch.setattr(sharded, "MERGE_LIMIT", 1)
    output = tmp_path / "large"
    assert merge_lmdb_shards([a, b], output).written == 2
    assert json.loads((output / "shards.json").read_text())["entries"] == 2
    monkeypatch.setattr(sharded, "MERGE_LIMIT", 5_000_000)
    small = tmp_path / "small"
    assert merge_lmdb_shards([a, b], small, map_size=1024**2).written == 2
    assert (small / "data.mdb").exists()
