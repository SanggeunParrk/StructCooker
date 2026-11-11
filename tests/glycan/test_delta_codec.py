"""The Re-Glyco delta codec: entry base + delta records, read back as ordinary BioMol values."""

from __future__ import annotations

import lmdb
import pytest
from glycan_fakes import PDB, SPEC, make_refs, same, value_for

from structcooker.build import build_verify, schemas
from structcooker.instructions.transforms.codecs import from_bytes, to_bytes
from structcooker.instructions.transforms.delta_codec import (
    DeltaReader,
    canonical_of,
    encode_entry,
    is_delta,
    parse_delta,
    read_record,
)


@pytest.fixture(scope="module")
def encoded():
    refs = make_refs(SPEC)
    return refs, encode_entry(refs, lambda k: value_for(k))


def test_layout_and_sizes(encoded):
    refs, values = encoded
    assert set(values) == {r.key for r in refs}
    first = min(values)
    assert not is_delta(values[first])
    assert same(
        from_bytes(values[first]), value_for(first)
    )  # an entry base is an ordinary to_bytes value
    deltas = {k: v for k, v in values.items() if k != first}
    assert all(is_delta(v) for v in deltas.values())
    plain = len(to_bytes(value_for(first)))
    assert (
        max(len(v) for v in deltas.values()) < plain / 8
    )  # a delta holds the glycan only
    # the first conformer of another group is a delta against the entry base; its siblings against that group base
    gb = parse_delta(values[f"{PDB}_A293_GS00007_conf000"])
    assert (gb.kind, gb.ref) == (1, first)
    sib = parse_delta(values[f"{PDB}_A293_GS00007_conf002"])
    assert (sib.kind, sib.ref) == (2, f"{PDB}_A293_GS00007_conf000")
    own = parse_delta(
        values[f"{PDB}_A293_GS00001_conf002"]
    )  # the entry base's own group: against the entry base
    assert (own.kind, own.ref) == (1, first)


def test_every_record_reads_as_the_converted_value(encoded, tmp_path):
    refs, values = encoded
    env = _db(tmp_path, values)
    reader = DeltaReader(env)
    for r in refs:
        assert same(reader.read(r.key), value_for(r.key))
        assert canonical_of(read_record(env, r.key)) == canonical_of(value_for(r.key))


def _db(tmp_path, values):
    env = lmdb.open(str(tmp_path / "db"), map_size=1 << 26)
    with env.begin(write=True) as txn:
        for k, v in values.items():
            txn.put(k.encode(), v)
    return env


def test_plain_from_bytes_refuses_a_delta_by_name(encoded):
    _, values = encoded
    delta = next(v for v in values.values() if is_delta(v))
    with pytest.raises(ValueError, match="read_record"):
        from_bytes(delta)


def test_base_cache(encoded, tmp_path):
    refs, values = encoded
    env = _db(tmp_path, values)
    keys = sorted(r.key for r in refs)
    warm = DeltaReader(env)
    for k in keys:
        warm.read(k)
    # entry base + the group bases: 3 groups, one of which is the entry base's own
    assert warm.base_decodes == 3
    cold = DeltaReader(env, cache_bytes=0)
    for k in keys:
        cold.read(k)
    assert cold.base_decodes > warm.base_decodes


def test_missing_key_and_corrupt_value(encoded, tmp_path):
    _, values = encoded
    env = _db(tmp_path, values)
    with pytest.raises(KeyError):
        read_record(env, "nope")
    k = next(k for k, v in values.items() if is_delta(v))
    bad = bytearray(values[k])
    bad[-5] ^= 0xFF
    with env.begin(write=True) as txn:
        txn.put(k.encode(), bytes(bad))
    with pytest.raises(Exception, match=r"."):
        read_record(env, k)


def test_a_wrong_delta_is_caught_when_encoding(monkeypatch):
    import structcooker.instructions.transforms.delta_codec as dc

    monkeypatch.setattr(dc._Base, "decode", lambda _self, _d: b"x")
    with pytest.raises(ValueError, match="does not decode"):
        encode_entry(make_refs({"A1_GS00001": 2}), lambda k: value_for(k, 200))


def test_canonical_bytes_ignore_the_random_array_names():
    a, b = (
        value_for("x_A1_GS00001_conf000", 100),
        value_for("x_A1_GS00001_conf000", 100),
    )
    assert canonical_of(a) == canonical_of(b)


def test_schema_k_and_l_validators():
    ok = value_for("x_A1_GS00001_conf000", 10)
    assert schemas.get("K").validate(ok) == []
    assert schemas.get("K").codec == schemas.Codec.biomol_zstd_delta
    assert schemas.get("K").validate(
        {"assembly_dict": {}, "metadata_dict": {"id": "x"}}
    )
    assert schemas.get("K").validate(
        {"assembly_dict": {"1": {}}, "metadata_dict": {"id": "x"}}
    ) == ["metadata_dict has no deposition_date"]
    rec = {
        "key": "k",
        "shard": "s",
        "offset": 0,
        "length": 1,
        "sha256": "h",
        "deposition_date": "2020-01-01",
    }
    assert schemas.get("L").validate({"records": {"a": rec}}) == []
    assert schemas.get("L").validate({"records": {"a": {"key": "k"}}})
    assert schemas.get("L").validate({"records": {}})


def test_build_verify_reads_a_delta_db_through_the_reader(
    encoded, tmp_path, monkeypatch
):
    refs, values = encoded
    path = tmp_path / "reglyco.lmdb"
    path.mkdir()
    env = lmdb.open(str(path), map_size=1 << 26)
    with env.begin(write=True) as txn:
        for k, v in values.items():
            txn.put(k.encode(), v)
    env.close()
    monkeypatch.setattr(
        build_verify, "verify_structure", lambda _value, _schema: None
    )  # the values are synthetic
    keys = [r.key for r in refs]
    assert build_verify._check_records(str(path), keys, "K", None) == len(keys)
