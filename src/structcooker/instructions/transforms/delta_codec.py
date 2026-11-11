"""The delta codec of the Re-Glyco DB (``glycan/reglyco``, schema K) and the one function a loader calls to read it.

Re-Glyco holds hundreds of conformers of the same protein: two records of one entry differ in a few dozen glycan
atoms. The DB stores each entry as

* an **entry base**: the first record of the entry (by key) as an ordinary ``to_bytes`` value;
* per (site, glycan) **group**, the lowest conformer as a delta against the entry base (the group base);
* every other conformer as a delta against its group base (conformers of the entry base's own group: against the
  entry base).

A delta is zstd of the record's *canonical serialised bytes* with the base's canonical bytes as a raw-content
dictionary. The canonical bytes are what ``codecs.to_raw`` writes with the array names renamed ``a0, a1, ...`` in
header order (``to_bytes`` names them with random UUIDs, which would make every record differ from its base).

Every key is a real record of the DB and decodes to an ordinary BioMol dict equal to the converted one::

    from datacooker.lmdb.sharded import open_env
    from structcooker.instructions.transforms.delta_codec import DeltaReader

    with open_env(path, readonly=True, lock=False) as env:
        reader = DeltaReader(env)                      # keeps the entry bases it has decoded
        value = reader.read("2pmv_A293_GS00178_conf012")   # {"assembly_dict": ..., "metadata_dict": ...}

``codecs.from_bytes`` works on the entry bases only and refuses a delta value by name.
"""

from __future__ import annotations

import json
import struct
from collections import OrderedDict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import zstandard as zstd

from structcooker.instructions.transforms.codecs import (
    DELTA_MAGIC,
    from_bytes,
    from_raw,
    to_raw,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

KIND_ENTRY = 1  # the reference is the entry base (a plain frame)
KIND_GROUP = 2  # the reference is a group base (itself a delta against the entry base)
LEVEL = 19  # BioMol arrays shift when the glycan changes; levels below ~15 lose most of the gain
MAX_WINDOW = 27  # log2 of the largest window a default zstd decoder accepts
_HEAD = struct.Struct(
    "<BBI"
)  # kind, len(ref key), length of the canonical bytes of this record


def is_delta(value: bytes) -> bool:
    """Tell whether a stored value is a delta record (not a plain ``to_bytes`` frame)."""
    return bytes(value[: len(DELTA_MAGIC)]) == DELTA_MAGIC


# ---------------------------------------------------------------------------------------------- canonical bytes
def _remap(t: object, names: dict[str, str]) -> object:
    if isinstance(t, dict):
        return {k: _remap(v, names) for k, v in t.items()}
    if isinstance(t, str) and t in names:
        return names[t]
    return t


def canonical(raw: bytes) -> bytes:
    """Rename the arrays of a serialised value ``a0, a1, ...`` in header order (the same dict when decoded)."""
    hlen = int.from_bytes(raw[:8], "little")
    header = json.loads(raw[8 : 8 + hlen])
    names = {key: f"a{i}" for i, key in enumerate(header["arrays"])}
    head = json.dumps(
        {
            "template": _remap(header["template"], names),
            "arrays": {names[k]: n for k, n in header["arrays"].items()},
        }
    ).encode()
    return len(head).to_bytes(8, "little") + head + raw[8 + hlen :]


def canonical_of(value: dict) -> bytes:
    """Canonical bytes of a BioMol value."""
    return canonical(to_raw(value))


def _window_log(n: int) -> int:
    return min(MAX_WINDOW, max(10, n.bit_length() + 1))


def _pack(kind: int, ref: str, n: int, frame: bytes) -> bytes:
    rb = ref.encode()
    return DELTA_MAGIC + _HEAD.pack(kind, len(rb), n) + rb + frame


@dataclass(frozen=True)
class Delta:
    """A parsed delta value."""

    kind: int
    ref: str
    n: int  # length of the record's canonical bytes
    frame: bytes


def parse_delta(value: bytes) -> Delta:
    """Split a stored delta value into its parts."""
    if not is_delta(value):
        msg = "not a Re-Glyco delta value"
        raise ValueError(msg)
    off = len(DELTA_MAGIC)
    kind, rlen, n = _HEAD.unpack_from(value, off)
    off += _HEAD.size
    ref = bytes(value[off : off + rlen]).decode()
    return Delta(kind, ref, n, bytes(value[off + rlen :]))


class _Base:
    """A base record's canonical bytes as a zstd raw-content dictionary, for encoding and decoding."""

    def __init__(self, canon: bytes, *, level: int | None = None) -> None:
        self.canon = canon
        # A record is about as long as its base; the window must cover dictionary plus record (the decoder reads the
        # window from the frame, so it needs no such sizing)
        window_log = _window_log(2 * len(canon) + (1 << 20))
        self.dict = zstd.ZstdCompressionDict(canon, dict_type=zstd.DICT_TYPE_RAWCONTENT)
        self.cctx = None
        if level is not None:  # a reader only decompresses: no match tables to build
            params = zstd.ZstdCompressionParameters.from_level(
                level, window_log=window_log
            )
            self.dict.precompute_compress(compression_params=params)
            self.cctx = zstd.ZstdCompressor(
                dict_data=self.dict, compression_params=params
            )
        self.dctx = zstd.ZstdDecompressor(
            dict_data=self.dict, max_window_size=1 << MAX_WINDOW
        )

    def encode(self, kind: int, ref: str, canon: bytes) -> bytes:
        assert self.cctx is not None  # noqa: S101  built with a level
        return _pack(kind, ref, len(canon), self.cctx.compress(canon))

    def decode(self, d: Delta) -> bytes:
        return self.dctx.decompress(d.frame, max_output_size=d.n + 1)


# ---------------------------------------------------------------------------------------------------- encoding
@dataclass(frozen=True)
class RecordRef:
    """One record of an entry: its key, its (site, glycan) group and its conformer number."""

    key: str
    group: str
    conf: int


def encode_entry(
    refs: Iterable[RecordRef],
    convert: Callable[[str], dict[str, Any]],
    *,
    level: int = LEVEL,
    verify: bool = True,
) -> dict[str, bytes]:
    """Encode one entry: the stored value of every record as ``{key: value}``.

    ``convert(key)`` returns the BioMol value of a record (called once per record, in an order that needs only the
    entry base and the current group base in memory). With ``verify`` each delta is decoded again and compared with
    the canonical bytes it encodes, so a wrong value never reaches the DB.
    """
    ordered = sorted(refs, key=lambda r: r.key)
    if not ordered:
        return {}
    first = ordered[0]
    raw0 = to_raw(convert(first.key))
    out: dict[str, bytes] = {first.key: zstd.ZstdCompressor(level=level).compress(raw0)}
    base = _Base(canonical(raw0), level=level)
    del raw0
    groups: dict[str, list[RecordRef]] = {}
    for r in ordered[1:]:
        groups.setdefault(r.group, []).append(r)
    for gid, members in groups.items():
        members.sort(key=lambda r: (r.conf, r.key))
        if (
            gid == first.group
        ):  # the entry base's own group: its conformers are deltas against the entry base
            gbase, gref, kind, rest = base, first.key, KIND_ENTRY, members
        else:
            head = members[0]
            canon = canonical_of(convert(head.key))
            value = base.encode(KIND_ENTRY, first.key, canon)
            _check(base, value, canon, head.key, verify=verify)
            out[head.key] = value
            gbase = _Base(canon, level=level)
            gref, kind, rest = head.key, KIND_GROUP, members[1:]
        for r in rest:
            canon = canonical_of(convert(r.key))
            value = gbase.encode(kind, gref, canon)
            _check(gbase, value, canon, r.key, verify=verify)
            out[r.key] = value
    return out


def _check(base: _Base, value: bytes, canon: bytes, key: str, *, verify: bool) -> None:
    if verify and base.decode(parse_delta(value)) != canon:
        msg = f"delta of {key} does not decode to the record it encodes"
        raise ValueError(msg)


# ---------------------------------------------------------------------------------------------------- reading
class DeltaReader:
    """Read the records of a Re-Glyco DB as ordinary BioMol dicts.

    Keeps the dictionaries of the bases it has used (an entry base is decoded once for all of its records; reading
    the records of one entry in a row is the cheap case). ``cache_bytes`` bounds that cache (the base records are
    1-14 MB each).
    """

    def __init__(self, env: Any, cache_bytes: int = 512 << 20) -> None:
        self.env = env
        self.cache_bytes = cache_bytes
        self._bases: OrderedDict[str, _Base] = OrderedDict()
        self._size = 0
        self.base_decodes = (
            0  # how many bases were decoded (a cache miss), for tests and tuning
        )

    def _get(self, txn: Any, key: str) -> bytes:
        value = txn.get(key.encode())
        if value is None:
            raise KeyError(key)
        return bytes(value)

    def _base(self, txn: Any, key: str) -> _Base:
        hit = self._bases.get(key)
        if hit is not None:
            self._bases.move_to_end(key)
            return hit
        value = self._get(txn, key)
        if is_delta(value):
            d = parse_delta(value)
            canon = self._base(txn, d.ref).decode(d)
        else:
            canon = canonical(
                zstd.ZstdDecompressor().decompress(value, max_output_size=1 << 31)
            )
        self.base_decodes += 1
        base = _Base(canon)
        self._bases[key] = base
        self._size += len(canon)
        while self._size > self.cache_bytes and len(self._bases) > 1:
            _, old = self._bases.popitem(last=False)
            self._size -= len(old.canon)
        return base

    def read(self, key: str) -> dict[str, Any]:
        """Return the BioMol value stored under ``key``."""
        with self.env.begin() as txn:
            value = self._get(txn, key)
            if not is_delta(value):
                return from_bytes(value)
            d = parse_delta(value)
            return from_raw(self._base(txn, d.ref).decode(d))


def read_record(env: Any, key: str) -> dict[str, Any]:
    """One record, decoded without keeping any cache (use ``DeltaReader`` to read many)."""
    return DeltaReader(env, cache_bytes=0).read(key)
