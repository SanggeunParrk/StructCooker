#!/usr/bin/env python3
"""Prototype and measure compact, random-access layouts for the Re-Glyco records (see docs/reglyco-compact.md).

A Re-Glyco record is one annotated mmCIF (about 1-2 MB) for one (PDB entry, glycosylation site, glycan, conformer),
named ``<xx>/<pdb>/<pdb>_<site>_<glycan>_conf<NNN>.cif``. Records of one entry share the whole protein context, and two
conformers of one site differ in about 80 glycan atom lines. This script reads one shipped shard plus the Hub index and
measures, with byte-exact round trips, these layouts:

  per-record   as shipped (one zstd frame per record)                        [baseline]
  A1           entry base + every other record as a zstd delta against that base (raw-content dictionary)
  A2           entry base + group base (site, glycan) delta vs the entry base + conformer deltas vs the group base
  B            one solid zstd frame per (entry, site, glycan) group (random access at group granularity)
  A3           entry base + explicit line diff (copy / literal ops), zstd-compressed (sample of records)

Every layout is decoded again and compared with the original bytes (and with ``cif_sha256`` of the index).
Heavy: run through sbatch (``reglyco_compact.sbatch``).
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import statistics
import struct
import sys
import time
from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import zstandard as zstd

NAME = re.compile(
    r"^(?P<xx>\w\w)/(?P<pdb>\w{4})/(?P=pdb)_(?P<site>.+)_(?P<glycan>GS\d+)_conf(?P<conf>\d+)\.cif$"
)


def parse_path(path: str) -> tuple[str, str, str, int]:
    m = NAME.match(path)
    if not m:
        raise ValueError(f"unexpected record path {path!r}")
    return m["pdb"], m["site"], m["glycan"], int(m["conf"])


def window_log(n: int) -> int:
    """Smallest window that covers dictionary plus record (zstd max 27 for the default decoder)."""
    return min(27, max(10, n.bit_length() + 1))


class DeltaCodec:
    """zstd delta of a record against a base record, the base used as a raw-content dictionary.

    The decoder needs the base bytes only; a loader keeps the bases of the entries it is reading in a small cache.
    """

    def __init__(
        self, base: bytes, level: int, max_len: int, *, encode: bool = True
    ) -> None:
        self.base = base
        self.level = level
        self.wl = window_log(len(base) + max_len)
        self.dict = zstd.ZstdCompressionDict(base, dict_type=zstd.DICT_TYPE_RAWCONTENT)
        if encode:  # a reader only decompresses: no match tables to build
            params = zstd.ZstdCompressionParameters.from_level(
                level, window_log=self.wl
            )
            self.dict.precompute_compress(compression_params=params)
            self.cctx = zstd.ZstdCompressor(
                dict_data=self.dict, compression_params=params
            )
        self.dctx = zstd.ZstdDecompressor(dict_data=self.dict, max_window_size=1 << 27)

    def encode(self, data: bytes) -> bytes:
        return struct.pack("<I", len(data)) + self.cctx.compress(data)

    def decode(self, blob: bytes) -> bytes:
        (n,) = struct.unpack_from("<I", blob)
        return self.dctx.decompress(blob[4:], max_output_size=n + 1)


def plain_encode(data: bytes, level: int = 19) -> bytes:
    return struct.pack("<I", len(data)) + zstd.ZstdCompressor(level=level).compress(
        data
    )


def plain_decode(blob: bytes) -> bytes:
    (n,) = struct.unpack_from("<I", blob)
    return zstd.ZstdDecompressor().decompress(blob[4:], max_output_size=n + 1)


# ---- A3: explicit line diff ---------------------------------------------------------------------------------------


def line_diff_encode(base: bytes, data: bytes) -> bytes:
    """Ops against the base lines: ('c', start, count) copy, ('l', n) then n literal lines; zstd-compressed text."""
    a, b = base.split(b"\n"), data.split(b"\n")
    ids: dict[bytes, int] = {}
    ai = [ids.setdefault(x, len(ids)) for x in a]
    bi = [ids.setdefault(x, len(ids)) for x in b]
    ops: list[bytes] = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(
        None, ai, bi, autojunk=False
    ).get_opcodes():
        if tag == "equal":
            ops.append(b"c %d %d" % (i1, i2 - i1))
        elif tag in ("replace", "insert"):
            ops.append(b"l %d" % (j2 - j1))
            ops.extend(b[j1:j2])
    return zstd.ZstdCompressor(level=19).compress(b"\n".join(ops))


def line_diff_decode(base: bytes, blob: bytes) -> bytes:
    a = base.split(b"\n")
    ops = zstd.ZstdDecompressor().decompress(blob, max_output_size=1 << 31).split(b"\n")
    out: list[bytes] = []
    i = 0
    while i < len(ops):
        op = ops[i].split()
        if op[0] == b"c":
            out.extend(a[int(op[1]) : int(op[1]) + int(op[2])])
            i += 1
        else:
            n = int(op[1])
            out.extend(ops[i + 1 : i + 1 + n])
            i += 1 + n
    return b"\n".join(out)


# ---- shard reading ------------------------------------------------------------------------------------------------


def load_records(shard: Path, index: Path) -> list[dict]:
    t = pq.read_table(index)
    t = t.filter(pc.equal(t["shard"], pa.scalar(shard.name)))
    rows = t.to_pylist()
    buf = shard.read_bytes()
    dctx = zstd.ZstdDecompressor()
    out = []
    for r in rows:
        data = dctx.decompress(
            buf[r["offset"] : r["offset"] + r["length"]],
            max_output_size=r["uncompressed_length"] + 1,
        )
        assert hashlib.sha256(data).hexdigest() == r["cif_sha256"], r["path"]
        pdb, site, glycan, conf = parse_path(r["path"])
        out.append(
            {
                "path": r["path"],
                "pdb": pdb,
                "site": site,
                "glycan": glycan,
                "conf": conf,
                "data": data,
                "shipped": r["length"],
                "sha": r["cif_sha256"],
            }
        )
    out.sort(key=lambda x: x["path"])
    return out


def ok(rec: dict, data: bytes) -> bool:
    return data == rec["data"] and hashlib.sha256(data).hexdigest() == rec["sha"]


def mb(x: float) -> str:
    return f"{x / 1e6:8.2f} MB"


# ---- layouts ------------------------------------------------------------------------------------------------------


def eval_a1(recs: list[dict], level: int) -> dict:
    by_entry = defaultdict(list)
    for r in recs:
        by_entry[r["pdb"]].append(r)
    base_bytes = delta_bytes = 0
    enc_t = dec_cold = dec_warm = 0.0
    n_dec = 0
    per_entry = {}
    store = {}
    for pdb, rs in by_entry.items():
        base = rs[0]
        t0 = time.time()
        b_blob = plain_encode(base["data"])
        enc_t += time.time() - t0
        codec = DeltaCodec(base["data"], level, max(len(x["data"]) for x in rs))
        blobs = {}
        t0 = time.time()
        for r in rs[1:]:
            blobs[r["path"]] = codec.encode(r["data"])
        enc_t += time.time() - t0
        store[pdb] = (b_blob, blobs)
        bb, db = len(b_blob), sum(len(v) for v in blobs.values())
        base_bytes += bb
        delta_bytes += db
        per_entry[pdb] = {
            "records": len(rs),
            "base": bb,
            "delta_total": db,
            "delta_mean": db / max(1, len(rs) - 1),
            "base_shipped": base["shipped"],
        }
    # decode everything: byte-exact check, plus timings
    exact = True
    for pdb, rs in by_entry.items():
        b_blob, blobs = store[pdb]
        t0 = time.time()
        base = plain_decode(b_blob)
        codec = DeltaCodec(base, level, max(len(x["data"]) for x in rs), encode=False)
        cold = time.time() - t0
        exact &= ok(rs[0], base)
        for r in rs[1:]:
            t1 = time.time()
            data = codec.decode(blobs[r["path"]])
            dec_warm += time.time() - t1
            n_dec += 1
            exact &= ok(r, data)
        dec_cold += cold
    return {
        "layout": "A1",
        "level": level,
        "total": base_bytes + delta_bytes,
        "base_bytes": base_bytes,
        "delta_bytes": delta_bytes,
        "n_delta": len(recs) - len(by_entry),
        "n_entries": len(by_entry),
        "exact": exact,
        "encode_s": enc_t,
        "decode_warm_ms": 1e3 * dec_warm / max(1, n_dec),
        "decode_cold_ms": 1e3 * dec_cold / len(by_entry),
        "per_entry": per_entry,
    }


def eval_a2(recs: list[dict], level: int) -> dict:
    by_entry = defaultdict(lambda: defaultdict(list))
    for r in recs:
        by_entry[r["pdb"]][(r["site"], r["glycan"])].append(r)
    base_b = gbase_b = conf_b = 0
    n_groups = n_conf = 0
    exact = True
    dec_warm = 0.0
    n_dec = 0
    enc_t = 0.0
    per_entry = {}
    for pdb, groups in by_entry.items():
        keys = sorted(groups)
        allr = [r for k in keys for r in groups[k]]
        ent = groups[keys[0]][0]
        t0 = time.time()
        b_blob = plain_encode(ent["data"])
        ecodec = DeltaCodec(ent["data"], level, max(len(x["data"]) for x in allr))
        base_b += len(b_blob)
        e_gb = e_cb = e_nc = 0
        stored = {}
        for k in keys:
            rs = groups[k]
            g0 = rs[0]
            if g0 is ent:
                g_blob = None  # the entry base is also this group's base
            else:
                g_blob = ecodec.encode(g0["data"])
                gbase_b += len(g_blob)
                e_gb += len(g_blob)
            n_groups += 1
            gcodec = DeltaCodec(g0["data"], level, max(len(x["data"]) for x in rs))
            c_blobs = [gcodec.encode(r["data"]) for r in rs[1:]]
            conf_b += sum(len(x) for x in c_blobs)
            n_conf += len(c_blobs)
            e_cb += sum(len(x) for x in c_blobs)
            e_nc += len(c_blobs)
            stored[k] = (g_blob, c_blobs)
        enc_t += time.time() - t0
        per_entry[pdb] = {
            "groups": len(keys),
            "records": len(allr),
            "base": len(b_blob),
            "group_base": e_gb,
            "conf": e_cb,
            "n_conf": e_nc,
            "base_shipped": ent["shipped"],
        }
        # decode check
        base = plain_decode(b_blob)
        exact &= ok(ent, base)
        ecodec2 = DeltaCodec(
            base, level, max(len(x["data"]) for x in allr), encode=False
        )
        for k in keys:
            g_blob, c_blobs = stored[k]
            rs = groups[k]
            g0 = base if g_blob is None else ecodec2.decode(g_blob)
            exact &= ok(rs[0], g0)
            gc = DeltaCodec(g0, level, max(len(x["data"]) for x in rs), encode=False)
            for r, blob in zip(rs[1:], c_blobs, strict=True):
                t1 = time.time()
                data = gc.decode(blob)
                dec_warm += time.time() - t1
                n_dec += 1
                exact &= ok(r, data)
    return {
        "layout": "A2",
        "level": level,
        "total": base_b + gbase_b + conf_b,
        "base_bytes": base_b,
        "group_base_bytes": gbase_b,
        "conf_bytes": conf_b,
        "n_groups": n_groups,
        "n_conf": n_conf,
        "n_entries": len(by_entry),
        "exact": exact,
        "encode_s": enc_t,
        "decode_warm_ms": 1e3 * dec_warm / max(1, n_dec),
        "per_entry": per_entry,
    }


def eval_b(recs: list[dict], level: int) -> dict:
    groups = defaultdict(list)
    for r in recs:
        groups[(r["pdb"], r["site"], r["glycan"])].append(r)
    params = zstd.ZstdCompressionParameters.from_level(
        level, window_log=27, enable_ldm=1
    )
    cctx = zstd.ZstdCompressor(compression_params=params)
    dctx = zstd.ZstdDecompressor(max_window_size=1 << 27)
    total = 0
    exact = True
    enc_t = dec_t = 0.0
    sizes = []
    for k, rs in groups.items():
        raw = b"".join(r["data"] for r in rs)
        t0 = time.time()
        blob = cctx.compress(raw)
        enc_t += time.time() - t0
        total += len(blob) + 4 * len(rs)  # plus the per-record length table
        t0 = time.time()
        back = dctx.decompress(blob, max_output_size=len(raw) + 1)
        dec_t += time.time() - t0
        pos = 0
        for r in rs:
            exact &= ok(r, back[pos : pos + len(r["data"])])
            pos += len(r["data"])
        sizes.append((len(rs), len(raw)))
    return {
        "layout": "B",
        "level": level,
        "total": total,
        "n_groups": len(groups),
        "exact": exact,
        "encode_s": enc_t,
        "decode_group_ms": 1e3 * dec_t / len(groups),
        "mean_group_records": statistics.mean(s[0] for s in sizes),
        "mean_group_raw_mb": statistics.mean(s[1] for s in sizes) / 1e6,
    }


def eval_a3(recs: list[dict], every: int) -> dict:
    by_entry = defaultdict(list)
    for r in recs:
        by_entry[r["pdb"]].append(r)
    tot = n = 0
    exact = True
    t_enc = t_dec = 0.0
    changed = []
    for pdb, rs in by_entry.items():
        base = rs[0]["data"]
        for r in rs[1::every]:
            t0 = time.time()
            blob = line_diff_encode(base, r["data"])
            t_enc += time.time() - t0
            t0 = time.time()
            back = line_diff_decode(base, blob)
            t_dec += time.time() - t0
            exact &= ok(r, back)
            tot += len(blob)
            n += 1
            changed.append(len(blob))
    return {
        "layout": "A3",
        "n_sample": n,
        "delta_mean": tot / max(1, n),
        "exact": exact,
        "encode_ms": 1e3 * t_enc / max(1, n),
        "decode_ms": 1e3 * t_dec / max(1, n),
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("shard", type=Path)
    ap.add_argument("index", type=Path)
    ap.add_argument("--out", type=Path, required=True, help="JSON results")
    ap.add_argument(
        "--levels", default="3,9,19", help="zstd levels tried for the delta layouts"
    )
    ap.add_argument("--a3-every", type=int, default=25)
    a = ap.parse_args()

    recs = load_records(a.shard, a.index)
    shipped = sum(r["shipped"] for r in recs)
    raw = sum(len(r["data"]) for r in recs)
    ents = {r["pdb"] for r in recs}
    grp = {(r["pdb"], r["site"], r["glycan"]) for r in recs}
    res: dict = {
        "shard": a.shard.name,
        "records": len(recs),
        "entries": len(ents),
        "groups": len(grp),
        "shipped": shipped,
        "raw": raw,
    }
    print(
        f"{a.shard.name}: {len(recs)} records, {len(ents)} entries, {len(grp)} (site, glycan) groups; "
        f"shipped {mb(shipped)}, raw {mb(raw)}",
        flush=True,
    )
    res["layouts"] = []
    for lv in map(int, a.levels.split(",")):
        r1 = eval_a1(recs, lv)
        print(
            f"A1 level {lv:2d}: {mb(r1['total'])} (bases {mb(r1['base_bytes'])}, deltas {mb(r1['delta_bytes'])}, "
            f"{r1['delta_bytes'] / r1['n_delta'] / 1e3:.1f} KB/record)  exact={r1['exact']}  enc {r1['encode_s']:.0f}s  "
            f"dec {r1['decode_warm_ms']:.1f} ms warm / {r1['decode_cold_ms']:.0f} ms cold",
            flush=True,
        )
        res["layouts"].append(r1)
        r2 = eval_a2(recs, lv)
        print(
            f"A2 level {lv:2d}: {mb(r2['total'])} (bases {mb(r2['base_bytes'])}, group bases {mb(r2['group_base_bytes'])}, "
            f"conformers {mb(r2['conf_bytes'])}, {r2['conf_bytes'] / max(1, r2['n_conf']) / 1e3:.1f} KB/conf)  "
            f"exact={r2['exact']}  enc {r2['encode_s']:.0f}s  dec {r2['decode_warm_ms']:.1f} ms",
            flush=True,
        )
        res["layouts"].append(r2)
    for lv in (12,):
        rb = eval_b(recs, lv)
        print(
            f"B  level {lv:2d}: {mb(rb['total'])}  groups {rb['n_groups']} (mean {rb['mean_group_records']:.1f} records, "
            f"{rb['mean_group_raw_mb']:.1f} MB raw)  exact={rb['exact']}  enc {rb['encode_s']:.0f}s  "
            f"decode a whole group {rb['decode_group_ms']:.0f} ms",
            flush=True,
        )
        res["layouts"].append(rb)
    r3 = eval_a3(recs, a.a3_every)
    print(
        f"A3 line diff ({r3['n_sample']} records): {r3['delta_mean'] / 1e3:.1f} KB/record  exact={r3['exact']}  "
        f"enc {r3['encode_ms']:.0f} ms  dec {r3['decode_ms']:.0f} ms",
        flush=True,
    )
    res["layouts"].append(r3)
    a.out.write_text(json.dumps(res, indent=1, default=str))
    print("done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
