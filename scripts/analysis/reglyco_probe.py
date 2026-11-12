#!/usr/bin/env python3
"""How compressible is one RFD4-Proteina Re-Glyco shard, and how much of it is repeated protein context?

A shard (`blob/data/rs_*.bin`) is a run of independent zstd frames, one annotated mmCIF per record.
Without the index we walk the frames, group the records by the entry id on their `data_` line and report:

  - the shard as shipped (sum of frame sizes) and decompressed,
  - the same records re-compressed one by one at level 6 (our `to_bytes` level) and level 19,
  - the same records compressed together per entry (long window), which is what removing the repeated
    protein context between conformers could at best reach,
  - the share of atom_site rows (identical lines) that recur across the records of one entry.

Heavy: run through sbatch (`reglyco_probe.sbatch`).
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import zstandard as zstd


def frames(buf: bytes):
    """Yield (compressed_size, decompressed_bytes) for each zstd frame in buf."""
    dctx = zstd.ZstdDecompressor()
    pos = 0
    n = len(buf)
    while pos < n:
        obj = dctx.decompressobj()
        out = []
        view = memoryview(buf)[pos:]
        # feed in slices until the frame ends
        used = 0
        step = 1 << 20
        while not obj.eof and used < len(view):
            out.append(obj.decompress(view[used : used + step]))
            used += step
        left = len(obj.unused_data)
        size = min(used, len(view)) - left
        yield size, b"".join(out)
        pos += size


def entry_id(text: bytes) -> str:
    m = re.search(rb"^data_(\S+)", text[:2000], re.M)
    name = m.group(1).decode() if m else "?"
    # records are named like ".tmp_7r2x.cif": the PDB id follows the temp prefix
    m2 = re.search(r"(?:^|[._])([0-9][a-z0-9]{3})(?:[._]|$)", name.lower())
    return m2.group(1) if m2 else name


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("shard", type=Path)
    ap.add_argument("--max-records", type=int, default=4000)
    ap.add_argument("--l19-sample", type=int, default=150, help="records to also compress at level 19")
    a = ap.parse_args()

    buf = a.shard.read_bytes()
    recs, shipped = [], 0
    for size, text in frames(buf):
        shipped += size
        recs.append(text)
        if len(recs) >= a.max_records:
            break
    unc = sum(map(len, recs))
    print(f"records read: {len(recs)}  shipped {shipped / 1e6:.1f} MB  decompressed {unc / 1e6:.1f} MB  "
          f"ratio {unc / shipped:.2f}x  mean record {shipped / len(recs) / 1e3:.0f} KB compressed")

    c6 = zstd.ZstdCompressor(level=6)
    t6 = sum(len(c6.compress(r)) for r in recs)
    print(f"re-zstd level 6 one by one: {t6 / 1e6:.1f} MB ({t6 / shipped:.2f} of shipped)")

    sample = recs[:: max(1, len(recs) // a.l19_sample)][: a.l19_sample]
    c19 = zstd.ZstdCompressor(level=19)
    s_ship = sum(len(c6.compress(r)) for r in sample)
    s19 = sum(len(c19.compress(r)) for r in sample)
    print(f"level 19 vs 6 on {len(sample)} records: {s19 / s_ship:.2f}")

    groups = defaultdict(list)
    for r in recs:
        groups[entry_id(r)].append(r)
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    print(f"entries: {len(groups)}  with >1 record: {len(multi)}  mean records/entry {len(recs) / len(groups):.1f}")

    solid = indiv = rows_total = rows_dup = 0
    cl = zstd.ZstdCompressor(
        compression_params=zstd.ZstdCompressionParameters.from_level(19, window_log=27)
    )
    for k, v in list(multi.items())[:60]:
        indiv += sum(len(c6.compress(r)) for r in v)
        solid += len(cl.compress(b"".join(v)))
        seen = set()
        for i, r in enumerate(v):
            lines = [ln for ln in r.split(b"\n") if ln.startswith((b"ATOM", b"HETATM"))]
            rows_total += len(lines)
            rows_dup += sum(1 for ln in lines if ln in seen)
            seen.update(lines)
    if indiv:
        print(f"per-entry solid (level 19, 128 MB window) on {min(60, len(multi))} multi-record entries: "
              f"{solid / 1e6:.1f} MB vs {indiv / 1e6:.1f} MB one by one -> {solid / indiv:.2f}")
        print(f"identical atom lines already seen in an earlier record of the same entry: {rows_dup / max(1, rows_total):.1%}")


if __name__ == "__main__":
    sys.exit(main())
