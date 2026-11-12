#!/usr/bin/env python3
"""Convert a random sample of each RFD4-Proteina component with the rfd4 reader + the pdb/cif recipe; report failures.

Heavy-ish (parses CIFs, loads the CCD): run through sbatch (``probe.sbatch``). Used to find what the adapter lacks
before a full build; prints per-component ok/fail counts, the first distinct errors and the mean time per record.
"""

from __future__ import annotations

import argparse
import random
import time
import traceback
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq
from datacooker import parse_dict

from structcooker.instructions.readers.cif import dot_transform
from structcooker.instructions.readers.rfd4 import prepare_cif_inputs, read_blob
from structcooker.instructions.transforms.rfd4_worklist import (
    COMPONENTS,
    component_dir,
    default_root,
    make_keys,
    resolve_shards,
)
from structcooker.workflows.ingest.cif import RECIPE, TARGETS


def sample_rows(name: str, root: str, n: int, seed: int) -> tuple[Path, list[dict]]:
    comp = COMPONENTS[name]
    base = component_dir(root, name)
    f = pq.ParquetFile(base / "blob/index.parquet")
    rng = random.Random(seed)
    groups = rng.sample(
        range(f.metadata.num_row_groups), min(f.metadata.num_row_groups, 8)
    )
    rows: list[dict] = []
    for g in groups:
        t = f.read_row_group(g)
        sha = next((c for c in comp.sha_cols if c in t.column_names), None)
        keys = make_keys(t.column(comp.key_col), comp.key_kind)
        cols = [t.column(c).to_pylist() for c in ("shard", "offset", "length")]
        shas = t.column(sha).to_pylist() if sha else [None] * len(keys)
        for k, s, o, ln, h in zip(keys, *cols, shas, strict=True):
            rows.append(
                {
                    "key": k,
                    "shard": s,
                    "offset": o,
                    "length": ln,
                    "sha256": h,
                    "date": comp.default_date,
                }
            )
    rows = rng.sample(rows, min(n, len(rows)))
    rel = resolve_shards(base / "blob", [r["shard"] for r in rows])
    for r in rows:
        r["shard"] = rel[r["shard"]]
    return base / "blob", rows


def worklist_rows(name: str, root: str, keys: list[str]) -> tuple[Path, list[dict]]:
    """The worklist rows of ``keys`` (a failed key of a build, say)."""
    import os

    import lmdb

    from structcooker.instructions.transforms.codecs import from_bytes

    biomol = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    env = lmdb.open(
        f"{biomol}/materials/intermediate/rfd4/{name}/worklist.lmdb",
        readonly=True,
        lock=False,
    )
    with env.begin() as txn:
        return component_dir(root, name) / "blob", [
            from_bytes(bytes(txn.get(k.encode()))) for k in keys
        ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--only", nargs="*", default=list(COMPONENTS))
    ap.add_argument(
        "--keys",
        nargs="*",
        help="probe these record keys (read from the worklist LMDB of the one --only set)",
    )
    a = ap.parse_args()
    root = default_root()
    for name in a.only:
        if a.keys:
            blob, rows = worklist_rows(name, root, a.keys)
        else:
            blob, rows = sample_rows(name, root, a.n, a.seed)
        errs: Counter[str] = Counter()
        first: dict[str, str] = {}
        ok = 0
        t0 = time.time()
        for r in rows:
            try:
                text = read_blob(blob, r)
                d = prepare_cif_inputs(text, r["key"], r["date"])
                out = parse_dict(
                    RECIPE, d, key_transform=dot_transform, targets=TARGETS
                )
                assert out["assembly_dict"], "empty assembly_dict"
                ok += 1
            except Exception as e:  # noqa: BLE001
                root_e = e
                while root_e.__cause__ is not None:
                    root_e = root_e.__cause__
                sig = f"{type(root_e).__name__}: {str(root_e)[:160]}"
                errs[sig] += 1
                first.setdefault(sig, f"{r['key']}\n{traceback.format_exc()[-900:]}")
        dt = (time.time() - t0) / max(1, len(rows))
        print(f"== {name}: {ok}/{len(rows)} ok, {dt * 1000:.0f} ms/record", flush=True)
        for sig, c in errs.most_common(5):
            print(f"  {c}x {sig}\n    first: {first[sig][:1100]}", flush=True)


if __name__ == "__main__":
    main()
