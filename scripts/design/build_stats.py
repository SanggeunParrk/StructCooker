"""Per-item statistics + backbone tar archives of an LMDB of monomers (sbatch only).

usage: build_stats.py --db LMDB_OR_SHARD_COLLECTION --out DIR [--workers 32] [--limit N] [--seed 0]
Each worker takes a share of the keys and writes DIR/part_<w>.tsv (id, length, mean pLDDT, helix, strand, coil, max_coil_run)
and DIR/part_<w>.tar (one backbone PDB per item, input of ``foldseek createdb``).
"""

import argparse
import io
import json
import multiprocessing as mp
import os
import random
import re
import tarfile
import time
import traceback
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import lmdb  # noqa: E402
import numpy as np  # noqa: E402

from structcooker.instructions.design.export import (  # noqa: E402
    CHAIN_LETTERS,
    backbone_pdb_text,
    n_protein_chains,
)
from structcooker.instructions.design.ss import secondary_structure  # noqa: E402
from structcooker.instructions.mols import CIFMolAttached  # noqa: E402
from structcooker.instructions.transforms.codecs import from_bytes  # noqa: E402


def first_models(rec):
    """One assembly key per assembly id: the lowest model, then the first alt (NMR bundles and alt locations are not separate items)."""

    def order(k):
        a, m, alt = k.split("_", 2)
        return (a, int(m) if m.isdigit() else 10**9, alt)

    best = {}
    for k in sorted(rec, key=order):
        best.setdefault(k.split("_", 1)[0], k)
    return [best[a] for a in sorted(best)]


def work(task):
    w, keys, db, out, assemblies = task
    envs = [
        lmdb.open(
            sp,
            readonly=True,
            lock=False,
            readahead=False,
            subdir=Path(sp).is_dir(),
            max_readers=512,
        )
        for sp in shard_paths(db)
    ]

    def get(key):
        for env in envs:
            with env.begin() as txn:
                v = txn.get(key)
            if v is not None:
                return v
        raise KeyError(key)

    n = fail = 0
    skipped = {}
    t0 = time.time()

    def add(tar, tsv, name, cm, extra):
        _, s = secondary_structure(cm)
        data = backbone_pdb_text(cm).encode()
        ti = tarfile.TarInfo(f"{name}.pdb")
        ti.size = len(data)
        tar.addfile(ti, io.BytesIO(data))
        plddt = (
            "nan" if assemblies else round(float(np.mean(cm.atoms.b_factor.value)), 2)
        )
        tsv.write(
            "\t".join(
                map(
                    str,
                    [
                        name,
                        s["n_assessed"],
                        plddt,
                        round(s["helix"], 4),
                        round(s["strand"], 4),
                        round(s["coil"], 4),
                        s["max_coil_run"],
                        *extra,
                    ],
                )
            )
            + "\n"
        )

    with (
        open(Path(out) / f"part_{w:03d}.tsv", "w") as tsv,
        tarfile.open(Path(out) / f"part_{w:03d}.tar", "w") as tar,
    ):
        for key in keys:
            try:
                rec = from_bytes(get(key))
                if not assemblies:
                    cm = CIFMolAttached.from_dict(
                        rec[sorted(rec)[0]]["cifmol_attached_dict"]
                    )
                    add(tar, tsv, key.decode(), cm, [])
                    n += 1
                else:
                    for ak in first_models(rec):
                        cm = CIFMolAttached.from_dict(rec[ak]["cifmol_attached_dict"])
                        nc = n_protein_chains(cm)
                        if nc < 2 or nc > len(CHAIN_LETTERS):
                            reason = (
                                "fewer than 2 protein chains"
                                if nc < 2
                                else f"more than {len(CHAIN_LETTERS)} protein chains"
                            )
                            skipped[reason] = skipped.get(reason, 0) + 1
                            continue
                        add(
                            tar,
                            tsv,
                            re.sub(r"[^A-Za-z0-9_]", "-", f"{key.decode()}_{ak}"),
                            cm,
                            [key.decode(), ak],
                        )
                        n += 1
            except Exception:  # noqa: BLE001  (a failed record is counted)
                fail += 1
                if fail <= 3:
                    traceback.print_exc()
            if (n + fail) % 20000 == 0 and w == 0:
                print(
                    f"worker0 {n + fail}/{len(keys)} {time.time() - t0:.0f}s",
                    flush=True,
                )
    return n, fail, skipped


def shard_paths(db):
    """The LMDB environments behind ``db``: itself, or the shards of a shard collection (shards.json)."""
    meta = Path(db) / "shards.json"
    if meta.exists():
        return [
            str(Path(db) / x["path"]) for x in json.loads(meta.read_text())["shards"]
        ]
    return [db]


def list_keys(db):
    """All keys of ``db``: the ``<db>.index.tsv`` beside it when there is one (listing a 900 GB LMDB reads all of it), else a scan."""
    idx = Path(str(db).rstrip("/") + ".index.tsv")
    if not idx.exists():
        idx = Path(db).parent / (Path(db).name + ".index.tsv")
    if idx.exists():
        with open(idx) as f:
            return [
                ln.split("\t", 1)[0].encode()
                for ln in f
                if ln[0] != "#" and not ln.startswith("key\t")
            ]
    keys = []
    for sp in shard_paths(db):
        env = lmdb.open(sp, readonly=True, lock=False, subdir=Path(sp).is_dir())
        with env.begin() as txn:
            keys += list(txn.cursor().iternext(values=False))
        env.close()
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=0, help="number of keys")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--assemblies",
        action="store_true",
        help="complexes: one item per assembly with at least 2 protein chains (62 at most); columns entry, asm appended",
    )
    ap.add_argument(
        "--part",
        default="0/1",
        help="i/n: this job takes every n-th key (after the seed shuffle) starting at i; files are prefixed p<i>_",
    )
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    keys = list_keys(a.db)
    random.Random(a.seed).shuffle(keys)
    if a.limit:
        keys = keys[: a.limit]
    pi, pn = map(int, a.part.split("/"))
    keys = keys[pi::pn]
    print(f"{len(keys)} keys", flush=True)
    tasks = [
        (pi * a.workers + w, keys[w :: a.workers], a.db, a.out, a.assemblies)
        for w in range(a.workers)
    ]
    with mp.get_context("fork").Pool(a.workers) as p:
        res = p.map(work, tasks, chunksize=1)
    skipped = {}
    for r in res:
        for k, v in r[2].items():
            skipped[k] = skipped.get(k, 0) + v
    print(
        "done",
        sum(r[0] for r in res),
        "failed",
        sum(r[1] for r in res),
        "skipped",
        skipped,
    )


if __name__ == "__main__":
    main()
