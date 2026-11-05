"""Design training DB of monomers in BioMol format (sbatch only).

usage: build_db.py --db SRC --items items.tsv --out DIR --source NAME [--workers 56] [--limit N] [--part i/n]
       build_db.py --finalize --out DIR --name NAME --source NAME
Every item whose ``pass`` column in items.tsv (summarize.py) is 1 is copied from the source shard collection unchanged
(``{assembly key: {"cifmol_attached_dict": ...}}``) with these keys added to the CIFMolAttached ``metadata``:
``source``, ``length`` (residues with a complete backbone), ``mean_plddt``, ``ss_helix`` / ``ss_strand`` / ``ss_coil`` (fractions of
the assessed residues, structcooker.instructions.design.ss), ``ss_max_coil_run``, ``fs_cluster`` (id of the foldseek representative),
``fs_cluster_size`` (members of that cluster in the source, passing or not). Worker w writes DIR/<name>.lmdb/shards/shard_<w>.lmdb and
DIR/manifest_<w>.tsv; --finalize writes shards.json, the meta file and manifest.tsv.
"""

import argparse
import json
import os
import random
import sys
from collections import Counter
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

sys.path.insert(0, str(Path(__file__).parent))

import multiprocessing as mp  # noqa: E402

import lmdb  # noqa: E402
from build_stats import shard_paths  # noqa: E402

from structcooker.instructions.transforms.codecs import (  # noqa: E402
    from_bytes,
    to_bytes,
)

LINES: list[list[str]] = []
SIZE: dict[str, int] = {}
SOURCE = ""


def work(task):
    w, lo, hi, db, out, name = task
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
    shard = Path(out) / f"{name}.lmdb" / "shards" / f"shard_{w:03d}.lmdb"
    shard.parent.mkdir(parents=True, exist_ok=True)
    dst = lmdb.open(
        str(shard), map_size=200 << 30, subdir=True, writemap=False, max_readers=8
    )
    n = 0
    txn = dst.begin(write=True)
    with open(Path(out) / f"manifest_{w:03d}.tsv", "w") as man:
        for x in LINES[lo:hi]:
            key = (
                x[9] if len(x) > 9 else x[0]
            ).encode()  # complexes by assembly: the entry holds the assembly x[10]
            raw = None
            for env in envs:
                with env.begin() as t:
                    raw = t.get(key)
                if raw is not None:
                    break
            rec = from_bytes(bytes(raw))
            if len(x) > 9:
                rec = {x[10]: rec[x[10]]}
            for v in rec.values():
                md = v["cifmol_attached_dict"]["metadata"]
                if x[2] != "nan":
                    md["mean_plddt"] = float(x[2])
                md.update(
                    {
                        "source": SOURCE,
                        "length": int(x[1]),
                        "ss_helix": float(x[3]),
                        "ss_strand": float(x[4]),
                        "ss_coil": float(x[5]),
                        "ss_max_coil_run": int(x[6]),
                        "fs_cluster": x[7],
                        "fs_cluster_size": SIZE[x[7]],
                    }
                )
            txn.put(x[0].encode(), to_bytes(rec))
            man.write("\t".join([*x[:8], str(SIZE[x[7]]), str(w), *x[9:]]) + "\n")
            n += 1
            if n % 20000 == 0:
                txn.commit()
                txn = dst.begin(write=True)
                if w % 56 == 0:
                    print(f"worker {w}: {n}/{hi - lo}", flush=True)
    txn.commit()
    dst.sync()
    dst.close()
    return n


def finalize(out, name, source):
    out = Path(out)
    shards = []
    for p in sorted((out / f"{name}.lmdb" / "shards").glob("shard_*.lmdb")):
        env = lmdb.open(str(p), readonly=True, lock=False, subdir=True)
        shards.append({"path": f"shards/{p.name}", "entries": env.stat()["entries"]})
        env.close()
    total = sum(s["entries"] for s in shards)
    (out / f"{name}.lmdb" / "shards.json").write_text(
        json.dumps(
            {
                "format": "shard_lmdb",
                "version": 1,
                "entries": total,
                "merge_limit": 5000000,
                "shards": shards,
            },
            indent=2,
        )
    )
    (out / f"{name}.lmdb.meta.json").write_text(
        json.dumps(
            {
                "schema": "design-monomer v1",
                "key_convention": "source key",
                "codec": "biomol_zstd",
                "entries": total,
                "source": source,
                "recipe": "scripts/design/build_db.py",
                "extra": {"shards": len(shards)},
            },
            indent=2,
        )
    )
    first = next((ln for p in sorted(out.glob("manifest_*.tsv")) for ln in open(p)), "")
    with_asm = first.count("\t") > 9
    with open(out / "manifest.tsv", "w") as m:
        m.write(
            "id\tlength\tmean_plddt\tss_helix\tss_strand\tss_coil\tss_max_coil_run\tfs_cluster\tfs_cluster_size\tshard"
            + ("\tentry\tasm" if with_asm else "")
            + "\n"
        )
        for p in sorted(out.glob("manifest_*.tsv")):
            m.write(p.read_text())
    print("entries", total, "shards", len(shards))


def main():
    global SOURCE
    ap = argparse.ArgumentParser()
    ap.add_argument("--db")
    ap.add_argument("--items")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="design_monomer")
    ap.add_argument("--source", default="")
    ap.add_argument("--workers", type=int, default=56)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--part", default="0/1")
    ap.add_argument("--finalize", action="store_true")
    a = ap.parse_args()
    if a.finalize:
        return finalize(a.out, a.name, a.source)
    SOURCE = a.source
    Path(a.out).mkdir(parents=True, exist_ok=True)
    sizes = Counter()
    rows = []
    with open(a.items) as f:
        next(f)
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            sizes[x[7]] += 1
            if x[8] == "1":
                rows.append(x)
    SIZE.update(sizes)
    random.Random(0).shuffle(rows)
    if a.limit:
        rows = rows[: a.limit]
    pi, pn = map(int, a.part.split("/"))
    rows = rows[pi::pn]
    LINES.extend(rows)
    print(f"{len(rows)} items", flush=True)
    step = -(-len(rows) // a.workers)
    tasks = [
        (
            pi * a.workers + w,
            w * step,
            min((w + 1) * step, len(rows)),
            a.db,
            a.out,
            a.name,
        )
        for w in range(a.workers)
    ]
    with mp.get_context("fork").Pool(a.workers) as p:
        res = p.map(work, tasks, chunksize=1)
    print("written", sum(res))


if __name__ == "__main__":
    main()
