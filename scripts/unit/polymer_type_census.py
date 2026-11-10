#!/usr/bin/env python3
"""Census of polymer entity types in cif_pdb (assembly 1 of each entry): which chain types the unit DB's node builder
sees, and how many chains of each type it would drop as ``skipped_polymer``.

Per entity type: chains, entries, chains with >= --min-res residues, median residue count, share of residues whose
canonical one-letter code is a standard amino acid / NA base, example ids. Heavy: run it with ``polymer_type_census.sbatch``.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import defaultdict
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif/cif_pdb.lmdb"
STD = set("ACDEFGHIKLMNPQRSTVWY") | set("UT")
_S: dict = {}


def _init(db):
    _S["env"] = lmdb.open(
        str(db),
        readonly=True,
        lock=False,
        readahead=False,
        max_readers=512,
        subdir=Path(db).is_dir(),
    )


def _work(key):
    with _S["env"].begin() as txn:
        raw = txn.get(key.encode())
    try:
        a = from_bytes(raw)["assembly_dict"]
        if not a:
            return key, []
        g = a[
            sorted(
                a, key=lambda s: [int(x) if x.isdigit() else 0 for x in s.split("_")]
            )[0]
        ]
        et = g["chains"]["nodes"]["entity_type"]["value"].astype(str)
        one = g["residues"]["nodes"]["one_letter_code_can"]["value"].astype(str)
        ix = g["index_table"]
        ptr, res = (
            np.asarray(ix["chain_res_indptr"]),
            np.asarray(ix["chain_res_indices"]),
        )
        out = []
        for c, t in enumerate(et):
            if t in ("non-polymer", "water", "branched"):
                continue
            codes = one[res[ptr[c] : ptr[c + 1]]]
            out.append((t, len(codes), int(sum(x in STD for x in codes))))
    except Exception as e:  # noqa: BLE001
        return key, f"ERR {type(e).__name__}: {e}"[:150]
    else:
        return key, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=80)
    ap.add_argument("--min-res", type=int, default=4)
    a = ap.parse_args()
    keys = []
    for ln in Path(str(DB) + ".index.tsv").read_text().splitlines():
        if ln and not ln.startswith(("#", "key\t")):
            k, s = ln.split("\t")[:2]
            keys.append((k, int(s)))
    small = [k for k, s in keys if s <= 20_000_000]
    mid = [k for k, s in keys if 20_000_000 < s <= 200_000_000]
    big = [k for k, s in keys if s > 200_000_000]
    print(len(small), len(mid), len(big), flush=True)
    stat = defaultdict(
        lambda: {
            "chains": 0,
            "entries": set(),
            "ge": 0,
            "lens": [],
            "std": 0,
            "res": 0,
            "ex": [],
        }
    )
    errs = []
    for ks, w, ch in ((small, a.workers, 64), (mid, 8, 1), (big, 1, 1)):
        if not ks:
            continue
        with mp.Pool(w, _init, (DB,), maxtasksperchild=None if ch > 1 else 1) as pool:
            for key, rows in pool.imap_unordered(_work, ks, chunksize=ch):
                if isinstance(rows, str):
                    errs.append((key, rows))
                    continue
                for t, n, ns in rows:
                    s = stat[t]
                    s["chains"] += 1
                    s["entries"].add(key)
                    s["ge"] += n >= a.min_res
                    s["lens"].append(n)
                    s["std"] += ns
                    s["res"] += n
                    if len(s["ex"]) < 6:
                        s["ex"].append(key)
    with open(Path(a.out), "w") as f:
        f.write(
            f"entity_type\tchains\tentries\tchains_ge{a.min_res}\tmedian_res\tstd_code_share\texamples\n"
        )
        for t, s in sorted(stat.items(), key=lambda kv: -kv[1]["chains"]):
            f.write(
                f"{t}\t{s['chains']}\t{len(s['entries'])}\t{s['ge']}\t{int(np.median(s['lens']))}\t"
                f"{s['std'] / max(1, s['res']):.3f}\t{','.join(s['ex'])}\n"
            )
        for k, e in errs:
            f.write(f"#ERR\t{k}\t{e}\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
