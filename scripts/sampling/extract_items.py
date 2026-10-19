"""Extract the sampling pool of a PDB training DB: every chain and every chain-chain interface.

One pass over a cif_attached LMDB (schema B: one record per entry, one cifmol per bioassembly).
For each bioassembly it writes

* ``chains.<pid>.tsv``      entry, assembly, chain_id, entity_type, cluster_id, n_resolved_res, atoms, seq_id
* ``interfaces.<pid>.tsv``  entry, assembly, chain_id1, chain_id2, cluster1, cluster2, in5, seq_id1, seq_id2
                            -- the stored contact edges (any atom pair within 6 A, the pool the
                               legacy ``edge_node`` file was cut from), with ``in5 = 1`` when the
                               pair also has a heavy-atom pair within 5 A (AF3's interface rule)

Only chains with at least one finite coordinate are kept, as the legacy extraction did.
No sampling happens here; ``compare_sampling.py`` turns these tables into sampling weights.

    python extract_items.py --db <cif_attached.lmdb> --out <dir> [--workers 32] [--limit N]
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.geometry import chain_contacts_grid

_STATE: dict = {}


def _keys(index_tsv: Path) -> list[tuple[str, int]]:
    rows = []
    with index_tsv.open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("key\t"):
                continue
            key, size = line.rstrip("\n").split("\t")[:2]
            rows.append((key, int(size)))
    rows.sort(key=lambda r: -r[1])          # biggest first: the tail of the pool balances the workers
    return rows


def _init(db: str, out: str) -> None:
    pid = os.getpid()
    _STATE["env"] = lmdb.open(db, readonly=True, lock=False, readahead=False, max_readers=512)
    _STATE["chains"] = open(Path(out) / f"chains.{pid}.tsv", "a", buffering=1 << 20)       # noqa: SIM115
    _STATE["ifaces"] = open(Path(out) / f"interfaces.{pid}.tsv", "a", buffering=1 << 20)   # noqa: SIM115


def _assembly(entry: str, asm: str, d: dict, chains_out, ifaces_out) -> tuple[int, int]:
    nodes = d["chains"]["nodes"]
    chain_id = np.asarray(nodes["chain_id"]["value"]).astype(str)
    etype = np.asarray(nodes["entity_type"]["value"]).astype(str)
    cluster = np.asarray(nodes["cluster_id"]["value"]).astype(str)
    seq_id = np.asarray(nodes["seq_id"]["value"]).astype(str)
    n_chains = len(chain_id)
    xyz = np.asarray(d["atoms"]["nodes"]["xyz"]["value"])
    atom_to_res = np.asarray(d["index_table"]["atom_to_res"])
    res_to_chain = np.asarray(d["index_table"]["res_to_chain"])
    chain_idx = res_to_chain[atom_to_res]
    finite = np.isfinite(xyz).all(axis=1)
    atoms_per_chain = np.bincount(chain_idx[finite], minlength=n_chains)
    res_ok = np.zeros(len(res_to_chain), dtype=bool)
    res_ok[atom_to_res[finite]] = True
    res_per_chain = np.bincount(res_to_chain[res_ok], minlength=n_chains)
    valid = atoms_per_chain > 0
    for i in np.nonzero(valid)[0]:
        chains_out.write(f"{entry}\t{asm}\t{chain_id[i]}\t{etype[i]}\t{cluster[i]}\t{res_per_chain[i]}\t{atoms_per_chain[i]}\t{seq_id[i]}\n")

    contact = d["chains"]["edges"]["contact"]
    s6 = np.asarray(contact["src_indices"])
    d6 = np.asarray(contact["dst_indices"])
    if len(s6) == 0:
        return int(valid.sum()), 0
    s5, d5, _ = chain_contacts_grid(xyz, chain_idx, 5.0)
    near = set((np.minimum(s5, d5) * n_chains + np.maximum(s5, d5)).tolist())
    n_if = 0
    for s, t in zip(s6.tolist(), d6.tolist(), strict=True):
        if not (valid[s] and valid[t]):
            continue
        in5 = int((min(s, t) * n_chains + max(s, t)) in near)
        ifaces_out.write(f"{entry}\t{asm}\t{chain_id[s]}\t{chain_id[t]}\t{cluster[s]}\t{cluster[t]}\t{in5}\t{seq_id[s]}\t{seq_id[t]}\n")
        n_if += 1
    return int(valid.sum()), n_if


def _work(keys: list[str]) -> tuple[int, int, int]:
    env, chains_out, ifaces_out = _STATE["env"], _STATE["chains"], _STATE["ifaces"]
    n_asm = n_ch = n_if = 0
    for key in keys:
        with env.begin() as txn:
            raw = txn.get(key.encode())
        if raw is None:
            continue
        record = from_bytes(bytes(raw))
        for asm, wrapped in record.items():
            c, i = _assembly(key, asm, wrapped["cifmol_attached_dict"], chains_out, ifaces_out)
            n_asm += 1
            n_ch += c
            n_if += i
    chains_out.flush()
    ifaces_out.flush()
    return n_asm, n_ch, n_if


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=0, help="only the N smallest entries (smoke test)")
    ap.add_argument("--chunk", type=int, default=16)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = _keys(Path(str(args.db) + ".index.tsv"))
    if args.limit:
        rows = rows[-args.limit:]
    keys = [k for k, _ in rows]
    chunks = [keys[i:i + args.chunk] for i in range(0, len(keys), args.chunk)]
    print(f"{len(keys):,} entries in {len(chunks):,} chunks, {args.workers} workers", flush=True)
    tot = [0, 0, 0]
    ctx = mp.get_context("forkserver")
    with ctx.Pool(args.workers, initializer=_init, initargs=(str(args.db), str(args.out))) as pool:
        for n, (a, c, i) in enumerate(pool.imap_unordered(_work, chunks), 1):
            tot[0] += a
            tot[1] += c
            tot[2] += i
            if n % 500 == 0 or n == len(chunks):
                print(f"  {n}/{len(chunks)} chunks  assemblies {tot[0]:,}  chains {tot[1]:,}  interfaces {tot[2]:,}", flush=True)
    print(f"done: assemblies {tot[0]:,}  chains {tot[1]:,}  interfaces {tot[2]:,}")


if __name__ == "__main__":
    main()
