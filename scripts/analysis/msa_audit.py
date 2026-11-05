"""Do the stored MSAs hold what the keys promise?

For sampled keys of every MSA DB: array shapes agree, row 0 is the query (tokens of the query letters), every token is in the
alphabet of the molecule (protein ARNDCQEGHILKMFPSTWYV = 0..19, X/other = 20, gap 31; RNA A21 U22 G23 C24 N/other 25, gap 31;
docs/stages.md), the query equals the ``seq_id_map`` sequence of the key, a depth-capped DB has at most its cap rows and its
rows are the first rows of the full DB's record.

    python msa_audit.py [--keys 200]
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes

L = Path("/data/shared/cssb_data/BioMol_clean/lmdb")
SEQ_ID_MAP = "/data/shared/cssb_data/BioMol_clean/metadata/seq_id_map.tsv"
PROT = {c: i for i, c in enumerate("ARNDCQEGHILKMFPSTWYV")} | {
    "U": 4
}  # selenocysteine is encoded as Cys
RNA = {"A": 21, "U": 22, "G": 23, "C": 24}

# (name, path, molecule, cap, full DB it is a cap of)
DBS = [
    ("pdb a3m", "pdb/msa/a3m", "P", None, None),
    ("pdb a3m_d16k", "pdb/msa/a3m_d16k", "P", 16384, "pdb a3m"),
    ("pdb a3m_d2k", "pdb/msa/a3m_d2k", "P", 2048, "pdb a3m"),
    ("pdb a3m_d512", "pdb/msa/a3m_d512", "P", 512, "pdb a3m"),
    ("pdb rna", "pdb/msa/msa_rna", "R", None, None),
    ("pdb rna_d2k", "pdb/msa/msa_rna_d2k", "R", 2048, "pdb rna"),
    ("pdb rna_d512", "pdb/msa/msa_rna_d512", "R", 512, "pdb rna"),
    ("afm", "afdb_multimer/msa/msa_afm", "P", None, None),
    ("afm d2k", "afdb_multimer/msa/msa_afm_d2k", "P", 2048, "afm"),
    ("teddymer", "teddymer/msa/msa_teddymer", "P", None, None),
    ("teddymer d2k", "teddymer/msa/msa_teddymer_d2k", "P", 2048, "teddymer"),
    ("short", "distillation_short/msa/msa_short", "P", None, None),
    ("short d2k", "distillation_short/msa/msa_short_d2k", "P", 2048, "short"),
    ("short d512", "distillation_short/msa/msa_short_d512", "P", 512, "short"),
    ("long", "distillation_long/msa/msa_long", "P", None, None),
    ("long d2k", "distillation_long/msa/msa_long_d2k", "P", 2048, "long"),
    ("long d512", "distillation_long/msa/msa_long_d512", "P", 512, "long"),
    ("ofd rna", "distillation_rna/msa/msa_rna", "R", None, None),
    ("ofd rna d2k", "distillation_rna/msa/msa_rna_d2k", "R", 2048, "ofd rna"),
    ("ofd rna d512", "distillation_rna/msa/msa_rna_d512", "R", 512, "ofd rna"),
    ("disordered", "distillation_disordered/msa/msa_disordered", "P", None, None),
    (
        "disordered d2k",
        "distillation_disordered/msa/msa_disordered_d2k",
        "P",
        2048,
        "disordered",
    ),
    (
        "disordered d512",
        "distillation_disordered/msa/msa_disordered_d512",
        "P",
        512,
        "disordered",
    ),
]


def encode(seq, mol):
    table, other = (PROT, 20) if mol == "P" else (RNA, 25)
    return np.array([table.get(str(c), other) for c in seq])


def record(txn, key):
    raw = txn.get(key.encode())
    return None if raw is None else from_bytes(bytes(raw))["msa_dict"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", type=int, default=200)
    ap.add_argument("--only", nargs="*", help="names of the DBs to check (default all)")
    a = ap.parse_args()
    sample, queries = {}, {}
    for name, path, mol, cap, full in DBS:
        keys = [
            ln.split("\t", 1)[0]
            for ln in open(f"{L / path}.lmdb.index.tsv")
            if not ln.startswith(("#", "key\t"))
        ]
        sample[name] = random.Random("msa-" + name).sample(keys, min(a.keys, len(keys)))
    want = {k for ks in sample.values() for k in ks}
    with open(SEQ_ID_MAP) as handle:
        for line in handle:
            k, _, rest = line.partition("\t")
            if k in want:
                queries[k] = rest.split("\t", 1)[0].strip()
    envs = {}
    for name, path, mol, cap, full in DBS:
        envs[name] = open_env(
            str(L / f"{path}.lmdb"), readonly=True, lock=False, readahead=False
        )
    for name, path, mol, cap, full in DBS:
        if a.only and name not in a.only:
            continue
        bad = {"shape": 0, "row0": 0, "token": 0, "query": 0, "cap": 0, "prefix": 0}
        no_query = n = 0
        ex = []
        env = envs[name]
        with env.begin() as txn:
            fulltxn = envs[full].begin() if full else None
            for key in sample[name]:
                m = record(txn, key)
                if m is None:
                    bad["shape"] += 1
                    ex.append((key, "in index, not in DB"))
                    continue
                n += 1
                kmol = (
                    "R" if key.startswith("R") else "P"
                )  # the disordered set holds both molecules under one config
                s = m["sequences"]
                q, al, de, pr = (
                    s["query_sequence"],
                    np.asarray(s["aligned_sequences"]),
                    np.asarray(s["deletions"]),
                    s["profile"],
                )
                length = len(q)
                if not (
                    al.ndim == 2
                    and al.shape[1] == length
                    and de.shape == al.shape
                    and pr.shape[0] == length
                    and all(len(v) == al.shape[0] for v in m["headers"].values())
                ):
                    bad["shape"] += 1
                    ex.append(
                        (key, f"shapes {al.shape} {de.shape} {pr.shape} q{length}")
                    )
                    continue
                if not np.array_equal(al[0], encode(q, kmol)):
                    bad["row0"] += 1
                    ex.append(
                        (
                            key,
                            f"row 0 is not the query: {al[0][:12].tolist()} vs {''.join(map(str, q))[:12]}",
                        )
                    )
                allowed = (
                    set(range(20)) | {20, 31}
                    if kmol == "P"
                    else {21, 22, 23, 24, 25, 31}
                )
                if not set(np.unique(al).tolist()) <= allowed:
                    bad["token"] += 1
                    ex.append(
                        (
                            key,
                            f"tokens {sorted(set(np.unique(al).tolist()) - allowed)[:5]}",
                        )
                    )
                ref = queries.get(key)
                if ref is None:
                    no_query += 1
                elif "".join(str(c) for c in q) != ref and "".join(
                    str(c) for c in q
                ) != ref.replace("U", "U"):
                    bad["query"] += 1
                    ex.append((key, "query differs from seq_id_map"))
                if cap is not None and al.shape[0] > cap:
                    bad["cap"] += 1
                    ex.append((key, f"{al.shape[0]} rows > cap {cap}"))
                if fulltxn is not None:
                    f = record(fulltxn, key)
                    if f is not None:
                        fal = np.asarray(f["sequences"]["aligned_sequences"])
                        k = min(al.shape[0], fal.shape[0])
                        if not np.array_equal(al[:k], fal[:k]):
                            bad["prefix"] += 1
                            ex.append(
                                (key, "rows are not the first rows of the full record")
                            )
        print(
            f"{name:16} keys={n} bad={ {k: v for k, v in bad.items() if v} or 'none' } no_seq_id_map_entry={no_query} ex={[e for e in ex if 'query differs' not in e[1]][:4]}",
            flush=True,
        )


if __name__ == "__main__":
    main()
