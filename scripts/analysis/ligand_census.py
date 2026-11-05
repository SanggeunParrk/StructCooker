"""Which AF3 SI Table 9 / Table 10 components are present in a cif_attached train DB (sbatch only).

usage: ligand_census.py --db DB --out DIR --workers N [--every K]
Counts, per CCD code of every non-polymer chain's residues, the chains and entries (assemblies) that carry it, over
every K-th key of the DB index. Writes census.tsv (ccd, chains, assemblies, in_table9, in_table10) and cluster_ccd.tsv (cluster, ccd, chains).
"""

import argparse
import json
import multiprocessing as mp
from collections import Counter
from pathlib import Path

import lmdb
import numpy as np
import yaml

from structcooker.instructions.transforms.codecs import from_bytes

_S: dict = {}


def _init(db):
    sj = Path(db) / "shards.json"
    paths = (
        [str(Path(db) / x["path"]) for x in json.loads(sj.read_text())["shards"]]
        if sj.exists()
        else [db]
    )
    _S["envs"] = [
        lmdb.open(
            p,
            readonly=True,
            lock=False,
            readahead=False,
            max_readers=512,
            subdir=Path(p).is_dir(),
        )
        for p in paths
    ]


def _work(keys):
    chains, asms, n_asm, cl = Counter(), Counter(), 0, Counter()
    for key in keys:
        raw = None
        for env in _S["envs"]:
            with env.begin() as txn:
                raw = txn.get(key.encode())
            if raw is not None:
                break
        if raw is None:
            continue
        for asm, rec in from_bytes(raw).items():
            d = rec["cifmol_attached_dict"]
            n_asm += 1
            etype = np.asarray(d["chains"]["nodes"]["entity_type"]["value"]).astype(str)
            res_to_chain = np.asarray(d["index_table"]["res_to_chain"])
            ccd = np.asarray(d["residues"]["nodes"]["chem_comp_id"]["value"]).astype(
                str
            )
            cluster = np.asarray(d["chains"]["nodes"]["cluster_id"]["value"]).astype(
                str
            )
            seen = set()
            for c in np.nonzero(etype == "non-polymer")[0]:
                for code in set(ccd[res_to_chain == c].tolist()):
                    chains[code] += 1
                    cl[(cluster[c], code)] += 1
                    seen.add(code)
            for code in seen:
                asms[code] += 1
    return chains, asms, n_asm, cl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=56)
    ap.add_argument("--every", type=int, default=1)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    tab = yaml.safe_load(
        (root / "db_db_configs/interface/params/params/af3_si_tables.yaml").read_text()
    )
    t9 = {str(c).upper() for c in tab["table9_crystallization_aids"]}
    t10 = {str(c).upper() for c in tab["table10_ligand_exclusion_list"]}
    index = next(Path(a.db).glob("*.index.tsv"), None) or Path(str(a.db) + ".index.tsv")
    keys = []
    for p in [index] if index.exists() else Path(a.db).glob("*.index.tsv"):
        for ln in p.open():
            if ln.startswith(("#", "key\t")):
                continue
            keys.append(ln.split("\t")[0])
    keys = keys[:: a.every]
    chunks = [keys[i : i + 200] for i in range(0, len(keys), 200)]
    chains, asms, n_asm, cl = Counter(), Counter(), 0, Counter()
    with mp.Pool(a.workers, _init, (a.db,)) as pool:
        for c, s, n, k in pool.imap_unordered(_work, chunks):
            cl.update(k)
            chains.update(c)
            asms.update(s)
            n_asm += n
    Path(a.out).mkdir(parents=True, exist_ok=True)
    with open(Path(a.out) / "census.tsv", "w") as f:
        f.write(
            f"# keys {len(keys)} assemblies {n_asm}\nccd\tchains\tassemblies\tin_table9\tin_table10\n"
        )
        for code, n in chains.most_common():
            f.write(
                f"{code}\t{n}\t{asms[code]}\t{int(code.upper() in t9)}\t{int(code.upper() in t10)}\n"
            )
    with open(Path(a.out) / "cluster_ccd.tsv", "w") as f:
        f.write("cluster\tccd\tchains\n")
        for (k, code), n in sorted(cl.items()):
            f.write(f"{k}\t{code}\t{n}\n")
    print(f"keys {len(keys)} assemblies {n_asm} ccd {len(chains)}")


if __name__ == "__main__":
    main()
