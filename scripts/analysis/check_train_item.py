"""Spot-check train_item rows against the LMDBs MiniWorld reads (sbatch only).

usage: check_train_item.py train_item.tsv [--n 2000] [--seed 0]
Per source, N random rows. Hard failures (MiniWorld would raise): record missing in the cif DB, assembly key missing,
chain id missing. Soft misses (MiniWorld falls back to a query-only MSA / empty template): protein chain's seq_id missing in
the MSA DB, template key missing ({record}_{label}; the bare {record} for the monomer sources, as MiniWorld's chain_template_key).
"""

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes

B = Path("/data/shared/cssb_data/BioMol_clean/lmdb")
SRC = {
    "pdb": {
        "cif": B
        / "pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb",
        "msa": [B / "pdb/msa/a3m.lmdb"],
        "tpl": B / "chain/template/template.lmdb",
    },
    "protein_monomer": {
        "cif": B / "distillation_long/cif_attached/cif_long_attached.lmdb",
        "msa": [B / "distillation_long/msa/msa_long_d2k.lmdb"],
        "tpl": B / "distillation_long/template/template_topn_long.lmdb",
    },
    "short_protein_monomer": {
        "cif": B / "distillation_short/cif_attached/cif_short_attached.lmdb",
        "msa": [B / "distillation_short/msa/msa_short_d2k.lmdb"],
        "tpl": B / "distillation_short/template/template_topn_short.lmdb",
    },
    "disordered_pdb": {
        "cif": B / "distillation_disordered/cif_attached/cif_disordered_attached.lmdb",
        "msa": [B / "distillation_disordered/msa/msa_disordered_d2k.lmdb"],
        "tpl": B / "distillation_disordered/template/template_disordered.lmdb",
    },
}


def envs(p):
    p = Path(p)
    sj = p / "shards.json"
    paths = (
        [p / x["path"] for x in json.loads(sj.read_text())["shards"]]
        if sj.exists()
        else [p]
    )
    return [
        lmdb.open(
            str(q),
            readonly=True,
            lock=False,
            readahead=False,
            max_readers=512,
            subdir=q.is_dir(),
        )
        for q in paths
    ]


def get(es, key):
    for e in es:
        with e.begin() as t:
            v = t.get(key.encode())
        if v is not None:
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    by = defaultdict(list)
    with open(a.tsv) as f:
        next(f)
        for i, ln in enumerate(f):
            by[ln.split("\t", 1)[0]].append(i)
    rng = random.Random(a.seed)
    pick = {s: set(rng.sample(v, min(a.n, len(v)))) for s, v in by.items() if s in SRC}
    rows = defaultdict(list)
    with open(a.tsv) as f:
        next(f)
        for i, ln in enumerate(f):
            s = ln.split("\t", 1)[0]
            if s in pick and i in pick[s]:
                rows[s].append(ln.rstrip("\n").split("\t"))
    for s, rs in rows.items():
        cif, msa, tpl = (
            envs(SRC[s]["cif"]),
            [e for p in SRC[s]["msa"] for e in envs(p)],
            envs(SRC[s]["tpl"]),
        )
        c = Counter()
        bad = []
        for r in rs:
            _, k1, k2, pid, asm, mdl, alt, ch1, ch2 = r
            rec_id = pid.lower() if s == "pdb" else pid
            raw = get(cif, rec_id)
            if raw is None:
                c["no_record"] += 1
                bad.append(("no_record", rec_id))
                continue
            rec = from_bytes(bytes(raw))
            ak = f"{asm}_{mdl}_{alt}"
            if ak not in rec:
                c["no_assembly"] += 1
                bad.append(("no_assembly", rec_id, ak))
                continue
            d = rec[ak]["cifmol_attached_dict"]
            nodes = d["chains"]["nodes"]
            ids = np.asarray(nodes["chain_id"]["value"]).astype(str).tolist()
            seq = np.asarray(nodes["seq_id"]["value"]).astype(str).tolist()
            et = np.asarray(nodes["entity_type"]["value"]).astype(str).tolist()
            ok = True
            for ch in [ch1] if k2 == "None" else [ch1, ch2]:
                if ch not in ids:
                    c["no_chain"] += 1
                    bad.append(("no_chain", rec_id, ak, ch))
                    ok = False
                    continue
                j = ids.index(ch)
                if et[j] != "polypeptide(L)":
                    continue
                c["protein_chains"] += 1
                c["msa_hit"] += get(msa, seq[j]) is not None
                key = (
                    rec_id
                    if s in ("protein_monomer", "short_protein_monomer")
                    else f"{rec_id}_{ch.rsplit('_', 1)[0]}"
                )
                c["tpl_hit"] += get(tpl, key) is not None
            c["ok_rows"] += ok
        print(
            f"{s}: rows {len(rs)}  hard-ok {c['ok_rows']}  no_record {c['no_record']} no_assembly {c['no_assembly']} no_chain {c['no_chain']}  "
            f"protein chains {c['protein_chains']}  msa hit {c['msa_hit']}  template hit {c['tpl_hit']}",
            flush=True,
        )
        for b in bad[:5]:
            print("   bad", b)


if __name__ == "__main__":
    main()
