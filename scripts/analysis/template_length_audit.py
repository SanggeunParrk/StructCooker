"""Do the stored templates fit the chains they are looked up for?

For sampled structure records of every set: each protein chain's template record (by the key contracts' rule) holds mols;
a template mol must have as many residues as the query chain, or the consumer's ``template.residues[crop_indices]`` raises
IndexError / reads the wrong residues. Prints, per set, the chains checked, with a template, with all mols the right length,
with some mol the wrong length, and examples.

    python template_length_audit.py [--records 150]
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.release.key_contracts import CONTRACTS, RULES

ROOT = Path("/data/shared/cssb_data/BioMol_clean/lmdb")


def col(nodes, name):
    c = nodes[name]
    return np.asarray(c["value"] if isinstance(c, dict) else c).ravel()


def nres(mol) -> int:
    nodes = mol["residues"]["nodes"]
    return len(col(nodes, "one_letter_code_can"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", type=int, default=150)
    a = ap.parse_args()
    for c in CONTRACTS:
        struct, templ = ROOT / f"{c.structure}.lmdb", ROOT / f"{c.template}.lmdb"
        names = [ln.split("\t", 1)[0] for ln in open(f"{struct}.index.tsv") if not ln.startswith(("#", "key\t"))]
        chosen = random.Random("len-" + c.name).sample(names, min(a.records, len(names)))
        n = with_t = ok = bad = nmols = badmols = 0
        ex = []
        with open_env(str(struct), readonly=True, lock=False) as se, open_env(str(templ), readonly=True, lock=False) as te, \
                se.begin() as st, te.begin() as tt:
            for name in chosen:
                rec = from_bytes(bytes(st.get(name.encode())))
                for asm in rec.values():
                    v = asm["cifmol_attached_dict"]
                    cn = v["chains"]["nodes"]
                    ptr = np.asarray(v["index_table"]["chain_res_indptr"]).ravel()
                    for i, (cid, et, sid) in enumerate(zip(col(cn, "chain_id"), col(cn, "entity_type"), col(cn, "seq_id"))):
                        if "polypeptide" not in str(et):
                            continue
                        n += 1
                        raw = tt.get(RULES[c.rule](name, str(cid), str(sid)).encode())
                        if raw is None:
                            continue
                        mols = from_bytes(bytes(raw))["template_mols"]
                        if not mols:
                            continue
                        with_t += 1
                        qlen = int(ptr[i + 1] - ptr[i])
                        lens = [nres(m) for m in mols.values()]
                        wrong = [x for x in lens if x != qlen]
                        nmols += len(lens)
                        badmols += len(wrong)
                        if wrong:
                            bad += 1
                            if len(ex) < 4:
                                ex.append((name, str(cid), qlen, wrong[:3]))
                        else:
                            ok += 1
        print(f"{c.name:16} chains={n} with_template={with_t} all_mols_fit={ok} some_mol_wrong={bad} "
              f"mols={nmols} wrong_mols={badmols} ex={ex}", flush=True)


if __name__ == "__main__":
    main()
