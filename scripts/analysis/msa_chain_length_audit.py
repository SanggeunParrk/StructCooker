"""Does each chain's MSA have the chain's length?

The MSA of a chain is looked up by the chain's ``seq_id`` (protein: ``P...``, RNA: ``R...``); its query must have as many
residues as the chain in the structure record, or the consumer's MSA features do not line up with the tokens.

    python msa_chain_length_audit.py [--records 150]
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes

L = Path("/data/shared/cssb_data/BioMol_clean/lmdb")
# name, structure DB, [MSA DBs tried in order]
SETS = [
    (
        "pdb train_20210930",
        "pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300",
        ["pdb/msa/a3m_d2k", "pdb/msa/msa_rna_d2k"],
    ),
    (
        "pdb valid_1",
        "pdb/valid/cif_attached_valid_1",
        ["pdb/msa/a3m_d2k", "pdb/msa/msa_rna_d2k"],
    ),
    (
        "pdb all",
        "pdb/cif_attached/cif_pdb_attached",
        ["pdb/msa/a3m_d2k", "pdb/msa/msa_rna_d2k"],
    ),
    (
        "afm homodimer",
        "afdb_multimer/cif_attached/cif_afm_homodimer_attached",
        ["afdb_multimer/msa/msa_afm_d2k"],
    ),
    (
        "afm heterodimer",
        "afdb_multimer/cif_attached/cif_afm_heterodimer_attached",
        ["afdb_multimer/msa/msa_afm_d2k"],
    ),
    (
        "teddymer",
        "teddymer/cif_attached/cif_teddymer_attached",
        ["teddymer/msa/msa_teddymer_d2k"],
    ),
    (
        "short",
        "distillation_short/cif_attached/cif_short_attached",
        ["distillation_short/msa/msa_short_d2k"],
    ),
    (
        "long",
        "distillation_long/cif_attached/cif_long_attached",
        ["distillation_long/msa/msa_long_d2k"],
    ),
    (
        "ofd rna",
        "distillation_rna/cif_attached/cif_rna_attached",
        ["distillation_rna/msa/msa_rna_d2k"],
    ),
    (
        "disordered",
        "distillation_disordered/cif_attached/cif_disordered_attached",
        ["distillation_disordered/msa/msa_disordered_d2k"],
    ),
]


def col(nodes, name):
    c = nodes[name]
    return np.asarray(c["value"] if isinstance(c, dict) else c).ravel()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", type=int, default=150)
    a = ap.parse_args()
    for name, sdb, mdbs in SETS:
        names = [
            ln.split("\t", 1)[0]
            for ln in open(f"{L / sdb}.lmdb.index.tsv")
            if not ln.startswith(("#", "key\t"))
        ]
        chosen = random.Random("ml-" + name).sample(names, min(a.records, len(names)))
        chains = found = equal = 0
        short, long_, ex = [], [], []
        with (
            open_env(str(L / f"{sdb}.lmdb"), readonly=True, lock=False) as se,
            se.begin() as st,
        ):
            mts = [
                open_env(str(L / f"{m}.lmdb"), readonly=True, lock=False) for m in mdbs
            ]
            txns = [m.begin() for m in mts]
            for rec_name in chosen:
                rec = from_bytes(bytes(st.get(rec_name.encode())))
                for asm in rec.values():
                    v = asm["cifmol_attached_dict"]
                    cn = v["chains"]["nodes"]
                    ptr = np.asarray(v["index_table"]["chain_res_indptr"]).ravel()
                    for i, (et, sid) in enumerate(
                        zip(col(cn, "entity_type"), col(cn, "seq_id"), strict=False)
                    ):
                        if not any(
                            t in str(et) for t in ("polypeptide", "polyribonucleotide")
                        ):
                            continue
                        chains += 1
                        raw = next(
                            (
                                t.get(str(sid).encode())
                                for t in txns
                                if t.get(str(sid).encode()) is not None
                            ),
                            None,
                        )
                        if raw is None:
                            continue
                        found += 1
                        qlen = len(
                            from_bytes(bytes(raw))["msa_dict"]["sequences"][
                                "query_sequence"
                            ]
                        )
                        clen = int(ptr[i + 1] - ptr[i])
                        if qlen == clen:
                            equal += 1
                        else:
                            (short if qlen < clen else long_).append(clen - qlen)
                            if len(ex) < 3:
                                ex.append((rec_name, str(sid), clen, qlen))
        print(
            f"{name:20} chains={chains} with_msa={found} same_length={equal} msa_shorter={len(short)} msa_longer={len(long_)} "
            f"ex={ex}",
            flush=True,
        )


if __name__ == "__main__":
    main()
