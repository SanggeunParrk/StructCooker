"""Per-entry weights for the training set: 1 / (number of entries with the same cluster signature) (sbatch only).

usage: entry_weights.py <out.tsv>
An entry's signature is the set of its protein chains' 40% sequence clusters (metadata/seq_cluster40_PDB.tsv,
through seq_id_map.tsv) plus its nucleic-acid chain sequences, over all chains of cif_pdb.fasta. Entries are those
of the training DB (cif_pdb_attached_train_20210930_res9_chain300). A protein sequence with no seq_id is its own
cluster. Output columns: entry, n_same_signature, weight.
"""

import hashlib
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
train = {
    l.split("\t")[0]
    for l in open(
        ROOT
        / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb.index.tsv"
    )
    if not l.startswith(("#", "key"))
}
prot, na = defaultdict(set), defaultdict(set)
with open(ROOT / "materials/intermediate/fasta/cif_pdb.fasta") as f:
    for h in f:
        s = f.readline().strip()
        e = h[1:5].lower()
        if e not in train:
            continue
        t = h.split("|")[1].strip()
        if t.startswith("polypeptide"):
            prot[e].add(s)
        elif t.startswith(("polydeoxy", "polyribo")):
            na[e].add(s)
need = {s for v in prot.values() for s in v}
seq2id = {}
with open(ROOT / "metadata/seq_id_map.tsv") as f:
    for line in f:
        if line[0] != "P":
            continue
        i, _, s = line.rstrip("\n").partition("\t")
        if s in need:
            seq2id[s] = i
ids = set(seq2id.values())
id2c = {}
with open(ROOT / "metadata/seq_cluster40_PDB.tsv") as f:
    for line in f:
        c, _, members = line.rstrip("\n").partition("\t")
        for m in members.split(","):
            if m in ids:
                id2c[m] = c
print(
    f"{len(train)} train entries, {len(need)} protein sequences, {len(seq2id)} with seq_id",
    flush=True,
)
sig = {}
for e in train:
    parts = set()
    for s in prot.get(e, ()):
        i = seq2id.get(s)
        parts.add(id2c.get(i) or "S" + hashlib.md5(s.encode()).hexdigest()[:12])  # noqa: S324
    for s in na.get(e, ()):
        parts.add("N" + hashlib.md5(s.encode()).hexdigest()[:12])  # noqa: S324
    sig[e] = tuple(sorted(parts))
cnt = Counter(sig.values())
with open(sys.argv[1], "w") as fo:
    for e in sorted(train):
        n = cnt[sig[e]]
        fo.write(f"{e}\t{n}\t{1 / n:.6g}\n")
print(
    f"{len(cnt)} distinct signatures; largest group {cnt.most_common(1)[0][1]} entries; "
    f"effective entries {sum(1 / cnt[sig[e]] for e in train):.0f} of {len(train)}"
)
