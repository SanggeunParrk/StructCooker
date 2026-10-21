"""Date split of the interface DB, AF3 style (sbatch only).

usage: make_split.py <db_root>   ->  <db_root>/split.tsv, split_summary.md

train    release_date <  2021-09-30, resolution <= 9.0   (the rule of train_20210930 / AF3 training)
valid    2021-09-30 <= release_date <= 2023-01-12, resolution <= 4.5   (AF3's PDB evaluation window)
test     release_date >  2023-01-12, resolution <= 4.5   (newer than AF3's window)
none     no release date / no resolution / resolution above the cut (reason column)

The split is per PDB entry, so all units of an entry share it. valid / test units also carry leakage
flags against the training set, nothing is dropped:
  n_prot_c30_in_train   protein nodes whose 30% cluster has a node in a train unit
  novel_protein         no protein node of the unit has its c30 cluster in train
  novel_pair            the unit's cluster_key30 (set of c30 clusters) does not occur in train
  ligands_seen_in_train every ligand CCD id of the unit occurs in train
"""
import csv
import sys
from collections import Counter
from pathlib import Path

db = Path(sys.argv[1])
TRAIN_END, VALID_END = "2021-09-30", "2023-01-12"
rows = list(csv.DictReader(open(db / "manifest_full.tsv"), delimiter="\t"))


def classify(r):
    d, res = r["release_date"], r["resolution"]
    if not d or d == "None":
        return "none", "no_release_date"
    if res in ("", "None", "nan"):
        return "none", "no_resolution"
    res = float(res)
    if d < TRAIN_END:
        return ("train", "") if res <= 9.0 else ("none", "resolution>9")
    if res > 4.5:
        return "none", "resolution>4.5"
    return ("valid", "") if d <= VALID_END else ("test", "")


for r in rows:
    r["split"], r["reason"] = classify(r)

tr_c, tr_pair, tr_lig = set(), set(), set()
for r in rows:
    if r["split"] != "train":
        continue
    tr_c.update(c for c in r["node_c30"].split("|") if c)
    tr_pair.add(r["cluster_key30"])
    tr_lig.update(c for k, c in zip(r["node_kinds"].split("|"), r["node_ccd"].split("|"), strict=True) if k == "ligand")

cols = ["unit_id", "pdb_id", "split", "reason", "release_date", "resolution", "type", "n_prot_c30_in_train",
        "n_protein", "novel_protein", "novel_pair", "ligands_seen_in_train", "cluster_key30"]
cnt, flags = Counter(), Counter()
with open(db / "split.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        n_in, novel_p, novel_pair, lig_seen = "", "", "", ""
        if r["split"] in ("valid", "test"):
            cl = [c for c in r["node_c30"].split("|") if c]
            n_in = sum(c in tr_c for c in cl)
            novel_p = int(n_in == 0)
            novel_pair = int(r["cluster_key30"] not in tr_pair)
            ligs = [c for k, c in zip(r["node_kinds"].split("|"), r["node_ccd"].split("|"), strict=True) if k == "ligand"]
            lig_seen = int(all(c in tr_lig for c in ligs)) if ligs else ""
            flags[(r["split"], "novel_protein")] += novel_p
            flags[(r["split"], "novel_pair")] += novel_pair
        cnt[r["split"]] += 1
        if r["split"] == "none":
            cnt["none:" + r["reason"]] += 1
        out = dict(r, n_prot_c30_in_train=n_in, novel_protein=novel_p, novel_pair=novel_pair, ligands_seen_in_train=lig_seen)
        fo.write("\t".join(str(out[c]) for c in cols) + "\n")

ents = {s: len({r["pdb_id"] for r in rows if r["split"] == s}) for s in ("train", "valid", "test")}
lines = ["# Split summary", "", f"train: release < {TRAIN_END}, resolution <= 9.0; valid: {TRAIN_END} .. {VALID_END}, resolution <= 4.5; "
         f"test: after {VALID_END}, resolution <= 4.5", "", "| split | units | entries |", "|---|---:|---:|"]
for s in ("train", "valid", "test", "none"):
    lines.append(f"| {s} | {cnt[s]:,} | {ents.get(s, '') if s != 'none' else ''} |")
lines += ["", "excluded: " + ", ".join(f"{k.split(':')[1]} {v:,}" for k, v in cnt.items() if k.startswith("none:")), "",
          "| split | units with no train-seen protein cluster | units with a novel cluster pair |", "|---|---:|---:|"]
for s in ("valid", "test"):
    lines.append(f"| {s} | {flags[(s, 'novel_protein')]:,} | {flags[(s, 'novel_pair')]:,} |")
(db / "split_summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
