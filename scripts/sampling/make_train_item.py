"""Assemble train_item.tsv (the sample list MiniWorld reads) from extract_items.py outputs (sbatch only).

usage: make_train_item.py OUT.tsv [--exclude CLUSTERS.txt] SOURCE=DIR[,DIR...] ...
  SOURCE is the name MiniWorld's DataConfig uses (pdb, protein_monomer, short_protein_monomer, disordered_pdb, rna);
  each DIR holds chains.*.tsv / interfaces.*.tsv written by extract_items.py.

Columns: source cluster1 cluster2 pdb_id assembly_id model_id alt_id chain_id1 chain_id2
* an interface row per chain pair in contact (6 A), once per unordered pair, clusters sorted (chains follow their clusters);
* --exclude: a file of cluster ids (one per line) whose chains are dropped as if absent from the structure (no interface or
  monomer row carries them; a partner that only touched them falls back to the monomer rule);
* a monomer row (cluster2 = None) per chain whose cluster appears in no interface of the source -- the legacy edge_node rule.
"""

import glob
import sys
from pathlib import Path

out_path, specs = sys.argv[1], sys.argv[2:]
drop: set[str] = set()
if specs and specs[0] == "--exclude":
    drop = {ln.strip() for ln in open(specs[1]) if ln.strip()}
    specs = specs[2:]
HEADER = "source\tcluster1\tcluster2\tpdb_id\tassembly_id\tmodel_id\talt_id\tchain_id1\tchain_id2\n"
with open(out_path, "w") as out:
    out.write(HEADER)
    for spec in specs:
        source, dirs = spec.split("=", 1)
        dirs = dirs.split(",")
        in_iface: set[str] = set()
        n_if = n_mono = 0
        for d in dirs:
            for p in sorted(glob.glob(str(Path(d) / "interfaces.*.tsv"))):
                cur, seen = None, set()
                with open(p) as f:
                    for ln in f:
                        e, asm, c1, c2, k1, k2 = ln.rstrip("\n").split("\t")[:6]
                        if k1 in drop or k2 in drop:
                            continue
                        if (e, asm) != cur:
                            cur, seen = (e, asm), set()
                        pair = (c1, c2) if c1 <= c2 else (c2, c1)
                        if pair in seen:
                            continue
                        seen.add(pair)
                        if k1 > k2:
                            k1, k2, c1, c2 = k2, k1, c2, c1
                        in_iface.add(k1)
                        in_iface.add(k2)
                        a, m, alt = asm.split("_", 2)
                        out.write(
                            f"{source}\t{k1}\t{k2}\t{e}\t{a}\t{m}\t{alt}\t{c1}\t{c2}\n"
                        )
                        n_if += 1
        for d in dirs:
            for p in sorted(glob.glob(str(Path(d) / "chains.*.tsv"))):
                with open(p) as f:
                    for ln in f:
                        x = ln.split("\t", 5)
                        e, asm, c, k = x[0], x[1], x[2], x[4]
                        if k in in_iface or k in drop:
                            continue
                        a, m, alt = asm.split("_", 2)
                        out.write(
                            f"{source}\t{k}\tNone\t{e}\t{a}\t{m}\t{alt}\t{c}\tNone\n"
                        )
                        n_mono += 1
        print(
            f"{source}: interface rows {n_if:,}  monomer rows {n_mono:,}  interface clusters {len(in_iface):,}",
            flush=True,
        )
