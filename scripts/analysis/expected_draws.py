"""Expected draws per PDB edge cluster under MiniWorld's train_item weighting (sbatch only).

usage: expected_draws.py train_item.tsv [total_draws]
P(source) = source_weights normalised; P(type) = sampler weight / sum over present types;
uniform over the type's edge clusters (edge_id = cluster1 or cluster1_cluster2); uniform over rows inside a cluster.
"""

import re
import sys
from collections import Counter, defaultdict

SAMPLER = {
    "protein_protein": 25.0,
    "protein_ligand": 25,
    "protein_dna": 10,
    "protein_rna": 10,
    "antibody_protein": 15,
    "dna_dna": 5,
    "rna_rna": 5,
    "dna_rna": 0.5,
    "antibody_antibody": 0.5,
    "antibody_ligand": 1.0,
    "na_ligand": 1.0,
    "etc_interface": 1.0,
    "sole": 1.0,
}
P_SRC = (
    0.5 / 1.02
)  # pdb among pdb 0.5, disordered 0.02, long 0.495, short 0.005 (all four loaded)


def edge_type(edge_id):
    cl = re.findall(r"c(?:PDB_)?([A-Z])\d", edge_id)
    if len(cl) < 2:
        return "sole"
    s = set(cl)
    if s == {"P"}:
        return "protein_protein"
    if s <= {"P", "L", "B"} and "P" in s:
        return "protein_ligand"
    if s == {"P", "D"}:
        return "protein_dna"
    if s == {"P", "R"}:
        return "protein_rna"
    if s == {"P", "A"}:
        return "antibody_protein"
    if s == {"D"}:
        return "dna_dna"
    if s == {"R"}:
        return "rna_rna"
    if s == {"D", "R"}:
        return "dna_rna"
    if s == {"A"}:
        return "antibody_antibody"
    if s <= {"A", "L", "B"} and "A" in s:
        return "antibody_ligand"
    if s <= {"N", "L", "B"} and "N" in s:
        return "na_ligand"
    return "etc_interface"


def main():
    total = float(sys.argv[2]) if len(sys.argv) > 2 else 25600 * 800
    rows = Counter()
    with open(sys.argv[1]) as f:
        next(f)
        for ln in f:
            x = ln.split("\t", 3)
            rows[x[1] if x[2] == "None" else f"{x[1]}_{x[2]}"] += 1
    by_type = defaultdict(list)
    for e, n in rows.items():
        by_type[edge_type(e)].append(n)
    wsum = sum(SAMPLER[t] for t in by_type)
    print(
        f"total draws {total:,.0f}; pdb share {P_SRC:.4f}; types present weight sum {wsum}"
    )
    print("type\tclusters\trows\tP(type)\texp_draws_per_cluster\tper_row_draws")
    for t, ns in sorted(by_type.items(), key=lambda kv: -SAMPLER[kv[0]]):
        pt = SAMPLER[t] / wsum
        per = total * P_SRC * pt / len(ns)
        print(
            f"{t}\t{len(ns):,}\t{sum(ns):,}\t{pt:.4f}\t{per:.2f}\t{per / (sum(ns) / len(ns)):.3f}"
        )


if __name__ == "__main__":
    main()
