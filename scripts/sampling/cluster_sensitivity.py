"""How the sampling pool and the sampling probabilities move with the sequence clustering.

The same items (``extract_items.py`` output, with chain seq_ids) are re-labelled with the clusters of
each clustering table, and ``compare_sampling.pdb_section`` is run on each: how many interface
clusters each class has, what mass AF3's eq. 1 and MiniWorld give it, and how often one cluster is
drawn. Nucleic acids, peptides and ligands only cluster with identical sequences in every table here;
what moves is the protein side.

    python cluster_sensitivity.py --items <dir> --out report.md \
        --map 30pdb=<seq_cluster30_PDB.tsv> --map 40pdb=<seq_cluster40_PDB.tsv> [--map name=<table> ...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_sampling as cs  # noqa: E402

CLASSES = ["protein_protein", "protein_ligand", "protein_rna", "protein_dna", "antibody_protein", "rna_rna",
           "dna_dna", "dna_rna", "antibody_antibody", "antibody_ligand", "na_ligand", "etc_interface",
           "sole", "chain:prot", "chain:nuc", "chain:lig"]
EPOCH = 128_000


def letter(rep: str) -> str:
    """Type letter of a cluster id: ``c<DB>_<letter><n>`` (per-DB tables) or ``c<letter><n>`` (legacy global)."""
    return rep.split("_", 1)[1][0] if "_" in rep else rep[1]


def load_map(path: Path, seq_names: list[str]) -> tuple[np.ndarray, list[str], list[str], int]:
    """Cluster code of every seq_id in the pool; sequences the table lacks become singletons."""
    want = {s: i for i, s in enumerate(seq_names)}
    cl_of = np.full(len(seq_names), -1, dtype=np.int32)
    names: list[str] = []
    letters: list[str] = []
    with open(path) as handle:
        for line in handle:
            rep, members = line.rstrip("\n").split("\t")
            code = None
            for m in members.split(","):
                i = want.get(m)
                if i is None:
                    continue
                if code is None:
                    code = len(names)
                    names.append(rep)
                    letters.append(letter(rep))
                cl_of[i] = code
    missing = 0
    for i in np.nonzero(cl_of < 0)[0].tolist():
        cl_of[i] = len(names)
        names.append(f"c?_{seq_names[i]}")
        letters.append(seq_names[i][0])
        missing += 1
    return cl_of, names, letters, missing


def relabel(pool: dict, path: Path) -> tuple[dict, int]:
    cl_of, names, letters, missing = load_map(path, pool["seq_names"])
    out = dict(pool)
    out.update(n_clusters=len(names), letters=np.array(letters), names=names,
               ch_cl=cl_of[pool["ch_seq"]], if_c1=cl_of[pool["if_s1"]], if_c2=cl_of[pool["if_s2"]])
    return out, missing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--items", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--map", action="append", required=True, help="name=path of a cluster table")
    args = ap.parse_args()
    pool = cs.load_pool(args.items)
    if not len(pool["ch_seq"]):
        sys.exit("the item tables carry no seq_id columns: rerun extract_items.py")
    res: dict[str, dict] = {}
    missing: dict[str, int] = {}
    for spec in args.map:
        name, path = spec.split("=", 1)
        alt, miss = relabel(pool, Path(path))
        missing[name] = miss
        sink: list[str] = []
        res[name] = cs.pdb_section(alt, sink)
        print(f"{name}: done ({alt['n_clusters']:,} clusters, {miss:,} sequences not in the table)", flush=True)
    names = list(res)
    out = ["# Clustering sensitivity of the PDB sampling pool\n",
           f"{pool['n_entries']:,} entries, {len(pool['ch_seq']):,} chain rows, {len(pool['if_s1']):,} stored contact rows. "
           "Tables: " + ", ".join(f"`{n}`" for n in names) + ". Sequences missing from a table become singletons: "
           + ", ".join(f"{n} {missing[n]:,}" for n in names) + ".\n"]

    def table(title: str, key: str, fmt=lambda v: f"{v:,}") -> None:
        out.append(f"### {title}\n")
        out.append("| class | " + " | ".join(names) + " |\n|---|" + "---:|" * len(names))
        for c in CLASSES:
            row = [res[n][key].get(c, 0) for n in names]
            if not any(row):
                continue
            base = row[0]
            cells = [fmt(v) + (f" ({v / base:.2f}x)" if i and base else "") for i, v in enumerate(row)]
            out.append(f"| {c} | " + " | ".join(cells) + " |")
        out.append("")

    table("Distinct item clusters per class, MiniWorld pool (6 A contact edges; `sole` = chains never in an interface)", "mw_clusters")
    table("Distinct item clusters per class, AF3 pool (5 A interfaces; single chains as chain:*)", "af_clusters")
    out.append("(`x` is the ratio to the first table.)\n")

    out.append("## Protein-RNA interfaces in detail\n")
    out.append("| | " + " | ".join(names) + " |\n|---|" + "---:|" * len(names))

    def line(label: str, fn) -> None:
        out.append(f"| {label} | " + " | ".join(fn(res[n]) for n in names) + " |")

    line("MiniWorld pool: protein_rna rows", lambda r: f"{r['mw_rows'].get('protein_rna', 0):,}")
    line("MiniWorld pool: protein_rna clusters", lambda r: f"{r['mw_clusters'].get('protein_rna', 0):,}")
    line("MiniWorld: probability of the class", lambda r: f"{r['mw_mass'].get('protein_rna', 0):.2%}")
    line(f"MiniWorld: draws per cluster in {EPOCH:,} draws",
         lambda r: f"{r['mw_mass'].get('protein_rna', 0) * EPOCH / max(r['mw_clusters'].get('protein_rna', 1), 1):.2f}")
    line("AF3 pool: protein_rna rows", lambda r: f"{r['af_rows'].get('protein_rna', 0):,}")
    line("AF3 pool: protein_rna clusters", lambda r: f"{r['af_clusters'].get('protein_rna', 0):,}")
    line("AF3: probability of the class", lambda r: f"{r['af_mass'].get('protein_rna', 0):.2%}")
    line(f"AF3: draws per cluster in {EPOCH:,} draws",
         lambda r: f"{r['af_mass'].get('protein_rna', 0) * EPOCH / max(r['af_clusters'].get('protein_rna', 1), 1):.2f}")
    out.append("")

    out.append("## Whole pool\n")
    out.append("| | " + " | ".join(names) + " |\n|---|" + "---:|" * len(names))
    line("chain clusters", lambda r: f"{r['n_chain_clusters']:,}")
    line("AF3 pool rows", lambda r: f"{r['af_pool_rows']:,}")
    line("MiniWorld pool rows", lambda r: f"{r['mw_pool_rows']:,}")
    line("effective chain clusters, MiniWorld / AF3", lambda r: f"{r['mw_eff_clusters']:,.0f} / {r['af_eff_clusters']:,.0f}")
    line("effective entries, MiniWorld / AF3", lambda r: f"{r['mw_eff_entries']:,.0f} / {r['af_eff_entries']:,.0f}")
    out.append("")
    args.out.write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
