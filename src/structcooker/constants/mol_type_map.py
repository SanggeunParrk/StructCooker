mol_type_map: dict[str, str] = {
    "polypeptide(L)": "P",
    "polypeptide(D)": "Q",
    "polydeoxyribonucleotide": "D",
    "polyribonucleotide": "R",
    "polydeoxyribonucleotide/polyribonucleotide hybrid": "N",
    "branched": "B",
    "non-polymer": "L",
    "unknown": "X",
}

cluster_types: set[str] = set({"P", "Q", "D", "R", "N", "A", "B", "L", "X"})
cluster_maps: dict[str, str] = {
    "P": "Protein(L)",
    "Q": "Protein(D)",
    "A": "Antibody",
    "D": "DNA",
    "R": "RNA",
    "N": "DNA/RNA hybrid",
    "B": "Branched",
    "L": "Non-polymer",
    "X": "Unknown",
}
polymer_cluster_types: set[str] = set({"P", "Q", "D", "R", "N", "A"})
