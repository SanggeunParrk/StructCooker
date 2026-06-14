from pathlib import Path


def load_fasta(fasta_path: str | Path) -> dict[str, str]:
    """Load fasta file into a dictionary."""
    if isinstance(fasta_path, str):
        fasta_path = Path(fasta_path)
    with fasta_path.open("r") as f:
        lines = f.readlines()
    current_header = ""
    fasta_dict = {}
    for _line in lines:
        line = _line.strip()
        if line.startswith(">"):
            current_header = line[1:]
            fasta_dict[current_header] = ""
        else:
            fasta_dict[current_header] += line
    return fasta_dict


def load_seq_id_map(seq_id_map_path: str | Path) -> dict[str, str] | None:
    """Load sequence ID map from a given path."""
    if isinstance(seq_id_map_path, str):
        seq_id_map_path = Path(seq_id_map_path)
    if not seq_id_map_path.exists():
        return None
    seq_id_map = {}
    with seq_id_map_path.open("r") as f:
        for line in f:
            seq_id, sequence = line.strip().split()
            seq_id_map[f"{seq_id[0]}{sequence}"] = seq_id
    return seq_id_map


def parse_sabdab_data(
    sabdab_summary_path: str | Path,
    fasta_dict: dict[str, str],
) -> tuple[list[str], list[str]]:
    ab_H_chain_list, ab_L_chain_list = [], []
    with sabdab_summary_path.open("r") as f:
        lines = f.readlines()
    for _line in lines[1:]:
        line = _line.strip()
        items = line.split("\t")
        pdb_ID, ab_H_chain, ab_L_chain = items[0], items[1], items[2]
        ab_H_chain = f"{pdb_ID}_{ab_H_chain}".upper() if ab_H_chain != "NA" else None
        ab_L_chain = f"{pdb_ID}_{ab_L_chain}".upper() if ab_L_chain != "NA" else None
        if ab_H_chain in chain_ID_to_seq:
            seq_hash = seq_to_hash[chain_ID_to_seq[ab_H_chain]]
            ab_H_chain_list.append(seq_hash)
        if ab_L_chain in chain_ID_to_seq:
            seq_hash = seq_to_hash[chain_ID_to_seq[ab_L_chain]]
            ab_L_chain_list.append(seq_hash)
