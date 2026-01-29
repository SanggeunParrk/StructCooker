from pathlib import Path


def write_fasta(data: dict[str, dict[str, dict[str, str]]], output_path: Path) -> None:
    """Write fasta dictionary to a fasta file."""
    for cif_key, fasta_dict in data.items():
        for chain_id, fasta_entry in fasta_dict["fasta"].items():
            fasta_file_path = (
                output_path / f"{cif_key[1:3]}" / f"{cif_key}_{chain_id}.fasta"
            )
            fasta_file_path.parent.mkdir(parents=True, exist_ok=True)
            with fasta_file_path.open("w") as f:
                f.write(fasta_entry)
