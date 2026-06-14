from __future__ import annotations

from pathlib import Path


def load_fasta(fasta_path: str | Path) -> dict[str, str]:
    """Load a FASTA file into a header-to-sequence mapping."""
    path = Path(fasta_path) if isinstance(fasta_path, str) else fasta_path
    current_header = ""
    fasta_dict: dict[str, str] = {}
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_header = line[1:]
                fasta_dict[current_header] = ""
                continue
            fasta_dict[current_header] += line
    return fasta_dict


def load_seq_id_map(seq_id_map_path: str | Path) -> dict[str, str] | None:
    """Load a sequence-id map from a TSV file if it exists."""
    path = Path(seq_id_map_path) if isinstance(seq_id_map_path, str) else seq_id_map_path
    if not path.exists():
        return None
    seq_id_map: dict[str, str] = {}
    with path.open() as handle:
        for raw_line in handle:
            seq_id, sequence = raw_line.strip().split()
            seq_id_map[f"{seq_id[0]}{sequence}"] = seq_id
    return seq_id_map
