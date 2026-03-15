from pathlib import Path
from typing import Any

from Bio.PDB.MMCIF2Dict import MMCIF2Dict as mmcif2dict  # noqa: N813


def get_a3m_data(a3m_path: Path) -> dict[str, Any]:
    """Parse a a3m file and return its data as a dictionary."""
    with a3m_path.open("r") as f:
        raw_lines = f.readlines()
    raw_sequences = []
    headers = []
    for line in raw_lines:
        line = line.strip()
        if line.startswith(">"):
            headers.append(line[1:])
            raw_sequences.append("")
        else:
            if not raw_sequences:
                msg = f"Invalid a3m format: sequence data found before any header in {a3m_path}"
                raise ValueError(msg)
            raw_sequences[-1] += line

    return {
        "raw_sequences": raw_sequences,
        "headers": headers,
    }
