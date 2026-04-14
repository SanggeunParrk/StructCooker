from pathlib import Path
from typing import Any

from biomol.core import FeatureContainer


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


def convert_to_msa_container(
    value: dict,
) -> dict[str, dict[str, FeatureContainer]]:
    """Convert a dictionary containing CIFMol data into a dictionary of CIFMol objects."""
    value = value["msa_container"]
    residue_container, chain_container = (
        value["residue_container"],
        value["chain_container"],
    )
    residue_container = FeatureContainer.from_dict(residue_container)
    chain_container = FeatureContainer.from_dict(chain_container)
    return {
        "msa_container": {
            "_residue_container": residue_container,
            "_chain_container": chain_container,
        },
    }
