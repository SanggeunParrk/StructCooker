import io
from pathlib import Path
from typing import Any

import zstandard as zstd
from biomol.core import FeatureContainer


def get_a3m_data(a3m_path: Path) -> dict[str, Any]:
    r"""Parse an a3m file (plain, or zstd-compressed ``.zst``) into headers and sequences.

    ``#`` lines are skipped: ColabFold-style a3m -- the format the AFDB complex release ships
    its MSAs in -- opens with ``#<lengths>\t<cardinalities>``, which is metadata, not a
    record. HHblits and MMseqs2 ``result2msa`` output has no such line, so those parse as
    before.
    """
    raw_sequences = []
    headers = []
    with a3m_path.open("rb") as raw:
        stream = (
            zstd.ZstdDecompressor().stream_reader(raw) if a3m_path.suffix == ".zst" else raw
        )
        handle = io.TextIOWrapper(stream, encoding="utf-8")
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
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
    value = value["msa_dict"]
    sequences, headers = (
        value["sequences"],
        value["headers"],
    )
    return {
        "msa_dict": {
            "_sequences": sequences,
            "_headers": headers,
        },
    }
