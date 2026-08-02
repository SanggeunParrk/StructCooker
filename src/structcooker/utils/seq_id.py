from __future__ import annotations

import re
from pathlib import Path


SEQ_ID_RE = re.compile(r"^([A-Za-z])(\d{20})(?:\..*)?$")


def seq_id_from_name(name: str) -> str | None:
    match = SEQ_ID_RE.match(name)
    if match is None:
        return None
    return f"{match.group(1)}{match.group(2)}"


def seq_id_shard_path(base_dir: Path, seq_id: str) -> Path:
    if not re.match(r"^[A-Za-z]\d{20}$", seq_id):
        msg = f"Expected type-prefixed 20-digit sequence id, got {seq_id!r}"
        raise ValueError(msg)
    return base_dir / seq_id[0] / seq_id[-3:]
