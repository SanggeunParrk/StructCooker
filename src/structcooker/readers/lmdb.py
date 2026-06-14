from pathlib import Path
from typing import Any

from datacooker.lmdb import read_lmdb as _read_lmdb

from structcooker.codecs import from_bytes


def read_lmdb(env_path: Path, key: str) -> dict[str, Any]:
    """Read a StructCooker LMDB entry through the shared DataCooker utility."""
    return _read_lmdb(env_path, key, deserialize=from_bytes)
