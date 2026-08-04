from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from datacooker.lmdb import read_lmdb as _read_lmdb

from structcooker.instructions.transforms.codecs import from_bytes

if TYPE_CHECKING:
    from datacooker.protocols import DeserializeFunc


def read_lmdb(env_path: Path, key: str) -> dict[str, Any]:
    """Read a StructCooker LMDB entry through the shared DataCooker utility."""
    return _read_lmdb(env_path, key, deserializer=cast("DeserializeFunc", from_bytes))
