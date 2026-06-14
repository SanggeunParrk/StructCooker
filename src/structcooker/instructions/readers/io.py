from pathlib import Path

from biomol.cif import CIFMol
from biomol.core.types import BioMolDict
from biomol.core.utils import load_bytes
from datacooker.lmdb import (
    read_lmdb_raw,
)


def load_raw_data(key: str, env_path: Path) -> bytes | None:
    """Read a raw LMDB payload by key."""
    return read_lmdb_raw(env_path, key)


def load_cif(key: str, env_path: Path) -> dict[str, dict[str, CIFMol]]:
    """Read a value from the LMDB database by key.

    Args:
        env_path: Path to the LMDB environment.
        key: Key of the data to retrieve.

    Returns
    -------
        dict
            The data dictionary retrieved from the LMDB database.

    """
    raw_data = load_raw_data(key, env_path)
    if raw_data is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)
    value = load_bytes(raw_data)
    value, metadata = value["assembly_dict"], value["metadata_dict"]

    cifmol_dict: dict[str, dict[str, CIFMol]] = {}
    for cif_key, _item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")

        md = dict(metadata)
        md["assembly_id"] = assembly_id
        md["model_id"] = model_id
        md["alt_id"] = alt_id

        item = BioMolDict(_item)
        item["metadata"] = md

        cifmol_dict[cif_key] = {"cifmol": CIFMol.from_dict(item)}

    return cifmol_dict


__all__ = ["load_bytes", "load_cif", "load_raw_data"]
