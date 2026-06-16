import gzip
import io
import os
from pathlib import Path
from typing import Any

import zstandard as zstd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict as mmcif2dict  # noqa: N813

# Optional OOM guard: skip CIF files whose compressed size exceeds this many MB
# (a cheap proxy for parse memory). Set per-tier via the SLURM job env so a rare
# monster raises a catchable error and is skipped instead of OOM-killing the
# worker. Unset / <= 0 disables the guard.
_MAX_GZIP_MB_ENV = "DC_CIF_MAX_GZIP_MB"


def dot_transform(key: str) -> list[str]:
    """Transform a dot-separated key into a list of keys."""
    return key.split(".")


def get_cif_data(cif_path: Path) -> dict[str, Any]:
    """Parse a CIF file and return its data as a dictionary."""
    max_mb = float(os.environ.get(_MAX_GZIP_MB_ENV, "0") or "0")
    if max_mb > 0 and cif_path.stat().st_size > max_mb * 1024 * 1024:
        msg = f"{cif_path.name} exceeds {max_mb} MB gzip cap; skipped to avoid OOM."
        raise ValueError(msg)
    if cif_path.suffix == ".gz":
        with gzip.open(cif_path, "rt") as f:
            cif_raw_data = mmcif2dict(f)
    elif cif_path.suffix == ".zst":
        with cif_path.open("rb") as f:
            dctx = zstd.ZstdDecompressor()
            with dctx.stream_reader(f) as reader:
                text_reader = io.TextIOWrapper(reader, encoding="utf-8")
                cif_raw_data = mmcif2dict(text_reader)
    elif cif_path.suffix == ".cif":
        cif_raw_data = mmcif2dict(cif_path)
    else:
        msg = f"Unsupported file format: {cif_path}"
        raise ValueError(msg)
    # Reformat the mmcif2dict output to a more organized structure
    # into a nested dictionary: {key1 : {key2: [values]}}
    organized_dict = {}
    key_list = list(cif_raw_data.keys())
    for key in key_list:
        if "." not in key:
            organized_dict[key] = cif_raw_data[key]
            continue
        main_key, sub_key = key.split(".")
        if main_key not in organized_dict:
            organized_dict[main_key] = {}
        organized_dict[main_key][sub_key] = cif_raw_data[key]
    return organized_dict
