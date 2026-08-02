"""Lightweight metadata patch: inject the initial PDB *release* date into an
already-built cif LMDB's stored metadata, WITHOUT reconstructing any CIFMol.

Used by ``datacooker lmdb rebuild`` to turn the deposition-date-only metadata of
``cif_pdb.lmdb`` / ``cif_pdb_attached.lmdb`` into release-date-aware metadata.
The release date is looked up from a ``pdbid -> release_date`` table (built by
``extract_release_dates.py``) keyed on the ``id`` already present in the stored
metadata, so no LMDB key access is needed.

Adapters rename the raw stored keys so the recipe can re-emit the canonical
target names (a recipe step whose output name equals an input name is skipped by
DataCooker, leaving the value unpatched -- hence the rename).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


# --------------------------- metadata-recipe loader --------------------------
def load_release_map(release_table_path: str | Path) -> dict[str, str]:
    """Load ``pdbid -> release_date`` from the extraction TSV (header-aware)."""
    path = Path(release_table_path)
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        pid_i = col.get("pdbid", 0)
        rel_i = col.get("release_date", 1)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rel = parts[rel_i] if rel_i < len(parts) else ""
            if rel:
                out[parts[pid_i].lower()] = rel
    return out


def _release_for(metadata: dict[str, Any], release_map: dict[str, str]) -> str | None:
    ids = metadata.get("id")
    pdbid = (ids[0] if isinstance(ids, list) else ids)
    if pdbid is None:
        return None
    return release_map.get(str(pdbid).lower())


# ------------------------------- adapters ------------------------------------
def adapt_cif_raw(data: dict[str, Any]) -> dict[str, Any]:
    """cif_pdb.lmdb value -> renamed inputs (avoid output/input name collision)."""
    return {
        "_assembly_dict": data["assembly_dict"],
        "_metadata_dict": data["metadata_dict"],
    }


def adapt_cif_attached(data: dict[str, Any]) -> dict[str, Any]:
    """cif_pdb_attached.lmdb value {akey: {cifmol_attached_dict: X}} ->
    {akey: {_cad: X}} so split_entries feeds each inner as ``_cad``."""
    return {akey: {"_cad": inner["cifmol_attached_dict"]} for akey, inner in data.items()}


# ------------------------------ instructions ---------------------------------
def passthrough(value: Any) -> Any:
    return value


def inject_release_metadata_dict(
    metadata_dict: dict[str, Any], release_map: dict[str, str]
) -> dict[str, Any]:
    """cif_pdb.lmdb: add release_date to the top-level metadata_dict."""
    md = dict(metadata_dict)
    md["release_date"] = _release_for(md, release_map)
    return md


def inject_release_cad(
    cad: dict[str, Any], release_map: dict[str, str]
) -> dict[str, Any]:
    """cif_pdb_attached.lmdb: add release_date to one assembly's metadata."""
    out = dict(cad)
    meta = dict(out["metadata"])
    meta["release_date"] = _release_for(meta, release_map)
    out["metadata"] = meta
    return out


def keep_cad_if_release_in_range(
    cad: dict[str, Any],
    start_date: str | None = None,
    end_date: str | None = None,
) -> dict[str, Any] | None:
    """Cheap date pre-filter (NO CIFMol reconstruction): keep one assembly only
    if its stored release_date is in [start, end). Used to build a date-scoped
    subset so the expensive train filter never reconstructs out-of-range giants.
    Falls back to deposition_date if release_date is missing.
    """
    from datetime import date

    meta = cad.get("metadata", {})
    ds = meta.get("release_date") or meta.get("deposition_date")
    if not ds:
        return None
    try:
        d = date.fromisoformat(str(ds))
    except (ValueError, TypeError):
        return None
    lo = date.fromisoformat(start_date) if start_date else date(1900, 1, 1)
    hi = date.fromisoformat(end_date) if end_date else date(5099, 1, 1)
    return cad if (lo <= d < hi) else None
