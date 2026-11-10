"""One ``cif_pdb`` entry -> its unit records (the ``unit/unit_pdb`` rebuild).

The stage is an exploding rebuild (``split_entries`` + ``explode_entries``): the engine writes one
record per sub-entry of the entry the reader hands to the recipe, keyed ``<pdb id>_<sub-entry>``.
Two functions split the work around that hook:

* :func:`explode_units` is the reader adapter. It sees the whole entry, which the selection needs
  (first model of each assembly, units deduplicated across assemblies), finds every unit and
  returns ``{"<assembly>_<model>_<index>": {"unit_input": UnitInput}}``: the sub-entries.
* :func:`build_unit_record` is the recipe step, run once per sub-entry: ``UnitInput`` -> the stored
  record (:func:`structcooker.instructions.unit.record.unit_to_record`).

Settings are ``params/unit.yaml`` (:func:`default_config`); the ANARCI cache is the TSV named
by ``$UNIT_ANTIBODY_CACHE``, which the stage config declares in ``runtime_environment`` (so the
release pins it and fingerprints it). :func:`entry_records` runs both halves for tests and checks.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from structcooker.instructions.unit.antibody import AntibodyClassifier
from structcooker.instructions.unit.config import Config, load_config
from structcooker.instructions.unit.pipeline import extract_units
from structcooker.instructions.unit.record import asm_fields, unit_to_record
from structcooker.instructions.mols import CIFMol

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from biomol.core.types import BioMolDict

    from structcooker.instructions.unit.adapt import Assembly
    from structcooker.instructions.unit.detect import Unit

ANTIBODY_CACHE_ENV = "UNIT_ANTIBODY_CACHE"


@dataclass
class UnitInput:
    """Everything :func:`build_unit_record` needs to write one unit (a sub-entry of an entry)."""

    asm: Assembly
    fields: dict[str, Any]
    unit: Unit
    duplex: dict[tuple[int, int], int]
    meta: dict[str, Any]
    ab_class: list[str | None]
    index: int


@lru_cache(maxsize=1)
def default_config() -> Config:
    """Return the thresholds of ``params/unit.yaml`` (loaded once per process)."""
    return load_config()


@lru_cache(maxsize=1)
def default_classifier() -> AntibodyClassifier:
    """ANARCI classifier over the cache TSV in ``$UNIT_ANTIBODY_CACHE`` (misses fall back to ANARCI)."""
    return AntibodyClassifier(os.environ.get(ANTIBODY_CACHE_ENV))


def select_assembly_keys(keys: Iterable[str]) -> list[str]:
    """First model per assembly; the no-altloc (``.``) entry if present, else the first altloc.

    ``keys`` are ``<assembly>_<model>_<altloc>`` strings; the result is ordered by assembly id.
    """
    best: dict[str, tuple[str, str]] = {}
    for k in sorted(keys):
        a, m, alt = k.split("_")
        cur = best.get(a)
        if cur is None or int(m) < int(cur[1]) or (m == cur[1] and alt == "."):
            best[a] = (k, m)
    return [v[0] for _, v in sorted(best.items(), key=lambda kv: kv[0])]


def entry_pdb_id(metadata: dict[str, Any]) -> str:
    """Lower-case PDB id of an entry (``metadata_dict['id']`` is ``['1TSR']`` or ``'1TSR'``)."""
    value = metadata["id"]
    if not isinstance(value, str):
        value = value[0]
    return value.lower()


def _resolution(metadata: dict[str, Any]) -> float | None:
    res = metadata.get("resolution")
    if res is None or not np.isfinite(float(res)):
        return None
    return float(res)


def entry_units(
    entry: dict[str, Any],
    cfg: Config,
    ab_lookup: Callable[[str], str | None] | None,
) -> dict[str, dict[str, UnitInput]]:
    """Every unit of one decoded ``cif_pdb`` entry, as sub-entries ``{"<asm>_<model>_<index>": {...}}``.

    Per assembly: its first model (altloc ``.`` preferred); an assembly with more polymer chains than
    ``assembly.max_polymer_chains`` is skipped; a unit equal in node sequences and (bucketed) contact
    size to one already kept in this entry is dropped (symmetry copies, repeated assemblies).
    """
    md = entry["metadata_dict"]
    base_meta = {
        "pdb_id": entry_pdb_id(md),
        "resolution": _resolution(md),
        "release_date": md.get("release_date"),
        "deposition_date": md.get("deposition_date"),
    }
    seen: set[tuple[tuple[str, ...], int]] = set()
    out: dict[str, dict[str, UnitInput]] = {}
    for ck in select_assembly_keys(entry["assembly_dict"]):
        aid, mid, alt = ck.split("_")
        item = dict(entry["assembly_dict"][ck])
        item["metadata"] = {**md, "assembly_id": aid, "model_id": mid, "alt_id": alt}
        cm = CIFMol.from_dict(cast("BioMolDict", item))
        asm, units, duplex = extract_units(cm, cfg, ab_lookup)
        n_polymer = sum(k in ("protein", "na") for k in asm.kinds)
        if n_polymer > cfg.unit.max_polymer_chains_in_assembly or not units:
            continue
        fields = asm_fields(cm, asm)
        idx = 0
        for u in units:
            sig = (tuple(sorted(asm.keys[x] for x in u.nodes)), round(u.weight / 10))
            if sig in seen:
                continue
            seen.add(sig)
            meta = {**base_meta, "assembly_id": aid, "model_id": mid, "alt_id": alt}
            out[f"{aid}_{mid}_{idx}"] = {
                "unit_input": UnitInput(
                    asm, fields, u, duplex, meta, asm.stats["ab_class"], idx
                )
            }
            idx += 1
    return out


def explode_units(entry: dict[str, Any]) -> dict[str, dict[str, UnitInput]]:
    """Reader adapter of the rebuild: a decoded ``cif_pdb`` entry -> its unit sub-entries."""
    return entry_units(entry, default_config(), default_classifier())


def build_unit_record(
    unit_input: UnitInput, cfg: Config | None = None
) -> dict[str, Any]:
    """Recipe step: the stored record of one unit (``record.unit_to_record``)."""
    return unit_to_record(
        unit_input.asm,
        unit_input.fields,
        unit_input.unit,
        unit_input.duplex,
        cfg or default_config(),
        unit_input.meta,
        unit_input.ab_class,
        unit_input.index,
    )


def entry_records(
    entry: dict[str, Any],
    cfg: Config | None = None,
    ab_lookup: Callable[[str], str | None] | None = None,
) -> dict[str, dict[str, Any]]:
    """Return ``{unit id: record}`` of one entry: what the rebuild stores (adapter + recipe step)."""
    cfg = cfg or default_config()
    subs = entry_units(entry, cfg, ab_lookup or default_classifier())
    return {
        f"{u['unit_input'].meta['pdb_id']}_{key}": build_unit_record(
            u["unit_input"], cfg
        )
        for key, u in subs.items()
    }
