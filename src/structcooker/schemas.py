"""Schema catalogue for BioMol LMDB value types.

Every BioMol DB stores one of a small set of value schemas. This module names
them (A-F) and, for each, records the one place-of-truth for:

* the value STRUCTURE  -- a ``validate()`` that returns issues (empty == valid),
* the KEY convention   -- pdbid / pdbid_chain / seqid / mgyp / ...,
* the CODEC framing    -- standard vs streaming zstd vs raw bytes,
* the memory EXPANSION E (peak live-mem / stored-value-bytes) the planner needs.

db configs reference a schema by name; the structcooker build layer pulls E +
validator from here and hands E to the (domain-agnostic) datacooker planner, so
nothing is hand-tuned per config and every build is verifiable against its
declared schema.

Structures below are the ones observed across ``/data/shared/cssb_data/BioMol``:
  A  {assembly_dict:{<assm>:{atoms,residues,chains,index_table}}, metadata_dict}
  B  {<assm>: {cifmol_attached_dict:{atoms,residues,chains,index_table,metadata}}}
  C  {atoms,residues,chains,index_table,metadata}                (one chain/key)
  D  {template_mols:{<id>:{...}}, [template_ids]}
  E  {msa_dict:{sequences:{...}, headers:{...}}}
  F  raw scalar bytes/str (e.g. a sequence or a seqid)
  H  {atoms,residues,chains,index_table,metadata}   (one assembly-model-altloc/key)
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from typing import Any


class Codec(str, Enum):
    """How a value is framed on disk (drives which decoder ``from_bytes`` uses)."""

    biomol_zstd = "biomol_zstd"                      # standard frame (content-size present)
    biomol_zstd_streaming = "biomol_zstd_streaming"  # streaming frame (no content-size)
    raw = "raw"                                      # bare bytes


# ---- per-schema validators (deserialized value -> list of issues) -----------
def _is_dict(v: Any) -> bool:
    return isinstance(v, dict)


def _v_A(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v):
        return [f"not a dict: {type(v).__name__}"]
    issues = []
    if "assembly_dict" not in v:
        issues.append("missing assembly_dict")
    md = v.get("metadata_dict")
    if not _is_dict(md):
        issues.append("missing metadata_dict")
    elif "id" not in md:
        issues.append("metadata_dict has no id")
    return issues


def _v_B(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v) or not v:
        return ["not a non-empty assembly map"]
    first = next(iter(v.values()))
    if not _is_dict(first) or "cifmol_attached_dict" not in first:
        return ["assembly value lacks cifmol_attached_dict"]
    cad = first["cifmol_attached_dict"]
    return [] if _is_dict(cad) and "atoms" in cad and "metadata" in cad else \
        ["cifmol_attached_dict missing atoms/metadata"]


def _v_C(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v):
        return [f"not a dict: {type(v).__name__}"]
    return [f"missing {k}" for k in ("atoms", "residues", "chains", "index_table") if k not in v]


def _v_D(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v) or "template_mols" not in v:
        return ["missing template_mols"]
    return [] if _is_dict(v["template_mols"]) else ["template_mols not a dict"]


def _v_E(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v) or "msa_dict" not in v:
        return ["missing msa_dict"]
    md = v["msa_dict"]
    return [] if _is_dict(md) and ("sequences" in md or "aligned_sequences" in md) else \
        ["msa_dict missing sequences"]


def _v_F(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    return [] if isinstance(v, (str, bytes, bytearray)) else \
        [f"raw scalar expected, got {type(v).__name__}"]


def _v_G(v: Any) -> list[str]:  # noqa: N802 (schema letter, referenced by name)
    if not _is_dict(v) or "chem_comp_dict" not in v:
        return ["missing chem_comp_dict"]
    return [] if _is_dict(v["chem_comp_dict"]) else ["chem_comp_dict not a dict"]


@dataclass(frozen=True)
class Schema:
    """Place-of-truth record for one BioMol LMDB value schema."""

    name: str                 # "A".."F"
    title: str
    key_convention: str       # pdbid | pdbid_chain | seqid | mgyp | mgyp_topn | mixed
    codec: Codec
    expansion: float | None   # E = peak live mem / stored value bytes (None = unmeasured)
    validate: Callable[[Any], list[str]]


SCHEMAS: dict[str, Schema] = {
    "A": Schema("A", "assembly + metadata_dict", "pdbid",
                Codec.biomol_zstd, 100.0, _v_A),
    "B": Schema("B", "cifmol_attached (per assembly)", "pdbid",
                Codec.biomol_zstd, 100.0, _v_B),
    "C": Schema("C", "chain BioMol", "pdbid_chain",
                Codec.biomol_zstd, 100.0, _v_C),
    # D/E carry conservative E estimates, not a precise measurement like CIFMol's
    # 100x: their values are small (template ~KB, msa ~MB) with no giant tail, so
    # the planner's size-tiering keeps them memory-safe regardless -- refine only
    # if a specific DB shows memory pressure.
    "D": Schema("D", "template_mols", "pdbid_chain",
                Codec.biomol_zstd, 20.0, _v_D),
    # E was 40, but that under-sized MSA rebuilds: depth-capping materializes a
    # dense depth x length array whose LIVE memory (~1-2 GB/record, measured) is
    # unrelated to the compressed stored bytes, so at n_jobs=112 the mem the planner
    # requested (209 g) was crossed and shards OOM-Killed. Worst-case measured
    # E ~= 1.9 GB / 17.9 MB ~= 104; 100 makes the planner request ~326 g (safe).
    "E": Schema("E", "msa_dict", "seqid",
                Codec.biomol_zstd_streaming, 100.0, _v_E),  # codec: standardize (see below)
    "F": Schema("F", "raw scalar", "mixed",
                Codec.raw, 1.0, _v_F),
    # CCD chem-component reference DB: small per-component dicts (~1-2 KB), no
    # giant tail, so a low conservative E is safe for the planner.
    "G": Schema("G", "CCD chem component", "comp_id",
                Codec.biomol_zstd, 5.0, _v_G),
    # Model-input DBs: a CIFMol pruned to the features one model reads, exploded
    # so a key is one training sample. Same value structure as C (hence _v_C),
    # different key convention. Its index_table may be compact -- only the two
    # parent maps -- so read it through
    # structcooker.instructions.transforms.pruning.load_pruned_cifmol.
    # E stays the source's 100.0: a rebuild's memory is driven by the entries it
    # reads, and pruning happens after the source entry is already live.
    "H": Schema("H", "pruned CIFMol (per assembly-model-altloc)",
                "pdbid_assembly_model_altloc", Codec.biomol_zstd, 100.0, _v_C),
}


def get(name: str) -> Schema:
    """Return the schema registered under ``name`` (case-insensitive)."""
    try:
        return SCHEMAS[name.upper()]
    except KeyError as exc:
        msg = f"unknown schema {name!r}; known: {sorted(SCHEMAS)}"
        raise KeyError(msg) from exc


def expansion(name: str, default: float = 100.0) -> float:
    """Return the schema's memory-expansion factor, or ``default`` if unmeasured."""
    e = get(name).expansion
    return default if e is None else e
