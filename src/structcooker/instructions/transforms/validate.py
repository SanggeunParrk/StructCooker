"""Post-build validators for the ingest LMDBs (msa / cif / template).

Each validator takes one decoded LMDB record and returns a list of issue dicts
(empty == valid), mirroring :func:`validate_chem_comp`. They check field
presence, dtype kinds and shape consistency against the schema of the reference
BioMolDB LMDBs (``a3m.lmdb`` / ``cif.lmdb`` / ``template.lmdb``), so a built DB
can be QC'd record-by-record after the build.

dtype *kinds* (``np.dtype.kind``) are checked rather than exact dtypes, because
byte-string widths (``|S<n>``) and unicode widths (``<U<n>``) are data-dependent
and legitimately vary per record.
"""

from typing import Any

import numpy as np

# Expected per-level node fields (BioMol container schema).
_ATOM_NODE_FIELDS = (
    "id", "element", "aromatic", "stereo", "charge",
    "model_xyz", "xyz", "b_factor", "occupancy",
)
_RESIDUE_NODE_FIELDS = (
    "id", "formula", "cif_idx", "auth_idx", "chem_comp_id",
    "hetero", "one_letter_code_can", "one_letter_code",
)
_CHAIN_NODE_FIELDS = ("entity_id", "entity_type", "auth_asym_id", "chain_id")
# Template BioMols (TemplateMol) carry only coordinates + residue/chain labels,
# no bond edges or per-atom chemistry, plus richer provenance metadata.
_TEMPLATE_ATOM_NODE_FIELDS = ("id", "xyz", "b_factor", "occupancy")
_TEMPLATE_RESIDUE_NODE_FIELDS = (
    "one_letter_code_can", "one_letter_code", "cif_idx", "auth_idx", "chem_comp_id", "hetero",
)
_TEMPLATE_METADATA_KEYS = (
    "id", "deposition_date", "resolution", "assembly_id", "model_id", "alt_id",
)
_INDEX_TABLE_KEYS = (
    "atom_to_res", "res_to_chain", "res_atom_indptr",
    "res_atom_indices", "chain_res_indptr", "chain_res_indices",
)
# mmCIF entity_poly.type vocabulary used by the reference cif/template LMDBs.
_ENTITY_TYPES = {
    "polypeptide(L)", "polypeptide(D)", "polyribonucleotide", "polydeoxyribonucleotide",
    "polydeoxyribonucleotide/polyribonucleotide hybrid", "polysaccharide(D)",
    "polysaccharide(L)", "cyclic-pseudo-peptide", "peptide nucleic acid", "other",
    "non-polymer", "branched", "water", "macrolide",
}


def _issue(code: str, message: str, **detail: object) -> dict[str, Any]:
    """Build one issue record (empty list of these == valid)."""
    return {"code": code, "message": message, "detail": detail}


def _node_value(container: dict, field: str) -> np.ndarray | None:
    """Return ``container['nodes'][field]['value']`` or ``None`` if absent."""
    node = container.get("nodes", {}).get(field)
    if isinstance(node, dict):
        return node.get("value")
    return None


def _check_array(
    issues: list[dict],
    value: object,
    name: str,
    *,
    kind: str | None = None,
    ndim: int | None = None,
) -> np.ndarray | None:
    """Append issues if ``value`` is not an ndarray of the expected kind/ndim."""
    if value is None:
        issues.append(_issue("missing_field", f"{name} is missing"))
        return None
    if not isinstance(value, np.ndarray):
        issues.append(_issue("wrong_type", f"{name} is not an ndarray", got=type(value).__name__))
        return None
    if kind is not None and value.dtype.kind not in kind:
        issues.append(_issue("wrong_dtype", f"{name} dtype {value.dtype} (kind not in '{kind}')"))
    if ndim is not None and value.ndim != ndim:
        issues.append(_issue("wrong_ndim", f"{name} ndim {value.ndim} != {ndim}"))
    return value


def validate_msa_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate one ``msa_dict`` record against the a3m LMDB schema."""
    issues: list[dict[str, Any]] = []
    payload = record.get("msa_dict", record)
    seq = payload.get("sequences")
    hdr = payload.get("headers")
    if not isinstance(seq, dict) or not isinstance(hdr, dict):
        issues.append(_issue("missing_field", "msa_dict.sequences/headers missing"))
        return issues

    query = _check_array(issues, seq.get("query_sequence"), "sequences.query_sequence", kind="U", ndim=1)
    aligned = _check_array(issues, seq.get("aligned_sequences"), "sequences.aligned_sequences", kind="u", ndim=2)
    deletions = _check_array(issues, seq.get("deletions"), "sequences.deletions", kind="iu", ndim=2)
    _check_array(issues, seq.get("deletion_mean"), "sequences.deletion_mean", kind="f", ndim=1)
    profile = _check_array(issues, seq.get("profile"), "sequences.profile", kind="f", ndim=2)

    length = len(query) if query is not None else None
    depth = aligned.shape[0] if aligned is not None else None
    if aligned is not None and length is not None and aligned.shape[1] != length:
        issues.append(_issue("shape_mismatch", f"aligned width {aligned.shape[1]} != query length {length}"))
    if aligned is not None and deletions is not None and aligned.shape != deletions.shape:
        issues.append(_issue("shape_mismatch", f"deletions {deletions.shape} != aligned {aligned.shape}"))
    if profile is not None and profile.shape[-1] != 32:  # noqa: PLR2004 - profile is (L, 32)
        issues.append(_issue("shape_mismatch", f"profile last dim {profile.shape[-1]} != 32"))
    if depth is not None and depth < 1:
        issues.append(_issue("empty", "alignment has no rows"))

    for field in ("database", "database_id", "species", "rep_id"):
        col = _check_array(issues, hdr.get(field), f"headers.{field}", kind="S", ndim=1)
        if col is not None and depth is not None and len(col) != depth:
            issues.append(_issue("shape_mismatch", f"headers.{field} len {len(col)} != depth {depth}"))
    return issues


def _validate_biomol(  # noqa: PLR0913 - per-level schemas passed explicitly
    bm: dict[str, Any],
    issues: list[dict[str, Any]],
    *,
    prefix: str,
    atom_fields: tuple[str, ...],
    residue_fields: tuple[str, ...],
    chain_fields: tuple[str, ...] = _CHAIN_NODE_FIELDS,
    metadata_keys: tuple[str, ...] = (),
) -> None:
    """Validate one BioMol payload (shared by cif and template records)."""
    issues.extend(
        _issue("missing_field", f"{prefix}{container} missing")
        for container in ("atoms", "residues", "chains", "index_table")
        if container not in bm
    )

    atoms = bm.get("atoms", {})
    residues = bm.get("residues", {})
    chains = bm.get("chains", {})
    index_table = bm.get("index_table", {})

    for level, fields in (
        (atoms, atom_fields),
        (residues, residue_fields),
        (chains, chain_fields),
    ):
        issues.extend(
            _issue("missing_field", f"{prefix}node field '{field}' missing")
            for field in fields
            if _node_value(level, field) is None
        )

    issues.extend(
        _issue("missing_field", f"{prefix}index_table.{key} missing")
        for key in _INDEX_TABLE_KEYS
        if key not in index_table
    )

    # Shape consistency across the hierarchy (chem_comp_id labels residues in
    # both the cif and template schemas; cif also has residues.id).
    atom_ids = _node_value(atoms, "id")
    res_ids = _node_value(residues, "chem_comp_id")
    chain_ids = _node_value(chains, "chain_id")
    xyz = _node_value(atoms, "xyz")
    n_atoms = len(atom_ids) if atom_ids is not None else None
    n_res = len(res_ids) if res_ids is not None else None
    n_chain = len(chain_ids) if chain_ids is not None else None

    if xyz is not None and n_atoms is not None and getattr(xyz, "shape", (0,))[0] != n_atoms:
        issues.append(_issue("shape_mismatch", f"{prefix}xyz rows {xyz.shape[0]} != n_atoms {n_atoms}"))
    if xyz is not None and (xyz.ndim != 2 or xyz.shape[1] != 3):  # noqa: PLR2004 - xyz is (n, 3)
        issues.append(_issue("shape_mismatch", f"{prefix}xyz shape {xyz.shape} != (n, 3)"))

    a2r = index_table.get("atom_to_res")
    r2c = index_table.get("res_to_chain")
    if a2r is not None and n_atoms is not None and len(a2r) != n_atoms:
        issues.append(_issue("shape_mismatch", f"{prefix}atom_to_res len {len(a2r)} != n_atoms {n_atoms}"))
    if r2c is not None and n_res is not None and len(r2c) != n_res:
        issues.append(_issue("shape_mismatch", f"{prefix}res_to_chain len {len(r2c)} != n_res {n_res}"))
    if a2r is not None and n_res is not None and len(a2r) and max(a2r) >= n_res:
        issues.append(_issue("index_out_of_range", f"{prefix}atom_to_res max {max(a2r)} >= n_res {n_res}"))
    if r2c is not None and n_chain is not None and len(r2c) and max(r2c) >= n_chain:
        issues.append(_issue("index_out_of_range", f"{prefix}res_to_chain max {max(r2c)} >= n_chain {n_chain}"))

    entity_type = _node_value(chains, "entity_type")
    if entity_type is not None:
        bad = {str(t) for t in entity_type.tolist()} - _ENTITY_TYPES
        if bad:
            issues.append(_issue("bad_value", f"{prefix}entity_type has unexpected {sorted(bad)}"))

    metadata = bm.get("metadata")
    issues.extend(
        _issue("missing_field", f"{prefix}metadata.{key} missing")
        for key in metadata_keys
        if not isinstance(metadata, dict) or key not in metadata
    )


def validate_cif_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate one ``cifmol_dict`` (BioMol) record against the cif LMDB schema."""
    issues: list[dict[str, Any]] = []
    bm = record.get("cifmol_dict", record)
    if not isinstance(bm, dict):
        return [_issue("missing_field", "cifmol_dict missing")]
    _validate_biomol(
        bm, issues, prefix="",
        atom_fields=_ATOM_NODE_FIELDS,
        residue_fields=_RESIDUE_NODE_FIELDS,
    )
    return issues


def validate_template_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate one ``template_mols`` record (a mapping of hit id -> BioMol)."""
    issues: list[dict[str, Any]] = []
    mols = record.get("template_mols", record)
    if not isinstance(mols, dict):
        return [_issue("missing_field", "template_mols missing")]
    # An empty hit set is valid: a query may legitimately have no templates.
    for hit_id, bm in mols.items():
        if not isinstance(bm, dict):
            issues.append(_issue("wrong_type", f"hit '{hit_id}' is not a BioMol dict"))
            continue
        _validate_biomol(
            bm, issues, prefix=f"{hit_id}.",
            atom_fields=_TEMPLATE_ATOM_NODE_FIELDS,
            residue_fields=_TEMPLATE_RESIDUE_NODE_FIELDS,
            metadata_keys=_TEMPLATE_METADATA_KEYS,
        )
    return issues
