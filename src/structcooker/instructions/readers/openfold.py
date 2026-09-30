"""Readers for the OpenFold3 distillation datasets (monomer / RNA).

Each distillation entry lives in a per-id folder::

    <dataset>/<entry_id>/alignment.npz
    <dataset>/<entry_id>/structure.npz
    <dataset>/<entry_id>/template.npz   # protein monomers only

The entry id is therefore the *parent directory name* rather than the file
name (every entry shares the same ``alignment.npz`` / ``structure.npz`` /
``template.npz`` file names). ``openfold_entry_key`` exposes that id so it can
be wired into a build config via ``key_builder``.
"""

import io
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import zstandard

from structcooker.instructions.readers.teddymer import pdb_atom_table


def openfold_entry_key(path: Path) -> str:
    """Return the distillation entry id (the parent folder name)."""
    return path.parent.name


def openfold_chain_key(path: Path) -> str:
    """Return the file stem as the key.

    The disordered (PDB-derived) alignment / template arrays are stored flat as
    ``<pdbid>_<chain>.npz`` rather than one folder per entry, so the key is the
    file stem instead of the parent folder name.
    """
    return path.stem


_MONOMER_SEQID: dict[str, str] | None = None
_STRUCTURE_OVERRIDES: dict[str, str] | None = None
_MAP_IDENTITIES: dict[str, tuple] = {}


def _map_changed(name: str, path: Path) -> bool:
    """Invalidate local caches on path/version changes; releases pin immutable copies."""
    if os.environ.get("STRUCTCOOKER_IMMUTABLE_REFERENCES") == "1":
        identity = (str(path.absolute()),)
    else:
        stat = path.stat()
        identity = (str(path.resolve()), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    changed = _MAP_IDENTITIES.get(name) != identity
    _MAP_IDENTITIES[name] = identity
    return changed



def _structure_path(path: Path) -> Path:
    """Resolve explicitly recorded, separately recovered structures for this run."""
    global _STRUCTURE_OVERRIDES  # noqa: PLW0603 - process-local cache
    manifest = os.environ.get("OPENFOLD_STRUCTURE_OVERRIDES")
    if not manifest:
        _STRUCTURE_OVERRIDES = None
        return path
    changed = _map_changed("structure", Path(manifest))
    if _STRUCTURE_OVERRIDES is None or changed:
        _STRUCTURE_OVERRIDES = None
        mapping = json.loads(Path(manifest).read_text()) if manifest else {}
        if not isinstance(mapping, dict) or any(
            not isinstance(key, str) or not isinstance(value, str) for key, value in mapping.items()
        ):
            msg = "Structure overrides must map original paths to recovered paths"
            raise ValueError(msg)
        _STRUCTURE_OVERRIDES = mapping
    return Path(_STRUCTURE_OVERRIDES.get(str(path), str(path)))


def _monomer_seqid_map() -> dict[str, str]:
    """Load (and cache) the monomer entry -> seq_id map.

    The MSA depends only on the sequence, so the distillation MSA DBs are keyed by
    seq_id (like the PDB a3m DBs), not by the per-structure entry id. This maps each
    monomer entry (``path.parent.name``: MGYP id for short/long, accession for rna) to
    its seq_id via ``$OUTPUT_ROOT/metadata/distillation_monomer_seqid.tsv`` -- a build
    artifact (distillation fasta sequences resolved through seq_id_map). Structures
    sharing a sequence collapse onto one seq_id, matching production's dedup.
    """
    global _MONOMER_SEQID  # noqa: PLW0603 - process-local cache
    root = os.environ.get("OUTPUT_ROOT", "/data/shared/cssb_data/BioMol_clean")
    map_path = Path(os.environ.get(
        "MONOMER_SEQID_MAP", str(Path(root) / "metadata" / "distillation_monomer_seqid.tsv"),
    ))
    changed = _map_changed("monomer", map_path)
    if _MONOMER_SEQID is None or changed:
        _MONOMER_SEQID = None
        mapping: dict[str, str] = {}
        with map_path.open(encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2:
                    mapping[parts[0]] = parts[1]
        _MONOMER_SEQID = mapping
    return _MONOMER_SEQID


def openfold_seqid_key(path: Path) -> str:
    """Key a monomer (short/long/rna) alignment npz by its sequence's seq_id."""
    return _monomer_seqid_map()[path.parent.name]


_DISORDERED_SEQID: dict[str, str] | None = None


def _disordered_seqid_map() -> dict[str, str]:
    """Load (and cache) the precomputed ``{id}_{chain}`` -> seq_id map.

    The path comes from the ``DISORDERED_SEQID_MAP`` env var (a TSV built by
    ``scripts/maintenance/precompute_disordered_seqid.py``).
    """
    global _DISORDERED_SEQID  # noqa: PLW0603 - process-local cache
    map_path = Path(os.environ["DISORDERED_SEQID_MAP"])
    changed = _map_changed("disordered", map_path)
    if _DISORDERED_SEQID is None or changed:
        _DISORDERED_SEQID = None
        mapping: dict[str, str] = {}
        with map_path.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2:
                    mapping[parts[0]] = parts[1]
        _DISORDERED_SEQID = mapping
    return _DISORDERED_SEQID


def openfold_disordered_seqid_key(path: Path) -> str:
    """Key a disordered chain npz by its chain seq_id.

    Maps the ``{id}_{chain}`` file stem to the chain's seq_id so the disordered
    MSA/template DBs are keyed exactly like the PDB pipeline (per-chain seq_id).
    Homomer chains sharing a sequence collapse onto one seq_id, matching PDB's
    sequence-level deduplication.
    """
    return _disordered_seqid_map()[path.stem]


def get_openfold_msa_data(alignment_path: Path) -> dict[str, Any]:
    """Load an ``alignment.npz`` file into its per-source MSA payloads.

    Protein monomers expose a single MSA source (``mmseqs_colabfold`` or
    ``concat_cfdb_uniref100_filtered``); RNA monomers expose several
    (``rnacentral_hits``, ``nt_hits``, ``rfam_hits``). Every source shares the
    same inner schema (``msa`` / ``deletion_matrix`` / ``metadata``), so all of
    them are returned uniformly keyed by source name.
    """
    # A few distillation structures ship an empty / truncated alignment.npz; production
    # excludes them from the MSA DBs, so surface a clean error and let the build's
    # error_mode="skip" drop them (they never become records).
    if alignment_path.stat().st_size == 0:
        msg = f"empty alignment.npz (no MSA): {alignment_path}"
        raise ValueError(msg)
    with np.load(alignment_path, allow_pickle=True) as handle:
        msa_sources = {source: handle[source].item() for source in handle.files}
    if not msa_sources:
        msg = f"alignment.npz contains no MSA sources: {alignment_path}"
        raise ValueError(msg)
    return {"msa_sources": msa_sources}


def get_openfold_structure_data(structure_path: Path) -> dict[str, Any]:
    """Load a ``structure.npz`` file into its flat per-atom feature arrays.

    The arrays describe a single predicted monomer structure
    (``coord`` / ``atom_name`` / ``res_id`` / ``chain_id`` / ...) and are
    shared across protein and RNA distillation sets.
    """
    with np.load(_structure_path(structure_path), allow_pickle=True) as handle:
        atom_site_dict = {field: handle[field] for field in handle.files}
    return {"atom_site_dict": atom_site_dict, "entry_id": structure_path.parent.name}


def get_disordered_template_data(npz_path: Path) -> dict[str, Any]:
    """Load one disordered template chain (``templates/<id>/<id>_<chain>.npz``).

    Same atom-table schema as ``structure.npz`` (reuses the structure recipe),
    but the entry id is the file stem (``<id>_<chain>``). The per-folder
    ``chain_id_to_moltype.npz`` index is skipped (raises so the runner drops it).
    """
    if npz_path.name == "chain_id_to_moltype.npz":
        msg = "moltype index, not a template chain"
        raise ValueError(msg)
    out = get_openfold_structure_data(npz_path)
    out["entry_id"] = npz_path.stem
    return out


def get_disordered_template_group(moltype_path: Path) -> dict[str, Any]:
    """Load all template chains of one disordered query folder.

    Anchored on the per-folder ``chain_id_to_moltype.npz`` (one per query, so the
    build key is the query id), this reads every sibling ``<id>_<chain>.npz`` atom
    table and returns them keyed by file stem (the template hit id).
    """
    templates: dict[str, Any] = {}
    for npz_path in sorted(moltype_path.parent.glob("*.npz")):
        if npz_path.name == "chain_id_to_moltype.npz":
            continue
        with np.load(npz_path, allow_pickle=True) as handle:
            templates[npz_path.stem] = {field: handle[field] for field in handle.files}
    return {"templates_atom_site_dict": templates}


def get_disordered_template_chain(npz_path: Path) -> dict[str, Any]:
    """Load ONE disordered template chain atom table (per-chain, seq_id-keyed build).

    The disordered template DB is keyed per chain by seq_id (like the PDB
    template DB), so each ``<id>_<chain>.npz`` is read individually. Returns the
    single chain's atom table keyed by its stem so ``build_disordered_template_mols``
    yields ``{template_mols: {"{id}_{chain}": mol}}`` -- matching the PDB layout.
    """
    with np.load(npz_path, allow_pickle=True) as handle:
        arrays = {field: handle[field] for field in handle.files}
    return {"templates_atom_site_dict": {npz_path.stem: arrays}}


def get_openfold_template_data(template_path: Path) -> dict[str, Any]:
    """Load a ``template.npz`` file into its per-hit template payloads.

    Each top-level key is a ``<pdb_id>_<chain_id>`` template hit mapping to an
    ``index`` / ``release_date`` / ``idx_map`` payload. Only protein monomers
    carry templates (RNA monomers have no ``template.npz``). The sibling
    ``structure.npz`` provides the query length needed to align templates over
    the full monomer.
    """
    with np.load(template_path, allow_pickle=True) as handle:
        template_hits = {hit: handle[hit].item() for hit in handle.files}
    query_len = 0
    structure_path = _structure_path(template_path.parent / "structure.npz")
    if structure_path.exists():
        with np.load(structure_path, allow_pickle=True) as handle:
            keys = np.stack(
                [
                    handle["chain_id"].astype(str),
                    handle["res_id"].astype(str),
                    handle["ins_code"].astype(str),
                ],
                axis=1,
            )
            _, first_idx = np.unique(keys, axis=0, return_index=True)
            query_len = len(first_idx)
    return {"template_hits": template_hits, "query_len": query_len}


def get_openfold_pdb_structure_data(pdb_path: Path) -> dict[str, Any]:
    """Load ``raw/<id>/best_structure_relaxed.pdb[.zst]`` -- the model OpenFold released.

    Used instead of ``preprocessed/structure.npz`` where the release ships the raw model
    (short and long monomers): the PDB carries pLDDT in its B-factor column, which the
    preprocessed npz drops. Same heavy atoms and coordinates as the npz, plus the
    C-terminal OXT the npz leaves out; hydrogens (the model is Amber-relaxed) are dropped.
    No connectivity -- the recipe derives bonds from the CCD, as for teddymer.
    """
    path = Path(pdb_path)
    if path.suffix == ".zst":
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as reader:
            lines = io.TextIOWrapper(reader, encoding="utf-8")
            table = pdb_atom_table(lines, name=path.name, drop_hydrogens=True)
    else:
        with path.open() as handle:
            table = pdb_atom_table(handle, name=path.name, drop_hydrogens=True)
    return {"raw_atom_site_dict": table, "entry_id": path.parent.name}

