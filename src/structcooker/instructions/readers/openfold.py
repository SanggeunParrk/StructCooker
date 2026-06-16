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

from pathlib import Path
from typing import Any

import numpy as np


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


def get_openfold_msa_data(alignment_path: Path) -> dict[str, Any]:
    """Load an ``alignment.npz`` file into its per-source MSA payloads.

    Protein monomers expose a single MSA source (``mmseqs_colabfold`` or
    ``concat_cfdb_uniref100_filtered``); RNA monomers expose several
    (``rnacentral_hits``, ``nt_hits``, ``rfam_hits``). Every source shares the
    same inner schema (``msa`` / ``deletion_matrix`` / ``metadata``), so all of
    them are returned uniformly keyed by source name.
    """
    with np.load(alignment_path, allow_pickle=True) as handle:
        msa_sources = {source: handle[source].item() for source in handle.files}
    return {"msa_sources": msa_sources}


def get_openfold_structure_data(structure_path: Path) -> dict[str, Any]:
    """Load a ``structure.npz`` file into its flat per-atom feature arrays.

    The arrays describe a single predicted monomer structure
    (``coord`` / ``atom_name`` / ``res_id`` / ``chain_id`` / ...) and are
    shared across protein and RNA distillation sets.
    """
    with np.load(structure_path, allow_pickle=True) as handle:
        atom_arrays = {field: handle[field] for field in handle.files}
    return {"atom_arrays": atom_arrays}


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
    structure_path = template_path.parent / "structure.npz"
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
