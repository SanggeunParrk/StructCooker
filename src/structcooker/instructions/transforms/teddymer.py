"""Derive the connectivity a PDB does not carry, so the openfold structure steps apply.

``openfold_structure`` gets ``bonds`` (an (m,2) atom-index array) from the source
``structure.npz`` and consults the CCD only for each bond's *type*. A Teddymer PDB has
no CONECT records, so the index pairs must be built here -- from the same CCD the recipe
already loads -- before ``derive_bond_edges`` runs.

Keeping this in its own step means ``openfold_structure`` is reused verbatim: this turns
the reader's ``raw_atom_site_dict`` into the ``atom_site_dict`` those steps expect.
"""

from typing import Any

import numpy as np

from .openfold_structure import _group_ids


def attach_bonds(
    raw_atom_site_dict: dict[str, np.ndarray],
    ccd_cache: dict[str, Any],
) -> dict[str, np.ndarray]:
    """Return the atom table with a ``bonds`` (m,2) atom-index array added."""
    atom_name = raw_atom_site_dict["atom_name"].astype(str)
    res_name = raw_atom_site_dict["res_name"].astype(str)
    chain_id = raw_atom_site_dict["chain_id"].astype(str)
    res_id = raw_atom_site_dict["res_id"]

    # Residue boundaries: [starts[r], ends[r]) is residue r's slice of the atom table.
    atom_to_res = _group_ids(chain_id, res_id, raw_atom_site_dict["ins_code"].astype(str))
    starts = np.flatnonzero(np.r_[True, atom_to_res[1:] != atom_to_res[:-1]])
    ends = np.r_[starts[1:], len(atom_to_res)]

    pairs = _intra_residue_bonds(atom_name, res_name, starts, ends, ccd_cache)
    pairs += _peptide_bonds(atom_name, chain_id, res_id, starts, ends)

    out = dict(raw_atom_site_dict)
    out["bonds"] = (
        np.array(pairs, dtype=np.int64) if pairs else np.empty((0, 2), dtype=np.int64)
    )
    return out


def _intra_residue_bonds(
    atom_name: np.ndarray,
    res_name: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    ccd_cache: dict[str, Any],
) -> list[tuple[int, int]]:
    """Map each CCD component's own edge list onto the atoms of every such residue."""
    # Resolve CCD node indices to atom names once per component, not once per residue:
    # a chain has hundreds of residues but only ~20 distinct names.
    edges_by_res: dict[str, list[tuple[str, str]]] = {}
    for name, entry in ccd_cache.items():
        if entry is None:
            continue
        nodes = entry["atom"]["nodes"]["id"]["value"]
        bond_type = entry["atom"]["edges"]["bond_type"]
        edges_by_res[name] = [
            (str(nodes[bond_type["src_indices"][k]]), str(nodes[bond_type["dst_indices"][k]]))
            for k in range(len(bond_type["value"]))
        ]

    pairs: list[tuple[int, int]] = []
    for lo, hi in zip(starts, ends, strict=True):
        edges = edges_by_res.get(res_name[lo])
        if not edges:
            continue
        index_of = {atom_name[i]: i for i in range(lo, hi)}
        for src, dst in edges:
            a, b = index_of.get(src), index_of.get(dst)
            if a is not None and b is not None:  # atom absent from the model -> no bond
                pairs.append((a, b))
    return pairs


def _peptide_bonds(
    atom_name: np.ndarray,
    chain_id: np.ndarray,
    res_id: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
) -> list[tuple[int, int]]:
    """Link C(i)-N(i+1), but only where residue i+1 really follows residue i.

    Adjacency in the file is not adjacency in the protein. A TED domain can be
    sequence-discontinuous -- TED02 of AF-A0A3N5GBR8 is 129-170 plus 229-367, wrapping
    around TED03 -- and both of its segments are one chain here, so a chain check does
    not separate them. 17.8% of chains in a 400-file sample carry such a break.

    The residue NUMBER is the authority: AFDB predicts every residue of the chain, so
    consecutive numbers are always bonded and a jump is always a domain boundary. (Over
    that same sample all 104,773 number-adjacent pairs had a C-N of 1.313-1.349 A and
    none exceeded 2.0 A, so a distance test would reject nothing this already rejects.)
    """
    chain_of_res = chain_id[starts]
    res_id_of_res = res_id[starts]
    # Candidates: same chain and numerically consecutive.
    linkable = np.flatnonzero(
        (chain_of_res[:-1] == chain_of_res[1:])
        & (res_id_of_res[1:] == res_id_of_res[:-1] + 1),
    )

    pairs: list[tuple[int, int]] = []
    for r in linkable.tolist():
        c = _atom_in_residue(atom_name, starts[r], ends[r], "C")
        n = _atom_in_residue(atom_name, starts[r + 1], ends[r + 1], "N")
        if c is not None and n is not None:
            pairs.append((c, n))
    return pairs


def _atom_in_residue(atom_name: np.ndarray, lo: int, hi: int, name: str) -> int | None:
    """Index of ``name`` within one residue's atom range, or None if it is absent."""
    for i in range(lo, hi):
        if atom_name[i] == name:
            return i
    return None
