"""Geometry helpers: grid neighbor search and chain-contact graph extraction.

Pure spatial operations used by the CIF ingest; not CIF-format parsing.

The chain-contact graph is built with a uniform-grid (cell-list) neighbor
search. Cell side = ``d_thr`` so any pair within ``d_thr`` lands in cells that
differ by at most one along each axis -> scanning the 27 neighbor cells finds
every true contact (no false negatives). Inter-chain contacts are accumulated
per chain-pair on the fly, one cell-offset at a time, so neither a dense
``(n_atom, n_max)`` neighbor matrix nor a global atom-pair edge list is ever
materialized. Peak memory is therefore ``O(n_atom)`` for the grid plus a single
offset's candidate pairs, which keeps even multi-million-atom symmetry
assemblies (e.g. icosahedral viral capsids) within a single node's RAM.
"""

from collections import defaultdict
from collections.abc import Callable

import numpy as np
from biomol.core.feature import EdgeFeature

# Little-endian (x, y, z) int64 cell key; structured-dtype compare matches the
# lexsort order used to group atoms into cells.
_CELL_DTYPE = np.dtype([("x", "<i8"), ("y", "<i8"), ("z", "<i8")])
# 27 neighbor-cell offsets (-1, 0, 1)^3.
_OFFSETS = (
    np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"))
    .reshape(3, -1)
    .T
)


def _as_cell_keys(cells_int64x3: np.ndarray) -> np.ndarray:
    """View an ``(n, 3)`` int64 array as a 1-D structured (x, y, z) key array."""
    return np.ascontiguousarray(cells_int64x3).view(_CELL_DTYPE).ravel()


def chain_contacts_grid(
    xyz: np.ndarray,
    chain_idx: np.ndarray,
    d_thr: float,
    *,
    count_atom_pairs_once: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Count inter-chain atom contacts per chain-pair via a streaming grid scan.

    Parameters
    ----------
    xyz : (n_atom, 3) float
        Atom coordinates; any row containing NaN/inf is ignored.
    chain_idx : (n_atom,) int
        Chain index of each atom.
    d_thr : float
        L2 distance threshold defining an atom-atom contact.
    count_atom_pairs_once : bool
        Count each unordered atom pair once. False preserves the historical PDB
        directed-pair multiplicity (both a->b and b->a contribute).

    Returns
    -------
    src, dst : (n_edge,) int64
        Undirected chain-pair endpoints with ``src < dst``.
    counts : (n_edge,) int64
        Number of contacting atom pairs, with the selected directional convention.
    """
    empty = (
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.int64),
    )
    # 1) Compress to finite atoms only.
    valid = np.all(np.isfinite(xyz), axis=1)
    n_valid = int(np.count_nonzero(valid))
    if n_valid == 0:
        return empty
    valid_xyz = xyz[valid]
    chain_valid = chain_idx[valid]
    n_chains = int(chain_idx.max()) + 1 if chain_idx.size else 0

    # 2) Discretize into d_thr-sided cells and group atoms by cell via lexsort.
    cell = np.floor(valid_xyz / d_thr).astype(np.int64)
    order = np.lexsort((cell[:, 2], cell[:, 1], cell[:, 0]))
    cell_sorted = cell[order]
    if n_valid > 1:
        change = np.any(np.diff(cell_sorted, axis=0) != 0, axis=1)
        starts = np.concatenate(([0], np.nonzero(change)[0] + 1))
    else:
        starts = np.array([0], dtype=np.int64)
    ends = np.concatenate((starts[1:], [n_valid]))
    unique_cells = cell_sorted[starts]
    n_unique = unique_cells.shape[0]
    unique_keys = _as_cell_keys(unique_cells)

    d_thr_sq = d_thr * d_thr
    # Accumulate per chain-pair counts keyed by packed (lo * n_chains + hi).
    acc: dict[int, int] = defaultdict(int)

    # 3) One cell-offset at a time: pair up atoms in adjacent cells, keep
    #    inter-chain pairs within d_thr, tally per chain-pair, then discard.
    for off in _OFFSETS:
        nei_keys = _as_cell_keys(unique_cells + off)
        pos = np.searchsorted(unique_keys, nei_keys, side="left")
        in_bounds = pos < n_unique
        match = np.zeros(n_unique, dtype=bool)
        if np.any(in_bounds):
            match[in_bounds] = unique_keys[pos[in_bounds]] == nei_keys[in_bounds]
        if not np.any(match):
            continue
        src_cells = np.nonzero(match)[0]
        dst_cells = pos[match]

        # Cartesian product of atom indices for each matched (src_cell, dst_cell).
        src_chunks = []
        dst_chunks = []
        for sc, dc in zip(src_cells.tolist(), dst_cells.tolist(), strict=True):
            s_rng = np.arange(int(starts[sc]), int(ends[sc]), dtype=np.int64)
            d_rng = np.arange(int(starts[dc]), int(ends[dc]), dtype=np.int64)
            src_chunks.append(np.repeat(s_rng, d_rng.size))
            dst_chunks.append(np.tile(d_rng, s_rng.size))
        src_sorted = np.concatenate(src_chunks)
        dst_sorted = np.concatenate(dst_chunks)
        if src_sorted.size == 0:
            continue

        # Map sorted-space -> compressed atom indices, then to chains.
        si = order[src_sorted]
        di = order[dst_sorted]
        cs = chain_valid[si]
        cd = chain_valid[di]

        # Keep inter-chain pairs first (drops self-pairs and intra-chain bulk).
        inter = cs < cd if count_atom_pairs_once else cs != cd
        if not np.any(inter):
            continue
        si, di, cs, cd = si[inter], di[inter], cs[inter], cd[inter]

        # Distance filter on the (now much smaller) inter-chain candidates.
        dvec = valid_xyz[si] - valid_xyz[di]
        keep = np.einsum("ij,ij->i", dvec, dvec) <= d_thr_sq
        if not np.any(keep):
            continue
        cs, cd = cs[keep], cd[keep]

        lo = np.minimum(cs, cd)
        hi = np.maximum(cs, cd)
        packed = lo.astype(np.int64) * n_chains + hi.astype(np.int64)
        keys, pair_counts = np.unique(packed, return_counts=True)
        for key, cnt in zip(keys.tolist(), pair_counts.tolist(), strict=True):
            acc[key] += cnt

    if not acc:
        return empty
    packed_keys = np.fromiter(acc.keys(), dtype=np.int64, count=len(acc))
    counts = np.fromiter(acc.values(), dtype=np.int64, count=len(acc))
    src = packed_keys // n_chains
    dst = packed_keys % n_chains
    return src, dst, counts


def extract_contact_graph(
    d_thr: float = 6.0,
) -> Callable[..., dict]:
    """Return an instruction that attaches a chain-level contact graph.

    Each undirected edge connects two chains with at least one atom pair within
    ``d_thr``; its value is the number of contacting atom pairs.
    """

    def _function(container_dict: dict) -> dict:
        xyz = container_dict["atoms"]["xyz"].value  # (L, 3)
        chain_idx = container_dict["index_table"].atoms_to_chains(
            np.arange(xyz.shape[0]),
        )  # (L,)
        src, dst, counts = chain_contacts_grid(xyz, chain_idx, d_thr)
        contact_edges = EdgeFeature(
            value=counts.astype(np.int32),
            src_indices=src.astype(chain_idx.dtype),
            dst_indices=dst.astype(chain_idx.dtype),
        )
        chain_container = container_dict["chains"]
        chain_container = chain_container.update(contact=contact_edges.copy())
        container_dict["chains"] = chain_container
        return container_dict

    def _worker(assembly_dict: dict) -> dict:
        output = {}
        for key, container_dict in assembly_dict.items():
            output[key] = _function(container_dict)
        return output

    return _worker
