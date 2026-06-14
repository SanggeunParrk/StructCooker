"""Geometry helpers: grid neighbor search and chain-contact graph extraction.

Pure spatial operations used by the CIF ingest; not CIF-format parsing. Moved
verbatim from ``cif.py`` (behavior unchanged).
"""

from collections.abc import Callable

import numpy as np
from biomol.core.container import FeatureContainer
from biomol.core.feature import EdgeFeature


def neighbor_list_grid(  # noqa: PLR0915
    xyz: np.ndarray,
    d_thr: float,
    n_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute (n_atom, n_max) neighbor indices using a uniform grid of cell size d_thr.

    Vectorized over atoms; only a tiny fixed loop over 27 neighbor-cell offsets.

    Parameters
    ----------
    xyz : (n_atom, 3) float32/64
        Coordinates; any row containing NaN is ignored.
    d_thr : float
        L2 distance threshold for neighbor definition.
    n_max : int
        Maximum number of neighbors to keep per atom.

    Returns
    -------
    nbrs : (n_atom, n_max) int64
        Neighbor indices, padded with -1.
    counts : (n_atom,) int32
        Actual neighbor counts for each atom (self excluded).
    """
    n_atom = xyz.shape[0]
    nbrs = np.full((n_atom, n_max), -1, dtype=np.int64)
    counts = np.zeros(n_atom, dtype=np.int32)

    # 1) Mask invalid atoms (any NaN)
    valid = np.all(np.isfinite(xyz), axis=1)
    if not np.any(valid):
        return nbrs, counts

    # 2) Compressed array of valid points
    valid_xyz = xyz[valid]
    n_valid = valid_xyz.shape[0]

    # 3) Discretize into cells of side length d_thr (int64 coordinates)
    cell = np.floor(valid_xyz / d_thr).astype(np.int64)  # (n_valid, 3)

    # 4) Group points by cell via lexicographic sort on (x, y, z)
    order = np.lexsort((cell[:, 2], cell[:, 1], cell[:, 0]))
    cell_sorted = cell[order]

    # 5) Unique cells and their spans [start, end) in the sorted index space
    if n_valid > 1:
        change = np.any(np.diff(cell_sorted, axis=0) != 0, axis=1)
        starts = np.concatenate(([0], np.nonzero(change)[0] + 1))
    else:
        starts = np.array([0], dtype=np.int64)
    ends = np.concatenate((starts[1:], [n_valid]))
    unique_cells = cell_sorted[starts]  # (n_unique, 3)
    n_unique = unique_cells.shape[0]

    # 6) Helper: view (n,3) int64 as a structured dtype for consistent numeric lex compare
    #    (little-endian int64 x,y,z). This matches the lexsort order above.
    def as_struct3(a_int64x3: np.ndarray) -> np.ndarray:
        a = np.ascontiguousarray(a_int64x3)
        dt = np.dtype([("x", "<i8"), ("y", "<i8"), ("z", "<i8")])
        return a.view(dt).ravel()

    unique_struct = as_struct3(unique_cells)

    # 7) Precompute inverse map from compressed->sorted if needed later
    inv_order = np.empty(n_valid, dtype=np.int64)
    inv_order[order] = np.arange(n_valid)

    # 8) 27 neighbor-cell offsets (-1,0,1)^3
    offsets = (
        np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"))
        .reshape(3, -1)
        .T
    )  # (27, 3)

    # 9) Accumulate candidate (i,j) pairs in compressed indices
    pair_i = []
    pair_j = []

    # For each offset, find matching neighbor cells and generate Cartesian pairs
    for off in offsets:
        # neighbor cells for ALL existing unique cells under this offset
        nei_cells = unique_cells + off  # (n_unique, 3)
        nei_struct = as_struct3(nei_cells)

        # Search where these neighbor cells would be, under the SAME structured ordering
        pos = np.searchsorted(unique_struct, nei_struct, side="left")
        in_bounds = pos < n_unique

        ok = np.zeros_like(in_bounds, dtype=bool)
        if np.any(in_bounds):
            ok[in_bounds] = unique_struct[pos[in_bounds]] == nei_struct[in_bounds]

        if not np.any(ok):
            continue

        # Source cell ids are those indices where a neighbor cell exists
        src_c = np.nonzero(ok)[0]  # indices in [0..n_unique)
        dst_c = pos[ok]  # matching neighbor cell ids

        # Build index ranges for points in each cell's span (sorted space indices)
        src_ranges = [np.arange(starts[c], ends[c], dtype=np.int64) for c in src_c]
        dst_ranges = [np.arange(starts[c], ends[c], dtype=np.int64) for c in dst_c]

        if len(src_ranges) == 0:
            continue

        # Cartesian product per (src_cell, dst_cell) pair (vectorized at cell level)
        src_idx_sorted = np.concatenate(
            [
                np.repeat(r, len(d))
                for r, d in zip(src_ranges, dst_ranges, strict=False)
            ],
        )
        if src_idx_sorted.size == 0:
            continue
        dst_idx_sorted = np.concatenate(
            [np.tile(d, len(r)) for r, d in zip(src_ranges, dst_ranges, strict=False)],
        )

        # Map back from sorted space → compressed (unsorted) space
        src_idx = order[src_idx_sorted]
        dst_idx = order[dst_idx_sorted]

        # Drop self-pairs
        keep = src_idx != dst_idx
        if not np.any(keep):
            continue

        pair_i.append(src_idx[keep])
        pair_j.append(dst_idx[keep])

    if not pair_i:
        # No candidate pairs; return empty neighbor lists
        return nbrs, counts

    pair_i = np.concatenate(pair_i)
    pair_j = np.concatenate(pair_j)

    # 10) Distance filtering: keep pairs with ||valid_xyz[i]-valid_xyz[j]|| <= d_thr
    dvec = valid_xyz[pair_i] - valid_xyz[pair_j]
    dist2 = np.einsum("ij,ij->i", dvec, dvec)
    keep = dist2 <= (d_thr * d_thr)
    if not np.any(keep):
        return nbrs, counts

    pair_i = pair_i[keep]
    pair_j = pair_j[keep]

    # 11) Remove duplicate pairs (same (i,j) can appear via multiple offsets)
    ij = np.stack([pair_i, pair_j], axis=1).astype(np.int64)
    ij_packed = ij.view(np.dtype((np.void, ij.dtype.itemsize * 2))).ravel()
    uniq_idx = np.unique(ij_packed, return_index=True)[1]
    ij = ij[uniq_idx]
    pair_i, pair_j = ij[:, 0], ij[:, 1]

    # 12) Map compressed indices back to original n_atom-space
    valid_true_idx = np.flatnonzero(valid)
    gi = valid_true_idx[pair_i]
    gj = valid_true_idx[pair_j]

    # 13) Fill neighbor matrix: group by gi, keep up to n_max in order of (gi, gj)
    order_fill = np.lexsort((gj, gi))
    gi = gi[order_fill]
    gj = gj[order_fill]

    uniq_i, first_pos, counts_all = np.unique(gi, return_index=True, return_counts=True)
    take_counts = np.minimum(counts_all, n_max)

    if uniq_i.size > 0:
        gather_idx = np.concatenate(
            [
                np.arange(s, s + t, dtype=np.int64)
                for s, t in zip(first_pos, take_counts, strict=True)
            ],
        )
        gi_take = gi[gather_idx]
        gj_take = gj[gather_idx]

        # Relative slot [0..taken-1] within each group
        rel = np.arange(gj.size, dtype=np.int64) - np.repeat(first_pos, counts_all)
        rel = rel[gather_idx]

        nbrs[gi_take, rel] = gj_take
        counts[uniq_i] = take_counts.astype(np.int32)

    return nbrs, counts


def extract_contact_graph(
    d_thr: float = 6.0,
    n_max: int = 128,
) -> Callable[..., FeatureContainer]:
    """Return a configured instruction function that extracts contact chains.

    Edge values will be the number of atom-atom contacts (pairs) between chain pairs.
    """

    def _function(container_dict: dict) -> FeatureContainer:
        xyz = container_dict["atoms"]["xyz"].value  # (L,3)
        chain_idx = container_dict["index_table"].atoms_to_chains(
            np.arange(xyz.shape[0]),
        )  # (L,)
        # 2) Neighborhood via grid cells
        nbrs, counts = neighbor_list_grid(
            xyz,
            d_thr,
            n_max,
        )  # nbrs: (L, N_max), -1 padded

        # 3) Build (chain_i, chain_j) for every valid inter-chain atom pair
        #    - broadcast chain_i over neighbor slots
        chain_idx_i_matrix = np.broadcast_to(chain_idx[:, np.newaxis], nbrs.shape)

        #    - map neighbor indices (including -1) to chain_j using padding trick
        padded_chain_idx = np.append(chain_idx, -1)
        chain_idx_j_matrix = padded_chain_idx[nbrs]

        #    - keep only valid neighbors and inter-chain pairs
        valid_neighbor_mask = nbrs != -1
        inter_chain_mask = chain_idx_i_matrix != chain_idx_j_matrix
        final_mask = valid_neighbor_mask & inter_chain_mask

        if not np.any(final_mask):
            # No contacts at all: return empty edge set
            contact_edges = EdgeFeature(
                value=np.empty((0,), dtype=np.int32),
                src_indices=np.empty((0,), dtype=chain_idx.dtype),
                dst_indices=np.empty((0,), dtype=chain_idx.dtype),
            )
            chain_container = container_dict["chains"]
            chain_container = chain_container.update(contact=contact_edges.copy())
            container_dict["chains"] = chain_container
            return container_dict

        # 4) Collect chain pair list for all contacting atom pairs
        chain_pairs_i = chain_idx_i_matrix[final_mask]
        chain_pairs_j = chain_idx_j_matrix[final_mask]

        # 5) Make edges undirected by sorting endpoints within each pair
        chain_pairs = np.stack([chain_pairs_i, chain_pairs_j], axis=1)  # (E_atom, 2)
        sorted_edges = np.sort(chain_pairs, axis=1)  # ensure (min, max)

        # 6) Count how many atom-atom contacts per chain-pair
        #    Use unique by rows with counts
        contact_edges_unique, counts_per_edge = np.unique(
            sorted_edges,
            axis=0,
            return_counts=True,
        )  # contact_edges_unique: (E_chain, 2), counts_per_edge: (E_chain,)

        # 7) Build EdgeFeature with counts as values
        contact_src = contact_edges_unique[:, 0]
        contact_dst = contact_edges_unique[:, 1]

        contact_edges = EdgeFeature(
            value=counts_per_edge.astype(np.int32),  # number of atom-atom contacts
            src_indices=contact_src,
            dst_indices=contact_dst,
        )
        chain_container = container_dict["chains"]
        chain_container = chain_container.update(contact=contact_edges.copy())
        container_dict["chains"] = chain_container
        return container_dict

    def _worker(assembly_dict: dict) -> FeatureContainer:
        output = {}
        for key, container_dict in assembly_dict.items():
            output[key] = _function(container_dict)
        return output

    return _worker
