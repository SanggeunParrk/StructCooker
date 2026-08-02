"""Per-structure QC checks for cif LMDB items (fully vectorized).

Pure-geometry anomaly detection on an assembly biomoldict, used to audit a
built cif LMDB. All checks operate on the finite (resolved), non-water atoms
only -- unresolved atoms are stored with NaN coordinates by the full-length
scheme (not anomalies), and water (HOH/DOD/...) is excluded entirely because
solvent on crystallographic special positions is duplicated by biological-
assembly symmetry expansion and would otherwise dominate the duplicate/clash
counts. A uniform-grid neighbor search (cell side = query radius)
finds candidate atom pairs without an ``O(n^2)`` scan, and *all* per-pair work
is vectorized over numpy arrays (no Python loops over atoms/pairs), so even
multi-million-atom symmetry assemblies finish in seconds within ``O(n_atom)``
memory.

Checks:
  - clash             : non-bonded heavy-atom pairs closer than ``clash_dist``
                        (1-2 and 1-3 bonded neighbours excluded). A hard
                        distance floor is used rather than van der Waals overlap
                        because the data is heavy-atom only -- hydrogen bonds and
                        plain vdW contacts (>=2.4 A) are normal and must not be
                        flagged; only physically impossible separations are.
  - cross_chain_close : heavy-atom pairs in different *author* chains at a
                        bond-like separation (sum of covalent radii +/- bond_tol)
                        and NOT a recorded struct_conn link -- a likely missing
                        inter-chain bond/contact. Author chains are used (not the
                        label_asym_id) because one author chain is routinely split
                        into several cif chains (polymer + glycans/ligands); those
                        intra-author links would otherwise dominate. The covalent
                        band also excludes clashes (too short).
  - broken_bond       : recorded atom-atom bonds whose finite endpoints are
                        farther apart than ``broken_bond_thr`` (stretched/wrong).
  - duplicate_atom    : distinct finite atoms within ``dup_thr`` (overlapping
                        coordinates, e.g. un-collapsed alt-locs).
  - nan_fraction      : fraction of atoms with non-finite coordinates
                        (informational; high values flag sparse models).
"""

from __future__ import annotations

import numpy as np

_MIN_ATOMS_FOR_PAIRS = 2
_WATER_COMPS = ("HOH", "DOD", "WAT", "H2O")
# Cordero covalent radii (angstrom); fallback for unlisted symbols. Used to flag
# inter-chain pairs at a *bond-like* separation (sum of covalent radii +/- tol),
# which excludes both clashing overlaps (too short) and mere contacts (too long).
_COVALENT = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57,
    "CL": 1.02, "BR": 1.20, "I": 1.39, "B": 0.84, "SI": 1.11, "SE": 1.20,
    "ZN": 1.22, "FE": 1.32, "MG": 1.41, "CA": 1.76, "NA": 1.66, "K": 2.03,
    "MN": 1.39, "CU": 1.32, "NI": 1.24, "CO": 1.26, "CD": 1.44, "MO": 1.54,
}
_COVALENT_DEFAULT = 0.75

# Connectivity-consistency constants -----------------------------------------
_METALS = frozenset({
    "ZN", "MG", "CA", "FE", "MN", "CU", "NA", "K", "NI", "CO", "CD", "MO", "HG",
    "PB", "BA", "SR", "CS", "RB", "LI", "AL", "V", "W", "PT", "AU", "AG", "HF",
})
# Common saccharide chem_comp ids (N-/O-glycans); branched-entity residues are
# also treated as glycan regardless of comp.
_SUGAR_COMPS = frozenset({
    "NAG", "NDG", "BMA", "MAN", "FUC", "FUL", "GAL", "GLC", "BGC", "GLA", "SIA", "NGA", "A2G", "XYP", "RIB", "ARA", "GLP", "MAL", "SUC",
})
# Distance (angstrom) above which a recorded inter-residue linkage of a given
# type is "broken" (geometry inconsistent with a real bond).
_BROKEN_THR = {
    "peptide": 2.0, "nucleic": 2.2, "glycan": 1.9, "disulfide": 2.6,
}
_BROKEN_OTHER_MARGIN = 0.6  # other linkages: broken if dist > cov_sum + margin
# "should-be-bonded" windows for the missing-linkage scan.
_SS_MIN, _SS_MAX = 1.8, 2.4          # disulfide SG-SG
_GLY_MIN, _GLY_MAX = 1.2, 1.7        # glycosidic C-O / C-N
_PROT_BB = frozenset({"C", "N"})
_NUC_P = frozenset({"P"})
_NUC_O = frozenset({"O3'", "O5'", "O3*", "O5*"})
_OFFSETS = (
    np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"))
    .reshape(3, -1)
    .T
)


def _grouped_cartesian(
    a_start: np.ndarray,
    a_len: np.ndarray,
    b_start: np.ndarray,
    b_len: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized per-group Cartesian product of two index ranges.

    For each group ``g`` returns every ``(a_start[g] + x, b_start[g] + y)`` with
    ``x in [0, a_len[g])`` and ``y in [0, b_len[g])`` -- no Python loop over groups.
    """
    sizes = (a_len * b_len).astype(np.int64)
    total = int(sizes.sum())
    if total == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    b_len_rep = np.repeat(b_len, sizes)
    base = np.repeat(np.cumsum(sizes) - sizes, sizes)
    local = np.arange(total, dtype=np.int64) - base
    ai = np.repeat(a_start, sizes) + local // b_len_rep
    bi = np.repeat(b_start, sizes) + local % b_len_rep
    return ai, bi


def grid_pairs(  # noqa: PLR0915
    xyz: np.ndarray,
    d_thr: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find finite atom pairs within ``d_thr`` via a uniform-grid scan.

    Returns undirected pairs ``(i, j)`` with ``i < j`` and their distances, in
    the ORIGINAL atom space (NaN rows skipped, not reindexed). Cell side =
    ``d_thr`` so the 27-cell neighbourhood captures every pair within range.
    """
    empty = (np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, np.float64))
    valid = np.all(np.isfinite(xyz), axis=1)
    vidx = np.flatnonzero(valid)
    n = vidx.size
    if n < _MIN_ATOMS_FOR_PAIRS:
        return empty
    vxyz = xyz[vidx]
    # Discretize into d_thr cells, shifted non-negative with a 1-cell margin so a
    # +/-1 offset never underflows, then pack (x, y, z) into a single int64 key.
    # int64 packing keeps searchsorted on the scalar fast path (a structured
    # dtype falls into a ~100x slower element-wise comparison path).
    cell = np.floor(vxyz / d_thr).astype(np.int64)
    cmin = cell.min(axis=0) - 1
    dims = cell.max(axis=0) - cell.min(axis=0) + 3
    mult = np.array([dims[1] * dims[2], dims[2], 1], dtype=np.int64)
    shifted = cell - cmin
    acell = shifted @ mult
    order = np.argsort(acell, kind="stable")
    acell_sorted = acell[order]
    shifted_sorted = shifted[order]
    change = np.diff(acell_sorted) != 0
    starts = np.concatenate(([0], np.nonzero(change)[0] + 1))
    ends = np.concatenate((starts[1:], [n]))
    ukeys = acell_sorted[starts]
    ucoords = shifted_sorted[starts]
    n_uniq = ukeys.shape[0]
    orig_sorted = vidx[order]
    lens = ends - starts
    d_sq = d_thr * d_thr

    out_i, out_j, out_d = [], [], []
    for off in _OFFSETS:
        nkeys = (ucoords + off) @ mult
        pos = np.searchsorted(ukeys, nkeys, side="left")
        ib = pos < n_uniq
        match = np.zeros(n_uniq, dtype=bool)
        if np.any(ib):
            match[ib] = ukeys[pos[ib]] == nkeys[ib]
        if not np.any(match):
            continue
        src = np.nonzero(match)[0]
        dst = pos[match]
        si_sorted, di_sorted = _grouped_cartesian(starts[src], lens[src], starts[dst], lens[dst])
        gi = orig_sorted[si_sorted]
        gj = orig_sorted[di_sorted]
        keep = gi < gj  # dedupes the two directions
        if not np.any(keep):
            continue
        gi, gj = gi[keep], gj[keep]
        dvec = xyz[gi] - xyz[gj]
        dist_sq = np.einsum("ij,ij->i", dvec, dvec)
        within = dist_sq <= d_sq
        if not np.any(within):
            continue
        out_i.append(gi[within])
        out_j.append(gj[within])
        out_d.append(np.sqrt(dist_sq[within]))
    if not out_i:
        return empty
    return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_d)


def _edge_keys(atoms: dict, field: str, n_atoms: int) -> np.ndarray:
    """Return recorded edges of ``field`` as packed int64 ``min*n_atoms+max`` keys."""
    edge = atoms["edges"].get(field)
    if edge is None:
        return np.empty(0, np.int64)
    src = np.asarray(edge["src_indices"], dtype=np.int64)
    dst = np.asarray(edge["dst_indices"], dtype=np.int64)
    if src.size == 0:
        return np.empty(0, np.int64)
    lo = np.minimum(src, dst)
    hi = np.maximum(src, dst)
    return np.unique(lo * n_atoms + hi)


def _one_three_keys(atoms: dict, n_atoms: int) -> np.ndarray:
    """Return packed keys for all 1-3 (angle) atom pairs from the bond graph."""
    bt = atoms["edges"].get("bond_type")
    if bt is None:
        return np.empty(0, np.int64)
    src = np.asarray(bt["src_indices"], dtype=np.int64)
    dst = np.asarray(bt["dst_indices"], dtype=np.int64)
    if src.size == 0:
        return np.empty(0, np.int64)
    # Symmetric adjacency sorted by source; neighbours of each centre form 1-3 pairs.
    centre = np.concatenate([src, dst])
    nbr = np.concatenate([dst, src])
    order = np.argsort(centre, kind="stable")
    centre, nbr = centre[order], nbr[order]
    _, first, counts = np.unique(centre, return_index=True, return_counts=True)
    ai, bi = _grouped_cartesian(first, counts, first, counts)
    na, nb = nbr[ai], nbr[bi]
    keep = na < nb
    na, nb = na[keep], nb[keep]
    if na.size == 0:
        return np.empty(0, np.int64)
    return np.unique(na * n_atoms + nb)


def _water_atom_mask(biomol: dict) -> np.ndarray:
    """Boolean mask over atoms belonging to water residues (HOH/DOD/...)."""
    res_comp = np.asarray(biomol["residues"]["nodes"]["chem_comp_id"]["value"]).astype("U5")
    is_water_res = np.isin(np.char.upper(res_comp), _WATER_COMPS)
    a2r = np.asarray(biomol["index_table"]["atom_to_res"])
    return is_water_res[a2r]


def _polymer_chain_mask(biomol: dict) -> np.ndarray:
    """Boolean mask over chains that are polymers (protein / nucleic acid)."""
    et = np.asarray(biomol["chains"]["nodes"]["entity_type"]["value"]).astype("U40")
    lower = np.char.lower(et)
    return np.char.startswith(lower, "polypeptide") | (np.char.find(lower, "ribonucleotide") >= 0)


def af3_chain_clash(
    biomol: dict,
    *,
    clash_radius: float = 1.1,
    max_clashes: int = 100,
    max_ratio: float = 0.5,
) -> dict:
    """AF3-style inter-chain clash flag over polymer chains.

    ``clashes(A, B) = #{ i in A, j in B : d_ij < clash_radius }``. A structure
    ``has_clash`` if any polymer chain pair has ``clashes > max_clashes`` or
    ``clashes / min(N_A, N_B) > max_ratio`` where ``N`` is the chain's atom count
    (AF3 SI 5.9.3). Returns the flag plus the worst-offending chain pair.
    """
    atoms = biomol["atoms"]
    xyz = np.array(atoms["nodes"]["xyz"]["value"], dtype=np.float64)
    xyz[_water_atom_mask(biomol)] = np.nan  # water excluded from clash counting
    it = biomol["index_table"]
    a2c = np.asarray(it["res_to_chain"])[np.asarray(it["atom_to_res"])]
    is_poly = _polymer_chain_mask(biomol)
    n_chains = is_poly.size
    atoms_per_chain = np.bincount(a2c, minlength=n_chains)

    gi, gj, dist = grid_pairs(xyz, clash_radius)
    out = {"has_clash": False, "worst_clash_pair": None}
    if gi.size == 0:
        return out
    strict = dist < clash_radius
    ca, cb = a2c[gi[strict]], a2c[gj[strict]]
    keep = (ca != cb) & is_poly[ca] & is_poly[cb]
    ca, cb = ca[keep], cb[keep]
    if ca.size == 0:
        return out
    lo = np.minimum(ca, cb)
    hi = np.maximum(ca, cb)
    packed, counts = np.unique(lo.astype(np.int64) * n_chains + hi.astype(np.int64), return_counts=True)
    lo_u = packed // n_chains
    hi_u = packed % n_chains
    min_n = np.minimum(atoms_per_chain[lo_u], atoms_per_chain[hi_u])
    ratio = counts / np.maximum(min_n, 1)
    flag = (counts > max_clashes) | (ratio > max_ratio)
    if not flag.any():
        return out
    worst = int(np.argmax(ratio))
    out["has_clash"] = True
    out["worst_clash_pair"] = [
        int(lo_u[worst]), int(hi_u[worst]), int(counts[worst]),
        int(min_n[worst]), round(float(ratio[worst]), 3),
    ]
    return out


def analyze_assembly_qc(  # noqa: PLR0913, PLR0915
    biomol: dict,
    *,
    clash_dist: float = 2.0,
    bond_tol: float = 0.45,
    broken_bond_thr: float = 2.5,
    dup_thr: float = 0.30,
    max_examples: int = 5,
) -> dict:
    """Run all geometry QC checks on one assembly biomoldict (vectorized).

    Returns per-check counts plus a few worst-offending example pairs.
    """
    atoms = biomol["atoms"]
    xyz = np.array(atoms["nodes"]["xyz"]["value"], dtype=np.float64)
    elem = np.asarray(atoms["nodes"]["element"]["value"]).astype("U2")
    n_atoms = xyz.shape[0]  # total, kept for consistent index packing
    water = _water_atom_mask(biomol)
    coord_finite = np.all(np.isfinite(xyz), axis=1)
    xyz[water] = np.nan  # exclude water from every geometry check
    finite = coord_finite & ~water  # resolved, non-water
    n_nonwater = int((~water).sum())
    it = biomol["index_table"]
    a2c = np.asarray(it["res_to_chain"])[np.asarray(it["atom_to_res"])]
    # Author-chain id per atom: a single author chain is often split into several
    # label_asym_id (cif) chains (polymer + its glycans/ligands, branched sugars),
    # so cross-chain checks must compare auth_asym_id -- not the label-chain index.
    auth_atom = np.asarray(biomol["chains"]["nodes"]["auth_asym_id"]["value"])[a2c]
    elem_up = np.char.upper(elem)
    is_h = elem_up == "H"
    cov = np.array([_COVALENT.get(e, _COVALENT_DEFAULT) for e in elem_up.tolist()])

    bonded_keys = _edge_keys(atoms, "bond_type", n_atoms)
    sconn_keys = _edge_keys(atoms, "struct_conn", n_atoms)
    one3_keys = _one_three_keys(atoms, n_atoms)

    d_query = max(clash_dist, broken_bond_thr)
    gi, gj, dist = grid_pairs(xyz, d_query)

    if gi.size:
        keys = gi.astype(np.int64) * n_atoms + gj.astype(np.int64)
        heavy = ~is_h[gi] & ~is_h[gj]
        is_bonded = np.isin(keys, bonded_keys)
        is_sconn = np.isin(keys, sconn_keys)
        is_13 = np.isin(keys, one3_keys)
        # exclude recorded covalent links (bond_type + struct_conn: disulfide,
        # glycosidic, covalent ligand) so real bonds are not counted as clashes.
        clash_mask = heavy & ~is_bonded & ~is_sconn & ~is_13 & (dist < clash_dist)
        # cross-chain "missing bond": inter-chain heavy pair at a bond-like
        # separation (sum of covalent radii +/- bond_tol), not already a recorded
        # struct_conn. The covalent-radius band excludes clashing overlaps
        # (too short) and incidental contacts (too long).
        cov_sum = cov[gi] + cov[gj]
        bond_like = (dist >= cov_sum - bond_tol) & (dist <= cov_sum + bond_tol)
        cross_mask = (
            heavy
            & (auth_atom[gi] != auth_atom[gj])
            & bond_like
            & ~is_sconn
        )
        dup_mask = dist < dup_thr
    else:
        clash_mask = cross_mask = dup_mask = np.empty(0, dtype=bool)

    def _examples(mask: np.ndarray, *, with_chains: bool = False) -> list:
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            return []
        idx = idx[np.argsort(dist[idx])][:max_examples]
        if with_chains:
            return [
                [int(a2c[gi[k]]), int(a2c[gj[k]]), int(gi[k]), int(gj[k]), round(float(dist[k]), 3)]
                for k in idx
            ]
        return [[int(gi[k]), int(gj[k]), round(float(dist[k]), 3)] for k in idx]

    # Broken bonds: recorded bonds whose finite endpoints are too far apart.
    bt = atoms["edges"].get("bond_type")
    broken_ex: list = []
    n_broken = 0
    if bt is not None:
        bsrc = np.asarray(bt["src_indices"], dtype=np.int64)
        bdst = np.asarray(bt["dst_indices"], dtype=np.int64)
        if bsrc.size:
            fin = finite[bsrc] & finite[bdst]
            bd = np.full(bsrc.shape, -1.0)
            bd[fin] = np.linalg.norm(xyz[bsrc[fin]] - xyz[bdst[fin]], axis=1)
            bmask = fin & (bd > broken_bond_thr)
            n_broken = int(bmask.sum())
            order_b = np.nonzero(bmask)[0]
            order_b = order_b[np.argsort(-bd[order_b])][:max_examples]
            broken_ex = [[int(bsrc[k]), int(bdst[k]), round(float(bd[k]), 3)] for k in order_b]

    af3 = af3_chain_clash(biomol)
    return {
        "n_atoms": n_nonwater,
        "n_finite": int(finite.sum()),
        "nan_fraction": round(1.0 - finite.sum() / n_nonwater, 4) if n_nonwater else 0.0,
        "n_clash": int(clash_mask.sum()),
        "n_cross_chain_close": int(cross_mask.sum()),
        "n_broken_bond": n_broken,
        "n_duplicate_atom": int(dup_mask.sum()),
        "has_clash": af3["has_clash"],
        "worst_clash_pair": af3["worst_clash_pair"],
        "clash_examples": _examples(clash_mask),
        "cross_chain_examples": _examples(cross_mask, with_chains=True),
        "broken_bond_examples": broken_ex,
        "duplicate_examples": _examples(dup_mask),
    }


def _is_glycan_res(comp: np.ndarray, chain_et_lower: np.ndarray, res_chain: np.ndarray) -> np.ndarray:
    """Per-residue mask: residue is a saccharide (branched entity or sugar comp)."""
    branched = np.char.find(chain_et_lower[res_chain], "branched") >= 0
    return branched | np.isin(comp, tuple(_SUGAR_COMPS))


def _classify_linkage(  # noqa: PLR0913
    ai: str, aj: str, ei: str, ej: str, *, poly_i: bool, poly_j: bool,
    nuc_i: bool, nuc_j: bool, gly_i: bool, gly_j: bool,
) -> str:
    """Classify an inter-residue bond into a linkage chemistry type."""
    if ei in _METALS or ej in _METALS:
        return "metal"
    if ei == "S" and ej == "S":
        return "disulfide"
    names = {ai, aj}
    if poly_i and poly_j and names == _PROT_BB:
        return "peptide"
    if nuc_i and nuc_j and (names & _NUC_P) and (names & _NUC_O):
        return "nucleic"
    if gly_i or gly_j:
        return "glycan"
    return "other"


def connectivity_consistency(biomol: dict, *, max_examples: int = 6) -> dict:  # noqa: PLR0912, PLR0915
    """Audit inter-residue connectivity against chemistry + geometry.

    Classifies recorded inter-residue bonds by linkage type (peptide / nucleic
    backbone, glycan, disulfide, metal, other) and flags those whose endpoints
    are too far apart to be a real bond (``broken``). Separately scans for
    high-confidence linkages that *should* exist by geometry but are not recorded
    (``missing``: disulfide SG-SG, glycosidic C-O). Water excluded.
    """
    atoms = biomol["atoms"]
    nd = atoms["nodes"]
    xyz = np.array(nd["xyz"]["value"], dtype=np.float64)
    elem = np.char.upper(np.asarray(nd["element"]["value"]).astype("U2"))
    aid = np.asarray(nd["id"]["value"]).astype("U6")
    water = _water_atom_mask(biomol)
    xyz[water] = np.nan
    it = biomol["index_table"]
    a2r = np.asarray(it["atom_to_res"])
    r2c = np.asarray(it["res_to_chain"])
    comp = np.asarray(biomol["residues"]["nodes"]["chem_comp_id"]["value"]).astype("U6")
    rauth = np.asarray(biomol["residues"]["nodes"]["auth_idx"]["value"]).astype(str)
    et_lower = np.char.lower(np.asarray(biomol["chains"]["nodes"]["entity_type"]["value"]).astype("U40"))
    auth = np.asarray(biomol["chains"]["nodes"]["auth_asym_id"]["value"])
    poly = np.char.startswith(et_lower, "polypeptide")
    nuc = np.char.find(et_lower, "ribonucleotide") >= 0
    res_chain = r2c
    gly_res = _is_glycan_res(comp, et_lower, res_chain)

    broken: dict[str, int] = dict.fromkeys(
        ("peptide", "nucleic", "glycan", "disulfide", "metal", "other"), 0,
    )
    broken_ex: list = []
    recorded: set[tuple[int, int]] = set()

    def _both(edge_field: str) -> tuple[np.ndarray, np.ndarray]:
        edge = atoms["edges"].get(edge_field)
        if edge is None:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        return np.asarray(edge["src_indices"], np.int64), np.asarray(edge["dst_indices"], np.int64)

    # Recorded set (for suppressing missing-linkage false positives) spans both
    # the covalent graph and struct_conn; but broken detection uses bond_type
    # ONLY -- struct_conn also carries hydrogen bonds / base pairs (non-covalent),
    # which must not be judged as "broken bonds".
    for field in ("bond_type", "struct_conn"):
        bs, bd = _both(field)
        for s, d in zip(bs.tolist(), bd.tolist(), strict=True):
            recorded.add((min(s, d), max(s, d)))

    bs, bd = _both("bond_type")
    for s, d in zip(bs.tolist(), bd.tolist(), strict=True):
        ri, rj = int(a2r[s]), int(a2r[d])
        if ri == rj:  # intra-residue (CCD), trusted
            continue
        if not (np.isfinite(xyz[s]).all() and np.isfinite(xyz[d]).all()):
            continue
        ci, cj = r2c[ri], r2c[rj]
        kind = _classify_linkage(
            str(aid[s]), str(aid[d]), str(elem[s]), str(elem[d]),
            poly_i=bool(poly[ci]), poly_j=bool(poly[cj]),
            nuc_i=bool(nuc[ci]), nuc_j=bool(nuc[cj]),
            gly_i=bool(gly_res[ri]), gly_j=bool(gly_res[rj]),
        )
        if kind == "metal":
            continue
        dist = float(np.linalg.norm(xyz[s] - xyz[d]))
        thr = _BROKEN_THR.get(kind)
        if thr is None:  # other
            thr = (
                _COVALENT.get(str(elem[s]), _COVALENT_DEFAULT)
                + _COVALENT.get(str(elem[d]), _COVALENT_DEFAULT)
                + _BROKEN_OTHER_MARGIN
            )
        if dist > thr:
            broken[kind] += 1
            if len(broken_ex) < max_examples:
                broken_ex.append(
                    f"{kind}: {comp[ri]}{rauth[ri]}.{aid[s]}({auth[ci]}) -- "
                    f"{comp[rj]}{rauth[rj]}.{aid[d]}({auth[cj]}) d={dist:.2f}",
                )

    # Missing high-confidence linkages: should-bond by geometry, not recorded.
    miss = {"disulfide": 0, "glycan": 0}
    miss_ex: list = []
    gi, gj, dist = grid_pairs(xyz, max(_SS_MAX, _GLY_MAX))
    for k in range(gi.size):
        i, j = int(gi[k]), int(gj[k])
        if (i, j) in recorded:
            continue
        ri, rj = int(a2r[i]), int(a2r[j])
        if ri == rj:
            continue
        ei, ej, d = str(elem[i]), str(elem[j]), float(dist[k])
        kind = None
        if ei == "S" and ej == "S" and _SS_MIN <= d <= _SS_MAX:
            kind = "disulfide"
        elif _GLY_MIN <= d <= _GLY_MAX and {ei, ej} == {"C", "O"} and (gly_res[ri] or gly_res[rj]):
            kind = "glycan"
        if kind is None:
            continue
        miss[kind] += 1
        if len(miss_ex) < max_examples:
            miss_ex.append(
                f"{kind}: {comp[ri]}{rauth[ri]}.{aid[i]} -- {comp[rj]}{rauth[rj]}.{aid[j]} d={d:.2f}",
            )

    return {
        "n_broken_backbone": broken["peptide"] + broken["nucleic"],
        "n_broken_glycan": broken["glycan"],
        "n_broken_disulfide": broken["disulfide"],
        "n_broken_other": broken["other"],
        "n_missing_disulfide": miss["disulfide"],
        "n_missing_glycan": miss["glycan"],
        "broken_examples": broken_ex,
        "missing_examples": miss_ex,
    }
