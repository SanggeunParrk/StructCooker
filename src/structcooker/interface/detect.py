"""Interface-unit detection on plain arrays (no BioMol dependency).

An assembly is a set of *nodes*: polymer chains (``protein`` / ``na``), ligand instances
(a CCD molecule with its covalently attached parts already merged) and ``ion`` atoms.
Water is never a node (see ``select_waters``).

Every significant edge (node pair) seeds one unit: R = the contacting residues of both
sides; every other node with a significant edge that touches R joins (one hop only, so
a complex does not percolate into one blob). Units with the same node set are merged;
symmetry copies are removed by ``unit_signature``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

_OFFSETS = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij")).reshape(3, -1).T


@dataclass
class ContactCfg:
    cutoff: float = 5.0
    min_residues: int = 5          # polymer-polymer: contacting residues on BOTH sides
    min_ligand_atoms: int = 3      # polymer-ligand: ligand atoms in contact
    min_ligand_residues: int = 2   # polymer-ligand: polymer residues in contact
    ion_cutoff: float = 3.0        # coordination distance (MetalPDB / InterMetalDB); only N/O/S donors count
    min_ion_residues: int = 1


@dataclass
class UnitCfg:
    max_polymer_chains: int = 6
    max_ligands: int = 4
    max_ions: int = 8
    join_min_residues: int = 3     # joiner's residues touching R
    allowed_pairs: frozenset = frozenset({
        ("protein", "protein"), ("na", "protein"), ("ligand", "protein"),
        ("ion", "protein"), ("ion", "na"),
    })
    max_polymer_chains_in_assembly: int = 60
    max_units: int = 200
    weight_bucket: int = 5
    max_forced_polymer: int = 10   # antibody partners / duplex strands bypass max_polymer_chains up to this


@dataclass
class Edge:
    a: int
    b: int
    n_atom_pairs: int
    res_pairs: np.ndarray  # (n, 2) unique residue pairs, column 0 on side a
    atoms_a: np.ndarray    # contacting atom indices on side a
    atoms_b: np.ndarray
    min_dist: float

    @property
    def res_a(self) -> np.ndarray:
        return np.unique(self.res_pairs[:, 0])

    @property
    def res_b(self) -> np.ndarray:
        return np.unique(self.res_pairs[:, 1])

    @property
    def weight(self) -> int:
        return len(self.res_a) + len(self.res_b)


@dataclass
class Unit:
    nodes: list[int]
    edges: list[Edge] = field(default_factory=list)
    seeds: list[tuple[int, int]] = field(default_factory=list)
    truncated: bool = False

    @property
    def weight(self) -> int:
        return sum(e.weight for e in self.edges)


def inter_node_pairs(xyz: np.ndarray, node_idx: np.ndarray, cutoff: float):
    """Atom pairs (i < j) of different nodes within ``cutoff`` (cell list, exact)."""
    ok = np.flatnonzero(np.isfinite(xyz).all(axis=1))
    if ok.size < 2:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0)
    cell = np.floor(xyz[ok] / cutoff).astype(np.int64)
    order = np.lexsort((cell[:, 2], cell[:, 1], cell[:, 0]))
    cs = cell[order]
    chg = np.any(np.diff(cs, axis=0) != 0, axis=1)
    starts = np.concatenate(([0], np.flatnonzero(chg) + 1))
    counts = np.diff(np.concatenate((starts, [len(cs)])))
    ucell = cs[starts]
    # pack cell coordinates into one sortable key
    lo = ucell.min(axis=0) - 1
    span = ucell.max(axis=0) - lo + 2
    key = (ucell - lo) @ np.array([span[1] * span[2], span[2], 1])
    ks = np.argsort(key)  # ucell is lexsorted so already ordered, kept explicit for safety
    key_sorted = key[ks]
    I, J = [], []
    for off in _OFFSETS:
        nk = (ucell + off - lo) @ np.array([span[1] * span[2], span[2], 1])
        pos = np.searchsorted(key_sorted, nk)
        pos[pos >= len(key_sorted)] = 0
        hit = key_sorted[pos] == nk
        if not hit.any():
            continue
        ca = np.flatnonzero(hit)
        cb = ks[pos[hit]]
        na, nb = counts[ca], counts[cb]
        tot = na * nb
        pid = np.repeat(np.arange(len(ca)), tot)
        within = np.arange(tot.sum()) - np.repeat(np.cumsum(tot) - tot, tot)
        ia = starts[ca][pid] + within // nb[pid]
        ib = starts[cb][pid] + within % nb[pid]
        I.append(order[ia])
        J.append(order[ib])
    if not I:
        e = np.empty(0, dtype=np.int64)
        return e, e, np.empty(0)
    i = ok[np.concatenate(I)]
    j = ok[np.concatenate(J)]
    keep = (i < j) & (node_idx[i] != node_idx[j])
    i, j = i[keep], j[keep]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    near = d < cutoff
    return i[near], j[near], d[near]


def node_edges(xyz, node_idx, res_idx, kinds, ccfg: ContactCfg, donor=None) -> dict[tuple[int, int], Edge]:
    """All node-node contacts (no significance filter yet), keyed ``(a, b)`` with a < b.

    Ion edges use ``ion_cutoff`` and, if ``donor`` (bool per atom) is given, only count donor atoms.
    Everything else uses ``cutoff``.
    """
    i, j, d = inter_node_pairs(xyz, node_idx, max(ccfg.cutoff, ccfg.ion_cutoff))
    if i.size == 0:
        return {}
    is_ion = np.array([k == "ion" for k in kinds])
    ion_pair = is_ion[node_idx[i]] | is_ion[node_idx[j]]
    keep = np.where(ion_pair, d < ccfg.ion_cutoff, d < ccfg.cutoff)
    if donor is not None:
        other = np.where(is_ion[node_idx[i]], j, i)
        keep &= ~ion_pair | donor[other]
    i, j, d = i[keep], j[keep], d[keep]
    if i.size == 0:
        return {}
    na, nb = node_idx[i], node_idx[j]
    swap = na > nb
    i, j = np.where(swap, j, i), np.where(swap, i, j)
    na, nb = node_idx[i], node_idx[j]
    nn = int(node_idx.max()) + 1
    pair = na * nn + nb
    order = np.argsort(pair, kind="stable")
    pair_s = pair[order]
    cut = np.concatenate(([0], np.flatnonzero(np.diff(pair_s)) + 1, [len(pair_s)]))
    out = {}
    for lo, hi in zip(cut[:-1], cut[1:], strict=True):
        m = order[lo:hi]
        p = int(pair_s[lo])
        a, b = p // nn, p % nn
        rp = np.unique(np.stack([res_idx[i[m]], res_idx[j[m]]], axis=1), axis=0)
        out[(a, b)] = Edge(a, b, len(m), rp, np.unique(i[m]), np.unique(j[m]), float(d[m].min()))
    return out


def _significant(e: Edge, kinds: list[str], ccfg: ContactCfg, ucfg: UnitCfg) -> bool:
    ka, kb = kinds[e.a], kinds[e.b]
    if (ka, kb) not in ucfg.allowed_pairs and (kb, ka) not in ucfg.allowed_pairs:
        return False
    for small, k in (("ion", ccfg.min_ion_residues), ("ligand", ccfg.min_ligand_residues)):
        if small in (ka, kb):
            if ka == small:
                atoms, pol_res = e.atoms_a, e.res_b
            else:
                atoms, pol_res = e.atoms_b, e.res_a
            min_atoms = 1 if small == "ion" else ccfg.min_ligand_atoms
            return len(atoms) >= min_atoms and len(pol_res) >= k
    return min(len(e.res_a), len(e.res_b)) >= ccfg.min_residues


def _side(e: Edge, node: int) -> int:
    return 0 if e.a == node else 1


def _touch(c: int, R: dict[int, np.ndarray], sig: dict[tuple[int, int], Edge]) -> int:
    """Residues of node ``c`` that contact a residue of R (R maps node -> residue ids)."""
    mine = [np.empty(0, dtype=np.int64)]
    for m, rm in R.items():
        e = sig.get((min(c, m), max(c, m)))
        if e is None:
            continue
        sc, sm = _side(e, c), _side(e, m)
        sel = np.isin(e.res_pairs[:, sm], rm)
        mine.append(e.res_pairs[sel, sc])
    return len(np.unique(np.concatenate(mine)))


def find_units(
    xyz: np.ndarray,
    node_idx: np.ndarray,
    res_idx: np.ndarray,
    kinds: list[str],
    ccfg: ContactCfg | None = None,
    ucfg: UnitCfg | None = None,
    *,
    donor: np.ndarray | None = None,
    ab_class: list[str | None] | None = None,
    duplex: dict[tuple[int, int], int] | None = None,
) -> list[Unit]:
    """One unit per significant edge (merged by node set); see the module docstring.

    ``ab_class`` (ANARCI chain type per node, None for non-antibody) and ``duplex`` (NA strand pairs
    from ``find_duplex``) add *forced* partners: whenever a unit holds an antibody / TCR chain, its
    best-contacting complementary chain joins; whenever it holds an NA strand, every strand it
    base-pairs with joins. Forced partners bypass ``max_polymer_chains`` (up to ``max_forced_polymer``).
    """
    from structcooker.interface.antibody import PAIR_PARTNERS

    ccfg, ucfg = ccfg or ContactCfg(), ucfg or UnitCfg()
    ab_class = ab_class or [None] * len(kinds)
    duplex = duplex or {}
    if sum(k in ("protein", "na") for k in kinds) > ucfg.max_polymer_chains_in_assembly:
        return []
    edges = node_edges(xyz, node_idx, res_idx, kinds, ccfg, donor)
    sig = {k: e for k, e in edges.items() if _significant(e, kinds, ccfg, ucfg)}
    partners: dict[int, set[int]] = {}
    for a, b in sig:
        partners.setdefault(a, set()).add(b)
        partners.setdefault(b, set()).add(a)
    dup_partners: dict[int, set[int]] = {}
    for a, b in duplex:
        dup_partners.setdefault(a, set()).add(b)
        dup_partners.setdefault(b, set()).add(a)

    def forced(n: int) -> set[int]:
        out = set(dup_partners.get(n, ()))
        if ab_class[n]:
            best, best_w = None, 0
            for m in partners.get(n, ()):
                if ab_class[m] in PAIR_PARTNERS.get(ab_class[n], ()):
                    w = sig[(min(n, m), max(n, m))].weight
                    if w > best_w:
                        best, best_w = m, w
            if best is not None:
                out.add(best)
        return out

    def is_bundle_edge(key) -> bool:
        a, b = key
        return bool(ab_class[a] and ab_class[b] and ab_class[b] in PAIR_PARTNERS.get(ab_class[a], ()))

    # an ion alone is not an interface: it seeds only when it bridges >= 2 nodes
    seeds = [k for k in sig if all(kinds[n] != "ion" or len(partners[n]) >= 2 for n in k)]
    limit = {"protein": ucfg.max_polymer_chains, "na": ucfg.max_polymer_chains,
             "ligand": ucfg.max_ligands, "ion": ucfg.max_ions}
    poly = ("protein", "na")
    by_nodes: dict[tuple[int, ...], Unit] = {}
    bundle_only: dict[tuple[int, ...], bool] = {}

    def close(group, truncated):
        queue = list(group)
        while queue:
            for p in forced(queue.pop()):
                if p in group:
                    continue
                if sum(kinds[g] in poly for g in group) >= ucfg.max_forced_polymer:
                    truncated = True
                    continue
                group.append(p)
                queue.append(p)
        return truncated

    for key in seeds:
        seed = sig[key]
        R = {seed.a: seed.res_a, seed.b: seed.res_b}
        group = list(key)
        truncated = close(group, False)
        cands = []
        for c in set(partners[seed.a]) | set(partners[seed.b]):
            if c in group:
                continue
            s = _touch(c, R, sig)
            if kinds[c] in poly:
                # a polymer linked through an ion / ligand is judged by that edge's threshold
                small = [kinds[m] for m in key if kinds[m] in ("ion", "ligand")]
                need = (ccfg.min_ion_residues if "ion" in small else ccfg.min_ligand_residues
                        if "ligand" in small else ucfg.join_min_residues)
            else:
                need = 1
            if s >= need:
                cands.append((s, c))
        cands.sort(key=lambda t: (-t[0], t[1]))
        for _, c in cands:
            if c in group:
                continue
            n_kind = sum(
                (kinds[g] in poly) == (kinds[c] in poly) and (kinds[c] in poly or kinds[g] == kinds[c])
                for g in group
            )
            if n_kind >= limit[kinds[c]]:
                truncated = True
                continue
            group.append(c)
            truncated = close(group, truncated)
        gk = tuple(sorted(group))
        if gk in by_nodes:
            by_nodes[gk].seeds.append(key)
            by_nodes[gk].truncated |= truncated
            bundle_only[gk] &= is_bundle_edge(key)
            continue
        members = {(a, b) for a in gk for b in gk if a < b}
        by_nodes[gk] = Unit(list(gk), [sig[k] for k in sorted(members) if k in sig], [key], truncated)
        bundle_only[gk] = is_bundle_edge(key)
    # a bare H-L (or TCR alpha-beta) unit is redundant when a larger unit holds the pair plus a partner
    sets = {gk: frozenset(gk) for gk in by_nodes}
    for gk in list(by_nodes):
        if bundle_only[gk] and any(sets[gk] < sets[o] for o in by_nodes if o != gk):
            del by_nodes[gk]
    units = sorted(by_nodes.values(), key=lambda u: (-u.weight, u.nodes))
    return units[: ucfg.max_units]


def select_waters(
    unit_xyz: np.ndarray,
    unit_node: np.ndarray,
    water_xyz: np.ndarray,
    near: float = 4.0,
    bridge: float = 3.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Waters (O atoms) kept as context of a unit.

    Returns ``(indices into water_xyz, bridging flag)``; a water is kept when within
    ``near`` of any unit atom and flagged bridging when it is within ``bridge`` of atoms of
    two or more different unit nodes. Water never affects edges or residue sets.
    """
    if len(water_xyz) == 0 or len(unit_xyz) == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=bool)
    n_u = len(unit_xyz)
    xyz = np.vstack([unit_xyz, water_xyz])
    node = np.concatenate([unit_node, np.full(len(water_xyz), int(unit_node.max()) + 1)])
    i, j, d = inter_node_pairs(xyz, node, max(near, bridge))
    # inter_node_pairs orders i < j, so a unit-water pair has the unit atom as i and the water as j
    is_uw = j >= n_u
    i, j, d = i[is_uw], j[is_uw], d[is_uw]
    w = j - n_u
    keep_mask = d < near
    kept = np.unique(w[keep_mask])
    bridging = np.zeros(len(kept), dtype=bool)
    near_b = d < bridge
    pair = np.unique(np.stack([w[near_b], unit_node[i[near_b]]], axis=1), axis=0)
    if len(pair):
        ids, counts = np.unique(pair[:, 0], return_counts=True)
        two = set(ids[counts >= 2].tolist())
        bridging = np.array([k in two for k in kept.tolist()], dtype=bool)
    return kept, bridging


def unit_signature(unit: Unit, node_keys: list[str], bucket: int = 5) -> tuple:
    """Symmetry-invariant identity: node keys (seq id / CCD id) + bucketed contact sizes."""
    nodes = tuple(sorted(node_keys[n] for n in unit.nodes))
    es = tuple(sorted(
        (*sorted((node_keys[e.a], node_keys[e.b])), round(e.weight / bucket))
        for e in unit.edges
    ))
    return nodes, es


def drop_ion_subsets(units: list[Unit], kinds: list[str]) -> list[Unit]:
    """Drop a unit whose nodes are a strict subset of another unit's, the difference being only ions."""
    sets = [frozenset(u.nodes) for u in units]
    keep = []
    for n, u in enumerate(units):
        covered = any(
            sets[n] < sets[m] and all(kinds[x] == "ion" for x in sets[m] - sets[n])
            for m in range(len(units)) if m != n
        )
        if not covered:
            keep.append(u)
    return keep


def dedup(units: list[Unit], node_keys: list[str], bucket: int = 5) -> list[Unit]:
    seen, out = set(), []
    for u in units:
        s = unit_signature(u, node_keys, bucket)
        if s not in seen:
            seen.add(s)
            out.append(u)
    return out
