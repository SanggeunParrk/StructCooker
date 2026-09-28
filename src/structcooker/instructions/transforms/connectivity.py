"""Connectivity QC transforms (CCD-driven bond reference + per-entry audit).

Two stages, both wired as DataCooker recipes:

* ``build_ccd_linker`` -- per chemical component, derive the atoms that can form
  an *inter-residue* covalent bond (a non-leaving atom bonded to a leaving atom
  gains free valence when the leaving group departs) plus the intra-residue
  ideal bond lengths (from the CCD model coordinates). Materializing this per
  component yields, by composition, the bond reference for every CCD pair: any
  ``(comp1.atom_a, comp2.atom_b)`` is bondable iff both are linker atoms.

* ``connectivity_qc`` -- for one assembly, classify every inter-residue atom
  pair (recorded bonds + close non-bonded neighbours) on three axes -- bondable
  (CCD pair reference) x distance (too-close / bond-band / far) x recorded -- into
  a 12-cell scheme, and attach the downstream action (add / merge / keep / ignore
  / relax / clash / dedup / review). Two of the twelve cells are the normal
  states; the rest are the things to fix.
"""

from __future__ import annotations

import functools
import re
from typing import cast

import lmdb
import numpy as np
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.readers.io import load_bytes as _load_bytes
from structcooker.instructions.transforms.qc import (
    _COVALENT,
    _COVALENT_DEFAULT,
    _water_atom_mask,
    grid_pairs,
)

_MAX_EX = 6
_BOND_ORDER = {"SING": 1.0, "DOUB": 2.0, "TRIP": 3.0, "QUAD": 4.0, "AROM": 1.5,
               "1": 1.0, "2": 2.0, "3": 3.0, "": 1.0}


def _order(value: str) -> float:
    return _BOND_ORDER.get(str(value).upper(), 1.0)


def _linker_for_comp(atoms: dict, bonds: dict | None) -> dict:
    """Linker atoms + intra ideal bond lengths for one component (H kept here)."""
    aid = [str(a) for a in atoms["atom_id"]]
    elem = [str(e).upper() for e in atoms["type_symbol"]]
    leaving = [str(f).upper() == "Y" for f in atoms.get("pdbx_leaving_atom_flag", ["N"] * len(aid))]
    idx = {a: i for i, a in enumerate(aid)}
    xyz = None
    if "model_Cartn_x" in atoms:
        def _f(vals: list) -> np.ndarray:
            return np.array([float(v) if str(v) not in ("?", ".", "") else np.nan for v in vals])
        xyz = np.stack([_f(atoms["model_Cartn_x"]), _f(atoms["model_Cartn_y"]), _f(atoms["model_Cartn_z"])], axis=1)

    freed: dict[str, float] = {}
    intra: list[list] = []
    if bonds is not None:
        b1 = [str(a) for a in bonds["atom_id_1"]]
        b2 = [str(a) for a in bonds["atom_id_2"]]
        bo = [_order(v) for v in bonds["value_order"]]
        for a1, a2, o in zip(b1, b2, bo, strict=True):
            if a1 not in idx or a2 not in idx:
                continue
            i1, i2 = idx[a1], idx[a2]
            l1, l2 = leaving[i1], leaving[i2]
            h1, h2 = elem[i1] != "H", elem[i2] != "H"
            # non-leaving heavy atom bonded to a leaving atom -> gains free valence
            if l2 and not l1 and h1:
                freed[a1] = freed.get(a1, 0.0) + o
            if l1 and not l2 and h2:
                freed[a2] = freed.get(a2, 0.0) + o
            # intra ideal length for heavy-heavy, non-leaving bonds
            if h1 and h2 and not l1 and not l2 and xyz is not None:
                d = float(np.linalg.norm(xyz[i1] - xyz[i2]))
                if np.isfinite(d):
                    intra.append([a1, a2, o, round(d, 3)])

    linkers = {a: [elem[idx[a]], freed[a]] for a in freed}
    return {"linkers": linkers, "intra_bonds": intra}


def build_ccd_linker(
    chem_comp_atom_dict: dict | None,
    chem_comp_bond_dict: dict | None,
    *,
    unwrap: bool = True,
) -> dict:
    """Build per-component linker descriptors from grouped CCD tables."""
    atom_dict = chem_comp_atom_dict or {}
    bond_dict = chem_comp_bond_dict or {}
    out = {cid: _linker_for_comp(atoms, bond_dict.get(cid)) for cid, atoms in atom_dict.items()}
    if unwrap:
        if len(out) != 1:
            return out
        return next(iter(out.values()))
    return out


# --- per-entry connectivity QC (uses the ccd_linker reference) ---------------
_QUERY_R = 2.6          # neighbour search radius (covers bond + clash range)
_BOND_TOL = 0.45        # bond-length tolerance around covalent-radii sum
_BROKEN_FACTOR = 1.6    # recorded bond broken if dist > ideal * factor
_CLASH_FACTOR = 0.65    # non-bondable pair clashes if dist < cov_sum * factor
_METALS = frozenset({
    "ZN", "MG", "CA", "FE", "MN", "CU", "NA", "K", "NI", "CO", "CD", "MO", "HG",
    "PB", "BA", "SR", "CS", "RB", "LI", "AL", "V", "W", "PT", "AU", "AG", "HF",
})


@functools.lru_cache(maxsize=1)
def _linker_env(path: str) -> lmdb.Environment:
    return open_env(path, readonly=True, lock=False, subdir=True, max_dbs=0)


@functools.lru_cache(maxsize=200_000)
def _linker_caps(path: str, comp: str) -> dict:
    """Per-atom inter-residue bond capacity (freed valence) for a component."""
    with _linker_env(path).begin() as txn:
        raw = txn.get(comp.encode())
    if raw is None:
        return {}
    return {a: v[1] for a, v in _load_bytes(bytes(raw))["ccd_linker"]["linkers"].items()}


def _ideal(e1: str, e2: str) -> float:
    return _COVALENT.get(e1, _COVALENT_DEFAULT) + _COVALENT.get(e2, _COVALENT_DEFAULT)


# 12-cell scheme: (bondable, distance regime, recorded) -> cell number.
# Cells 1 (normal recorded bond) and 9 (unrelated far pair) are the normal
# states; cell 4 (free-valence atom far + unrecorded = a chain terminus) is also
# benign. Everything else is something to fix. Far+unrecorded pairs (4, 9) are
# never produced -- the neighbour grid only reaches bonding/clash distance.
_CELL = {
    (True, "CLOSE", True): 1, (True, "CLOSE", False): 2,
    (True, "FAR", True): 3, (True, "FAR", False): 4,
    (True, "TOOCLOSE", True): 5, (True, "TOOCLOSE", False): 6,
    (False, "CLOSE", True): 7, (False, "CLOSE", False): 8,
    (False, "FAR", True): 10, (False, "FAR", False): 9,
    (False, "TOOCLOSE", True): 11, (False, "TOOCLOSE", False): 12,
}
_EMIT = (2, 3, 5, 6, 7, 8, 10, 11, 12)  # the non-normal cells we report
_QC_KEYS = tuple(f"n_cell{c}" for c in _EMIT)
_CELL_REF_GAP = 8   # B-, bond-band, unrecorded: likely a bond the reference misses
_CELL_CLASH = 12    # B-, overlapping, unrecorded: steric clash

_COVALENT_CONN = frozenset({
    "covale", "disulf", "covale_base", "covale_phosphate", "covale_sugar", "modres",
})
_STRUCT_CONN_NDIM = 2  # struct_conn edge value is [conn_type, bond_order]
_EMPTY: frozenset[int] = frozenset()
_BB_PEP = frozenset({"C", "N"})
_BB_NUC_P = frozenset({"P"})
_BB_NUC_O = frozenset({"O3'", "O5'", "O3*", "O5*"})


def _auth_num(x: object) -> int | None:
    m = re.match(r"-?\d+", str(x))
    return int(m.group()) if m else None


def _regime(d: float, ideal: float, *, recorded: bool) -> str | None:
    """Distance band: too-close (clash) / close (bond band) / far / skip."""
    if d < ideal * _CLASH_FACTOR:
        return "TOOCLOSE"
    if ideal - _BOND_TOL <= d <= ideal + _BOND_TOL:
        return "CLOSE"
    if d > ideal * _BROKEN_FACTOR:
        return "FAR"
    # in-between: a recorded bond is just a slightly long/short normal bond;
    # an unrecorded pair here is a benign near-contact -> skip.
    return "CLOSE" if recorded else None


def _action(cell: int, *, cross_op: bool, same_chain: bool, backbone_adj: bool) -> str:
    """Downstream handling for a classified pair."""
    if cell in (2, 6):  # should-be bond, not recorded
        if cross_op:
            return "ignore"      # symmetry-operator copy coincidence
        if same_chain:
            return "add"         # unrecorded intra-chain bond (e.g. cyclic closure)
        return "merge" if backbone_adj else "keep"
    if cell in (3, 5):           # recorded bond, geometry off (long / compressed)
        return "relax"
    if cell in (7, 10, 11):      # recorded bond neither end can make per CCD
        return "review"
    if cell == _CELL_CLASH:      # overlapping non-bondable, unrecorded -> steric clash
        return "dedup" if cross_op else "clash"
    if cell == _CELL_REF_GAP:    # at bond distance but not bondable, unrecorded:
        return "dedup" if cross_op else "review"  # likely a bond the reference misses
    return "review"


def connectivity_qc(assembly: dict, ccd_linker_db_path: str) -> dict:
    """Classify one assembly's inter-residue pairs into the 12-cell scheme.

    Returns ``{"cells": {cell: count}, "actions": {action: count},
    "examples": {cell: [...]}}`` over the non-normal cells.
    """
    atoms = assembly["atoms"]
    nd = atoms["nodes"]
    xyz = np.array(nd["xyz"]["value"], dtype=np.float64)
    elem = np.char.upper(np.asarray(nd["element"]["value"]).astype("U2"))
    aid = np.asarray(nd["id"]["value"]).astype("U6")
    xyz[_water_atom_mask(assembly)] = np.nan
    it = assembly["index_table"]
    a2r = np.asarray(it["atom_to_res"])
    r2c = np.asarray(it["res_to_chain"])
    comp = np.asarray(assembly["residues"]["nodes"]["chem_comp_id"]["value"]).astype("U6")
    rauth = np.asarray(assembly["residues"]["nodes"]["auth_idx"]["value"]).astype(str)
    chnodes = assembly["chains"]["nodes"]
    clab = np.asarray(chnodes.get("chain_id", chnodes["auth_asym_id"])["value"]).astype(str)

    recorded: set[tuple[int, int]] = set()
    bt = atoms["edges"].get("bond_type")
    if bt is not None:
        for s, d in zip(np.asarray(bt["src_indices"]).tolist(),
                        np.asarray(bt["dst_indices"]).tolist(), strict=True):
            recorded.add((min(s, d), max(s, d)))
    # struct_conn carries non-covalent links too (its value is [conn_type, order]);
    # keep only covalent connections -- hydrogen bonds (base pairs), salt bridges
    # and metal coordination are not bonds.
    sc = atoms["edges"].get("struct_conn")
    if sc is not None:
        src = np.asarray(sc["src_indices"]).tolist()
        dst = np.asarray(sc["dst_indices"]).tolist()
        val = np.asarray(sc["value"])
        ctype = (val[:, 0] if val.ndim == _STRUCT_CONN_NDIM else val).astype(str)
        for s, d, ct in zip(src, dst, ctype.tolist(), strict=True):
            if str(ct).lower() in _COVALENT_CONN:
                recorded.add((min(s, d), max(s, d)))

    # Per-atom inter-residue bond capacity (CCD freed valence) and how much is
    # already consumed by recorded inter-residue bonds. A RECORDED bond is
    # bondable if both ends have capacity (cap>0); an UNRECORDED pair is bondable
    # only if both ends still have *free* capacity (cap - used > 0).
    n_at = xyz.shape[0]
    cap = np.zeros(n_at)
    comp_caps: dict[str, dict] = {}
    for a in range(n_at):
        cmp = str(comp[a2r[a]])
        cc = comp_caps.get(cmp)
        if cc is None:
            cc = _linker_caps(ccd_linker_db_path, cmp)
            comp_caps[cmp] = cc
        cap[a] = cc.get(str(aid[a]), 0.0)  # disulfides excluded (S not auto-linker)
    used = np.zeros(n_at)
    for s, d in recorded:
        if a2r[s] != a2r[d]:
            used[s] += 1
            used[d] += 1
    free = cap - used

    cells: dict[int, int] = dict.fromkeys(_EMIT, 0)
    actions: dict[str, int] = {}
    examples: dict[int, list] = {}

    def loc(i: int) -> str:
        r = a2r[i]
        return f"{comp[r]}{rauth[r]}.{aid[i]}({clab[r2c[r]]})"

    def classify(i: int, j: int, d: float, *, recorded_pair: bool) -> None:
        ei, ej = str(elem[i]), str(elem[j])
        if ei in _METALS or ej in _METALS:
            return  # metal coordination, not a covalent CCD bond
        is_ss = ei == "S" and ej == "S"
        if is_ss and not recorded_pair:
            return  # unrecorded disulfide: excluded by design
        regime = _regime(d, _ideal(ei, ej), recorded=recorded_pair)
        if regime is None:
            return
        if is_ss:
            bondable = True  # a recorded S-S is a normal disulfide (-> cell 1)
        elif recorded_pair:
            # a recorded bond is a real linkage if either end is a CCD linker;
            # the acceptor side (e.g. ASN.ND2 in N-glycosylation) need not be one.
            bondable = cap[i] > 0 or cap[j] > 0
        else:
            # proposing a NEW bond needs free valence on BOTH ends -- this is what
            # stops an already-peptide-bonded backbone C from pairing with a
            # coincidentally-close N of another residue.
            bondable = free[i] > 0 and free[j] > 0
        cell = _CELL[(bondable, regime, recorded_pair)]
        if cell not in cells:
            return  # normal cell (1) or benign far (4, 9)
        ri, rj = a2r[i], a2r[j]
        ci, cj = r2c[ri], r2c[rj]
        full_i, full_j = str(clab[ci]), str(clab[cj])
        cross_op = full_i != full_j and full_i.split("_")[0] == full_j.split("_")[0]
        same_chain = full_i == full_j
        names = {str(aid[i]), str(aid[j])}
        ni, nj = _auth_num(rauth[ri]), _auth_num(rauth[rj])
        adj = ni is not None and nj is not None and abs(ni - nj) == 1
        is_bb = names == _BB_PEP or ((names & _BB_NUC_P) and (names & _BB_NUC_O))
        backbone_adj = bool(is_bb and adj)
        act = _action(cell, cross_op=cross_op, same_chain=same_chain, backbone_adj=backbone_adj)
        cells[cell] += 1
        actions[act] = actions.get(act, 0) + 1
        ex = examples.setdefault(cell, [])
        if len(ex) < _MAX_EX:
            ex.append(f"{loc(i)} -- {loc(j)} d={d:.2f} [{act}]")

    # (1) recorded inter-residue bonds (found by iterating edges: a broken bond
    #     is long, beyond the neighbour radius).
    for s, d in recorded:
        if a2r[s] == a2r[d] or not (np.isfinite(xyz[s]).all() and np.isfinite(xyz[d]).all()):
            continue
        classify(s, d, float(np.linalg.norm(xyz[s] - xyz[d])), recorded_pair=True)

    # (2) close NON-recorded neighbour pairs. Exclude 1-2 (recorded) and 1-3
    #     (share a bonded neighbour, e.g. backbone O...N geminal across a peptide
    #     bond) -- those are not contacts, just bond geometry.
    nbr: dict[int, set[int]] = {}
    for s, d in recorded:
        nbr.setdefault(s, set()).add(d)
        nbr.setdefault(d, set()).add(s)
    gi, gj, dist = grid_pairs(xyz, _QUERY_R)
    for i, j, d in zip(gi.tolist(), gj.tolist(), dist.tolist(), strict=True):
        if a2r[i] == a2r[j] or (i, j) in recorded:
            continue
        if nbr.get(i, _EMPTY) & nbr.get(j, _EMPTY):
            continue  # 1-3 neighbour
        classify(i, j, d, recorded_pair=False)

    return {
        "cells": {c: v for c, v in cells.items() if v},
        "actions": actions,
        "examples": {c: v for c, v in examples.items() if v},
    }


def connectivity_qc_record(record: dict, ccd_linker_db_path: str) -> dict:
    """Aggregate :func:`connectivity_qc` over all assemblies of one cif record.

    Returns ``{}`` for a clean entry; otherwise per-cell counts (``n_cell<N>``),
    an ``actions`` tally and located example pairs -- a per-item "what is wrong
    here, and what to do about it" annotation.
    """
    cells = dict.fromkeys(_EMIT, 0)
    actions: dict[str, int] = {}
    examples: dict[int, list] = {}
    for assembly in record["assembly_dict"].values():
        res = connectivity_qc(assembly, ccd_linker_db_path)
        for c, v in res["cells"].items():
            cells[c] += v
        for a, v in res["actions"].items():
            actions[a] = actions.get(a, 0) + v
        for c, items in res["examples"].items():
            examples.setdefault(c, []).extend(items)
    if not any(cells.values()):
        return {}
    out: dict = {f"n_cell{c}": cells[c] for c in _EMIT if cells[c]}
    out["actions"] = {a: v for a, v in actions.items() if v}
    out["examples"] = {f"cell{c}": v[:_MAX_EX] for c, v in examples.items() if v}
    return out


# --- 3-tier bondability (free-valence) connectivity QC ------------------------
import json as _json  # noqa: E402

_HETERO = frozenset({"N", "O", "S", "P", "SE"})


@functools.cache
def _valence_env(path: str) -> lmdb.Environment:
    return open_env(path, readonly=True, lock=False, subdir=True, max_dbs=0)


@functools.cache
def _valence_comp(path: str, comp: str) -> tuple:
    with _valence_env(path).begin() as txn:
        raw = txn.get(comp.encode())
    if raw is None:
        return ()
    return tuple((a, v[0], float(v[1]), int(v[2])) for a, v in _json.loads(cast("bytes", raw)).items())


def connectivity_qc_3tier(assembly: dict, ccd_valence_db_path: str,
                          ccd_linker_db_path: str = "/data/psk6950/CCD/ccd_linker.lmdb") -> dict:
    """Run connectivity QC with 3-tier bondability (결합 불가능 / 가능성 / 필요).

    Per atom, from the deposited heavy bonds + CCD ideal valence:
      slack = V_ideal - sum(bond_order of present heavy bonds)
      req  = max(0, slack - nH)                 dangling heavy valence -> MUST bond
      opt  = min(nH, slack) if heteroatom else 0   displaceable H -> MAY bond
    Pair tiers (unrecorded, bond-band): required(one req + partner accepts) -> missing;
    possible(both opt) -> candidate; else -> review/clash. Recorded bonds are legit
    if either end had room (heme/glycan -> normal, not review).
    Returns {"actions": {...}, "tiers": {...}, "examples": {...}} or {} if clean.
    """
    atoms = assembly["atoms"]
    nd = atoms["nodes"]
    xyz = np.array(nd["xyz"]["value"], dtype=np.float64)
    elem = np.char.upper(np.asarray(nd["element"]["value"]).astype("U2"))
    aid = np.asarray(nd["id"]["value"]).astype("U6")
    xyz[_water_atom_mask(assembly)] = np.nan
    it = assembly["index_table"]
    a2r = np.asarray(it["atom_to_res"])
    r2c = np.asarray(it["res_to_chain"])
    comp = np.asarray(assembly["residues"]["nodes"]["chem_comp_id"]["value"]).astype("U6")
    rauth = np.asarray(assembly["residues"]["nodes"]["auth_idx"]["value"]).astype(str)
    chn = assembly["chains"]["nodes"]
    clab = np.asarray(chn.get("chain_id", chn["auth_asym_id"])["value"]).astype(str)
    n_at = xyz.shape[0]

    # all bonds (+order), deduped by (min,max) so bidirectional edges count once
    bond_order: dict = {}
    bt = atoms["edges"].get("bond_type")
    if bt is not None:
        s = np.asarray(bt["src_indices"]).tolist()
        d = np.asarray(bt["dst_indices"]).tolist()
        vv = np.asarray(bt["value"]).astype(str).tolist()
        for si, di, vo in zip(s, d, vv, strict=True):
            bond_order[(min(si, di), max(si, di))] = _order(vo)
    sc = atoms["edges"].get("struct_conn")
    if sc is not None:
        s = np.asarray(sc["src_indices"]).tolist()
        d = np.asarray(sc["dst_indices"]).tolist()
        val = np.asarray(sc["value"])
        ct = (val[:, 0] if val.ndim == _STRUCT_CONN_NDIM else val).astype(str).tolist()
        for si, di, c in zip(s, d, ct, strict=True):
            if str(c).lower() in _COVALENT_CONN:
                bond_order.setdefault((min(si, di), max(si, di)), 1.0)
    bonded = np.zeros(n_at)
    recorded: dict = {}
    for (si, di), o in bond_order.items():
        bonded[si] += o
        bonded[di] += o
        if a2r[si] != a2r[di]:
            recorded[(si, di)] = o

    # used = recorded inter-residue bonds per atom
    used = np.zeros(n_at)
    for (si, di) in recorded:
        used[si] += 1
        used[di] += 1
    # req capacity: leaving-flag linker (CCD), like the binary scheme -> lfree = cap - used.
    #   captures backbone/link atoms (C via OXT, N via H2, O3'/P, glycosidic C1) that CAN
    #   form an inter-residue bond even while their leaving group is still modelled.
    # opt capacity: displaceable heteroatom-H (N/O/S/P with H in CCD) -> disulfide, glycan
    #   acceptor etc. that the leaving-flag misses.
    lfree = np.zeros(n_at)
    hetavail = np.zeros(n_at, dtype=bool)
    caps: dict = {}
    vcache: dict = {}
    for a in range(n_at):
        cmp = str(comp[a2r[a]])
        cc = caps.get(cmp)
        if cc is None:
            cc = _linker_caps(ccd_linker_db_path, cmp)
            caps[cmp] = cc
        lfree[a] = cc.get(str(aid[a]), 0.0) - used[a]
        vc = vcache.get(cmp)
        if vc is None:
            vc = {t[0]: t[3] for t in _valence_comp(ccd_valence_db_path, cmp)}  # atom -> nH
            vcache[cmp] = vc
        hetavail[a] = str(elem[a]) in _HETERO and vc.get(str(aid[a]), 0) > 0

    actions: dict = {}
    tiers: dict = {}
    examples: dict = {}

    def loc(i: int) -> str:
        r = a2r[i]
        return f"{comp[r]}{rauth[r]}.{aid[i]}({clab[r2c[r]]})"

    cells: dict = {}

    def emit(cell: str, act: str, tier: str, i: int, j: int, d: float) -> None:
        cells[cell] = cells.get(cell, 0) + 1
        actions[act] = actions.get(act, 0) + 1
        tiers[tier] = tiers.get(tier, 0) + 1
        ex = examples.setdefault(cell, [])
        if len(ex) < _MAX_EX:
            ex.append(f"{loc(i)} -- {loc(j)} d={d:.2f} [{act}]")

    def cross_op(i: int, j: int) -> bool:
        ci, cj = str(clab[r2c[a2r[i]]]), str(clab[r2c[a2r[j]]])
        return ci != cj and ci.split("_")[0] == cj.split("_")[0]

    def chain_action(i: int, j: int) -> str:
        if cross_op(i, j):
            return "ignore"
        ci, cj = str(clab[r2c[a2r[i]]]), str(clab[r2c[a2r[j]]])
        if ci == cj:
            return "add"
        ni, nj = _auth_num(rauth[a2r[i]]), _auth_num(rauth[a2r[j]])
        adj = ni is not None and nj is not None and abs(ni - nj) == 1
        names = {str(aid[i]), str(aid[j])}
        is_bb = names == _BB_PEP or ((names & _BB_NUC_P) and (names & _BB_NUC_O))
        return "merge" if (is_bb and adj) else "keep"

    def classify(i: int, j: int, d: float, rec: bool) -> None:
        # cell = {tier}_{O|X}_{regime}. req = leaving-flag linker (backbone/link, incl bb-break),
        # opt = displaceable heteroatom-H (disulfide/glycan), imp = neither.
        ei, ej = str(elem[i]), str(elem[j])
        if ei in _METALS or ej in _METALS:
            return
        regime = _regime(d, _ideal(ei, ej), recorded=rec)
        if regime is None:
            return
        is_ss = ei == "S" and ej == "S"
        if rec:
            cap_i, cap_j = lfree[i] + used[i], lfree[j] + used[j]   # cap = lfree + used
            bondable = (is_ss or cap_i > 0 or cap_j > 0
                        or hetavail[i] or hetavail[j])              # heme/glycan acceptor 포함
            if not bondable:
                emit("imp_O_" + regime, "review", "imp", i, j, d)  # 불가능인데 기록됨
                return
            if regime == "CLOSE":
                return                                              # 정상 결합
            tier = "req" if (lfree[i] > 0 or lfree[j] > 0) else "opt"
            emit(f"{tier}_O_{regime}", "relax", tier, i, j, d)      # broken(FAR)/과압축(TOOCLOSE)
            return
        if regime == "FAR":
            return                                                 # 미기록·먼 접촉 = 정상
        if is_ss:                                                  # 미기록 이황화 → candidate
            emit("opt_X_" + regime, "ignore" if cross_op(i, j) else "candidate", "opt", i, j, d)
            return
        required = lfree[i] > 0 and lfree[j] > 0                    # 둘 다 leaving-linker → missing
        possible = (not required) and hetavail[i] and hetavail[j]   # 둘 다 헤테로-H → candidate
        if required:
            emit("req_X_" + regime, chain_action(i, j), "req", i, j, d)   # merge/add/keep/ignore
        elif possible:
            emit("opt_X_" + regime, "ignore" if cross_op(i, j) else "candidate", "opt", i, j, d)
        elif regime == "TOOCLOSE":
            emit("imp_X_TOOCLOSE", "dedup" if cross_op(i, j) else "clash", "imp", i, j, d)
        else:
            emit("imp_X_CLOSE", "dedup" if cross_op(i, j) else "review", "imp", i, j, d)

    for s, d0 in recorded:
        if a2r[s] == a2r[d0] or not (np.isfinite(xyz[s]).all() and np.isfinite(xyz[d0]).all()):
            continue
        classify(s, d0, float(np.linalg.norm(xyz[s] - xyz[d0])), True)
    nbr: dict = {}
    for s, d0 in recorded:
        nbr.setdefault(s, set()).add(d0)
        nbr.setdefault(d0, set()).add(s)
    gi, gj, dist = grid_pairs(xyz, _QUERY_R)
    for i, j, d in zip(gi.tolist(), gj.tolist(), dist.tolist(), strict=True):
        if a2r[i] == a2r[j] or (i, j) in recorded:
            continue
        if nbr.get(i, _EMPTY) & nbr.get(j, _EMPTY):
            continue
        classify(i, j, d, False)

    if not cells:
        return {}
    return {"cells": cells, "actions": actions, "tiers": tiers,
            "examples": {c: v for c, v in examples.items() if v}}


def _slack_plus(
    vcache: dict,
    comp: np.ndarray,
    a2r: np.ndarray,
    aid: np.ndarray,
    bonded: np.ndarray,
    i: int,
    o: float,
) -> float:
    """Return the CCD ideal-valence slack for atom ``i`` if bond order ``o`` were absent."""
    vc = vcache.get(str(comp[a2r[i]]), {})
    info = vc.get(str(aid[i]))
    if info is None:
        return 0.0
    return info[1] - (bonded[i] - o)


def _valV(*a: object) -> float:  # noqa: N802, ARG001 -- placeholder (unused legit path)
    return 0.0
