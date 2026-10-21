"""Double-stranded nucleic-acid detection from base pairing (DNA, RNA, hybrids).

Two residues of different NA chains pair when their C1' atoms are 8.0-12.0 A apart (a base pair
spans ~10.5 A; stacked neighbours and backbone contacts are ~4-7 A) and at least ``min_hbonds``
base N/O atom pairs lie within ``hb_cutoff``. Each residue keeps its best partner. Two chains are
a duplex when they share >= ``min_bp`` such pairs; mismatches, wobble and Hoogsteen pairs count.
This is a geometric stand-in for DSSR / MC-Annotate (no external tool or licence); it is checked
against known duplexes in ``scripts/interface/run_sample.py``.
"""

from __future__ import annotations

import numpy as np

from structcooker.interface.detect import inter_node_pairs

BASE_ATOMS = {"N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"}
PYRIMIDINES = {"C", "U", "T"}


def find_duplex(
    xyz: np.ndarray,          # (n, 3) atoms of the assembly's nodes
    node_idx: np.ndarray,
    res_idx: np.ndarray,
    atom_name: np.ndarray,
    res_code: np.ndarray,     # parent one-letter code per atom (A/C/G/U/T/N)
    is_na: np.ndarray,        # (n,) bool, atom belongs to an NA node
    min_bp: int = 4,
    min_hbonds: int = 2,
    hb_cutoff: float = 3.5,
    c1_range: tuple[float, float] = (8.0, 12.0),
) -> dict[tuple[int, int], int]:
    """``{(node_a, node_b): n_base_pairs}`` for NA node pairs with >= ``min_bp`` pairs (a < b)."""
    c1 = np.flatnonzero(is_na & (atom_name == "C1'"))
    if len(c1) < 2:
        return {}
    i, j, d = inter_node_pairs(xyz[c1], node_idx[c1], c1_range[1])
    ok = d >= c1_range[0]
    i, j, d = c1[i[ok]], c1[j[ok]], d[ok]
    if len(i) == 0:
        return {}
    base = is_na & np.isin(atom_name, list(BASE_ATOMS))
    base &= ~((atom_name == "N1") & np.isin(res_code, list(PYRIMIDINES)))  # N1 is glycosidic there
    by_res: dict[int, list] = {}
    for a in np.flatnonzero(base):
        by_res.setdefault(int(res_idx[a]), []).append(a)
    by_res = {r: np.array(v) for r, v in by_res.items()}
    scored = []
    for a, b, dist in zip(i.tolist(), j.tolist(), d.tolist(), strict=True):
        ra, rb = int(res_idx[a]), int(res_idx[b])
        if ra not in by_res or rb not in by_res:
            continue
        A, B = by_res[ra], by_res[rb]
        dd = np.linalg.norm(xyz[A][:, None] - xyz[B][None], axis=-1)
        n_hb = int((dd <= hb_cutoff).sum())
        if n_hb >= min_hbonds:
            scored.append((-n_hb, abs(dist - 10.5), ra, rb, int(node_idx[a]), int(node_idx[b])))
    scored.sort()
    used: set[int] = set()
    counts: dict[tuple[int, int], int] = {}
    for _, _, ra, rb, na, nb in scored:
        if ra in used or rb in used:
            continue
        used.update((ra, rb))
        key = (min(na, nb), max(na, nb))
        counts[key] = counts.get(key, 0) + 1
    return {k: v for k, v in counts.items() if v >= min_bp}
