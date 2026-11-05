"""CIFMol -> interface nodes (arrays for ``detect``), plus unit extraction.

Node rules
- polymer chain -> node ``protein`` (polypeptide) or ``na`` (DNA / RNA / hybrid); other polymer types are skipped
- branched chain (glycan) / non-polymer residue -> node ``ligand``; a single-atom non-polymer -> ``ion``
- water (HOH / DOD) is never a node; it is returned separately
- anything in the exclusion set, or a ligand below ``min_heavy_atoms``, is dropped (atoms ignored)
- a ligand / glycan covalently bonded (struct_conn ``covale``) to a polymer is merged into that
  polymer's node: one party. ``metalc`` / ``hydrog`` / ``disulf`` never merge.
Missing (NaN) atoms are ignored in contacts. Modified residues are mapped to their parent via the
CCD canonical one-letter code (``one_letter_code_can``); unknown -> ``X`` / ``N``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from biomol.exceptions import FeatureKeyError

NA_TYPES = {"polydeoxyribonucleotide", "polyribonucleotide", "polydeoxyribonucleotide/polyribonucleotide hybrid"}
WATER = {"HOH", "DOD"}


def pick_alt_keys(keys) -> list[str]:
    """One key per (assembly, model): the no-altloc ``.`` entry if present, else the first altloc.

    Altloc entries in cif_pdb are full copies of the assembly (same atoms, alternate conformers),
    so one is enough; the ``.`` entry is the reference conformer.
    """
    best: dict[tuple[str, str], str] = {}
    for k in sorted(keys):
        a, m, alt = k.split("_")
        if (a, m) not in best or alt == ".":
            best[(a, m)] = k
    return list(best.values())


@dataclass
class Assembly:
    xyz: np.ndarray           # (n, 3) atoms that belong to a node (may hold NaN)
    node_idx: np.ndarray      # (n,) compact node index
    res_idx: np.ndarray       # (n,) global residue index
    atom_ids: np.ndarray      # (n,) index into the CIFMol atom axis
    kinds: list[str]
    keys: list[str]           # dedup key: "P:<seq>" / "N:<seq>" / "L:<ccd[-ccd]>" / "I:<ion>"
    labels: list[str]         # human label, e.g. chain id(s)
    atom_name: np.ndarray     # (n,) atom ids (N, CA, C1' ...)
    res_code: np.ndarray      # (n,) parent one-letter code of the atom's residue
    donor: np.ndarray         # (n,) bool: N / O / S (coordination donors for ions)
    water_xyz: np.ndarray
    water_atom_ids: np.ndarray
    stats: dict = field(default_factory=dict)


def build_assembly(cm, excluded: frozenset[str], min_heavy: int = 6) -> Assembly:
    ch, rs, at, ix = cm.chains, cm.residues, cm.atoms, cm.index_table
    n_atom = len(at.xyz.value)
    all_atoms = np.arange(n_atom)
    atom_res = ix.atoms_to_residues(all_atoms)
    xyz_all = at.xyz.value
    elem = at.element.value
    comp = rs.chem_comp_id.value
    one = rs.one_letter_code_can.value
    etype = ch.entity_type.value
    cid = ch.chain_id.value

    kinds, keys, labels, node_atoms = [], [], [], []
    dropped = {"excluded": 0, "small": 0, "skipped_polymer": 0}
    excluded_names: dict[str, int] = {}
    water_atoms: list[np.ndarray] = []
    for c in range(len(cid)):
        et = str(etype[c])
        ridx = ix.chains_to_residues(np.array([c]))
        if et.startswith("polypeptide") or et in NA_TYPES:
            kind = "protein" if et.startswith("polypeptide") else "na"
            seq = "".join(one[ridx].tolist())
            kinds.append(kind); keys.append(("P:" if kind == "protein" else "N:") + seq)
            labels.append(str(cid[c])); node_atoms.append(ix.chains_to_atoms(np.array([c])))
        elif et == "branched":
            names = sorted(set(comp[ridx].tolist()))
            atoms = ix.chains_to_atoms(np.array([c]))
            if len(atoms) < min_heavy:
                dropped["small"] += 1; continue
            kinds.append("ligand"); keys.append("L:" + "-".join(names)); labels.append(str(cid[c])); node_atoms.append(atoms)
        elif et in ("non-polymer", "water"):
            for r in ridx:
                name = str(comp[r])
                atoms = ix.residues_to_atoms(np.array([r]))
                if name in WATER:
                    water_atoms.append(atoms); continue
                if name in excluded:
                    dropped["excluded"] += 1
                    excluded_names[name] = excluded_names.get(name, 0) + 1
                    continue
                if len(atoms) == 1:
                    kinds.append("ion"); keys.append("I:" + name)
                else:
                    heavy = int((elem[atoms] != "H").sum())
                    if heavy < min_heavy:
                        dropped["small"] += 1; continue
                    kinds.append("ligand"); keys.append("L:" + name)
                labels.append(f"{cid[c]}:{name}"); node_atoms.append(atoms)
        else:
            dropped["skipped_polymer"] += 1

    # covalent merge of ligands into polymers
    node_of_atom = np.full(n_atom, -1, dtype=np.int64)
    for n, a in enumerate(node_atoms):
        node_of_atom[a] = n
    parent = list(range(len(kinds)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    try:
        sc = at.struct_conn
    except FeatureKeyError:  # entries without any struct_conn row lack the feature
        sc = None
    n_merged = 0
    if sc is not None and len(sc.src_indices):
        types = sc.value[:, 0]
        for s, d, t in zip(sc.src_indices, sc.dst_indices, types, strict=True):
            if t != "covale":
                continue
            a, b = node_of_atom[s], node_of_atom[d]
            if a < 0 or b < 0 or a == b:
                continue
            ka, kb = kinds[a], kinds[b]
            poly = {"protein", "na"}
            if ka in poly and kb in poly:
                continue
            if "ion" in (ka, kb):
                continue
            # polymer wins as root; otherwise the lower index
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            n_merged += 1
            if kinds[ra] in poly:
                parent[rb] = ra
            elif kinds[rb] in poly:
                parent[ra] = rb
            else:
                parent[max(ra, rb)] = min(ra, rb)
    roots = sorted({find(n) for n in range(len(kinds))})
    remap = {r: i for i, r in enumerate(roots)}
    final_of_node = np.array([remap[find(n)] for n in range(len(kinds))], dtype=np.int64)
    if len(final_of_node):
        node_idx_all = np.where(node_of_atom >= 0, final_of_node[np.maximum(node_of_atom, 0)], -1)
    else:  # every non-polymer was excluded and there is no polymer
        node_idx_all = np.full(n_atom, -1, dtype=np.int64)
    sel = np.flatnonzero(node_idx_all >= 0)
    w_atoms = np.concatenate(water_atoms) if water_atoms else np.empty(0, dtype=np.int64)
    w_atoms = w_atoms[elem[w_atoms] == "O"] if len(w_atoms) else w_atoms
    return Assembly(
        xyz_all[sel], node_idx_all[sel], atom_res[sel], sel,
        [kinds[r] for r in roots], [keys[r] for r in roots],
        [labels[r] for r in roots], at.id.value[sel], one[atom_res[sel]],
        np.isin(elem[sel], ["N", "O", "S"]), xyz_all[w_atoms], w_atoms,
        {"n_nodes_raw": len(kinds), "dropped": dropped, "excluded_names": excluded_names, "n_water": len(w_atoms), "n_merged": n_merged},
    )
