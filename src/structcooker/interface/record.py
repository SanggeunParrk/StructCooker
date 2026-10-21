"""Interface unit -> self-contained record (a dict of numpy arrays + JSON-able scalars).

Schema v1, one record per unit (``to_bytes`` / ``from_bytes`` of ``instructions.transforms.codecs``)::

    meta       pdb_id, assembly_id, model_id, alt_id, unit_id, resolution, release_date,
               deposition_date, signature, truncated, n_seed_edges, schema
    nodes      per node (a polymer chain, a ligand with its covalent parts, or an ion):
               kind (protein|na|ligand|ion), key ("P:<seq>" | "N:<seq>" | "L:<ccd[-ccd]>" | "I:<ion>"),
               hash (md5-16 of key), chain_id (asym_oper), auth_asym_id, entity_id, ab_class (ANARCI or "")
    atoms      xyz f4 (n,3) original PDB frame, element, atom_name, node i2, res i4 (index into residues),
               b_factor f4, occupancy f2, contact bool (within the contact cutoff of another node)
    residues   chem_comp_id, parent (canonical 1-letter), cif_idx, auth_idx (author number, may carry an insertion code: '95.A'), node i2, interface bool
               (only residues with >= 1 resolved atom; cif_idx / the node key sequence give the gaps)
    edges      a, b (node indices), n_atom_pairs, n_res_a, n_res_b, min_dist, center (3), res_pair_ptr and
               res_pairs (m, 2) of residue indices (a-side, b-side)
    duplex     a, b, n_bp  (base-paired strand pairs inside the unit; also bundled when not in contact)
    water      xyz f4 (w,3), chain_id, auth_idx (str), bridging bool: every water within ``water.near`` of the unit
               (does NOT depend on resolution; none modelled below ~2.5 A is normal)

Only resolved atoms (finite xyz) are stored. Hydrogens are absent in cif_pdb.
"""

from __future__ import annotations

import hashlib

import numpy as np

from structcooker.interface.detect import select_waters

SCHEMA = 1


def _md5(s: str) -> str:
    return hashlib.md5(s.encode()).hexdigest()[:16]  # noqa: S324


def asm_fields(cm, asm):
    """Per-assembly arrays shared by all units of the assembly."""
    at, rs, ch, ix = cm.atoms, cm.residues, cm.chains, cm.index_table
    f = {
        "elem": at.element.value, "bf": at.b_factor.value, "occ": at.occupancy.value,
        "comp": rs.chem_comp_id.value, "one": rs.one_letter_code_can.value,
        "cif_idx": rs.cif_idx.value, "auth_idx": rs.auth_idx.value,
        "chain_id": ch.chain_id.value, "auth_asym": ch.auth_asym_id.value, "entity": ch.entity_id.value,
        "ix": ix,
    }
    w = asm.water_atom_ids
    ok = np.isfinite(asm.water_xyz).all(axis=1) if len(w) else np.zeros(0, bool)
    f["w_ok"] = ok
    if len(w):
        f["w_chain"] = ix.atoms_to_chains(w)
        f["w_res"] = ix.atoms_to_residues(w)
    return f


def unit_to_record(asm, f, unit, duplex, cfg, meta: dict, ab_class: list, unit_idx: int) -> dict:
    ix = f["ix"]
    k = len(unit.nodes)
    local = {n: i for i, n in enumerate(unit.nodes)}
    node_of_atom = np.full(int(asm.node_idx.max()) + 1, -1, dtype=np.int64)
    for n, i in local.items():
        node_of_atom[n] = i
    atom_local_node = node_of_atom[asm.node_idx]
    m = (atom_local_node >= 0) & np.isfinite(asm.xyz).all(axis=1)
    sel = np.flatnonzero(m)
    pos = np.full(len(asm.xyz), -1, dtype=np.int64)
    pos[sel] = np.arange(len(sel))
    gres = asm.res_idx[sel]
    ures, rinv = np.unique(gres, return_inverse=True)
    r_first = np.unique(rinv, return_index=True)[1]
    atoms_global = asm.atom_ids[sel]

    contact = np.zeros(len(sel), dtype=bool)
    interface = np.zeros(len(ures), dtype=bool)
    e_a, e_b, e_pairs, e_ra, e_rb, e_min, e_cen, ptr, rp = [], [], [], [], [], [], [], [0], []
    for e in unit.edges:
        la, lb = local[e.a], local[e.b]
        ca, cb = pos[e.atoms_a], pos[e.atoms_b]
        ca, cb = ca[ca >= 0], cb[cb >= 0]
        contact[ca] = contact[cb] = True
        pairs = np.stack([np.searchsorted(ures, e.res_pairs[:, 0]), np.searchsorted(ures, e.res_pairs[:, 1])], axis=1)
        interface[pairs.ravel()] = True
        rp.append(pairs.astype(np.int32))
        ptr.append(ptr[-1] + len(pairs))
        e_a.append(la); e_b.append(lb); e_pairs.append(e.n_atom_pairs)
        e_ra.append(len(e.res_a)); e_rb.append(len(e.res_b)); e_min.append(e.min_dist)
        e_cen.append(asm.xyz[np.concatenate([e.atoms_a, e.atoms_b])].mean(axis=0))
    d_a, d_b, d_n = [], [], []
    for (a, b), n_bp in duplex.items():
        if a in local and b in local:
            d_a.append(local[a]); d_b.append(local[b]); d_n.append(n_bp)

    # nodes
    nodes = {"kind": [], "key": [], "hash": [], "chain_id": [], "auth_asym_id": [], "entity_id": [], "ab_class": []}
    first_atom_of_node = {}
    for ai, nl in zip(atoms_global, atom_local_node[sel], strict=True):
        first_atom_of_node.setdefault(int(nl), int(ai))
    for n in unit.nodes:
        c = int(ix.atoms_to_chains(np.array([first_atom_of_node[local[n]]]))[0])
        key = asm.keys[n]
        nodes["kind"].append(asm.kinds[n]); nodes["key"].append(key); nodes["hash"].append(_md5(key))
        nodes["chain_id"].append(str(f["chain_id"][c])); nodes["auth_asym_id"].append(str(f["auth_asym"][c]))
        nodes["entity_id"].append(str(f["entity"][c])); nodes["ab_class"].append(ab_class[n] or "")

    # water context
    w_all = asm.water_xyz[f["w_ok"]] if len(asm.water_xyz) else asm.water_xyz
    w_idx = np.flatnonzero(f["w_ok"]) if len(asm.water_xyz) else np.zeros(0, dtype=np.int64)
    keep, bridging = select_waters(asm.xyz[sel], atom_local_node[sel], w_all, cfg.water_near, cfg.water_bridge)
    wk = w_idx[keep]
    sig = _md5(repr((sorted(nodes["hash"]), sorted(zip(e_a, e_b, [round((x + y) / cfg.unit.weight_bucket) for x, y in zip(e_ra, e_rb, strict=True)], strict=True)))))
    meta = dict(meta, unit_id=f"{meta['pdb_id']}_{meta['assembly_id']}_{meta['model_id']}_{unit_idx}",
                signature=sig, truncated=bool(unit.truncated), n_seed_edges=len(unit.seeds), schema=SCHEMA,
                n_nodes=k)
    return {
        "meta": meta,
        "nodes": nodes,
        "atoms": {
            "xyz": asm.xyz[sel].astype(np.float32), "element": f["elem"][atoms_global].astype("U2"),
            "atom_name": asm.atom_name[sel].astype("U4"), "node": atom_local_node[sel].astype(np.int16),
            "res": rinv.astype(np.int32), "b_factor": f["bf"][atoms_global].astype(np.float32),
            "occupancy": f["occ"][atoms_global].astype(np.float16), "contact": contact,
        },
        "residues": {
            "chem_comp_id": f["comp"][ures].astype("U5"), "parent": f["one"][ures].astype("U1"),
            "cif_idx": f["cif_idx"][ures].astype(np.int32), "auth_idx": f["auth_idx"][ures].astype("U8"),
            "node": atom_local_node[sel][r_first].astype(np.int16), "interface": interface,
        },
        "edges": {
            "a": np.array(e_a, np.int16), "b": np.array(e_b, np.int16), "n_atom_pairs": np.array(e_pairs, np.int32),
            "n_res_a": np.array(e_ra, np.int32), "n_res_b": np.array(e_rb, np.int32),
            "min_dist": np.array(e_min, np.float32), "center": np.array(e_cen, np.float32).reshape(-1, 3),
            "res_pair_ptr": np.array(ptr, np.int64),
            "res_pairs": np.concatenate(rp).astype(np.int32) if rp else np.zeros((0, 2), np.int32),
        },
        "duplex": {"a": np.array(d_a, np.int16), "b": np.array(d_b, np.int16), "n_bp": np.array(d_n, np.int32)},
        "water": {
            "xyz": asm.water_xyz[wk].astype(np.float32) if len(wk) else np.zeros((0, 3), np.float32),
            "chain_id": f["chain_id"][f["w_chain"][wk]].astype("U8") if len(wk) else np.zeros(0, "U8"),
            "auth_idx": f["auth_idx"][f["w_res"][wk]].astype("U8") if len(wk) else np.zeros(0, "U8"),
            "bridging": bridging,
        },
    }


def manifest_row(rec: dict) -> dict:
    """Flat summary of a record (one manifest row)."""
    meta, nodes, at, rs = rec["meta"], rec["nodes"], rec["atoms"], rec["residues"]
    kinds = nodes["kind"]
    ed = rec["edges"]
    row = {
        **{k: meta.get(k) for k in ("unit_id", "pdb_id", "assembly_id", "model_id", "resolution", "release_date",
                                   "deposition_date", "signature", "truncated", "n_seed_edges")},
        "type": "+".join(f"{kinds.count(x)}{x}" for x in ("protein", "na", "ligand", "ion") if x in kinds),
        "n_nodes": len(kinds), "n_protein": kinds.count("protein"), "n_na": kinds.count("na"),
        "n_ligand": kinds.count("ligand"), "n_ion": kinds.count("ion"),
        "n_atoms": len(at["xyz"]), "n_residues": len(rs["chem_comp_id"]),
        "n_interface_res": int(rs["interface"].sum()), "n_contact_atoms": int(at["contact"].sum()),
        "n_atom_pairs": int(ed["n_atom_pairs"].sum()),
        "min_dist": float(ed["min_dist"].min()) if len(ed["min_dist"]) else "",
        "n_water": len(rec["water"]["xyz"]), "n_bridging_water": int(rec["water"]["bridging"].sum()),
        "has_ab": int(any(nodes["ab_class"])), "n_duplex_pairs": len(rec["duplex"]["a"]),
        "node_kinds": "|".join(kinds), "node_hashes": "|".join(nodes["hash"]),
        "node_len": "|".join(str(len(k) - 2) if kd in ("protein", "na") else "0" for k, kd in zip(nodes["key"], kinds, strict=True)),
        "node_ab": "|".join(nodes["ab_class"]),
        "node_ccd": "|".join(k[2:] if kd in ("ligand", "ion") else "" for k, kd in zip(nodes["key"], kinds, strict=True)),
        "node_chains": "|".join(nodes["auth_asym_id"]),
    }
    return row
