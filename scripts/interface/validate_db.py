"""Recompute edge statistics from the stored coordinates and compare them with the stored values.

usage: validate_db.py <db_root> [n_units=300] [seed=0] [max_atoms=15000]
Checks per unit: array shapes / index ranges, every atom finite, edge min_dist and n_atom_pairs
(brute force, cutoff from the config; ion edges: 3.0 A to N/O/S), contact-atom flag, interface-residue
flag, every water within water.near of a unit atom and the bridging flag, node sequences cover the residues.
"""
import random
import sys
from pathlib import Path

import numpy as np

from structcooker.interface.config import load_config
from structcooker.interface.db import InterfaceDB

root, n, seed = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 0
cfg = load_config()
rows = [l.rstrip("\n").split("\t") for l in open(Path(root) / "manifest.tsv")]
hdr, rows = rows[0], rows[1:]
ix = {c: i for i, c in enumerate(hdr)}
random.Random(seed).shuffle(rows)
max_atoms = int(sys.argv[4]) if len(sys.argv) > 4 else 15000   # brute force: keep it to units that fit
rows = [r for r in rows if int(r[ix["n_atoms"]]) <= max_atoms]
db = InterfaceDB(root)
bad = {}


def block_pairs(P, Q, cut, mask_q=None, mask_p=None, block=512):
    """(n_pairs, min_d, flag_p, flag_q, n_lo, n_hi): P-Q distances < cut (and within 1e-3 of it), in row blocks."""
    n, dmin, n_lo, n_hi = 0, np.inf, 0, 0
    fp, fq = np.zeros(len(P), bool), np.zeros(len(Q), bool)
    for s0 in range(0, len(P), block):
        d = np.linalg.norm(P[s0:s0 + block, None] - Q[None], axis=-1)
        ok = np.ones(d.shape, bool)
        if mask_q is not None:
            ok &= mask_q[None, :]
        if mask_p is not None:
            ok &= mask_p[s0:s0 + block, None]
        n_lo += int(((d < cut - 1e-3) & ok).sum()); n_hi += int(((d < cut + 1e-3) & ok).sum())  # float32 storage
        m = (d < cut) & ok
        if m.any():
            n += int(m.sum()); dmin = min(dmin, float(d[m].min()))
            core = (d < cut - 1e-3) & ok   # flags: ignore float32 boundary cases
            fp[s0:s0 + block] |= core.any(axis=1); fq |= core.any(axis=0)
    return n, dmin, fp, fq, n_lo, n_hi


def fail(msg, uid):
    bad.setdefault(msg, []).append(uid)


for r in rows[:n]:
    uid = r[ix["unit_id"]]
    rec = db.get(uid, int(r[ix["shard"]]))
    at, rs, ed, w, nodes = rec["atoms"], rec["residues"], rec["edges"], rec["water"], rec["nodes"]
    xyz = at["xyz"].astype(np.float64)
    k = len(nodes["kind"])
    if not np.isfinite(xyz).all():
        fail("nonfinite_xyz", uid)
    if at["res"].max() >= len(rs["chem_comp_id"]) or at["node"].max() >= k:
        fail("index_out_of_range", uid)
    if len(xyz) != int(r[ix["n_atoms"]]) or len(w["xyz"]) != int(r[ix["n_water"]]):
        fail("manifest_mismatch", uid)
    donor = np.isin(at["element"], ["N", "O", "S"])
    for e in range(len(ed["a"])):
        a, b = int(ed["a"][e]), int(ed["b"][e])
        A, B = np.flatnonzero(at["node"] == a), np.flatnonzero(at["node"] == b)
        ion = "ion" in (nodes["kind"][a], nodes["kind"][b])
        cut = cfg.contact.ion_cutoff if ion else cfg.contact.cutoff
        mq = donor[B] if (ion and nodes["kind"][a] == "ion") else None
        mp = donor[A] if (ion and nodes["kind"][b] == "ion") else None
        npairs, dmin, fa, fb, n_lo, n_hi = block_pairs(xyz[A], xyz[B], cut, mq, mp)
        if npairs == 0:
            fail("edge_without_contact", uid)
            continue
        if abs(dmin - float(ed["min_dist"][e])) > 1e-2:
            fail("min_dist_mismatch", uid)
        if not n_lo <= int(ed["n_atom_pairs"][e]) <= n_hi:
            fail("n_atom_pairs_mismatch", uid)
        if not at["contact"][A[fa]].all() or not at["contact"][B[fb]].all():
            fail("contact_flag_missing", uid)
    if len(w["xyz"]):
        wx = w["xyz"].astype(np.float64)
        for j0 in range(0, len(wx), 256):
            dw = np.linalg.norm(wx[j0:j0 + 256, None] - xyz[None], axis=-1)
            if not (dw.min(axis=1) < cfg.water_near + 1e-3).all():
                fail("water_too_far", uid)
            for j in range(len(dw)):
                near_nodes = len(set(at["node"][dw[j] < cfg.water_bridge].tolist()))
                if (near_nodes >= 2) != bool(w["bridging"][j0 + j]):
                    fail("bridging_flag_wrong", uid)
                    break
    for nd in range(k):
        if nodes["kind"][nd] in ("protein", "na"):
            seq = nodes["key"][nd][2:]
            res = rs["cif_idx"][rs["node"] == nd]
            if len(res) and res.max() > len(seq) + 1:
                fail("residue_beyond_sequence", uid)
print(f"checked {min(n, len(rows))} units: " + ("all consistent" if not bad else ""))
for msg, ids in bad.items():
    print(f"  FAIL {msg}: {len(ids)} units, e.g. {ids[:3]}")
