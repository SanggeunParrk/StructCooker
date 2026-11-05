"""A duck-typed stand-in for ``CIFMol`` built from plain Python specs (importable from every test folder)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from biomol.exceptions import FeatureKeyError


def _f(value):
    return SimpleNamespace(value=np.asarray(value))


class FakeIndexTable:
    """The subset of ``IndexTable`` the extractors use, over contiguous chain/residue/atom blocks."""

    def __init__(self, atom_res, res_chain):
        self.atom_res = np.asarray(atom_res, dtype=np.int64)
        self.res_chain = np.asarray(res_chain, dtype=np.int64)

    def atoms_to_residues(self, a):
        """Mirror of the IndexTable method of the same name."""
        return self.atom_res[np.asarray(a)]

    def residues_to_chains(self, r):
        """Mirror of the IndexTable method of the same name."""
        return self.res_chain[np.asarray(r)]

    def atoms_to_chains(self, a):
        """Mirror of the IndexTable method of the same name."""
        return self.res_chain[self.atom_res[np.asarray(a)]]

    def residues_to_atoms(self, r):
        """Mirror of the IndexTable method of the same name."""
        return np.flatnonzero(np.isin(self.atom_res, np.asarray(r)))

    def chains_to_residues(self, c):
        """Mirror of the IndexTable method of the same name."""
        return np.flatnonzero(np.isin(self.res_chain, np.asarray(c)))

    def chains_to_atoms(self, c):
        """Mirror of the IndexTable method of the same name."""
        return np.flatnonzero(np.isin(self.res_chain[self.atom_res], np.asarray(c)))


class _Atoms:
    def __init__(self, names, elems, xyz, bf, occ, conn):
        self.id, self.element = _f(names), _f(elems)
        self.xyz, self.b_factor, self.occupancy = _f(xyz), _f(bf), _f(occ)
        self._conn = conn

    @property
    def struct_conn(self):
        if self._conn is None:
            msg = "struct_conn"
            raise FeatureKeyError(msg)
        return self._conn


def residue(comp, one, atoms, auth=None):
    """``atoms`` is a list of ``(name, element, (x, y, z))``."""
    return {"comp": comp, "one": one, "auth": auth, "atoms": atoms}


def backbone_residue(i, one="A", comp="ALA", z=0.0, sign=1.0, cb=True, x0=0.0):
    """Build an extended-strand residue i (3.8 A rise along x); the side chain points to ``-sign * z``."""
    x = x0 + 3.8 * i
    s = sign
    at = [
        ("N", "N", (x - 1.3, 0.5, z)),
        ("CA", "C", (x, 0.0, z)),
        ("C", "C", (x + 1.3, 0.5, z)),
        ("O", "O", (x + 1.3, 1.7, z)),
    ]
    if cb:
        at.append(("CB", "C", (x, -0.5, z + s * 1.5)))
    return residue(comp, one, at, auth=str(i + 1))


def chain(chain_id, residues, etype="polypeptide(L)", cluster=None, entity="1"):
    return {
        "chain_id": chain_id,
        "etype": etype,
        "residues": residues,
        "cluster": cluster or f"c_{chain_id}",
        "entity": entity,
    }


def protein_chain(chain_id, n, z=0.0, sign=1.0, seq=None, cluster=None, x0=0.0):
    seq = seq or "A" * n
    return chain(
        chain_id,
        [
            backbone_residue(
                i,
                seq[i],
                {"A": "ALA", "G": "GLY"}.get(seq[i], "ALA"),
                z=z,
                sign=sign,
                cb=seq[i] != "G",
                x0=x0,
            )
            for i in range(n)
        ],
        cluster=cluster,
    )


def make_mol(chains, struct_conn=None):
    """Build the fake ``CIFMol``; ``struct_conn`` = list of ``(src_atom, dst_atom, type)`` or None (feature absent)."""
    names, elems, xyz, atom_res = [], [], [], []
    comp, one, cif, auth, res_chain = [], [], [], [], []
    for c, ch in enumerate(chains):
        for k, r in enumerate(ch["residues"]):
            rid = len(comp)
            comp.append(r["comp"])
            one.append(r["one"])
            cif.append(k)
            auth.append(r["auth"] or str(k + 1))
            res_chain.append(c)
            for n, e, p in r["atoms"]:
                names.append(n)
                elems.append(e)
                xyz.append(p)
                atom_res.append(rid)
    na = len(names)
    conn = None
    if struct_conn is not None:
        conn = SimpleNamespace(
            src_indices=np.array([s for s, _, _ in struct_conn], dtype=np.int64),
            dst_indices=np.array([d for _, d, _ in struct_conn], dtype=np.int64),
            value=np.array([[t] for _, _, t in struct_conn]),
        )
    return SimpleNamespace(
        atoms=_Atoms(
            np.array(names, dtype="U6"),
            np.array(elems, dtype="U2"),
            np.array(xyz, dtype=float).reshape(-1, 3),
            np.linspace(10, 20, na).astype(np.float32),
            np.ones(na, dtype=np.float32),
            conn,
        ),
        residues=SimpleNamespace(
            chem_comp_id=_f(np.array(comp, dtype="U5")),
            one_letter_code_can=_f(np.array(one, dtype="U1")),
            cif_idx=_f(np.array(cif, dtype=np.int32)),
            auth_idx=_f(np.array(auth, dtype="U8")),
        ),
        chains=SimpleNamespace(
            entity_type=_f(np.array([c["etype"] for c in chains], dtype="U64")),
            chain_id=_f(np.array([c["chain_id"] for c in chains], dtype="U6")),
            auth_asym_id=_f(
                np.array([c["chain_id"] + "a" for c in chains], dtype="U6")
            ),
            entity_id=_f(np.array([c["entity"] for c in chains], dtype="U6")),
            cluster_id=_f(np.array([c["cluster"] for c in chains], dtype="U20")),
        ),
        index_table=FakeIndexTable(atom_res, res_chain),
    )
