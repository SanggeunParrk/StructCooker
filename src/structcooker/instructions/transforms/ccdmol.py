"""The CCD as inference reads it: one BioMol graph per component (atoms / residues / chains).

The release CCD (``ccd/biomol_CCD.lmdb``, schema G) stores ``{"chem_comp_dict": {atom, residue}}``.
Inference (FoldBench and the model's data loader) reads the same component as a CCDMol graph:
``{atoms, residues, chains, index_table, metadata}``. This is a pure reshape of the CCD's own primary
fields -- id, element, aromatic, stereo, charge, model_xyz and the three bond features -- with no
derived chemistry (no RDKit), so single-atom and bare-ion components pass through with an empty bond
section. Built by scripts/staging/build_inference_ccd.py.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def _node(value: Any) -> dict:
    return {"value": value}


def _edge(value: Any, src: np.ndarray, dst: np.ndarray) -> dict:
    return {"value": value, "src_indices": src, "dst_indices": dst}


def convert_chem_comp_to_ccdmol(record: dict, ccd_code: str) -> dict:
    """``{"chem_comp_dict": ...}`` (release CCD record) -> CCDMol graph dict for ``ccd_code``."""
    cc = record["chem_comp_dict"]
    if cc is None:
        msg = f"{ccd_code}: the CCD defines no atoms"
        raise ValueError(msg)
    an, ae = cc["atom"]["nodes"], cc["atom"]["edges"]
    ids = np.asarray(an["id"]["value"])
    n = len(ids)
    src = np.asarray(ae["bond_type"]["src_indices"], dtype=np.int64)
    dst = np.asarray(ae["bond_type"]["dst_indices"], dtype=np.int64)
    atoms = {
        "nodes": {key: _node(np.asarray(an[key]["value"]))
                  for key in ("id", "element", "aromatic", "stereo", "charge", "model_xyz")},
        "edges": {key: _edge(np.asarray(ae[key]["value"]), src, dst)
                  for key in ("bond_type", "bond_aromatic", "bond_stereo")},
    }
    rn = cc["residue"]["nodes"]
    residues = {"nodes": {"id": _node(np.asarray(rn["id"]["value"])),
                          "formula": _node(np.asarray(rn["formula"]["value"]))}, "edges": {}}
    chains = {"nodes": {"id": _node(np.array([ccd_code]))}, "edges": {}}
    index_table = {
        "atom_to_res": np.zeros(n, dtype=np.int64),
        "res_to_chain": np.zeros(1, dtype=np.int64),
        "res_atom_indptr": np.array([0, n], dtype=np.int64),
        "res_atom_indices": np.arange(n, dtype=np.int64),
        "chain_res_indptr": np.array([0, 1], dtype=np.int64),
        "chain_res_indices": np.zeros(1, dtype=np.int64),
    }
    return {"atoms": atoms, "residues": residues, "chains": chains, "index_table": index_table, "metadata": {}}
