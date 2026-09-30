"""Read one Teddymer dimer PDB into the flat atom table the structure recipe consumes.

Teddymer ships 10M ``*.pdb`` files, each a pair of TED domains of ONE AFDB v4 model
split into chains A/B (they are domains of the same protein, not two proteins). The
files carry ``ATOM``/``TER``/``END`` only -- no header, no SEQRES, no HETATM, no
CRYST1 -- so the entry id comes from the filename and the sequence from the ATOM
records. The B-factor column holds pLDDT.

The table is the same shape ``openfold_structure`` builds CIFMol from, with one
difference: a PDB has no connectivity, so ``bonds`` is absent here and
``structcooker.instructions.transforms.teddymer.attach_bonds`` derives it from the
CCD (which the recipe already loads) before the openfold steps run. That is why this
returns ``raw_atom_site_dict`` rather than ``atom_site_dict``.
"""

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np

# Teddymer is protein-only; 0 is polypeptide(L) in the openfold entity-type map.
_MOL_TYPE_PROTEIN = 0


def teddymer_entry_key(pdb_path: Path) -> str:
    """Entry key = the filename stem (``<DimerIndex>DI_<accession>_<TED pair>``).

    The DimerIndex prefix is the join key against ``cluster.tsv`` and
    ``nonsingletonrep_metadata.tsv``; neither of those carries the TED suffix.
    """
    return Path(pdb_path).stem


def get_teddymer_structure_data(pdb_path: Path) -> dict[str, Any]:
    """Parse one dimer PDB into per-atom columns (no connectivity)."""
    path = Path(pdb_path)
    with path.open() as handle:
        table = pdb_atom_table(handle, name=path.name)
    return {"raw_atom_site_dict": table, "entry_id": teddymer_entry_key(path)}


def pdb_atom_table(lines: Iterable[str], *, name: str, drop_hydrogens: bool = False) -> dict[str, np.ndarray]:
    """Parse the ATOM/HETATM records of a predicted-model PDB into the flat atom table.

    Shared by the readers of PDB-format predicted models (teddymer, AFDB homodimers, the
    OpenFold distillation relaxed models). ``drop_hydrogens`` is for relaxed models, which
    carry explicit hydrogens the heavy-atom CIFMol does not keep. The B-factor column
    (pLDDT for these models) is kept as ``b_factor``.
    """
    atom_name: list[str] = []
    res_name: list[str] = []
    chain_id: list[str] = []
    res_id: list[int] = []
    ins_code: list[str] = []
    element: list[str] = []
    hetero: list[bool] = []
    occupancy: list[float] = []
    b_factor: list[float] = []
    coord: list[tuple[float, float, float]] = []

    for line in lines:
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        # Keep one conformer. Predicted models have no altlocs, but a stray one must not
        # silently double the residue's atoms.
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        # Predicted models always write the element column; fall back to the atom
        # name's leading letter rather than emitting an empty element.
        el = line[76:78].strip() or line[12:16].strip()[:1]
        if drop_hydrogens and el in ("H", "D"):
            continue
        atom_name.append(line[12:16].strip())
        res_name.append(line[17:20].strip())
        chain_id.append(line[21])
        res_id.append(int(line[22:26]))
        ins_code.append(line[26].strip())
        coord.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        occupancy.append(float(line[54:60] or 1.0))
        b_factor.append(float(line[60:66] or 0.0))
        element.append(el)
        hetero.append(record == "HETATM")

    if not atom_name:
        msg = f"{name} has no ATOM records"
        raise ValueError(msg)

    chain_arr = np.array(chain_id, dtype="<U4")
    res_name_arr = np.array(res_name, dtype="<U5")
    return {
        "atom_name": np.array(atom_name, dtype="<U4"),
        "res_name": res_name_arr,
        "chain_id": chain_arr,
        "res_id": np.array(res_id, dtype=np.int64),
        "ins_code": np.array(ins_code, dtype="<U1"),
        "element": np.array(element, dtype="<U2"),
        "hetero": np.array(hetero, dtype=bool),
        "occupancy": np.array(occupancy, dtype=np.float64),
        "b_factor": np.array(b_factor, dtype=np.float64),
        "coord": np.array(coord, dtype=np.float64),
        "molecule_type_id": np.full(len(atom_name), _MOL_TYPE_PROTEIN, dtype=np.int64),
        "entity_id": _entity_ids(chain_arr, res_name_arr, np.array(res_id, dtype=np.int64)),
    }


def _entity_ids(
    chain_id: np.ndarray,
    res_name: np.ndarray,
    res_id: np.ndarray,
) -> np.ndarray:
    """Group chains by sequence: identical sequences share one entity id.

    The two TED domains normally differ, so a dimer has two entities -- but a
    domain repeat would make them one, and ``build_hierarchy`` splits chains on
    (chain_id, entity_id), so this must reflect the sequence, not the chain label.
    """
    seq_of: dict[str, str] = {}
    for c in dict.fromkeys(chain_id.tolist()):
        mask = chain_id == c
        # One residue per (res_id) run; res_name of each distinct residue in order.
        rid = res_id[mask]
        rnm = res_name[mask]
        keep = np.empty(len(rid), dtype=bool)
        keep[0] = True
        keep[1:] = rid[1:] != rid[:-1]
        seq_of[c] = "-".join(rnm[keep].tolist())

    entity_of: dict[str, str] = {}
    by_seq: dict[str, str] = {}
    for c, seq in seq_of.items():
        entity_of[c] = by_seq.setdefault(seq, str(len(by_seq) + 1))
    return np.array([entity_of[c] for c in chain_id.tolist()], dtype="<U8")
