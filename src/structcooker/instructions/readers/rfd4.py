"""Reader of the RFD4-Proteina structure sets (``rfd4/*`` stages, docs/rfd4-sets.md).

Seven components share one record format: an annotated mmCIF (written by AtomWorks) in a zstd frame of a shard
(``blob/**/*.bin``), located by ``blob/index.parquet`` (shard, offset, length, sha256 of the decompressed CIF). The CIFs
are written by prediction / simulation pipelines, so they are thin: some carry only ``_atom_site`` (+ ``_chem_comp_bond``),
some have no ``_entity`` at all. This module

* ``synthesize_tables``       derives the entity tables the ``pdb/cif`` recipe reads and a record lacks, from its
                              ``_atom_site`` (``_entity``, ``_entity_poly``, ``_entity_poly_seq``, ``_chem_comp``); a table
                              the record has is kept;
* ``prepare_cif_inputs``      one record's CIF text -> the dict the ``pdb/cif`` recipe takes (after
                              ``readers.reglyco.adapt_annotated_cif`` has made the scheme tables);
* ``rfd4_record_to_cif_inputs``  reader adapter of every ``rfd4/<set>`` stage: a worklist row -> that dict;
* ``read_blob``               the shard read (sha256 checked).

The recipe is the PDB's own (``workflows/ingest/cif.py``), so a record is an ordinary schema-A BioMol value
``{assembly_dict, metadata_dict}``. Nothing is invented: a table is derived only from what ``_atom_site`` states, and a
canonical sequence letter that cannot be known is ``X`` (``N`` for nucleic acids), as for any unknown component.
"""

from __future__ import annotations

import hashlib
import io
import os
import re
from pathlib import Path
from typing import Any

import zstandard as zstd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict as mmcif2dict  # noqa: N813

from structcooker.instructions.readers.reglyco import adapt_annotated_cif

BLOB_ENV = "RFD4_BLOB_DIR"  # optional override of the component's blob/ directory (default: from BIOMOL_ROOT and the set)
CCD_ENV = "RFD4_CCD_DB"
CCD_REF_ENV = "RFD4_CCD_REFERENCE"
MAX_CIF_BYTES = 1 << 28

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "UNK": "X",
}  # fmt: skip
DNA = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I", "DN": "N"}
RNA = {"A": "A", "C": "C", "G": "G", "U": "U", "I": "I", "N": "N"}
WATER = {"HOH", "DOD"}
_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


# ------------------------------------------------------------------------------------------------- the CIF text
def parse_cif_text(text: str) -> dict[str, Any]:
    """Return a CIF's tables as ``{category: {item: [values]}}`` (what ``readers.cif.get_cif_data`` returns)."""
    raw = mmcif2dict(io.StringIO(text))
    out: dict[str, Any] = {}
    for key, value in raw.items():
        if "." not in key:
            out[key] = value
            continue
        cat, item = key.split(".", 1)
        out.setdefault(cat, {})[item] = value
    return out


def _residues(d: dict) -> list[dict[str, str]]:
    """Residues of the first model in ``_atom_site`` order, one dict per (asym, label seq, auth seq)."""
    site = d["_atom_site"]
    model = site.get("pdbx_PDB_model_num")
    seen: set[tuple[str, str, str]] = set()
    out: list[dict[str, str]] = []
    for i in range(len(site["label_asym_id"])):
        if model is not None and model[i] != model[0]:
            continue
        key = (
            site["label_asym_id"][i],
            site["label_seq_id"][i],
            site["auth_seq_id"][i],
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "asym": key[0],
                "seq": key[1],
                "comp": site["label_comp_id"][i],
                "entity": site["label_entity_id"][i]
                if "label_entity_id" in site
                else ".",
            }
        )
    return out


def _is_polymer(res: list[dict[str, str]]) -> bool:
    """Tell whether a chain is a polymer: it holds a standard residue or several residues numbered in sequence."""
    comps = {r["comp"] for r in res}
    if comps & (AA3.keys() | DNA.keys()) or (
        len(res) > 1 and comps <= (RNA.keys() | DNA.keys())
    ):
        return True
    return len(res) > 1 and all(r["seq"] not in {".", "?"} for r in res)


def _poly_type(res: list[dict[str, str]]) -> str:
    comps = [r["comp"] for r in res]
    n_dna = sum(c in DNA for c in comps)
    n_rna = sum(c in RNA and c not in AA3 for c in comps)
    n_aa = sum(c in AA3 for c in comps)
    if n_dna > max(n_aa, n_rna) or (n_dna and n_dna * 2 >= len(comps)):
        return "polydeoxyribonucleotide"
    if n_rna > n_aa and n_rna * 2 >= len(comps):
        return "polyribonucleotide"
    return "polypeptide(L)"


def _letter(poly_type: str, comp: str) -> str:
    table = (
        AA3 if poly_type == "polypeptide(L)" else DNA if "deoxy" in poly_type else RNA
    )
    return table.get(comp, "X" if poly_type == "polypeptide(L)" else "N")


def synthesize_tables(d: dict, sequence: str | None = None) -> dict[str, int]:
    """Add the entity tables a record lacks, derived from ``_atom_site``; return small diagnostics.

    ``sequence`` (a canonical one-letter sequence the dataset states for a record with one polymer chain) is used for
    ``pdbx_seq_one_letter_code_can`` when its length equals the chain's residue count; otherwise the letter of each
    residue comes from its component name (``X`` / ``N`` when the component is not a standard one).
    """
    res = _residues(d)
    by_asym: dict[str, list[dict[str, str]]] = {}
    for r in res:
        by_asym.setdefault(r["asym"], []).append(r)
    # entity of each asym: the id the atom_site states, or one new id per distinct (kind, sequence)
    ent_of: dict[str, str] = {}
    for asym, rs in by_asym.items():
        e = rs[0]["entity"]
        if e in {".", "?", ""}:
            e = ""
        ent_of[asym] = e
    if not all(ent_of.values()):
        sig: dict[tuple, str] = {}
        for asym, rs in by_asym.items():
            if not ent_of[asym]:
                key = tuple(r["comp"] for r in rs)
                ent_of[asym] = sig.setdefault(key, f"e{len(sig) + 1}")
    first_asym: dict[str, str] = {}
    for asym, e in ent_of.items():
        first_asym.setdefault(e, asym)

    ent = d.get("_entity")
    made_entity = 0
    if not ent:
        ids = list(first_asym)
        types = []
        for e in ids:
            rs = by_asym[first_asym[e]]
            if {r["comp"] for r in rs} <= WATER:
                types.append("water")
            else:
                types.append("polymer" if _is_polymer(rs) else "non-polymer")
        d["_entity"] = {"id": ids, "type": types}
        made_entity = len(ids)
    ent_type = dict(zip(d["_entity"]["id"], d["_entity"]["type"], strict=True))

    poly = d.get("_entity_poly") or {}
    have_poly = set(poly.get("entity_id", []))
    seq = d.get("_entity_poly_seq") or {}
    have_seq = set(seq.get("entity_id", []))
    n_poly_entities = sum(t == "polymer" for t in ent_type.values())
    made_poly = 0
    for e, t in ent_type.items():
        if t != "polymer" or e not in first_asym:
            continue
        rs = by_asym[first_asym[e]]
        ptype = _poly_type(rs)
        if e not in have_poly:
            letters = [_letter(ptype, r["comp"]) for r in rs]
            if sequence and n_poly_entities == 1 and len(sequence) == len(letters):
                letters = list(sequence)
            code = "".join(letters)
            row = {
                "entity_id": e,
                "type": ptype,
                "nstd_linkage": "no",
                "pdbx_seq_one_letter_code": code,
                "pdbx_seq_one_letter_code_can": code,
            }
            for col in set(poly) | set(row):
                poly.setdefault(col, [])
                poly[col].append(row.get(col, "?"))
            made_poly += 1
        if e not in have_seq:
            for r in rs:
                num = r["seq"] if r["seq"] not in {".", "?"} else str(rs.index(r) + 1)
                for col, val in (
                    ("entity_id", e),
                    ("num", num),
                    ("mon_id", r["comp"]),
                    ("hetero", "n"),
                ):
                    seq.setdefault(col, []).append(val)
    if made_poly:
        d["_entity_poly"] = poly
    if seq:
        d["_entity_poly_seq"] = seq
    if "_chem_comp" not in d:
        comps = list(dict.fromkeys(r["comp"] for r in res))
        d["_chem_comp"] = {
            "id": comps,
            "name": comps,
            "type": ["?"] * len(comps),
            "formula": ["?"] * len(comps),
        }
    return {"entities_made": made_entity, "polymers_made": made_poly}


def _ensure_poly_codes(d: dict) -> None:
    """Derive the sequence-code columns a record's own ``_entity_poly`` lacks.

    The recipe reads them; they come from the record's
    ``_entity_poly_seq`` (the residue names of the entity, in order).
    """
    poly, seq = d.get("_entity_poly"), d.get("_entity_poly_seq")
    if not poly or not seq:
        return
    by_ent: dict[str, list[str]] = {}
    for e, mon in zip(seq["entity_id"], seq["mon_id"], strict=True):
        by_ent.setdefault(e, []).append(mon)
    codes = [
        "".join(_letter(t, m) for m in by_ent.get(e, []))
        for e, t in zip(poly["entity_id"], poly["type"], strict=True)
    ]
    for col in ("pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"):
        poly.setdefault(col, list(codes))
    poly.setdefault("nstd_linkage", ["no"] * len(codes))


_CCD_ENVS: dict[int, Any] = {}


def _ccd_has(comp: str) -> bool:
    """Whether the release's CCD defines component ``comp`` (the environment is opened once per process)."""
    import lmdb

    pid = os.getpid()
    env = _CCD_ENVS.get(pid)
    if env is None:
        root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
        env = lmdb.open(
            os.environ.get(CCD_ENV, f"{root}/lmdb/ccd/biomol_CCD.lmdb"),
            readonly=True,
            lock=False,
            readahead=False,
        )
        _CCD_ENVS[pid] = env
    with env.begin() as txn:
        return txn.get(comp.encode()) is not None


def _keep_custom_components(d: dict) -> None:
    """Keep ``_chem_comp_atom`` / ``_chem_comp_bond`` rows only of components the CCD does not define.

    A record that carries these tables lists every component it uses (even ``ALA``) without the CCD's stereo flags;
    the recipe reads such a table over the CCD, so a CCD component is dropped from it and comes from the CCD, whole.
    What stays are the record's own components (an OMol molecule ``A:1``), which nothing else defines; a stereo flag
    they lack is ``N`` (no stereo centre stated), as in a CCD entry of an achiral component.
    """
    atom = d.get("_chem_comp_atom")
    if not atom:
        d.pop("_chem_comp_bond", None)
        return
    custom = {c for c in dict.fromkeys(atom["comp_id"]) if not _ccd_has(c)}
    for tab in ("_chem_comp_atom", "_chem_comp_bond"):
        t = d.get(tab)
        if not t:
            continue
        keep = [i for i, c in enumerate(t["comp_id"]) if c in custom]
        for col in t:
            t[col] = [t[col][i] for i in keep]
        t.setdefault("pdbx_stereo_config", ["N"] * len(keep))
    if not custom:
        d.pop("_chem_comp_atom")
        d.pop("_chem_comp_bond", None)


def _require_definitions(d: dict) -> None:
    """Every component of the record must be defined by the CCD or by the record's own ``_chem_comp_atom``.

    A component with neither (an OMol molecule written without its atom table) has atoms but no stated bonds or
    elements beyond the coordinates' own; the recipe would fail on it obscurely, and inventing its chemistry is not
    an option, so the record fails here by name.
    """
    own = set((d.get("_chem_comp_atom") or {}).get("comp_id", []))
    used = dict.fromkeys(d["_atom_site"]["label_comp_id"])
    undefined = sorted(
        c for c in used if c not in own and c != "UNL" and not _ccd_has(c)
    )
    if undefined:
        msg = f"component(s) {undefined} are in neither the CCD nor the record's _chem_comp_atom"
        raise ValueError(msg)


def _fill_atom_site(site: dict) -> None:
    """Columns the recipe reads that a thin record may lack: a prediction has no occupancy or B factor of its own."""
    n = len(site["id"])
    site.setdefault("occupancy", ["1.0"] * n)
    site.setdefault("B_iso_or_equiv", ["0.0"] * n)
    site.setdefault("label_alt_id", ["."] * n)
    site.setdefault("pdbx_PDB_ins_code", ["?"] * n)
    site.setdefault("pdbx_PDB_model_num", ["1"] * n)
    if "auth_asym_id" not in site:
        site["auth_asym_id"] = list(site["label_asym_id"])


def entry_date(d: dict, default: str) -> str:
    """Return the record's own ``_entry.date`` when it is a date, else the component's release date."""
    got = (d.get("_entry") or {}).get("date")
    got = got[0] if isinstance(got, list) and got else got
    return got if isinstance(got, str) and _DATE.match(got) else default


def prepare_cif_inputs(
    text: str, key: str, default_date: str, sequence: str | None = None
) -> dict[str, Any]:
    """One record's CIF text -> the dict the ``pdb/cif`` recipe takes (``readers.cif.dot_transform`` keys)."""
    d = parse_cif_text(text)
    _fill_atom_site(d["_atom_site"])
    date = entry_date(d, default_date)
    d["_entry"] = {"id": [key]}
    # A thin record's _chem_comp_bond covers only some components (and has no atoms): the recipe then takes each
    # component from the CCD, as for any CIF without _chem_comp_atom. A record with _chem_comp_atom defines its own
    # components (OMol molecules): those stay, the CCD's come from the CCD.
    _keep_custom_components(d)
    _require_definitions(d)
    synthesize_tables(d, sequence)
    _ensure_poly_codes(d)
    adapt_annotated_cif(d)
    d["_pdbx_database_status"] = {"recvd_initial_deposition_date": [date]}
    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    d["ccd_db_path"] = Path(os.environ.get(CCD_ENV, f"{root}/lmdb/ccd/biomol_CCD.lmdb"))
    d["ccd_reference_path"] = Path(
        os.environ.get(CCD_REF_ENV, f"{root}/materials/reference/reference_ccd.lmdb")
    )
    return d


# ------------------------------------------------------------------------------------------------- the shard read
def read_blob(blob_dir: Path, rec: dict[str, Any]) -> str:
    """Read one record's CIF through its worklist row and check it against the index's sha256."""
    with (Path(blob_dir) / rec["shard"]).open("rb") as f:
        f.seek(int(rec["offset"]))
        frame = f.read(int(rec["length"]))
    if len(frame) != int(rec["length"]):
        msg = f"{rec['key']}: shard {rec['shard']} is shorter than the index says (incomplete download?)"
        raise ValueError(msg)
    cif = zstd.ZstdDecompressor().decompress(frame, max_output_size=MAX_CIF_BYTES)
    if rec.get("sha256") and hashlib.sha256(cif).hexdigest() != rec["sha256"]:
        msg = f"{rec['key']}: sha256 of the decompressed CIF does not match the index"
        raise ValueError(msg)
    return cif.decode()


def blob_dir(rec: dict[str, Any]) -> Path:
    """Return the ``blob/`` directory of the record's component: ``$BIOMOL_ROOT/materials/raw/rfd4_proteina/...`` by its set.

    Directories cannot be declared as a stage's pinned reader inputs, so the path is derived from the root (a sampled
    end-to-end run sets ``BIOMOL_ROOT`` to its sample) and the row's ``set``; ``RFD4_BLOB_DIR`` overrides it.
    """
    if os.environ.get(BLOB_ENV):
        return Path(os.environ[BLOB_ENV])
    from structcooker.instructions.transforms.rfd4_worklist import COMPONENTS, SYNTHETIC

    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    return (
        Path(root)
        / "materials/raw/rfd4_proteina"
        / SYNTHETIC
        / COMPONENTS[rec["set"]].subdir
        / "blob"
    )


def rfd4_record_to_cif_inputs(rec: dict[str, Any]) -> dict[str, Any]:
    """Reader adapter of every ``rfd4/<set>`` stage: a worklist row -> the ``pdb/cif`` recipe's inputs."""
    text = read_blob(blob_dir(rec), rec)
    return prepare_cif_inputs(text, rec["key"], rec["date"], rec.get("sequence"))
