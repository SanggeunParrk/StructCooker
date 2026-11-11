"""Reader of the RFD4-Proteina Re-Glyco set (``glycan/reglyco*`` stages, docs/reglyco-compact.md).

A Re-Glyco record is an annotated mmCIF (written by AtomWorks) named ``<xx>/<pdb>/<pdb>_<site>_<glycan>_conf<NNN>.cif``,
kept in a shard (``blob/data/*.bin``, one zstd frame per record) and located by ``blob/index.parquet``. This module

* ``adapt_annotated_cif``   derives the tables the ``pdb/cif`` recipe reads and an AtomWorks CIF lacks, from what the
                            CIF does carry (nothing is invented);
* ``convert_record``        runs the ``pdb/cif`` recipe on one record -> the BioMol value ``{assembly_dict,
                            metadata_dict}`` (schema A shape; the record has the deposition date and resolution the CIF
                            lacks, taken from ``examples.parquet``);
* ``read_worklist``         loader of ``glycan/reglyco_entries``: one entry's worklist file;
* ``explode_reglyco_entry`` reader adapter of ``glycan/reglyco``: one entry's worklist -> the stored value of each of
                            its records (``transforms.delta_codec.encode_entry``).

The raw shards and the CCD are found through environment variables the stage declares (``runtime_environment``).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import zstandard as zstd

from structcooker.instructions.transforms.delta_codec import RecordRef, encode_entry

BLOB_ENV = "REGLYCO_BLOB_DIR"  # the directory holding the shards (blob/data of the downloaded component)
CCD_ENV = "REGLYCO_CCD_DB"
CCD_REF_ENV = "REGLYCO_CCD_REFERENCE"

KEY = re.compile(r"^(?P<pdb>\w{4})_(?P<site>.+)_(?P<glycan>GS\d+)_conf(?P<conf>\d+)$")
MAX_CIF_BYTES = 1 << 28


def record_key(path: str) -> str:
    """Return the key of a record: the Hub file name without ``.cif``."""
    return path.rsplit("/", 1)[-1].removesuffix(".cif")


def split_key(key: str) -> RecordRef:
    """Return the (site, glycan) group and the conformer number of a record key."""
    m = KEY.match(key)
    if not m:
        msg = f"unexpected Re-Glyco record key {key!r}"
        raise ValueError(msg)
    return RecordRef(key, f"{m['site']}_{m['glycan']}", int(m["conf"]))


# -------------------------------------------------------------------------------------------------------- the CIF
def adapt_annotated_cif(d: dict) -> dict:
    """Derive the mmCIF tables ``pdb/cif`` reads and an AtomWorks-written CIF lacks, from what the CIF carries.

    Missing: ``_pdbx_poly_seq_scheme`` / ``_pdbx_nonpoly_scheme`` / ``_pdbx_branch_scheme`` / ``_pdbx_entity_nonpoly``
    / ``_pdbx_entity_branch_list``, ``hetero`` of ``_entity_poly_seq``, ``ptnr*_auth_seq_id`` of ``_struct_conn``
    (the CIF has ``ptnr*_label_seq_id``), and the glycan's components in ``_chem_comp``. The sugar-sugar links of a
    branched chain move from ``_struct_conn`` to ``_pdbx_entity_branch_link``, as in a PDB mmCIF; their leaving atoms
    are not in the CIF and are written ``?``. Everything is read off ``_atom_site`` (label / auth ids, residue names)
    and ``_entity.type``. Returns small diagnostics.
    """
    site = d["_atom_site"]
    ent_type = dict(zip(d["_entity"]["id"], d["_entity"]["type"], strict=True))
    seen: dict[tuple[str, str, str], dict] = {}
    order: list[dict] = []
    n = len(site["label_asym_id"])
    model = site.get("pdbx_PDB_model_num")
    for i in range(n):
        if model is not None and model[i] != model[0]:
            continue
        key = (
            site["label_asym_id"][i],
            site["label_seq_id"][i],
            site["auth_seq_id"][i],
        )
        if key in seen:
            continue
        res = {
            "asym": site["label_asym_id"][i],
            "entity": site["label_entity_id"][i],
            "label_seq": site["label_seq_id"][i],
            "comp": site["label_comp_id"][i],
            "auth_seq": site["auth_seq_id"][i],
            "ins": site["pdbx_PDB_ins_code"][i],
        }
        seen[key] = res
        order.append(res)
    # _chem_comp lists the components of the deposited entry only; the rebuilt glycan's (BMA, MAN, FUC, ...) are missing.
    # Add them by id (name = id, the rest unknown); their atoms and bonds come from the CCD like any component without
    # a _chem_comp_atom table.
    cc = d["_chem_comp"]
    have = set(cc["id"])
    added = sorted({r["comp"] for r in order} - have)
    for comp in added:
        for col, vals in cc.items():
            vals.append(comp if col in {"id", "name"} else "?")
    poly, nonpoly, branch = [], [], []
    for r in order:
        t = ent_type.get(r["entity"], "")
        (poly if t == "polymer" else branch if t == "branched" else nonpoly).append(r)

    def cols(rows: list[dict], names: dict[str, str | None]) -> dict[str, list[str]]:
        return {c: [(r[k] if k else "n") for r in rows] for c, k in names.items()}

    d["_pdbx_poly_seq_scheme"] = cols(
        poly,
        {
            "asym_id": "asym",
            "entity_id": "entity",
            "seq_id": "label_seq",
            "mon_id": "comp",
            "pdb_seq_num": "auth_seq",
            "pdb_ins_code": "ins",
            "hetero": None,
        },
    )
    d["_pdbx_nonpoly_scheme"] = cols(
        nonpoly,
        {
            "asym_id": "asym",
            "entity_id": "entity",
            "mon_id": "comp",
            "pdb_seq_num": "auth_seq",
            "pdb_ins_code": "ins",
        },
    )
    # a branched chain is numbered 1..n in residue order (the entity's branch list position); the CIF only has auth ids
    pos: dict[str, int] = defaultdict(int)
    for r in branch:
        pos[r["asym"]] += 1
        r["num"] = str(pos[r["asym"]])
    d["_pdbx_branch_scheme"] = cols(
        branch,
        {
            "asym_id": "asym",
            "entity_id": "entity",
            "num": "num",
            "mon_id": "comp",
            "pdb_seq_num": "auth_seq",
            "hetero": None,
        },
    )
    bl: dict[
        tuple[str, str, str], None
    ] = {}  # (entity, num, comp) once, over every asym of the entity
    for r in branch:
        bl.setdefault((r["entity"], r["num"], r["comp"]), None)
    d["_pdbx_entity_branch_list"] = {
        "entity_id": [k[0] for k in bl],
        "num": [k[1] for k in bl],
        "comp_id": [k[2] for k in bl],
        "hetero": ["n"] * len(bl),
    }
    eps0 = d.get("_entity_poly_seq")
    if eps0 and "hetero" not in eps0:
        eps0["hetero"] = ["n"] * len(eps0["num"])
    seen_ent: dict[str, str] = {}
    for r in nonpoly:
        seen_ent.setdefault(r["entity"], r["comp"])
    d["_pdbx_entity_nonpoly"] = {
        "entity_id": list(seen_ent),
        "comp_id": list(seen_ent.values()),
    }
    sc = d.get("_struct_conn")
    n_conn = n_intra = 0
    if sc:
        num_of = {(r["asym"], r["label_seq"]): r["num"] for r in branch}
        ent_of_asym = {r["asym"]: r["entity"] for r in branch}
        links: dict[tuple, None] = {}
        keep = []
        for j in range(len(sc["id"])):
            a1, a2 = sc["ptnr1_label_asym_id"][j], sc["ptnr2_label_asym_id"][j]
            if a1 == a2 and a1 in ent_of_asym:
                links.setdefault(
                    (
                        ent_of_asym[a1],
                        sc["ptnr1_label_comp_id"][j],
                        sc["ptnr2_label_comp_id"][j],
                        sc["ptnr1_label_atom_id"][j],
                        sc["ptnr2_label_atom_id"][j],
                        sc["pdbx_value_order"][j],
                        num_of[(a1, sc["ptnr1_label_seq_id"][j])],
                        num_of[(a2, sc["ptnr2_label_seq_id"][j])],
                    ),
                    None,
                )
            else:
                keep.append(j)
        n_intra = len(sc["id"]) - len(keep)
        if links:
            names = [
                "entity_id",
                "comp_id_1",
                "comp_id_2",
                "atom_id_1",
                "atom_id_2",
                "value_order",
                "entity_branch_list_num_1",
                "entity_branch_list_num_2",
            ]
            d["_pdbx_entity_branch_link"] = {
                name: [k[i] for k in links] for i, name in enumerate(names)
            }
            d["_pdbx_entity_branch_link"]["leaving_atom_id_1"] = ["?"] * len(links)
            d["_pdbx_entity_branch_link"]["leaving_atom_id_2"] = ["?"] * len(links)
        for col in sc:
            sc[col] = [sc[col][j] for j in keep]
        auth = {(r["asym"], r["label_seq"]): r["auth_seq"] for r in order}
        for p in ("ptnr1", "ptnr2"):
            if f"{p}_auth_seq_id" not in sc:
                sc[f"{p}_auth_seq_id"] = [
                    auth.get((a, s), s)
                    for a, s in zip(
                        sc[f"{p}_label_asym_id"], sc[f"{p}_label_seq_id"], strict=True
                    )
                ]
            sc.setdefault(f"pdbx_{p}_label_alt_id", ["?"] * len(sc["id"]))
        n_conn = len(sc["id"])
        if not sc["id"]:
            d.pop("_struct_conn")
    return {
        "chem_comp_added": len(added),
        "poly_res": len(poly),
        "nonpoly_res": len(nonpoly),
        "branch_res": len(branch),
        "struct_conn": n_conn,
        "branch_links_moved": n_intra,
    }


def convert_record(
    cif: bytes, key: str, meta: dict[str, Any], tmp_dir: str | None = None
) -> dict[str, Any]:
    """Convert one record's CIF to its BioMol value ``{assembly_dict, metadata_dict}`` (the ``pdb/cif`` recipe).

    ``meta`` carries ``deposition_date`` (``YYYY-MM-DD``) and ``resolution`` from ``examples.parquet``; the CIF has
    neither, and its ``_entry.id`` is a temp file name (the record key replaces it).
    """
    from datacooker import parse_dict

    from structcooker.instructions.readers.cif import dot_transform, get_cif_data
    from structcooker.workflows.ingest.cif import RECIPE, TARGETS

    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    with tempfile.TemporaryDirectory(
        prefix="reglyco_", dir=tmp_dir or os.environ.get("TMPDIR")
    ) as tmp:
        f = Path(tmp) / f"{key}.cif"
        f.write_bytes(cif)
        d = get_cif_data(f)
    d["_entry"]["id"] = [key]
    # AtomWorks writes a subset of the CCD columns (no pdbx_stereo_config). The recipe then takes the component from the
    # CCD itself, as it does for any CIF without _chem_comp_atom; filling the column would invent chemistry.
    d.pop("_chem_comp_atom", None)
    d.pop("_chem_comp_bond", None)
    adapt_annotated_cif(d)
    d["_pdbx_database_status"] = {
        "recvd_initial_deposition_date": [meta["deposition_date"]]
    }
    if meta.get("resolution") is not None:
        d["_refine"] = {"ls_d_res_high": [str(meta["resolution"])]}
    d["ccd_db_path"] = Path(os.environ.get(CCD_ENV, f"{root}/lmdb/ccd/biomol_CCD.lmdb"))
    d["ccd_reference_path"] = Path(
        os.environ.get(CCD_REF_ENV, f"{root}/materials/reference/reference_ccd.lmdb")
    )
    try:
        out = parse_dict(RECIPE, d, key_transform=dot_transform, targets=TARGETS)
    except Exception as e:
        root_error = e
        while (
            root_error.__cause__ is not None
        ):  # the recipe engine wraps the real error
            root_error = root_error.__cause__
        msg = f"{key}: {type(root_error).__name__}: {root_error}"
        raise RuntimeError(msg) from e
    value = {k: out[k] for k in TARGETS}
    if not value["assembly_dict"]:
        msg = f"{key}: the recipe returned an empty assembly_dict"
        raise ValueError(msg)
    return value


# ------------------------------------------------------------------------------------------------- the work list
def read_worklist(path: Path) -> dict[str, Any]:
    """Loader of ``glycan/reglyco_entries``: the JSON worklist of one entry -> ``{"records": {sub key: record}}``."""
    return {"records": json.loads(Path(path).read_text())}


def read_blob(blob_dir: Path, rec: dict[str, Any]) -> bytes:
    """Read the CIF of one record through its index row and check it against the index's sha256."""
    with (Path(blob_dir) / rec["shard"]).open("rb") as f:
        f.seek(int(rec["offset"]))
        frame = f.read(int(rec["length"]))
    if len(frame) != int(rec["length"]):
        msg = f"{rec['key']}: shard {rec['shard']} is shorter than the index says (incomplete download?)"
        raise ValueError(msg)
    cif = zstd.ZstdDecompressor().decompress(frame, max_output_size=MAX_CIF_BYTES)
    if hashlib.sha256(cif).hexdigest() != rec["sha256"]:
        msg = f"{rec['key']}: sha256 of the decompressed CIF does not match the index"
        raise ValueError(msg)
    return cif


def explode_reglyco_entry(entry: dict[str, Any]) -> dict[str, dict[str, bytes]]:
    """Reader adapter of ``glycan/reglyco``: an entry's worklist -> the stored value of each of its records.

    The result is ``{sub key: {"reglyco_record": value}}``; the rebuild writes each under ``<pdb>_<sub key>``, which is
    the record's key. Each record is converted once and encoded as a delta (``delta_codec.encode_entry``; every delta
    is decoded again and compared).
    """
    recs = entry["records"]
    blob_dir = Path(os.environ[BLOB_ENV])
    by_key = {r["key"]: r for r in recs.values()}
    refs = [split_key(k) for k in by_key]

    def convert(key: str) -> dict[str, Any]:
        rec = by_key[key]
        meta = {
            "deposition_date": rec["deposition_date"],
            "resolution": rec.get("resolution"),
        }
        return convert_record(read_blob(blob_dir, rec), key, meta)

    values = encode_entry(refs, convert)
    prefix = next(iter(by_key)).split("_", 1)[0] + "_"
    return {
        key.removeprefix(prefix): {"reglyco_record": value}
        for key, value in values.items()
    }
