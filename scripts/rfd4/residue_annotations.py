#!/usr/bin/env python3
"""Per-residue ``ig_region`` / ``cdr`` of the antibody and TCR sets, read from the shards' CIFs (they are not BioMol features).

One row per residue: ``key, label_asym_id, auth_asym_id, label_seq_id, auth_seq_id, comp, ig_region[, cdr]`` (the value of the
residue's first atom; the CIF writes one value per atom, constant over a residue). Written to
``metadata/rfd4_<set>_residues.parquet``. Heavy: run through sbatch (``residue_annotations.sbatch``).
"""

from __future__ import annotations

import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from structcooker.instructions.readers.rfd4 import parse_cif_text, read_blob
from structcooker.instructions.transforms import rfd4_worklist as wl

COLUMNS = {"antibody_distill": ("ig_region", "cdr"), "tcr_distill": ("ig_region",)}
_BLOB: Path


def residues(rec: dict) -> list[tuple]:
    """Return the residue rows of one record."""
    cols = COLUMNS[rec["set"]]
    site = parse_cif_text(read_blob(_BLOB, rec))["_atom_site"]
    seen: set[tuple] = set()
    out = []
    for i in range(len(site["id"])):
        k = (site["label_asym_id"][i], site["label_seq_id"][i], site["auth_seq_id"][i])
        if k in seen:
            continue
        seen.add(k)
        out.append(
            (
                rec["key"],
                k[0],
                site["auth_asym_id"][i],
                k[1],
                k[2],
                site["label_comp_id"][i],
                *(site[c][i] for c in cols),
            )
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("name", choices=list(COLUMNS))
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    global _BLOB
    root = wl.default_root()
    _BLOB = wl.component_dir(root, a.name) / "blob"
    recs = [r for batch in wl.iter_worklist_rows(a.name, root) for r in batch]
    names = [
        "key",
        "label_asym_id",
        "auth_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "comp",
        *COLUMNS[a.name],
    ]
    rows: list[tuple] = []
    with Pool(a.workers) as pool:
        for part in pool.imap(residues, recs, chunksize=64):
            rows += part
    table = pa.table({n: [r[i] for r in rows] for i, n in enumerate(names)})
    out = (
        Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
        / f"metadata/rfd4_{a.name}_residues.parquet"
    )
    pq.write_table(table, out, compression="zstd")
    print(f"{a.name}: {len(recs)} records, {len(rows)} residues -> {out}")


if __name__ == "__main__":
    main()
