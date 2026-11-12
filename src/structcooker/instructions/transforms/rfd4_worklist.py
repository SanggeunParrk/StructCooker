"""Materialize step of every ``rfd4/<set>_worklist`` stage: a component's Hub tables -> a worklist LMDB and a sidecar.

Reads ``blob/index.parquet`` (which shard, offset and length hold each record, and the sha256 of its CIF) and
``examples.parquet`` (everything else the Hub says about a record) of one RFD4-Proteina component, then writes

* ``<worklist>`` (LMDB): key = the record key, value = ``to_bytes`` of ``{key, set, shard, offset, length, sha256, date[, sequence]}``;
  ``shard`` is the path of the shard relative to the component's ``blob/`` directory (the index's own for most
  components; the shards of the 2026-10 reannotation sit next to ``blob/data``);
* ``<sidecar>`` (parquet): ``examples.parquet`` with a ``key`` column, joined on the same column as the worklist (a record
  may have several rows: the TCR set has one per receptor/pMHC unit pair). The
  BioMol value holds the structure only, so every per-record score, cluster and flag stays here.

Record keys (stable, unique per component): ``example_id`` for most; the file stem of ``path`` for the TCR set; the
zero-padded ``lmdb_idx`` for OMol, whose ``example_id`` is a path-like string. ``COMPONENTS`` is the one place that
says which is which.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import lmdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from structcooker.instructions.transforms.codecs import to_bytes

SYNTHETIC = "train_datasets/streaming-latest/synthetic"


@dataclass(frozen=True)
class Component:
    """One RFD4-Proteina component: where it lives, how its records are keyed, what date its records carry."""

    subdir: str  # under synthetic/
    key_col: (
        str  # column of index.parquet and examples.parquet the record key is made from
    )
    key_kind: str  # "id" | "stem" | "idx"
    default_date: (
        str  # YYYY-MM-DD of the component's release (a record's own _entry.date wins)
    )
    sha_cols: tuple[str, ...] = ("sha256", "cif_sha256", "cleaned_cif_sha256")
    sequence_col: str | None = (
        None  # canonical one-letter sequence of a record with one polymer chain
    )


COMPONENTS: dict[str, Component] = {
    "omol": Component("omol_annotated", "lmdb_idx", "idx", "2026-10-02"),
    "rna_monomer": Component("2026_09_21_rna_monomers", "example_id", "id", "2026-09-21"),
    "na_complex": Component("2026_09_22_protein_dna_esmfold2", "example_id", "id", "2026-09-22"),
    "antibody_distill": Component("2026_09_18_antibody_antigen_esmfold2", "example_id", "id", "2026-09-18"),
    "tcr_distill": Component("2026_09_09_tcr_af2", "path", "stem", "2026-09-09"),
    "synthetic_ncaa": Component("ncaa/2026_09_21_ptm_esmfold2", "example_id", "id", "2026-09-21", sequence_col="sequence"),
    "ncaa_scaffolds": Component("ncaa/2026_09_26_ncaa_esmfold2_scaffolds", "example_id", "id", "2026-09-26", sequence_col="sequence"),
}  # fmt: skip


def component_dir(root: str | Path, name: str) -> Path:
    """Directory of component ``name`` under the downloaded dataset root (``materials/raw/rfd4_proteina``)."""
    return Path(root) / SYNTHETIC / COMPONENTS[name].subdir


def make_keys(column: pa.ChunkedArray | pa.Array, kind: str) -> list[str]:
    """Record keys from the key column of an index or examples table."""
    if kind == "id":
        return [str(v) for v in column.to_pylist()]
    if kind == "stem":
        return [
            str(v).rsplit("/", 1)[-1].removesuffix(".cif") for v in column.to_pylist()
        ]
    if kind == "idx":
        return [f"{int(v):09d}" for v in column.to_pylist()]
    msg = f"unknown key kind {kind!r}"
    raise ValueError(msg)


def resolve_shards(blob: Path, shards: list[str]) -> dict[str, str]:
    """Map each distinct index ``shard`` to its path relative to ``blob`` (it must exist)."""
    out = {}
    for s in dict.fromkeys(shards):
        for cand in (s, f"data/{s}", f"data/{s.rsplit('/', 1)[-1]}"):
            if (blob / cand).is_file():
                out[s] = cand
                break
        else:
            msg = f"shard {s!r} not found under {blob} (incomplete download?)"
            raise FileNotFoundError(msg)
    return out


def iter_worklist_rows(
    name: str, root: str | Path, batch: int = 500_000
) -> Iterator[list[dict]]:
    """Yield the worklist rows of component ``name`` in index order, ``batch`` at a time (OMol has 53.7M records)."""
    comp = COMPONENTS[name]
    base = component_dir(root, name)
    f = pq.ParquetFile(base / "blob/index.parquet")
    sha_col = next((c for c in comp.sha_cols if c in f.schema_arrow.names), None)
    seq_of: dict[str, str] = {}
    if comp.sequence_col:
        ex = pq.read_table(
            base / "examples.parquet", columns=[comp.key_col, comp.sequence_col]
        )
        seq_of = dict(
            zip(
                make_keys(ex.column(comp.key_col), comp.key_kind),
                ex.column(comp.sequence_col).to_pylist(),
                strict=True,
            )
        )
    rel: dict[str, str] = {}
    cols = [comp.key_col, "shard", "offset", "length", *([sha_col] if sha_col else [])]
    for b in f.iter_batches(batch_size=batch, columns=cols):
        keys = make_keys(b.column(comp.key_col), comp.key_kind)
        shards = b.column("shard").to_pylist()
        rel.update(
            resolve_shards(
                base / "blob", [s for s in dict.fromkeys(shards) if s not in rel]
            )
        )
        offs, lens = b.column("offset").to_pylist(), b.column("length").to_pylist()
        shas = b.column(sha_col).to_pylist() if sha_col else [None] * len(keys)
        rows = []
        for k, s, o, n, h in zip(keys, shards, offs, lens, shas, strict=True):
            row = {
                "key": k,
                "set": name,
                "shard": rel[s],
                "offset": int(o),
                "length": int(n),
                "sha256": h,
                "date": comp.default_date,
            }
            if seq_of.get(k):
                row["sequence"] = seq_of[k]
            rows.append(row)
        yield rows


def write_worklist(
    name: str, root: str, out_lmdb: str, out_sidecar: str, out_list: str
) -> str:
    """Write the worklist LMDB, the sidecar and the one-line list file; return a one-line summary."""
    comp = COMPONENTS[name]
    base = component_dir(root, name)
    path = Path(out_lmdb)
    path.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(path), map_size=1 << 40)
    n_rows = 0
    try:
        for rows in iter_worklist_rows(name, root):
            with env.begin(write=True) as txn:
                for r in rows:
                    txn.put(r["key"].encode(), to_bytes(r))
            n_rows += len(rows)
        n_keys = env.stat()["entries"]
    finally:
        env.close()
    if n_keys != n_rows:
        msg = f"{name}: {n_rows} index rows but {n_keys} distinct record keys (the key column is not unique)"
        raise ValueError(msg)
    examples = pq.read_table(base / "examples.parquet")
    ex_keys = make_keys(examples.column(comp.key_col), comp.key_kind)
    table = examples.add_column(0, "key", pa.array(ex_keys))
    in_blob = pc.is_in(table.column("key"), value_set=pa.array(_index_keys(name, root)))
    n_no_blob = len(table) - pc.sum(in_blob.cast(pa.int64())).as_py()
    Path(out_sidecar).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table.sort_by("key"), out_sidecar, compression="zstd")
    Path(out_list).parent.mkdir(parents=True, exist_ok=True)
    Path(out_list).write_text(f"{out_lmdb}\n")
    return f"{name}: {n_rows} records; {n_no_blob} examples rows have no blob record -> {out_lmdb}"


def _index_keys(name: str, root: str | Path) -> pa.Array:
    comp = COMPONENTS[name]
    t = pq.read_table(
        component_dir(root, name) / "blob/index.parquet", columns=[comp.key_col]
    )
    return pa.array(make_keys(t.column(comp.key_col), comp.key_kind))


def default_root() -> str:
    """Return the downloaded dataset root (``$BIOMOL_ROOT/materials/raw/rfd4_proteina``)."""
    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    return f"{root}/materials/raw/rfd4_proteina"
