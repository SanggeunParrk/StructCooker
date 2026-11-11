"""Materialize step of ``glycan/reglyco_worklist``: the Hub's index and example tables -> one worklist file per entry.

Reads ``blob/index.parquet`` (which shard, offset and length hold each record, and the sha256 of its CIF) and
``examples.parquet`` (deposition date, resolution and 88 more columns per record), then writes

* ``<out_dir>/<pdb>.json``: ``{sub key: {key, shard, offset, length, sha256, deposition_date, resolution, example_id}}``
  for the records of one entry (``sub key`` is the key without ``<pdb>_``);
* ``out_list``: the worklist files, one per line, relative to ``out_dir`` (what ``glycan/reglyco_entries`` builds from);
* ``out_sidecar``: ``examples.parquet`` with a ``key`` column (the record key; the Hub's ``path`` stays), plus the
  columns of ``examples_glycan.parquet`` that differ from ``examples.parquet`` as ``gv_<name>`` (all of them when new)
  when it is row-aligned with it. The record values hold the
  structure only; everything else the Hub says about a record is here, joined on ``key``.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq


def _key(path: str) -> str:
    return path.rsplit("/", 1)[-1].removesuffix(".cif")


def build_worklists(
    index: pa.Table, examples: pa.Table, pdbs: set[str] | None = None
) -> dict[str, dict[str, Any]]:
    """Return ``{pdb id: {sub key: record}}`` from the Hub's index and example tables (``pdbs``: only these entries)."""
    if (
        "path" not in examples.column_names
        or "deposition_date" not in examples.column_names
    ):
        msg = "examples.parquet lacks path / deposition_date"
        raise ValueError(msg)
    paths = index.column("path").to_pylist()
    row_of = {p: i for i, p in enumerate(examples.column("path").to_pylist())}
    missing = [p for p in paths if p not in row_of]
    if missing:
        msg = f"{len(missing)} indexed records have no examples row, e.g. {missing[:3]}"
        raise ValueError(msg)
    dates = examples.column("deposition_date").to_pylist()
    resolution = examples.column("resolution").to_pylist()
    ids = (
        examples.column("example_id").to_pylist()
        if "example_id" in examples.column_names
        else [None] * len(dates)
    )
    shard, offset, length, sha = (
        index.column(c).to_pylist() for c in ("shard", "offset", "length", "cif_sha256")
    )
    entries: dict[str, dict[str, Any]] = defaultdict(dict)
    for i, p in enumerate(paths):
        key = _key(p)
        pdb = key.split("_", 1)[0]
        if pdbs is not None and pdb not in pdbs:
            continue
        j = row_of[p]
        if dates[j] is None:
            msg = f"{key}: no deposition_date in examples.parquet"
            raise ValueError(msg)
        d = dates[j]
        entries[pdb][key.removeprefix(pdb + "_")] = {
            "key": key,
            "shard": str(shard[i]),
            "offset": int(offset[i]),
            "length": int(length[i]),
            "sha256": sha[i],
            "deposition_date": d.strftime("%Y-%m-%d")
            if hasattr(d, "strftime")
            else str(d)[:10],
            "resolution": resolution[j],
            "example_id": ids[j],
        }
    return dict(entries)


def write_reglyco_worklists(
    index_path: str,
    examples_path: str,
    examples_glycan_path: str,
    out_dir: str,
    out_list: str,
    out_sidecar: str,
) -> str:
    """Write the worklist files, their list and the sidecar table; returns a one-line summary."""
    index = pq.read_table(
        index_path, columns=["path", "shard", "offset", "length", "cif_sha256"]
    )
    examples = pq.read_table(examples_path)
    entries = build_worklists(index, examples)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = []
    for pdb in sorted(entries):
        (out / f"{pdb}.json").write_text(json.dumps(entries[pdb], sort_keys=True))
        names.append(f"{pdb}.json")
    Path(out_list).parent.mkdir(parents=True, exist_ok=True)
    Path(out_list).write_text("\n".join(names) + "\n")
    _write_sidecar(examples, examples_glycan_path, out_sidecar)
    return f"{sum(map(len, entries.values()))} records of {len(entries)} entries -> {out_dir}"


def _write_sidecar(examples: pa.Table, glycan_path: str, sidecar_path: str) -> None:
    keys = pa.array([_key(p) for p in examples.column("path").to_pylist()])
    table = examples.add_column(0, "key", keys)
    gp = Path(glycan_path)
    if gp.exists():
        gly = pq.read_table(gp)
        if gly.column("path").to_pylist() == examples.column("path").to_pylist():
            have = set(table.column_names)
            for name in gly.column_names:
                if name == "path":
                    continue
                if name not in have:
                    table = table.append_column(name, gly.column(name))
                elif not gly.column(name).equals(
                    examples.column(name)
                ):  # a glycan-view value (e.g. size fields)
                    table = table.append_column(f"gv_{name}", gly.column(name))
        # not row-aligned: examples_glycan is left out, said in the summary of the stage
    Path(sidecar_path).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table.sort_by("key"), sidecar_path, compression="zstd")


def keep_records(records: dict[str, Any]) -> dict[str, Any]:
    """Recipe step of ``glycan/reglyco_entries``: the worklist is the stored value."""
    return records


def keep_value(reglyco_record: bytes) -> bytes:
    """Recipe step of ``glycan/reglyco``: the adapter already made the stored value."""
    return reglyco_record
