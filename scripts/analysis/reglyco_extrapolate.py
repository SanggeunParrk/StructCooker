#!/usr/bin/env python3
"""Extrapolate the measured compact Re-Glyco layouts (reglyco_compact.py) to the whole release and the training view.

Population numbers (entries, (site, glycan) groups, records, shipped bytes) are exact: they come from the Hub's
``blob/index.parquet`` and ``examples.parquet``. Sizes are a model on top of a few-entry sample:

  base of an entry     = r_base * shipped length of the entry's first record   (r_base measured on the sample)
  A1 delta per record  = measured mean of the sample                         (records - entries of them)
  A2 group base / conformer delta = measured means of the sample             (groups - entries, records - groups)

The error bars are the smallest and largest per-entry means in the sample. With five entries they are indicative only.
Heavy (reads 3M-row tables): run through sbatch (``reglyco_extrapolate.sbatch``).
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

PAT = (
    r"^\w\w/(?P<pdb>\w{4})/\w{4}_(?P<site>.+)_(?P<glycan>GS\d+)_conf(?P<conf>\d+)\.cif$"
)


def population(index: pa.Table, mask: pa.Array | None) -> dict:
    t = index if mask is None else index.filter(mask)
    parts = pc.extract_regex(t["path"].combine_chunks(), PAT)
    t = (
        t.append_column("pdb", parts.field("pdb"))
        .append_column("site", parts.field("site"))
        .append_column("glycan", parts.field("glycan"))
    )
    t = t.append_column(
        "grp", pc.binary_join_element_wise(t["pdb"], t["site"], t["glycan"], "|")
    )
    entries = t.group_by("pdb").aggregate(
        [("path", "min"), ("shard", "count_distinct"), ("path", "count")]
    )
    # the first record (by path) of every entry: its shipped length models the base
    first = t.join(
        entries.select(["pdb", "path_min"]).rename_columns(["pdb", "path"]),
        keys=["pdb", "path"],
        join_type="inner",
    )
    return {
        "records": t.num_rows,
        "entries": entries.num_rows,
        "groups": pc.count_distinct(t["grp"]).as_py(),
        "shipped_bytes": pc.sum(t["length"]).as_py(),
        "raw_bytes": pc.sum(t["uncompressed_length"]).as_py(),
        "first_shipped_bytes": pc.sum(first["length"]).as_py(),
        "entries_in_one_shard": int(
            pc.sum(pc.equal(entries["shard_count_distinct"], 1)).as_py()
        ),
        "records_per_entry_mean": t.num_rows / entries.num_rows,
        "records_per_entry_max": pc.max(entries["path_count"]).as_py(),
    }


def sample_stats(js: dict, layout: str, level: int) -> dict | None:
    for r in js["layouts"]:
        if r["layout"] == layout and r.get("level") == level:
            return r
    return None


def model(pop: dict, a1: dict, a2: dict) -> dict:
    out = {}
    pe = a1["per_entry"]
    r_base = [v["base"] / v["base_shipped"] for v in pe.values()]
    rb = (
        sum(v["base"] for v in pe.values())
        / sum(v["base_shipped"] for v in pe.values()),
        min(r_base),
        max(r_base),
    )
    d1 = [v["delta_mean"] for v in pe.values() if v["records"] > 1]
    d1_mean = a1["delta_bytes"] / a1["n_delta"]
    pe2 = a2["per_entry"]
    gb = [v["group_base"] / (v["groups"] - 1) for v in pe2.values() if v["groups"] > 1]
    gb_mean = a2["group_base_bytes"] / max(
        1, sum(v["groups"] - 1 for v in pe2.values())
    )
    cb = [v["conf"] / v["n_conf"] for v in pe2.values() if v["n_conf"]]
    cb_mean = a2["conf_bytes"] / a2["n_conf"]
    E, G, R = pop["entries"], pop["groups"], pop["records"]

    def base_b(r: float) -> float:
        return r * pop["first_shipped_bytes"]

    for name, lo, mid, hi in (
        (
            "A1",
            base_b(rb[1]) + (R - E) * min(d1),
            base_b(rb[0]) + (R - E) * d1_mean,
            base_b(rb[2]) + (R - E) * max(d1),
        ),
        (
            "A2",
            base_b(rb[1]) + (G - E) * min(gb) + (R - G) * min(cb),
            base_b(rb[0]) + (G - E) * gb_mean + (R - G) * cb_mean,
            base_b(rb[2]) + (G - E) * max(gb) + (R - G) * max(cb),
        ),
    ):
        out[name] = {"low_gb": lo / 1e9, "mid_gb": mid / 1e9, "high_gb": hi / 1e9}
    out["params"] = {
        "r_base": rb,
        "A1_delta_kb": (min(d1) / 1e3, d1_mean / 1e3, max(d1) / 1e3),
        "A2_group_base_kb": (min(gb) / 1e3, gb_mean / 1e3, max(gb) / 1e3),
        "A2_conf_kb": (min(cb) / 1e3, cb_mean / 1e3, max(cb) / 1e3),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--examples", type=Path, required=True)
    ap.add_argument(
        "--sample", type=Path, required=True, help="JSON of reglyco_compact.py"
    )
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    index = pq.read_table(
        a.index, columns=["path", "shard", "length", "uncompressed_length"]
    )
    ex = pq.read_table(
        a.examples,
        columns=[
            "path",
            "deposition_date",
            "glycan_experiment_split",
            "has_in_vocabulary_glycan",
        ],
    )
    # training view of the configs: deposited before 2024-01-01, split train, in-vocabulary glycan
    cutoff = pa.scalar(datetime.datetime(2024, 1, 1))  # noqa: DTZ001  (the column is tz-naive)
    m = pc.and_(
        pc.and_(
            pc.less(ex["deposition_date"], cutoff),
            pc.equal(ex["glycan_experiment_split"], "train"),
        ),
        pc.fill_null(ex["has_in_vocabulary_glycan"], False),
    )
    train_paths = ex.filter(m)["path"]
    train_mask = pc.is_in(index["path"], value_set=train_paths)
    res = {"full": population(index, None), "train_view": population(index, train_mask)}
    print("train view rows:", res["train_view"]["records"], "(README: 2,307,945)")
    js = json.loads(a.sample.read_text())
    res["models"] = {}
    levels = sorted({r["level"] for r in js["layouts"] if r["layout"] == "A1"})
    for lv in levels:
        a1, a2 = sample_stats(js, "A1", lv), sample_stats(js, "A2", lv)
        res["models"][lv] = {
            k: model(v, a1, a2)
            for k, v in (("full", res["full"]), ("train_view", res["train_view"]))
        }
    a.out.write_text(json.dumps(res, indent=1))
    for k in ("full", "train_view"):
        p = res[k]
        print(
            f"\n[{k}] records {p['records']:,}  entries {p['entries']:,}  groups {p['groups']:,}  "
            f"records/entry {p['records_per_entry_mean']:.0f} (max {p['records_per_entry_max']:,})  "
            f"entries in one shard {p['entries_in_one_shard']:,}  shipped {p['shipped_bytes'] / 1e9:.1f} GB  "
            f"raw {p['raw_bytes'] / 1e12:.2f} TB"
        )
        for lv, mdl in res["models"].items():
            for lay in ("A1", "A2"):
                x = mdl[k][lay]
                print(
                    f"   level {lv:2d} {lay}: {x['mid_gb']:8.1f} GB  (range {x['low_gb']:.1f} - {x['high_gb']:.1f})"
                )
    print("\nparams:", json.dumps(next(iter(res["models"].values()))["full"]["params"]))


if __name__ == "__main__":
    os.environ.setdefault("ARROW_DEFAULT_MEMORY_POOL", "system")
    main()
