#!/usr/bin/env python3
"""Extrapolate the BioMol-format Re-Glyco measurements (reglyco_biomol.py) to the whole release and the training view.

Population numbers (entries, (site, glycan) groups, records, shipped bytes of the first record of each entry) are the
exact ones of ``reglyco_extrapolate.py`` (``extrapolation.json``). Sizes are a model on top of the few-entry sample:

  plain BioMol (one ``to_bytes`` value per record) = measured ratio to the shipped CIF frame x shipped bytes
  entry base           = r_base x shipped length of the entry's first record      (r_base measured, per entry)
  A1 delta per record  = measured mean of the sample                              (records - entries of them)
  A2 group base        = measured mean of the sample                              (groups - entries of them)
  A2 conformer delta   = measured mean of the sample                              (records - groups of them)

Low / high are the smallest / largest per-entry means of the sample (indicative with five entries, not a confidence
interval). Heavy (reads the 3M-row Hub index): run through sbatch (``reglyco_biomol_extrapolate.sbatch``).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "workdir",
        type=Path,
        help="reglyco_biomol/<shard> with the converted *.bin and evaluate_L*.json",
    )
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument(
        "--population",
        required=True,
        type=Path,
        help="reglyco_compact/extrapolation.json",
    )
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    keys = {f.stem: f.stat().st_size for f in a.workdir.glob("*.bin")}
    t = pq.read_table(a.index, columns=["path", "length"])
    names = pc.replace_substring_regex(
        t["path"].combine_chunks(), r"^.*/(.*)\.cif$", r"\1"
    )
    shipped = {}
    want = set(keys)
    for k, n in zip(names.to_pylist(), t["length"].to_pylist(), strict=True):
        if k in want:
            shipped[k] = n
    pdb_of = {k: k[:4] for k in keys}
    ent: dict[str, dict] = {}
    for k in sorted(keys):
        e = ent.setdefault(pdb_of[k], {"plain": 0, "shipped": 0, "first": k})
        e["plain"] += keys[k]
        e["shipped"] += shipped[k]
    plain_ratio_entries = [e["plain"] / e["shipped"] for e in ent.values()]
    plain_ratio = sum(e["plain"] for e in ent.values()) / sum(
        e["shipped"] for e in ent.values()
    )

    pop = json.loads(a.population.read_text())
    out: dict = {
        "sample": {
            "entries": len(ent),
            "records": len(keys),
            "plain_to_shipped_ratio": [
                plain_ratio,
                min(plain_ratio_entries),
                max(plain_ratio_entries),
            ],
        }
    }
    for view in ("full", "train_view"):
        p = pop[view]
        out[view] = {
            "plain_biomol_gb": [
                plain_ratio * p["shipped_bytes"] / 1e9,
                min(plain_ratio_entries) * p["shipped_bytes"] / 1e9,
                max(plain_ratio_entries) * p["shipped_bytes"] / 1e9,
            ]
        }
    for ev in sorted(a.workdir.glob("evaluate_L*.json")):
        level = ev.stem.removeprefix("evaluate_L")
        e = json.loads(ev.read_text())
        for variant in ("A1_canonical", "A2_canonical"):
            v = e[variant]
            r_base = [
                pe["base_bytes"]
                / shipped[sorted(k for k in keys if pdb_of[k] == pe["pdb"])[0]]
                for pe in v["per_entry"]
            ]
            conf = [
                pe["conf_mean"] for pe in v["per_entry"] if pe["conf_mean"] is not None
            ]
            gb = [
                pe["group_base_mean"]
                for pe in v["per_entry"]
                if pe["group_base_mean"] is not None
            ]
            tot_r = v["entry_bases_bytes"] / sum(
                shipped[sorted(k for k in keys if pdb_of[k] == pe["pdb"])[0]]
                for pe in v["per_entry"]
            )
            conf_mid = v["conformer_bytes"] / v["n_conformers"]
            gb_mid = (
                v["group_bases_bytes"] / v["n_group_bases"]
                if v["n_group_bases"]
                else 0.0
            )
            res = {
                "r_base": [tot_r, min(r_base), max(r_base)],
                "conf_mean_b": [conf_mid, min(conf), max(conf)],
                "group_base_mean_b": [gb_mid, min(gb), max(gb)] if gb else None,
            }
            for view in ("full", "train_view"):
                p = pop[view]
                n_gb = p["groups"] - p["entries"] if gb else 0
                n_conf = (
                    p["records"] - p["groups"] if gb else p["records"] - p["entries"]
                )
                tot = []
                for i in range(3):
                    tot.append(
                        (
                            res["r_base"][i] * p["first_shipped_bytes"]
                            + n_gb * (res["group_base_mean_b"][i] if gb else 0)
                            + n_conf * res["conf_mean_b"][i]
                        )
                        / 1e9
                    )
                res[f"{view}_gb"] = tot  # mid, low, high
            res["encode_ms_mean"] = v.get("encode_ms_mean")
            res["decode_ms_mean"] = v["decode_ms_mean"]
            out.setdefault("models", {}).setdefault(level, {})[variant] = res
    a.out.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
