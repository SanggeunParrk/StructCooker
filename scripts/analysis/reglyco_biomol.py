#!/usr/bin/env python3
"""Re-Glyco as ordinary BioMol records: conversion and a transparent delta codec (see docs/reglyco-compact.md).

A Re-Glyco record is an annotated mmCIF. The release keeps ONE value format (BioMol, schema A: ``{assembly_dict,
metadata_dict}``), so the compaction has to work on the BioMol records, not on CIF text. Two steps:

  convert    run the same ``workflows/ingest/cif.py`` recipe that builds ``pdb/cif`` on every record of one shipped
             shard (CCD from ``lmdb/ccd``), serialise with ``to_bytes`` (zstd 6), write one file per record. The
             record's key is the Hub file name without ``.cif``; the deposition date and the resolution come from
             ``examples.parquet`` (the CIF has none), so the recipe's metadata_dict is complete. Everything else in
             ``examples.parquet`` stays in the sidecar table, keyed by ``path``.
  evaluate   measure, on the converted records, (a) the per-record BioMol bytes, (b) a delta codec: the serialised bytes
             of a record compressed with zstd against the serialised bytes of its entry base (A1) or its group base (A2)
             as a raw-content dictionary, and check that every record decodes to an ORDINARY BioMol dict equal to the
             converted one. The array keys of ``to_bytes`` are random UUIDs; the codec renames them to ``a<i>`` in header
             order before the delta (the decoded dict is unchanged, the keys never reach the consumer).

Heavy: run through sbatch (``reglyco_biomol.sbatch``).
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import shutil
import statistics
import sys
import tempfile
import time
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import zstandard as zstd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reglyco_compact as rc

from structcooker.instructions.readers.reglyco import adapt_annotated_cif

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))


# ---------------------------------------------------------------------------------------------------------------- convert
def _key(path: str) -> str:
    return path.rsplit("/", 1)[-1].removesuffix(".cif")


def _count_ops(expr: str) -> int | None:
    """Number of operators in an oper_expression such as ``1``, ``1,2`` or ``(1-3)``; None for a product ``(1-2)(3-4)``."""
    e = expr.strip()
    if ")(" in e:
        return None
    n = 0
    for tok in e.strip("()").split(","):
        if "-" in tok:
            a, b = tok.split("-")
            n += int(b) - int(a) + 1
        else:
            n += 1
    return n


def check_conversion(d: dict, value: dict) -> dict:
    """Compare the BioMol assemblies with the CIF's own _atom_site.

    BioMol records carry every atom of a residue's CCD definition, the unobserved ones with NaN coordinates, so the
    comparison is over atoms with finite coordinates: their count, and the sums of x, y, z, occupancy and B factor.
    """
    site = d["_atom_site"]
    keep = [i for i, t in enumerate(site["type_symbol"]) if t not in ("H", "D")]

    def col(name: str) -> np.ndarray:
        return np.array([float(site[name][i]) for i in keep])

    xyz = np.stack([col("Cartn_x"), col("Cartn_y"), col("Cartn_z")], axis=1)
    occ, bfac = col("occupancy"), col("B_iso_or_equiv")
    asym = np.array([site["label_asym_id"][i] for i in keep])
    gen = d["_pdbx_struct_assembly_gen"]
    res = {
        "asm": 0,
        "finite_count_ok": 0,
        "sums_ok": 0,
        "skipped": 0,
        "bad": [],
        "nan_out": 0,
        "nan_site": 0,
    }
    for aid, expr, asyms in zip(
        gen["assembly_id"], gen["oper_expression"], gen["asym_id_list"], strict=True
    ):
        keys = [k for k in value["assembly_dict"] if k.split("_")[0] == aid]
        n_ops = _count_ops(expr)
        if (
            n_ops != 1 or not keys
        ):  # only identity assemblies are compared; others are counted as skipped
            res["skipped"] += 1
            continue
        res["asm"] += 1
        sel = np.isin(asym, asyms.split(","))
        fin_s = sel & np.isfinite(xyz).all(axis=1)
        o_xyz = np.concatenate(
            [
                np.asarray(value["assembly_dict"][k]["atoms"]["xyz"].value, dtype=float)
                for k in keys
            ]
        )
        o_occ = np.concatenate(
            [
                np.asarray(
                    value["assembly_dict"][k]["atoms"]["occupancy"].value, dtype=float
                )
                for k in keys
            ]
        )
        o_b = np.concatenate(
            [
                np.asarray(
                    value["assembly_dict"][k]["atoms"]["b_factor"].value, dtype=float
                )
                for k in keys
            ]
        )
        fin_o = np.isfinite(o_xyz).all(axis=1)
        res["nan_out"] += int((~fin_o).sum())
        res["nan_site"] += int((sel & ~np.isfinite(xyz).all(axis=1)).sum())
        cnt_ok = int(fin_o.sum()) == int(fin_s.sum())
        tol = 1e-3 * max(1, int(fin_s.sum()))
        sums_ok = (
            np.allclose(o_xyz[fin_o].sum(axis=0), xyz[fin_s].sum(axis=0), atol=tol)
            and abs(np.nansum(o_occ[fin_o]) - np.nansum(occ[fin_s])) <= tol
            and abs(np.nansum(o_b[fin_o]) - np.nansum(bfac[fin_s])) <= tol
        )
        res["finite_count_ok"] += cnt_ok
        res["sums_ok"] += bool(sums_ok)
        if not (cnt_ok and sums_ok) and len(res["bad"]) < 3:
            res["bad"].append(
                [aid, expr, int(fin_o.sum()), int(fin_s.sum()), bool(sums_ok)]
            )
    return res


def _convert_one(args: tuple[str, str, bytes, dict, str, bool]) -> dict:
    """Convert one record; returns stats (and the written file's size) or the error."""
    path, out_dir, data, meta, tmp_dir, inspect = args
    from datacooker import parse_dict

    from structcooker.instructions.readers.cif import dot_transform, get_cif_data
    from structcooker.instructions.transforms.codecs import to_bytes
    from structcooker.workflows.ingest.cif import RECIPE, TARGETS

    key = _key(path)
    t0 = time.time()
    f = Path(tmp_dir) / f"{key}.cif"
    try:
        f.write_bytes(data)
        d = get_cif_data(f)
        t_read = time.time() - t0
        d["_entry"]["id"] = [
            key
        ]  # the record key, not the temp file name the Hub CIF carries
        # AtomWorks writes a subset of the CCD columns (no pdbx_stereo_config). Our recipe then takes the component from
        # the CCD itself, as it does for any CIF without _chem_comp_atom; filling the column would invent chemistry.
        d.pop("_chem_comp_atom", None)
        d.pop("_chem_comp_bond", None)
        diag = adapt_annotated_cif(d)
        d["_pdbx_database_status"] = {
            "recvd_initial_deposition_date": [meta["deposition_date"]]
        }
        if meta.get("resolution") is not None:
            d["_refine"] = {"ls_d_res_high": [str(meta["resolution"])]}
        d["ccd_db_path"] = ROOT / "lmdb/ccd/biomol_CCD.lmdb"
        d["ccd_reference_path"] = ROOT / "materials/reference/reference_ccd.lmdb"
        out = parse_dict(RECIPE, d, key_transform=dot_transform, targets=TARGETS)
        t_recipe = time.time() - t0 - t_read
        value = {k: out[k] for k in TARGETS}
        stats = {
            "atom_site_rows": len(d["_atom_site"]["id"]),
            "atoms_out": 0,
            "occ0_out": 0,
            "nan_b_out": 0,
            "bonds_out": 0,
        }
        for asm in value["assembly_dict"].values():
            at = asm["atoms"]
            occ = np.asarray(at["occupancy"].value)
            stats["atoms_out"] += len(occ)
            stats["occ0_out"] += int((occ == 0).sum())
            stats["nan_b_out"] += int(np.isnan(np.asarray(at["b_factor"].value)).sum())
            stats["bonds_out"] += len(at["bond_type"].value)
        stats["assemblies"] = len(value["assembly_dict"])
        chk = check_conversion(d, value)
        blob = to_bytes(value, level=6)
        (Path(out_dir) / f"{key}.bin").write_bytes(blob)
        res = {
            "key": key,
            "ok": True,
            "bytes": len(blob),
            "t_read": t_read,
            "t_recipe": t_recipe,
            "diag": diag,
            "stats": stats,
            "check": chk,
            "cif_categories": sorted(k for k in d if k.startswith("_")),
        }
        if inspect:
            res["tree"] = _tree(value)
        return res  # noqa: TRY300
    except Exception as e:  # noqa: BLE001  one bad record must not stop the shard; it is counted
        root = e
        while root.__cause__ is not None:  # the recipe engine wraps the real error
            root = root.__cause__
        return {
            "key": key,
            "ok": False,
            "error": f"{type(root).__name__}: {root}"[:300],
            "tb": "".join(traceback.format_exception(root))[-900:],
        }
    finally:
        f.unlink(missing_ok=True)


def _tree(v: object, depth: int = 0) -> object:
    from biomol.core.container import FeatureContainer
    from biomol.core.index import IndexTable

    if isinstance(v, np.ndarray):
        return f"ndarray{v.shape} {v.dtype}"
    if isinstance(v, dict):
        return {k: _tree(x, depth + 1) for k, x in v.items()}
    if isinstance(v, FeatureContainer):
        return {"<FeatureContainer>": _tree(v.to_dict(), depth + 1)}
    if isinstance(v, IndexTable):
        return "<IndexTable>"
    return repr(v)[:80]


def cmd_convert(a: argparse.Namespace) -> None:
    recs = rc.load_records(Path(a.shard), Path(a.index))
    if a.limit:
        recs = recs[:: max(1, len(recs) // a.limit)][: a.limit]
    paths = {r["path"] for r in recs}
    ex = pq.read_table(
        a.examples,
        columns=["path", "deposition_date", "resolution", "example_id"],
        filters=[("path", "in", list(paths))],
    ).to_pylist()
    meta = {}
    for r in ex:
        dd = r["deposition_date"]
        meta[r["path"]] = {
            "deposition_date": dd.strftime("%Y-%m-%d"),
            "resolution": r["resolution"],
            "example_id": r["example_id"],
        }
    missing = paths - set(meta)
    print(
        f"records {len(recs)}, with examples row {len(meta)}, missing {len(missing)}",
        flush=True,
    )
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tmp_root = os.environ.get("TMPDIR", "/tmp")  # noqa: S108  node-local scratch, removed afterwards
    tmp = tempfile.mkdtemp(prefix="reglyco_cif_", dir=tmp_root)
    jobs = [
        (
            r["path"],
            str(out),
            r["data"],
            meta.get(r["path"], {"deposition_date": "1970-01-01", "resolution": None}),
            tmp,
            i == 0,
        )
        for i, r in enumerate(recs)
    ]
    t0 = time.time()
    results = []
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(_convert_one, jobs, chunksize=1)):
            results.append(res)
            if (i + 1) % 25 == 0 or i + 1 == len(jobs):
                print(f"  {i + 1}/{len(jobs)}  {time.time() - t0:.0f}s", flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    okr = [r for r in results if r["ok"]]
    bad = [r for r in results if not r["ok"]]
    errs: dict[str, int] = defaultdict(int)
    for r in bad:
        errs[r["error"]] += 1
    first = next((r for r in okr if "tree" in r), None)
    summary = {
        "records": len(recs),
        "ok": len(okr),
        "failed": len(bad),
        "errors": dict(errs),
        "example_tb": bad[0]["tb"] if bad else None,
        "mean_bytes": statistics.mean(r["bytes"] for r in okr) if okr else None,
        "sum_bytes": sum(r["bytes"] for r in okr),
        "check_sums": {
            k: sum(r["check"][k] for r in okr)
            for k in (
                "asm",
                "finite_count_ok",
                "sums_ok",
                "skipped",
                "nan_out",
                "nan_site",
            )
        },
        "check_bad_examples": [
            (r["key"], r["check"]["bad"]) for r in okr if r["check"]["bad"]
        ][:5],
        "stats_sums": {
            k: sum(r["stats"][k] for r in okr) for k in (okr[0]["stats"] if okr else {})
        },
        "diag_sums": {
            k: sum(r["diag"][k] for r in okr) for k in (okr[0]["diag"] if okr else {})
        },
        "mean_t_read_s": statistics.mean(r["t_read"] for r in okr) if okr else None,
        "mean_t_recipe_s": statistics.mean(r["t_recipe"] for r in okr) if okr else None,
        "wall_s": time.time() - t0,
        "cif_categories_seen": sorted({c for r in okr for c in r["cif_categories"]}),
        "tree_of_first_ok": first["tree"] if first else None,
        "meta": {_key(p): m for p, m in meta.items()},
    }
    (out / "convert_summary.json").write_text(
        json.dumps(summary, indent=1, default=str)
    )
    print(
        json.dumps(
            {
                k: v
                for k, v in summary.items()
                if k not in ("tree_of_first_ok", "meta", "cif_categories_seen")
            },
            indent=1,
            default=str,
        )
    )


# --------------------------------------------------------------------------------------------------------------- evaluate
def split_raw(raw: bytes) -> tuple[dict, bytes]:
    hlen = int.from_bytes(raw[:8], "little")
    return json.loads(raw[8 : 8 + hlen]), raw[8 + hlen :]


def _remap(t: object, m: dict) -> object:
    if isinstance(t, dict):
        return {k: _remap(v, m) for k, v in t.items()}
    if isinstance(t, str) and t in m:
        return m[t]
    return t


def canon(raw: bytes) -> bytes:
    """Same serialised value, array keys renamed a0, a1, ... in header order (to_bytes uses random UUIDs)."""
    header, payload = split_raw(raw)
    m = {k: f"a{i}" for i, k in enumerate(header["arrays"])}
    h = json.dumps(
        {
            "template": _remap(header["template"], m),
            "arrays": {m[k]: n for k, n in header["arrays"].items()},
        }
    ).encode()
    return len(h).to_bytes(8, "little") + h + payload


def from_raw(raw: bytes) -> dict:
    """What ``codecs.from_bytes`` does after the zstd step."""
    from structcooker.instructions.transforms.codecs import reconstruct_data

    header, payload = split_raw(raw)
    flat, off = {}, 0
    for k, n in header["arrays"].items():
        flat[k] = payload[off : off + n]
        off += n
    return reconstruct_data(header["template"], flat)


def same(a: object, b: object) -> bool:
    if isinstance(a, dict):
        return (
            isinstance(b, dict)
            and a.keys() == b.keys()
            and all(same(a[k], b[k]) for k in a)
        )
    if isinstance(a, np.ndarray):
        if not isinstance(b, np.ndarray) or a.dtype != b.dtype or a.shape != b.shape:
            return False
        if a.dtype.kind == "f":
            return bool(np.array_equal(a, b, equal_nan=True))
        return bool(np.array_equal(a, b))
    return bool(a == b)


def load_converted(d: Path) -> list[dict]:
    from structcooker.instructions.transforms.codecs import (
        from_bytes,  # noqa: F401  (import check)
    )

    recs = []
    for f in sorted(d.glob("*.bin")):
        blob = f.read_bytes()
        raw = zstd.ZstdDecompressor().decompress(blob, max_output_size=1 << 31)
        pdb, site, glycan, conf = rc.parse_path(f"xx/{f.stem[:4]}/{f.stem}.cif")
        recs.append(
            {
                "key": f.stem,
                "pdb": pdb,
                "site": site,
                "glycan": glycan,
                "conf": conf,
                "blob": blob,
                "raw": raw,
                "craw": canon(raw),
            }
        )
    return recs


def cmd_evaluate(a: argparse.Namespace) -> None:
    from structcooker.instructions.transforms.codecs import from_bytes

    recs = load_converted(Path(a.dir))
    n = len(recs)
    ref_s = json.loads((Path(a.dir) / "convert_summary.json").read_text())
    print(f"converted records: {n}", flush=True)
    shipped_total = (
        a.shipped_bytes
    )  # bytes the same records occupy as shipped CIF frames (from the index)
    out: dict = {"records": n, "per_record_to_bytes": sum(len(r["blob"]) for r in recs)}
    out["raw_serialised"] = sum(len(r["raw"]) for r in recs)
    out["shipped_cif_zstd"] = shipped_total

    c9 = zstd.ZstdCompressor(level=a.level)
    by_entry: dict[str, list[dict]] = defaultdict(list)
    for r in recs:
        by_entry[r["pdb"]].append(r)
    for v in by_entry.values():
        v.sort(key=lambda r: r["key"])

    def evaluate(two_level: bool, canonical: bool) -> dict:
        field = "craw" if canonical else "raw"
        base_b = gbase_b = conf_b = 0
        t_enc = [0.0, 0]

        def enc_with(codec: object, data: bytes) -> bytes:
            t0 = time.time()
            blob = codec.encode(data)
            t_enc[0] += time.time() - t0
            t_enc[1] += 1
            return blob

        n_gb = n_conf = 0
        per_entry = []
        ok = 0
        t_dec = []
        for pdb, rs in by_entry.items():
            eb = rs[0]
            enc = len(c9.compress(eb[field]))
            max_len = max(len(x[field]) for x in rs)
            ecodec = rc.DeltaCodec(eb[field], a.level, max_len)
            groups: dict[tuple, list[dict]] = defaultdict(list)
            for r in rs[1:]:
                groups[(r["site"], r["glycan"])].append(r)
            # the entry base's own group: its other conformers are deltas against it in A2
            gb_sum = cf_sum = 0
            ngb = ncf = 0
            stored: dict[str, tuple[str, bytes]] = {eb["key"]: ("plain", b"")}
            own = (eb["site"], eb["glycan"])
            for g, members in groups.items():
                members = sorted(members, key=lambda r: r["conf"])
                if two_level and g != own:
                    gbase = members[0]
                    blob = enc_with(ecodec, gbase[field])
                    gb_sum += len(blob)
                    ngb += 1
                    stored[gbase["key"]] = ("e", blob)
                    gcodec = rc.DeltaCodec(gbase[field], a.level, max_len)
                    for r in members[1:]:
                        blob = enc_with(gcodec, r[field])
                        cf_sum += len(blob)
                        ncf += 1
                        stored[r["key"]] = (gbase["key"], blob)
                elif two_level and g == own:
                    gcodec = rc.DeltaCodec(eb[field], a.level, max_len)
                    for r in members:
                        blob = enc_with(gcodec, r[field])
                        cf_sum += len(blob)
                        ncf += 1
                        stored[r["key"]] = (eb["key"], blob)
                else:
                    for r in members:
                        blob = enc_with(ecodec, r[field])
                        cf_sum += len(blob)
                        ncf += 1
                        stored[r["key"]] = (eb["key"], blob)
            by_key = {r["key"]: r for r in rs}
            for r in rs:  # decode every record through its chain and compare with the converted dict
                t0 = time.time()
                kind, blob = stored[r["key"]]
                if kind == "plain":
                    raw = eb[field]
                elif kind == "e":
                    raw = ecodec.decode(blob)
                else:
                    gb = by_key[kind]
                    gcodec = rc.DeltaCodec(gb[field], a.level, max_len, encode=False)
                    raw = gcodec.decode(blob)
                got = from_raw(raw)
                t_dec.append(time.time() - t0)
                want = from_bytes(r["blob"])
                ok += same(got, want)
            base_b += enc
            gbase_b += gb_sum
            conf_b += cf_sum
            n_gb += ngb
            n_conf += ncf
            per_entry.append(
                {
                    "pdb": pdb,
                    "records": len(rs),
                    "base_bytes": enc,
                    "base_raw_bytes": len(eb[field]),
                    "group_bases": ngb,
                    "group_base_mean": gb_sum / ngb if ngb else None,
                    "conformers": ncf,
                    "conf_mean": cf_sum / ncf if ncf else None,
                }
            )
        return {
            "total_bytes": base_b + gbase_b + conf_b,
            "entry_bases_bytes": base_b,
            "group_bases_bytes": gbase_b,
            "conformer_bytes": conf_b,
            "n_group_bases": n_gb,
            "n_conformers": n_conf,
            "roundtrip_ok": ok,
            "roundtrip_total": n,
            "decode_ms_mean": 1000 * statistics.mean(t_dec),
            "encode_ms_mean": 1000 * t_enc[0] / max(1, t_enc[1]),
            "per_entry": per_entry,
        }

    out["A1_canonical"] = evaluate(False, True)
    print("A1 done", flush=True)
    out["A2_canonical"] = evaluate(True, True)
    print("A2 done", flush=True)
    out["A1_uuid_keys"] = evaluate(False, False)
    print("A1 uuid done", flush=True)
    out["convert_ok"] = ref_s["ok"]
    out["convert_failed"] = ref_s["failed"]
    Path(a.out).write_text(json.dumps(out, indent=1))
    for k in ("A1_canonical", "A2_canonical", "A1_uuid_keys"):
        v = out[k]
        print(
            f"{k}: {v['total_bytes'] / 1e6:.2f} MB  roundtrip {v['roundtrip_ok']}/{v['roundtrip_total']}  "
            f"decode {v['decode_ms_mean']:.1f} ms"
        )
    print(
        f"per-record to_bytes {out['per_record_to_bytes'] / 1e6:.1f} MB, shipped CIF {shipped_total / 1e6:.1f} MB"
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sp = ap.add_subparsers(required=True)
    c = sp.add_parser("convert")
    c.add_argument("shard")
    c.add_argument("index")
    c.add_argument("--examples", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--workers", type=int, default=8)
    c.add_argument(
        "--limit",
        type=int,
        default=0,
        help="convert only about this many records (spread over the shard)",
    )
    c.set_defaults(fn=cmd_convert)
    e = sp.add_parser("evaluate")
    e.add_argument("dir")
    e.add_argument("--out", required=True)
    e.add_argument("--level", type=int, default=9)
    e.add_argument("--shipped-bytes", type=int, default=0)
    e.set_defaults(fn=cmd_evaluate)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
