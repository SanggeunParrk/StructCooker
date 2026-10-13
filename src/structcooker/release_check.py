"""Deployability check of a release tree (OUTPUT_ROOT): every stage's output, self-contained.

For each stage in the release DAG:

* an LMDB output must sit at its configured path with every physical file inside its own
  directory (a shard collection's shards under ``<db>.lmdb/``, not in build or release
  directories), its entry count must equal ``.meta.json`` and the ``.index.tsv`` rows, and
  a seeded sample of records must decode and pass the schema / structure / MSA checks the
  build verification applies;
* a file or directory output must exist and be non-empty.

Run through ``structcooker release-check`` (a SLURM job); the report goes to
``<OUTPUT_ROOT>/release_check.json``.
"""
from __future__ import annotations

import json
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

SAMPLE = 200


def _files_outside(target: Path) -> list[str]:
    """Physical LMDB files the output reads from outside its own directory."""
    from datacooker.lmdb.sharded import MANIFEST

    root = os.path.realpath(target)
    if (target / MANIFEST).is_file():
        shards = json.loads((target / MANIFEST).read_text())["shards"]
        paths = [os.path.realpath(target / s["path"]) for s in shards]
    else:
        paths = [os.path.realpath(target / "data.mdb")]
    return [p for p in paths if not p.startswith(root + os.sep)]


def check_lmdb(name: str, target: str, schema: str, max_depth: int | None) -> dict[str, Any]:
    from datacooker.lmdb.sharded import open_env

    from structcooker.build_verify import _check_records

    path = Path(target)
    report: dict[str, Any] = {"stage": name, "path": target, "schema": schema, "problems": []}
    problems = report["problems"]
    if not path.exists():
        problems.append("missing")
        return report
    outside = _files_outside(path)
    if outside:
        problems.append(f"{len(outside)} shard(s) outside its directory, e.g. {outside[0]}")
    with open_env(target, readonly=True, lock=False, readahead=False) as env, env.begin() as txn:
        report["entries"] = txn.stat()["entries"]
    meta_path, index_path = Path(target + ".meta.json"), Path(target + ".index.tsv")
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    keys, sizes = [], []
    if index_path.exists():
        with index_path.open() as handle:
            for line in handle:
                if not line.startswith(("#", "key\t")):
                    key, size = line.split("\t", 2)[:2]
                    keys.append(key)
                    sizes.append(int(size))
    if meta.get("entries") != report["entries"]:
        problems.append(f"meta.json entries {meta.get('entries')} != {report['entries']}")
    if sizes and (meta.get("total_bytes"), meta.get("max_bytes")) != (sum(sizes), max(sizes)):
        problems.append(f"meta.json bytes {meta.get('total_bytes')}/{meta.get('max_bytes')} "
                        f"!= index {sum(sizes)}/{max(sizes)}")
    if len(keys) != report["entries"]:
        problems.append(f"index rows {len(keys)} != {report['entries']}")
    if meta.get("schema") and meta["schema"] != schema:
        problems.append(f"meta.json schema {meta['schema']} != config {schema}")
    sample = random.Random(name).sample(keys, min(SAMPLE, len(keys)))
    try:
        report["sampled"] = _check_records(target, sample, schema, max_depth)
    except Exception as exc:  # noqa: BLE001 - a bad record is a finding, not a crash
        problems.append(f"sample check: {type(exc).__name__}: {str(exc)[:200]}")
    return report


def check_files(name: str, paths: list[str]) -> dict[str, Any]:
    problems = []
    for p in map(Path, paths):
        if p.is_dir():
            if not any(c.is_file() and c.stat().st_size for c in p.rglob("*")):
                problems.append(f"empty directory {p}")
        elif not (p.is_file() and p.stat().st_size):
            problems.append(f"missing or empty {p}")
    return {"stage": name, "paths": paths, "problems": problems}


def run(manifest: Path, out: Path) -> dict[str, Any]:
    import yaml
    from datacooker.config import load_config

    from structcooker import preflight
    from structcooker.build_verify import _verify_workers

    repo = Path(__file__).resolve().parents[2]
    jobs = []
    with ProcessPoolExecutor(max_workers=min(16, _verify_workers())) as pool:
        for name in yaml.safe_load(manifest.read_text()):
            cfg = load_config(repo / "db" / f"{name}.yaml")
            target = cfg.get("new_env_path") or cfg.get("env_path")
            if target:
                depth = (cfg.get("parameters") or {}).get("max_depth")
                jobs.append(pool.submit(check_lmdb, name, str(target), str(cfg["schema"]), depth))
            else:
                jobs.append(pool.submit(check_files, name, [str(p) for p in preflight.output_paths(cfg)]))
        reports = [j.result() for j in jobs]
    result = {"stages": len(reports), "with_problems": sum(bool(r["problems"]) for r in reports),
              "reports": sorted(reports, key=lambda r: r["stage"])}
    out.write_text(json.dumps(result, indent=1))
    return result


def main() -> None:
    if not os.environ.get("SLURM_JOB_ID"):
        msg = "release-check reads every output; submit it with `structcooker release-check`"
        raise RuntimeError(msg)
    manifest, out = Path(sys.argv[1]), Path(sys.argv[2])
    result = run(manifest, out)
    for r in result["reports"]:
        mark = "OK  " if not r["problems"] else "FAIL"
        detail = f"{r.get('entries', '')} entries, {r.get('sampled', 0)} sampled" if "entries" in r else ""
        print(f"{mark} {r['stage']:42} {detail}  {'; '.join(r['problems'])}")
    print(f"{result['stages']} stages, {result['with_problems']} with problems -> {out}")
    if result["with_problems"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
