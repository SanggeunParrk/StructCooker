"""Write a release tree's own documentation: ``<BIOMOL_ROOT>/docs``.

Hand-written pages (``release_docs/README.md``, ``provenance.md``) are copied; the rest is
generated from this tree, so it cannot drift from what is on disk:

* ``databases.md`` -- every database of the release DAG (``db/MANIFEST_all.yaml``): path,
  schema, entries (``.meta.json``), size on disk (its ``data.mdb`` files) and what it is
  built from;
* ``scenarios.md`` -- ``release_docs/scenarios.yaml`` with the sizes of this tree;
* ``check.md`` -- the latest ``release_check.json``.

Run ``structcooker release-docs`` after a build, a patch or a move.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SOURCES = REPO / "release_docs"


def disk_bytes(path: Path) -> int | None:
    """Bytes an LMDB (its data.mdb files, shards included) or a file takes; None for other dirs."""
    if path.is_file():
        return path.stat().st_size
    if (path / "data.mdb").is_file():
        return (path / "data.mdb").stat().st_size
    if (path / "shards.json").is_file():
        shards = json.loads((path / "shards.json").read_text())["shards"]
        return sum((path / s["path"] / "data.mdb").stat().st_size for s in shards)
    return None


def human(n: int | None) -> str:
    if n is None:
        return "—"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1000 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1000
    return ""


def _rel(path: Any, root: Path) -> str:
    try:
        return str(Path(str(path)).relative_to(root))
    except ValueError:
        return str(path)


def databases(root: Path, manifest: Path) -> tuple[str, int]:
    import yaml
    from datacooker.config import load_config

    from structcooker import schemas
    from structcooker.workflow_ops import op_of

    rows = []
    for name in yaml.safe_load(manifest.read_text()):
        cfg = load_config(REPO / "db" / f"{name}.yaml")
        target = cfg.get("new_env_path") or cfg.get("env_path")
        if not target:
            continue
        target = Path(str(target))
        meta_path = Path(str(target) + ".meta.json")
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        op = op_of(cfg)
        source = cfg.get("old_env_path") or cfg.get("file_list") or cfg.get("data_dir") or ""
        rows.append((_rel(target, root), name, str(cfg.get("schema", "")), meta.get("entries"),
                     disk_bytes(target) if target.exists() else None, op, _rel(source, root)))
    # Databases made outside the build DAG (FoldBench assets, the inference CCD): any other
    # <db>.lmdb with a .meta.json, outside hidden and superseded directories.
    known = {r[0] for r in rows}
    for meta_path in sorted((root / "lmdb").rglob("*.lmdb.meta.json")):
        rel = _rel(Path(str(meta_path)[: -len(".meta.json")]), root)
        if rel in known or any(part.startswith((".", "_")) for part in Path(rel).parts):
            continue
        meta = json.loads(meta_path.read_text())
        rows.append((rel, "(outside the DAG)", str(meta.get("schema", "")), meta.get("entries"),
                     disk_bytes(root / rel), "staging", str(meta.get("source") or "")))
    rows.sort()
    total = sum(r[4] or 0 for r in rows)
    lines = ["# Databases", "",
             f"{len(rows)} databases, {human(total)} on disk. Sizes are the LMDB files; entries come",
             "from each database's `.meta.json`. Schemas are described in the [README](README.md).", ""]
    group = None
    for path, name, schema, entries, size, op, source in rows:
        top = "/".join(path.split("/")[:2])
        if top != group:
            group = top
            lines += ["", f"## `{top}`", "", "| database | schema | entries | size | built by | from |",
                      "|---|---|---|---|---|---|"]
        title = f"{schema} — {schemas.SCHEMAS[schema].title}" if schema in schemas.SCHEMAS else schema
        count = f"{entries:,}" if isinstance(entries, int) else "—"
        lines.append(f"| `{path}` | {title} | {count} | {human(size)} | `{name}` ({op}) | {f"`{source}`" if source else ""} |")
    return "\n".join(lines) + "\n", total


def scenarios(root: Path, spec: Path) -> str:
    import yaml

    out = ["# What to copy", "",
           "Sizes are this tree's. MSA databases come in depth-capped versions (`_d512`, `_d2k`,",
           "`_d5k`, `_d16k`, or uncapped); one per molecule type is enough.", ""]
    for item in yaml.safe_load(spec.read_text()):
        out += [f"## {item['title']}", "", f"*{item['when']}*", ""]
        for part, label in (("take", "needed"), ("optional", "optional")):
            entries = item.get(part) or []
            if not entries:
                continue
            sizes = [disk_bytes(root / e["path"]) if (root / e["path"]).exists() else None for e in entries]
            total = sum(s or 0 for s in sizes)
            known = all(s is not None for s in sizes)
            out += [f"**{label}**{f' — {human(total)}' if known else ''}", "", "| path | size | |", "|---|---|---|"]
            for e, s in zip(entries, sizes, strict=True):
                missing = "" if (root / e["path"]).exists() else " (missing)"
                out.append(f"| `{e['path']}`{missing} | {human(s)} | {e.get('note', '')} |")
            out.append("")
        if item.get("notes"):
            out += [item["notes"], ""]
    return "\n".join(out)


def check(root: Path) -> str:
    report_path = root / "release_check.json"
    if not report_path.exists():
        return "# Release check\n\nNo `release_check.json` yet: run `structcooker release-check`.\n"
    report = json.loads(report_path.read_text())
    mtime = __import__("datetime").datetime.fromtimestamp(report_path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    lines = ["# Release check", "",
             f"`structcooker release-check`, {mtime}: **{report['stages']} stages, "
             f"{report['with_problems']} with problems.**", "",
             "Each database: at its configured path, shards inside its own directory, entries and bytes",
             "agreeing with `.meta.json` and `.index.tsv`, a seeded sample of 200 records decoding and",
             "passing the build verification's checks. File outputs: present and non-empty.", ""]
    bad = [r for r in report["reports"] if r["problems"]]
    if bad:
        lines += ["| stage | problems |", "|---|---|"] + [f"| `{r['stage']}` | {'; '.join(r['problems'])} |" for r in bad]
    return "\n".join(lines) + "\n"


def write(root: Path, manifest: Path) -> list[Path]:
    out = root / "docs"
    out.mkdir(exist_ok=True)
    written = []
    for name in ("README.md", "provenance.md"):
        shutil.copyfile(SOURCES / name, out / name)
        written.append(out / name)
    text, _total = databases(root, manifest)
    for name, body in (("databases.md", text), ("scenarios.md", scenarios(root, SOURCES / "scenarios.yaml")),
                       ("check.md", check(root))):
        (out / name).write_text(body)
        written.append(out / name)
    return written
