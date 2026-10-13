"""Write a release tree's own documentation: ``<BIOMOL_ROOT>/docs``.

Hand-written pages (``release_docs/README.md``, ``provenance.md``, ``licenses.md``) are copied;
the rest is generated from this tree, so it cannot drift from what is on disk:

* ``databases.md`` -- every database (the release DAG's, ``db/MANIFEST_all.yaml``, and any
  other ``<db>.lmdb`` with a ``.meta.json``): path, schema, entries, size on disk (its
  ``data.mdb`` files) and what it is built from;
* ``scenarios.md`` -- ``release_docs/scenarios.yaml`` with the sizes of this tree;
* ``check.md`` -- the latest ``release_check.json``;
* ``index.html`` -- all of the above on one page (overview, sets, what to copy, every DB).

The tables are computed once and rendered twice, so the Markdown and the HTML agree.
Run ``structcooker release-docs`` after a build, a patch or a move.
"""
from __future__ import annotations

import datetime
import html
import json
import shutil
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SOURCES = REPO / "release_docs"
PAGES = ("README.md", "provenance.md", "licenses.md")


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


def human(n: float | None) -> str:
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


# --------------------------------------------------------------------------- data
def database_rows(root: Path, manifest: Path) -> list[dict[str, Any]]:
    """One row per database: path, stage, schema (+ title), entries, size, op, source."""
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
        source = cfg.get("old_env_path") or cfg.get("file_list") or cfg.get("data_dir") or ""
        rows.append({"path": _rel(target, root), "stage": name, "schema": str(cfg.get("schema", "")),
                     "entries": meta.get("entries"), "size": disk_bytes(target) if target.exists() else None,
                     "op": op_of(cfg), "source": _rel(source, root)})
    # Databases made outside the build DAG (FoldBench assets, the inference CCD): any other
    # <db>.lmdb with a .meta.json, outside hidden and superseded directories.
    known = {r["path"] for r in rows}
    for meta_path in sorted((root / "lmdb").rglob("*.lmdb.meta.json")):
        rel = _rel(Path(str(meta_path)[: -len(".meta.json")]), root)
        if rel in known or any(part.startswith((".", "_")) for part in Path(rel).parts):
            continue
        meta = json.loads(meta_path.read_text())
        rows.append({"path": rel, "stage": "", "schema": str(meta.get("schema", "")), "entries": meta.get("entries"),
                     "size": disk_bytes(root / rel), "op": "staging", "source": str(meta.get("source") or "")})
    for r in rows:
        r["group"] = "/".join(r["path"].split("/")[:2])
        r["title"] = schemas.SCHEMAS[r["schema"]].title if r["schema"] in schemas.SCHEMAS else ""
    return sorted(rows, key=lambda r: r["path"])


def scenario_data(root: Path, spec: Path) -> list[dict[str, Any]]:
    """scenarios.yaml with each path's size in this tree and a total per part."""
    import yaml

    items = yaml.safe_load(spec.read_text())
    for item in items:
        for part in ("take", "optional"):
            entries = item.get(part) or []
            for e in entries:
                e["exists"] = (root / e["path"]).exists()
                e["size"] = disk_bytes(root / e["path"]) if e["exists"] else None
            item[f"{part}_total"] = (sum(e["size"] for e in entries)
                                     if entries and all(e["size"] is not None for e in entries) else None)
    return items


def check_data(root: Path) -> dict[str, Any] | None:
    path = root / "release_check.json"
    if not path.exists():
        return None
    report = json.loads(path.read_text())
    report["when"] = datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return report


# --------------------------------------------------------------------------- markdown
def databases_md(rows: list[dict[str, Any]]) -> str:
    total = sum(r["size"] or 0 for r in rows)
    lines = ["# Databases", "",
             f"{len(rows)} databases, {human(total)} on disk. Sizes are the LMDB files; entries come",
             "from each database's `.meta.json`. Schemas are described in the [README](README.md).", ""]
    group = None
    for r in rows:
        if r["group"] != group:
            group = r["group"]
            lines += ["", f"## `{group}`", "", "| database | schema | entries | size | built by | from |",
                      "|---|---|---|---|---|---|"]
        schema = f"{r['schema']} — {r['title']}" if r["title"] else r["schema"]
        count = f"{r['entries']:,}" if isinstance(r["entries"], int) else "—"
        built = f"`{r['stage']}` ({r['op']})" if r["stage"] else r["op"]
        source = f"`{r['source']}`" if r["source"] else ""
        lines.append(f"| `{r['path']}` | {schema} | {count} | {human(r['size'])} | {built} | {source} |")
    return "\n".join(lines) + "\n"


def scenarios_md(items: list[dict[str, Any]]) -> str:
    out = ["# What to copy", "",
           "Sizes are this tree's. MSA databases come in depth-capped versions (`_d512`, `_d2k`,",
           "`_d5k`, `_d16k`, or uncapped); one per molecule type is enough.", ""]
    for item in items:
        out += [f"## {item['title']}", "", f"*{item['when']}*", ""]
        for part, label in (("take", "needed"), ("optional", "optional")):
            entries = item.get(part) or []
            if not entries:
                continue
            total = item[f"{part}_total"]
            out += [f"**{label}**{f' — {human(total)}' if total is not None else ''}", "", "| path | size | |", "|---|---|---|"]
            out += [f"| `{e['path']}`{'' if e['exists'] else ' (missing)'} | {human(e['size'])} | {e.get('note', '')} |"
                    for e in entries]
            out.append("")
        if item.get("notes"):
            out += [item["notes"], ""]
    return "\n".join(out)


def check_md(report: dict[str, Any] | None) -> str:
    if report is None:
        return "# Release check\n\nNo `release_check.json` yet: run `structcooker release-check`.\n"
    lines = ["# Release check", "",
             f"`structcooker release-check`, {report['when']}: **{report['stages']} stages, "
             f"{report['with_problems']} with problems.**", "",
             "Each database: at its configured path, shards inside its own directory, entries and bytes",
             "agreeing with `.meta.json` and `.index.tsv`, a seeded sample of 200 records decoding and",
             "passing the build verification's checks. File outputs: present and non-empty.", ""]
    bad = [r for r in report["reports"] if r["problems"]]
    if bad:
        lines += ["| stage | problems |", "|---|---|"] + [f"| `{r['stage']}` | {'; '.join(r['problems'])} |" for r in bad]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- html
_CSS = """
:root{--bg:#fbfbfa;--card:#fff;--fg:#1d1f23;--muted:#5d6470;--line:#e3e5e8;--accent:#2d5bd7;--ok:#1f8a4c;--bad:#c2410c;--code:#f2f3f5}
@media (prefers-color-scheme:dark){:root{--bg:#121417;--card:#1a1d21;--fg:#e6e8eb;--muted:#9aa3ae;--line:#2b3036;--accent:#7aa2ff;--ok:#4ade80;--bad:#fb923c;--code:#22262b}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Noto Sans KR",sans-serif}
main{max-width:1180px;margin:0 auto;padding:32px 16px 64px}h1{font-size:28px;margin:0 0 4px}h2{font-size:19px;margin:40px 0 12px}
.sub{color:var(--muted);margin:0 0 24px}code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12.5px}
code{background:var(--code);padding:1px 5px;border-radius:4px}pre{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px;overflow-x:auto;line-height:1.45}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .k{color:var(--muted);font-size:12.5px}.card .v{font-size:24px;font-weight:600;margin-top:2px}.ok{color:var(--ok)}.bad{color:var(--bad)}
table{width:100%;border-collapse:collapse;background:var(--card);border:1px solid var(--line);border-radius:8px;overflow:hidden;font-size:13.5px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}th{color:var(--muted);font-weight:600;font-size:12.5px;background:var(--code)}
td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.grp td{background:var(--code);font-weight:600}
.wrap{overflow-x:auto}.scen{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:12px}
.scen .card h3{margin:0 0 2px;font-size:15.5px}.scen .when{color:var(--muted);font-size:13px;margin-bottom:8px}.scen ul{margin:6px 0 0;padding-left:18px;font-size:13px}
.scen li{margin:2px 0}.tot{font-weight:600}.note{color:var(--muted)}input{width:100%;padding:8px 10px;margin:0 0 10px;border:1px solid var(--line);border-radius:8px;background:var(--card);color:var(--fg);font:inherit}
nav a{color:var(--accent);text-decoration:none;margin-right:14px}a{color:var(--accent)}
"""

_TREE = """<root>/
├── docs/             this page and the Markdown pages beside it
├── lmdb/             the databases (catalogue below)
├── metadata/         released tables: seq_id_map, per-DB clusters, cif_metadata, ...
├── materials/
│   ├── raw/          downloads exactly as fetched
│   ├── reference/    supplied inputs the downloads cannot reproduce
│   └── intermediate/ anything derived on the way
└── release_check.json"""


def _e(text: Any) -> str:
    return html.escape(str(text))


def render_html(root: Path, rows: list[dict[str, Any]], items: list[dict[str, Any]],
                report: dict[str, Any] | None) -> str:
    total = sum(r["size"] or 0 for r in rows)
    records = sum(r["entries"] or 0 for r in rows)
    groups: dict[str, dict[str, Any]] = {}
    for r in rows:
        g = groups.setdefault(r["group"], {"dbs": 0, "size": 0})
        g["dbs"] += 1
        g["size"] += r["size"] or 0
    if report is None:
        check = '<span class="bad">not run</span>'
    elif report["with_problems"]:
        check = f'<span class="bad">{report["with_problems"]} problem(s)</span>'
    else:
        check = f'<span class="ok">{report["stages"]}/{report["stages"]} OK</span>'
    when = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out = [f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>BioMol release</title><style>{_CSS}</style></head>
<body><main><h1>BioMol release</h1>
<p class="sub"><code>{_e(root)}</code> · generated {when} by <code>structcooker release-docs</code></p>
<nav><a href="README.md">README</a><a href="databases.md">databases</a><a href="scenarios.md">what to copy</a>
<a href="provenance.md">provenance</a><a href="licenses.md">licenses</a><a href="check.md">check</a></nav>
<h2>At a glance</h2><div class="cards">
<div class="card"><div class="k">databases</div><div class="v">{len(rows)}</div></div>
<div class="card"><div class="k">on disk</div><div class="v">{human(total)}</div></div>
<div class="card"><div class="k">records</div><div class="v">{records / 1e6:,.1f} M</div></div>
<div class="card"><div class="k">data sets</div><div class="v">{len(groups)}</div></div>
<div class="card"><div class="k">release-check{f" · {_e(report['when'])}" if report else ""}</div><div class="v">{check}</div></div>
</div>
<h2>Layout</h2><pre>{_e(_TREE)}</pre>
<h2>Data sets</h2><div class="wrap"><table><tr><th>lmdb/…</th><th>databases</th><th>size</th></tr>"""]
    out += [f'<tr><td><code>{_e(g)}</code></td><td class="n">{v["dbs"]}</td><td class="n">{human(v["size"])}</td></tr>'
            for g, v in sorted(groups.items())]
    out.append('</table></div><h2>What to copy</h2><div class="scen">')
    for item in items:
        need = item.get("take_total")
        out.append(f'<div class="card"><h3>{_e(item["title"])}</h3><div class="when">{_e(item["when"])}</div>'
                   f'<div class="tot">needed: {human(need) if need is not None else "see list"}</div><ul>')
        out += [f'<li><code>{_e(e["path"])}</code> <span class="note">{human(e["size"])}'
                f'{" · " + _e(e["note"]) if e.get("note") else ""}</span></li>' for e in item.get("take") or []]
        out.append("</ul>")
        if item.get("optional"):
            opt = item.get("optional_total")
            out.append(f'<div class="note" style="margin-top:8px">optional{f" ({human(opt)})" if opt is not None else ""}:</div><ul>')
            out += [f'<li><code>{_e(e["path"])}</code> <span class="note">{human(e["size"])}'
                    f'{" · " + _e(e["note"]) if e.get("note") else ""}</span></li>' for e in item["optional"]]
            out.append("</ul>")
        if item.get("notes"):
            out.append(f'<div class="note" style="margin-top:8px">{_e(item["notes"])}</div>')
        out.append("</div>")
    out.append('</div><h2>Every database</h2><input id="q" placeholder="filter: path, stage, schema…" '
               'oninput="f(this.value)"><div class="wrap"><table id="t"><tr><th>database</th><th>schema</th>'
               '<th>entries</th><th>size</th><th>built by</th><th>from</th></tr>')
    group = None
    for r in rows:
        if r["group"] != group:
            group = r["group"]
            out.append(f'<tr class="grp"><td colspan="6"><code>{_e(group)}</code></td></tr>')
        count = f"{r['entries']:,}" if isinstance(r["entries"], int) else "—"
        schema = f"{r['schema']} · {r['title']}" if r["title"] else r["schema"]
        built = f"{r['stage']} ({r['op']})" if r["stage"] else r["op"]
        out.append(f'<tr><td><code>{_e(r["path"])}</code></td><td>{_e(schema)}</td><td class="n">{count}</td>'
                   f'<td class="n">{human(r["size"])}</td><td>{_e(built)}</td><td><code>{_e(r["source"])}</code></td></tr>')
    out.append("""</table></div></main><script>
function f(q){q=q.toLowerCase();for(const tr of document.querySelectorAll('#t tr')){if(tr.querySelector('th')||tr.classList.contains('grp'))continue;
tr.style.display=tr.textContent.toLowerCase().includes(q)?'':'none'}}</script></body></html>""")
    return "\n".join(out)


def write(root: Path, manifest: Path) -> list[Path]:
    out = root / "docs"
    out.mkdir(exist_ok=True)
    written = []
    for name in PAGES:
        shutil.copyfile(SOURCES / name, out / name)
        written.append(out / name)
    rows = database_rows(root, manifest)
    items = scenario_data(root, SOURCES / "scenarios.yaml")
    report = check_data(root)
    for name, body in (("databases.md", databases_md(rows)), ("scenarios.md", scenarios_md(items)),
                       ("check.md", check_md(report)), ("index.html", render_html(root, rows, items, report))):
        (out / name).write_text(body)
        written.append(out / name)
    return written
