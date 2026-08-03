"""structcooker CLI — declarative, planning-first BioMol DB builds.

One BioMol DB = one ``db/**/<name>.yaml`` = a datacooker *engine config*
(recipe / reader / writer + source + target) plus a single ``schema:`` tag.

``structcooker build <name>`` hands the config to ``datacooker pipeline``, which
sizes the items, snake-balances them across nodes, and submits tier-arrays → merge →
index, all afterok-chained. Within-node memory is **Ray's** admission control (the
engine's only fan-out backend), not a predicted budget, so no human tunes n_jobs /
mem / chunk / shards. The ``schema:`` tag records each DB's value type and its (now
advisory) expansion factor E from :mod:`structcooker.schemas`.

The ``schema:`` key is a structcooker concept, not a datacooker one, so build/rebuild
must not see it: we materialize a schema-stripped ``engine.yaml`` in the workdir and
hand *that* to the pipeline, passing the schema + E along as ``--schema/--expansion``.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import click
from omegaconf import OmegaConf

from datacooker.config import load_config

from structcooker import schemas

REPO = Path(__file__).resolve().parents[2]          # StructCooker/
DB_ROOT = REPO / "db"
PIPELINE_CLI = [sys.executable, "-u", "-m", "datacooker.cli.lmdb", "pipeline"]
_SCHEMA_LINE = re.compile(r"^\s*schema\s*:.*$", re.MULTILINE)


def _resolve_db(name: str) -> Path:
    """Find a db config by file path, ``sub/dir/name``, or bare ``name``."""
    p = Path(name)
    if p.suffix in {".yaml", ".yml"} and p.exists():
        return p
    direct = DB_ROOT / f"{name}.yaml"
    if direct.exists():
        return direct
    hits = sorted(DB_ROOT.rglob(f"{Path(name).name}.yaml"))
    if not hits:
        msg = f"no db config for {name!r} under {DB_ROOT}"
        raise click.ClickException(msg)
    if len(hits) > 1:
        joined = ", ".join(str(h.relative_to(DB_ROOT)) for h in hits)
        msg = f"{name!r} is ambiguous — qualify it: {joined}"
        raise click.ClickException(msg)
    return hits[0]


def _schema_of(cfg: dict, db_path: Path) -> str:
    schema = cfg.get("schema")
    if schema is None:
        msg = f"{db_path} has no `schema:` (need one of {sorted(schemas.SCHEMAS)})"
        raise click.ClickException(msg)
    schemas.get(str(schema))            # validates the name
    return str(schema)


def _write_engine_yaml(db_path: Path, workdir: Path) -> Path:
    """Copy the db yaml minus its ``schema:`` line so datacooker never sees it."""
    workdir.mkdir(parents=True, exist_ok=True)
    stripped = _SCHEMA_LINE.sub("", db_path.read_text(encoding="utf-8"))
    engine = workdir / "engine.yaml"
    engine.write_text(stripped, encoding="utf-8")
    return engine


@click.group()
def cli() -> None:
    """Declarative, planning-first BioMol DB builds."""


@cli.command("list")
def list_dbs() -> None:
    """List every db/*.yaml with its inferred op, schema, and target."""
    rows: list[tuple[str, str, str, str]] = []
    for f in sorted(DB_ROOT.rglob("*.yaml")):
        if f.name == "MANIFEST.yaml":
            continue
        try:
            cfg = load_config(f)
        except Exception as exc:                       # noqa: BLE001
            rows.append((str(f.relative_to(DB_ROOT)), "ERR", str(exc)[:32], ""))
            continue
        op = "rebuild" if cfg.get("old_env_path") else "build"
        target = cfg.get("new_env_path") or cfg.get("env_path") or "?"
        rows.append((str(f.relative_to(DB_ROOT)), op, str(cfg.get("schema", "-")), str(target)))
    if not rows:
        click.echo(f"(no db configs under {DB_ROOT})")
        return
    w = max(len(r[0]) for r in rows)
    for name, op, schema, target in rows:
        click.echo(f"{name:<{w}}  {op:<7}  schema={schema:<3}  {target}")


_TERMINAL_RE = re.compile(r"^PIPELINE_TERMINAL_JOB=(.*)$", re.MULTILINE)


def _submit_pipeline(
    name: str,
    workdir: Path | None,
    *,
    dry_run: bool,
    depends_on: tuple[str, ...] = (),
    condition: str | None = None,
) -> tuple[int, str | None]:
    """Materialize the engine config for ``name`` and submit its pipeline.

    Returns ``(returncode, terminal_job_id)`` -- the terminal id (the pipeline's
    final index stage) lets ``build-all`` afterok-chain the next config onto it. When a
    ``condition`` (dotted predicate) evaluates False in the engine, the pipeline no-ops
    and the terminal id comes back ``None`` (the caller treats that as "condition said
    skip, nothing to chain onto"). Streams the pipeline output through while capturing it.
    """
    db_path = _resolve_db(name)
    cfg = load_config(db_path)
    schema = _schema_of(cfg, db_path)
    expansion = schemas.expansion(schema)
    op = "rebuild" if cfg.get("old_env_path") else "build"
    wd = Path(workdir) if workdir else REPO / "logs" / "pipeline" / Path(name).name
    click.echo(f"[structcooker] {db_path.relative_to(REPO)}  op={op}  "
               f"schema={schema}  E={expansion}"
               + (f"  afterok={','.join(depends_on)}" if depends_on else ""))
    engine = _write_engine_yaml(db_path, wd)
    argv = [*PIPELINE_CLI, str(engine), "--workdir", str(wd),
            "--schema", schema, "--expansion", str(expansion), "--repo", str(REPO)]
    if depends_on:
        argv += ["--depends-on", ",".join(depends_on)]
    if condition:
        argv += ["--condition", condition]
    if dry_run:
        argv.append("--dry-run")
    proc = subprocess.run(argv, capture_output=True, text=True, check=False)  # noqa: S603
    click.echo(proc.stdout, nl=False)
    if proc.stderr:
        click.echo(proc.stderr, nl=False, err=True)
    m = _TERMINAL_RE.search(proc.stdout or "")
    job_id = (m.group(1).strip() or None) if m else None
    return proc.returncode, job_id


@cli.command("build")
@click.argument("name")
@click.option("--workdir", type=click.Path(path_type=Path), default=None,
              help="Scratch for tier item-lists + sbatch scripts "
                   "(default: logs/pipeline/<name>).")
@click.option("--dry-run", is_flag=True, help="Plan + write scripts, do not sbatch.")
@click.option("--depends-on", "depends_on", default="",
              help="Comma-separated upstream job ids to afterok-wait on.")
@click.option("--show", is_flag=True,
              help="Print the resolved schema/E + datacooker invocation and exit.")
def build(name: str, workdir: Path | None, dry_run: bool,
          depends_on: str, show: bool) -> None:
    """Plan + build (or rebuild) the DB named by its db/*.yaml (op auto-inferred)."""
    if show:
        db_path = _resolve_db(name)
        cfg = load_config(db_path)
        schema = _schema_of(cfg, db_path)
        wd = Path(workdir) if workdir else REPO / "logs" / "pipeline" / Path(name).name
        argv = [*PIPELINE_CLI, str(wd / "engine.yaml"), "--workdir", str(wd),
                "--schema", schema, "--expansion", str(schemas.expansion(schema)),
                "--repo", str(REPO)]
        click.echo("  (engine.yaml materialized on real build)\n  " + " ".join(argv))
        return
    deps = tuple(d for d in depends_on.split(",") if d.strip())
    rc, _ = _submit_pipeline(name, workdir, dry_run=dry_run, depends_on=deps)
    raise SystemExit(rc)


_INCREMENTAL_CONDITION = "datacooker.conditions.output_absent"


@cli.command("build-all")
@click.option("--manifest", type=click.Path(path_type=Path), default=None,
              help="DB dependency manifest (default: db/MANIFEST.yaml).")
@click.option("--workdir", type=click.Path(path_type=Path), default=None,
              help="Parent scratch dir (each DB gets a <workdir>/<name> subdir).")
@click.option("--dry-run", is_flag=True, help="Plan + write scripts, do not sbatch.")
@click.option("--force", is_flag=True,
              help="Rebuild every node, even ones already built (drops the condition).")
def build_all(manifest: Path | None, workdir: Path | None,
              dry_run: bool, force: bool) -> None:
    """Reproduce the whole DB set incrementally — a DAG of pipelines (each itself a DAG).

    Reads ``db/MANIFEST.yaml`` (``name: [upstream, ...]``), topologically sorts it, and
    submits each DB's pipeline with ``--depends-on`` the terminal job ids of its upstream
    DBs, so SLURM enforces the order. Each pipeline is gated by the datacooker
    ``condition`` ``output_absent``: a node already built no-ops (returns no job), and
    downstream nodes afterok-wait only on upstream actually submitted this run — so
    re-running fills in just what's missing. ``--force`` drops the condition (rebuild all).
    """
    man_path = Path(manifest) if manifest else DB_ROOT / "MANIFEST.yaml"
    if not man_path.exists():
        msg = f"no manifest at {man_path}"
        raise click.ClickException(msg)
    deps_map = {str(k): [str(d) for d in (v or [])]
                for k, v in OmegaConf.to_container(OmegaConf.load(man_path)).items()}
    order = _toposort(deps_map)
    click.echo(f"[build-all] {len(order)} DBs in dependency order:\n  "
               + " -> ".join(order))
    condition = None if force else _INCREMENTAL_CONDITION
    terminal: dict[str, str | None] = {}   # only nodes actually submitted this run
    built: list[str] = []
    skipped: list[str] = []
    for name in order:
        # afterok only on upstream nodes we actually submitted; a condition-skipped
        # upstream already has its output on disk, so no dependency is needed.
        upstream = tuple(j for d in deps_map[name] if (j := terminal.get(d)))
        rc, job_id = _submit_pipeline(name, (Path(workdir) / name) if workdir else None,
                                      dry_run=dry_run, depends_on=upstream,
                                      condition=condition)
        if rc != 0:
            msg = f"{name} pipeline submit failed (rc={rc}); aborting chain."
            raise click.ClickException(msg)
        if job_id:
            terminal[name] = job_id
            built.append(name)
        else:                       # condition returned False -> pipeline no-op'd
            skipped.append(name)
    click.echo(f"[build-all] submitted {len(built)}, skipped {len(skipped)} "
               f"(condition) ({len(order)} total).")


def _toposort(deps_map: dict[str, list[str]]) -> list[str]:
    """Topological sort (deterministic); raises on unknown deps or cycles."""
    for name, ds in deps_map.items():
        for d in ds:
            if d not in deps_map:
                msg = f"{name!r} depends on unknown {d!r}"
                raise click.ClickException(msg)
    order: list[str] = []
    remaining = dict(deps_map)
    while remaining:
        ready = sorted(n for n, ds in remaining.items() if all(d in order for d in ds))
        if not ready:
            msg = f"dependency cycle among: {sorted(remaining)}"
            raise click.ClickException(msg)
        order.extend(ready)
        for n in ready:
            del remaining[n]
    return order


def main() -> None:
    """Entry point for ``python -m structcooker`` / the ``structcooker`` script."""
    cli()


if __name__ == "__main__":
    main()
