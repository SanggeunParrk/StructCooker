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
from typing import TYPE_CHECKING

import click
from datacooker.api import execute
from datacooker.conditions import output_absent
from datacooker.config import load_config
from datacooker.errors import StepExecutionError
from datacooker.recipe import RecipeBook
from omegaconf import OmegaConf

from structcooker import schemas

if TYPE_CHECKING:
    from collections.abc import Callable

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
    """Build BioMol DBs declaratively (planning-first)."""


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
) -> tuple[int, str | None]:
    """Materialize the engine config for ``name`` and submit its pipeline.

    Returns ``(returncode, terminal_job_id)`` -- the terminal id (the pipeline's
    final index stage) lets ``build-all`` afterok-chain the next config onto it.
    Streams the pipeline output through while capturing it.
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


def _meta_book(deps_map: dict[str, list[str]], *,
               workdir: Path | None, dry_run: bool, force: bool) -> RecipeBook:
    """Turn the MANIFEST dep-map into a datacooker ``RecipeBook`` — one step per DB.

    Each DB is a step ``<name> <- submit(*upstream_job_ids)``: the instruction submits
    the DB's SLURM pipeline afterok-chained on whatever upstream ids it receives and
    returns this DB's terminal job id — or ``None`` when the DB is already built. That
    ``None`` **is** the "if": datacooker's only conditional is an instruction returning
    None, and the engine propagates it downstream (a skipped upstream arrives as a None
    arg, i.e. no afterok). ``--force`` drops the skip. Topological order, cycle detection,
    the unknown-dependency check, and the afterok data-flow are all the engine's job
    (``execution_order`` / ``resolve`` / ``validate``), not ours — build-all is just a
    recipe, cooked in-process by :func:`datacooker.api.execute`.

    Typing: an **input** dep is ``str`` (a real upstream job id is mandatory, so an
    unknown MANIFEST dep — not a declared target — trips ``MissingDependencyError``),
    while an **output** is ``str | None`` (this DB's id, or None when skipped). A skipped
    upstream still flows because ``resolve`` reads it from the context by name before any
    type check.
    """
    book = RecipeBook()

    def make_submit(name: str) -> Callable[..., str | None]:
        def submit(*upstream: str | None) -> str | None:
            if not force and not output_absent(load_config(_resolve_db(name))):
                click.echo(f"[build-all] {name}: output present -> skip")
                return None
            depends_on = tuple(j for j in upstream if j)
            wd = (Path(workdir) / name) if workdir else None
            rc, job_id = _submit_pipeline(name, wd, dry_run=dry_run, depends_on=depends_on)
            if rc != 0:
                msg = f"{name} pipeline submit failed (rc={rc})."
                raise click.ClickException(msg)
            return job_id
        submit.__name__ = f"submit[{name}]"
        return submit

    for name, deps in deps_map.items():
        book.step(
            outputs=(name, str | None),
            instruction=make_submit(name),
            args=[(d, str) for d in deps],
        )
    return book


@cli.command("build-all")
@click.option("--manifest", type=click.Path(path_type=Path), default=None,
              help="DB dependency manifest (default: db/MANIFEST.yaml).")
@click.option("--workdir", type=click.Path(path_type=Path), default=None,
              help="Parent scratch dir (each DB gets a <workdir>/<name> subdir).")
@click.option("--dry-run", is_flag=True, help="Plan + write scripts, do not sbatch.")
@click.option("--force", is_flag=True,
              help="Rebuild every node, even ones already built (drops the skip).")
def build_all(manifest: Path | None, workdir: Path | None,
              dry_run: bool, force: bool) -> None:
    """Reproduce the whole DB set incrementally — a DAG of pipelines (each itself a DAG).

    ``db/MANIFEST.yaml`` (``name: [upstream, ...]``) becomes a datacooker ``RecipeBook``
    (one step per DB), which the engine cooks in-process: it topologically orders the
    nodes, threads each DB's terminal job id into its dependents as ``--depends-on`` (so
    SLURM enforces the order), and skips any DB already built — a step whose instruction
    returns ``None``. Re-running therefore fills in only what's missing. ``--force``
    rebuilds every node. No hand-rolled toposort/afterok/condition: it's all data-flow.
    """
    man_path = Path(manifest) if manifest else DB_ROOT / "MANIFEST.yaml"
    if not man_path.exists():
        msg = f"no manifest at {man_path}"
        raise click.ClickException(msg)
    raw_manifest = OmegaConf.to_container(OmegaConf.load(man_path))
    if not isinstance(raw_manifest, dict):
        msg = f"manifest {man_path} must be a mapping of name -> [deps]"
        raise click.ClickException(msg)
    deps_map = {str(k): [str(d) for d in (v or [])] for k, v in raw_manifest.items()}
    book = _meta_book(deps_map, workdir=workdir, dry_run=dry_run, force=force)
    order = [r.target_names[0] for r in book.execution_order()]
    click.echo(f"[build-all] {len(order)} DBs in dependency order:\n  "
               + " -> ".join(order))
    try:
        results = execute(book, {})            # cook the meta-graph (engine orders + runs)
    except StepExecutionError as exc:
        raise click.ClickException(str(exc.original_exception or exc)) from exc
    built = [n for n, job_id in results.items() if job_id]
    skipped = [n for n in results if n not in built]
    click.echo(f"[build-all] submitted {len(built)}, skipped {len(skipped)} "
               f"(already built) ({len(order)} total).")


def main() -> None:
    """Entry point for ``python -m structcooker`` / the ``structcooker`` script."""
    cli()


if __name__ == "__main__":
    main()
