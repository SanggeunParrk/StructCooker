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

import os
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

from structcooker import preflight, schemas

if TYPE_CHECKING:
    from collections.abc import Callable

REPO = Path(__file__).resolve().parents[2]          # StructCooker/
_DATA_ROOT_DEFAULT = "/data/shared/cssb_data"
_OUTPUT_ROOT_DEFAULT = "/data/shared/cssb_data/BioMol_clean"
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


def _op_of(cfg: dict) -> str:
    """Infer the op from a db config's shape.

    ``old_env_path`` -> rebuild; ``output_data_path`` -> a projection op
    (``extract`` when it reads a source DB, else ``materialize``); ``env_path``
    -> build. Projection ops produce a file (TSV/fasta), not an LMDB, so they run
    through ``datacooker.cli.workflow`` instead of the planning-first pipeline.
    """
    if cfg.get("old_env_path"):
        return "rebuild"
    if cfg.get("split_recipe") or cfg.get("split_recipe_path"):
        return "parallel"                       # split items, process each (e.g. hmmsearch)
    if cfg.get("output_data_path"):
        reads_db = any(cfg.get(k) for k in ("db_path", "extract_recipe", "extract_recipe_path"))
        return "extract" if reads_db else "materialize"
    if cfg.get("env_path"):
        return "build"
    msg = ("config needs old_env_path (rebuild) / env_path (build) / "
           "output_data_path (project) / split_recipe (parallel)")
    raise click.ClickException(msg)


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
        try:
            op = _op_of(cfg)
        except click.ClickException:
            op = "?"
        target = cfg.get("new_env_path") or cfg.get("env_path") or cfg.get("output_data_path") or "?"
        rows.append((str(f.relative_to(DB_ROOT)), op, str(cfg.get("schema", "-")), str(target)))
    if not rows:
        click.echo(f"(no db configs under {DB_ROOT})")
        return
    w = max(len(r[0]) for r in rows)
    for name, op, schema, target in rows:
        click.echo(f"{name:<{w}}  {op:<7}  schema={schema:<3}  {target}")


def _manifest_names(manifest: Path | None) -> list[str]:
    """Return the node names declared in the manifest (default db/MANIFEST.yaml)."""
    man_path = Path(manifest) if manifest else DB_ROOT / "MANIFEST.yaml"
    if not man_path.exists():
        msg = f"no manifest at {man_path}"
        raise click.ClickException(msg)
    raw = OmegaConf.to_container(OmegaConf.load(man_path))
    if not isinstance(raw, dict):
        msg = f"manifest {man_path} must be a mapping of name -> [deps]"
        raise click.ClickException(msg)
    return [str(k) for k in raw]


@cli.command("inspect")
@click.option("--manifest", type=click.Path(path_type=Path), default=None,
              help="DB dependency manifest (default: db/MANIFEST.yaml).")
@click.argument("name", required=False)
def inspect_cmd(manifest: Path | None, name: str | None) -> None:
    """Preflight: check every build-all node's external inputs are already in place.

    Reports the env roots, the external tools the recipes need, and per-node readiness
    -- which raw inputs (mmCIF, CCD, OpenFold, SabDab, ...) are on disk vs missing, and
    for the missing ones the ``structcooker download`` target that fetches them. Inputs
    built by another node (upstream) are noted, not flagged -- build-all builds those.
    """
    data_root = os.environ.get("DATA_ROOT", _DATA_ROOT_DEFAULT)
    names = [name] if name else _manifest_names(manifest)
    configs = {n: load_config(_resolve_db(n)) for n in names}
    reports = preflight.inspect(configs)

    env = preflight.check_env()
    click.echo("[env]")
    click.echo(f"  DATA_ROOT   = {env['DATA_ROOT'] or f'(unset -> {data_root})'}")
    click.echo(f"  OUTPUT_ROOT = {env['OUTPUT_ROOT'] or '(unset -> BioMol_clean default)'}")
    click.echo(f"  SEQID_SEED  = {env['SEQID_SEED'] or '(unset -> fresh seq_id space; set to match production)'}")

    click.echo("[tools]")
    for tool, found, optional in preflight.check_tools():
        mark = "OK     " if found else ("--     " if optional else "MISSING")
        note = "  (optional, licensed)" if optional and not found else ""
        click.echo(f"  {mark} {tool}{note}")

    click.echo("[nodes]")
    for r in sorted(reports, key=lambda r: (r.ready, r.name)):
        detail = f"ext {len(r.present)} ok"
        if r.missing:
            detail += f", {len(r.missing)} MISSING"
        if r.upstream:
            detail += f"; upstream {len(set(r.upstream))}"
        click.echo(f"  {'READY  ' if r.ready else 'BLOCKED'}  {r.name:<34} {detail}")

    _report_manual_fixes(data_root)

    missing = preflight.missing_externals(reports, data_root)
    ready = sum(1 for r in reports if r.ready)
    click.echo(f"[summary] {ready}/{len(reports)} nodes ready; "
               f"{len(missing)} distinct external inputs missing")
    if missing:
        click.echo("[missing external inputs -- provide before build-all]")
        for path, hint in sorted(missing.items()):
            action = f"structcooker download {hint}" if hint else "provide (tool output / lab-supplied)"
            click.echo(f"  {path}\n      -> {action}")


_MANUAL_FIXES_PATH = REPO / "db" / "pdb" / "manual_cif_fixes.txt"
_MMCIF_SUBPATH = Path("mmcif_files_latest") / "mmcif_files"


def _mmcif_dir(data_root: str | Path) -> Path:
    return Path(data_root) / _MMCIF_SUBPATH


def _report_manual_fixes(data_root: str) -> None:
    """Print how many of the manual mmCIF substitutions are in place (for inspect)."""
    from structcooker import cif_fixes

    ids = cif_fixes.load_fix_ids(_MANUAL_FIXES_PATH)
    applied = cif_fixes.applied_ids(_mmcif_dir(data_root))
    n_applied = len(applied & set(ids))
    click.echo("[manual cif fixes]")
    if n_applied == len(ids):
        click.echo(f"  OK      all {len(ids)} substitutions applied")
    else:
        click.echo(f"  PENDING {n_applied}/{len(ids)} applied "
                   f"-> structcooker fix-cif --source <corrected-cif-dir> "
                   f"(see docs/manual-cif-fixes.md)")


@cli.command("fix-cif")
@click.option("--source", "source", type=click.Path(exists=True, path_type=Path),
              required=True, help="Corrected-cif snapshot (divided <id[1:3]>/ or flat; "
                                  ".cif or .cif.gz). Provided input; see the docs.")
@click.option("--dry-run", is_flag=True, help="Report what would be substituted, write nothing.")
def fix_cif_cmd(source: Path, dry_run: bool) -> None:
    """Overlay the manually-fixed mmCIFs the pdb/cif build needs, from a corrected snapshot.

    A set of PDB entries (``db/pdb/manual_cif_fixes.txt``) error out from the current
    wwPDB mmCIF -- mostly NMR ensembles with per-model-renumbered ligands -- so the
    production build substituted an older, known-good cif for each before ingest. This
    ports that step: it copies each id's corrected cif from ``--source`` into the mmCIF
    input dir. Run it before ``structcooker build pdb/cif``. See docs/manual-cif-fixes.md.
    """
    from structcooker import cif_fixes

    data_root = os.environ.get("DATA_ROOT", _DATA_ROOT_DEFAULT)
    mmcif_dir = _mmcif_dir(data_root)
    if not mmcif_dir.is_dir():
        msg = f"mmCIF input dir not found: {mmcif_dir} (download mmcif first)."
        raise click.ClickException(msg)
    ids = cif_fixes.load_fix_ids(_MANUAL_FIXES_PATH)
    applied, missing = cif_fixes.apply_fixes(ids, source, mmcif_dir, dry_run=dry_run)
    verb = "would substitute" if dry_run else "substituted"
    click.echo(f"{verb} {len(applied)}/{len(ids)} mmCIFs from {source} -> {mmcif_dir}")
    if missing:
        click.echo(f"  {len(missing)} not found in source (provide them): "
                   f"{', '.join(missing)}")


@cli.command("download")
@click.argument("target", type=click.Choice(["ccd", "sabdab", "mmcif", "openfold"]))
@click.option("--yes", is_flag=True, help="Confirm large downloads (mmcif is ~90 GB+).")
def download_cmd(target: str, yes: bool) -> None:
    """Fetch a raw external input into the DATA_ROOT / OUTPUT_ROOT layout.

    ``ccd`` / ``sabdab`` are small and download directly; ``mmcif`` (~90 GB+) needs
    ``--yes``; ``openfold`` (TB-scale, portal-hosted) only prints instructions. The
    ``seq_id_map`` seed is NOT here -- it is the Hugging Face dataset ``biomol/seq-id-map``
    (see the README); omit it for a fresh id space.
    """
    from structcooker import downloads

    data_root = Path(os.environ.get("DATA_ROOT", _DATA_ROOT_DEFAULT))
    try:
        if target == "ccd":
            dest = downloads.download_ccd(data_root)
        elif target == "sabdab":
            dest = downloads.download_sabdab(data_root)
        elif target == "mmcif":
            dest = downloads.download_mmcif(data_root, confirmed=yes)
        else:
            dest = downloads.download_openfold(data_root, confirmed=yes)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"downloaded {target} -> {dest}")


_TERMINAL_RE = re.compile(r"^PIPELINE_TERMINAL_JOB=(.*)$", re.MULTILINE)


def _submit_pipeline(
    name: str,
    workdir: Path | None,
    *,
    dry_run: bool,
    depends_on: tuple[str, ...] = (),
    cfg: dict | None = None,
) -> str | None:
    """Materialize the engine config for ``name`` and submit its pipeline.

    Returns the terminal job id (the pipeline's final index stage) so ``build-all`` can
    afterok-chain the next config onto it. Raises ``ClickException`` on any failure -- a
    non-zero exit, or a clean exit with no terminal id (output format drift) -- so a
    caller never mistakes a failed submit for a skip. Pass ``cfg`` to reuse an
    already-loaded config.
    """
    db_path = _resolve_db(name)
    if cfg is None:
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
    if proc.returncode != 0:
        msg = f"{name} pipeline submit failed (rc={proc.returncode})."
        raise click.ClickException(msg)
    m = _TERMINAL_RE.search(proc.stdout or "")
    job_id = (m.group(1).strip() or None) if m else None
    if job_id is None:
        # A clean pipeline always emits PIPELINE_TERMINAL_JOB. Missing it means the
        # output format drifted -- fail loud rather than let build-all treat this as
        # "skipped" and submit dependents with no afterok (silent out-of-order run).
        msg = f"{name}: pipeline exited 0 but emitted no terminal job id (format drift?)"
        raise click.ClickException(msg)
    return job_id


_WORKFLOW_CMD = {"materialize": "run", "extract": "extract-lmdb", "parallel": "parallel-run"}
# Single-node resource sizing per projection op. Config keys can't carry hints (the
# workflow CLI re-reads the yaml and rejects unknown keys), so size by op: parallel
# fans out (full node); extract runs joblib over a whole DB; materialize is mostly
# single-process (a couple stream a multi-GB fasta, hence the mem headroom).
_WORKFLOW_RESOURCES = {"parallel": (490, 112), "extract": (200, 32), "materialize": (128, 8)}
WORKFLOW_CLI = [sys.executable, "-u", "-m", "datacooker.cli.workflow"]


def _submit_workflow(
    name: str,
    workdir: Path | None,
    *,
    op: str,
    dry_run: bool,
    depends_on: tuple[str, ...] = (),
    cfg: dict | None = None,
) -> str | None:
    """Submit a projection op (materialize / extract / parallel) as a single SLURM job.

    These configs produce a file (a TSV/fasta projection), not an LMDB, so they run
    through ``datacooker.cli.workflow`` (single-process ``run`` / ``extract-lmdb`` /
    ``parallel-run``), not the planning-first tier pipeline. Returns the ``run_once`` job
    id (afterok-chainable, so build-all threads these nodes into the DAG like any other);
    ``run_once`` raises on submit failure, so a returned id is always real. Pass ``cfg``
    to reuse an already-loaded config.
    """
    from datacooker.executors.slurm import SlurmExecutor

    db_path = _resolve_db(name)
    if cfg is None:
        cfg = load_config(db_path)
    target = cfg.get("output_data_path")
    wd = Path(workdir) if workdir else REPO / "logs" / "workflow" / Path(name).name
    click.echo(f"[structcooker] {db_path.relative_to(REPO)}  op={op}  -> {target}"
               + (f"  afterok={','.join(depends_on)}" if depends_on else ""))
    argv = [*WORKFLOW_CLI, _WORKFLOW_CMD[op], str(db_path)]
    mem_gb, cores = _WORKFLOW_RESOURCES[op]
    execu = SlurmExecutor(workdir=wd, repo=REPO, submit=not dry_run)
    handle = execu.run_once(name=Path(name).name, argv=argv,
                            mem_gb=mem_gb, cores=cores, depends_on=depends_on)
    return handle.job_id


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
    """Plan + build the DB named by its db/*.yaml (op auto-inferred).

    build / rebuild run the planning-first LMDB pipeline; materialize / extract
    (projection ops that write a TSV/fasta, not an LMDB) run through
    ``datacooker.cli.workflow`` as a single job.
    """
    db_path = _resolve_db(name)
    cfg = load_config(db_path)
    op = _op_of(cfg)
    deps = tuple(d for d in depends_on.split(",") if d.strip())
    if op in _WORKFLOW_CMD:
        if show:
            click.echo(f"  op={op} -> {' '.join([*WORKFLOW_CLI, _WORKFLOW_CMD[op], str(db_path)])}")
            return
        _submit_workflow(name, workdir, op=op, dry_run=dry_run, depends_on=deps, cfg=cfg)
        return
    if show:
        schema = _schema_of(cfg, db_path)
        wd = Path(workdir) if workdir else REPO / "logs" / "pipeline" / Path(name).name
        argv = [*PIPELINE_CLI, str(wd / "engine.yaml"), "--workdir", str(wd),
                "--schema", schema, "--expansion", str(schemas.expansion(schema)),
                "--repo", str(REPO)]
        click.echo("  (engine.yaml materialized on real build)\n  " + " ".join(argv))
        return
    _submit_pipeline(name, workdir, dry_run=dry_run, depends_on=deps, cfg=cfg)


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
            cfg = load_config(_resolve_db(name))
            if not force and not output_absent(cfg):
                click.echo(f"[build-all] {name}: output present -> skip")
                return None
            depends_on = tuple(j for j in upstream if j)
            wd = (Path(workdir) / name) if workdir else None
            # dispatch by op: projection ops (materialize/extract/parallel) go through
            # the workflow runner, build/rebuild through the planning-first pipeline.
            # Both submitters raise on failure, so a returned id is always real.
            op = _op_of(cfg)
            if op in _WORKFLOW_CMD:
                return _submit_workflow(name, wd, op=op, dry_run=dry_run,
                                        depends_on=depends_on, cfg=cfg)
            return _submit_pipeline(name, wd, dry_run=dry_run,
                                    depends_on=depends_on, cfg=cfg)
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
