"""structcooker CLI — declarative, planning-first BioMol DB builds.

One BioMol DB = one ``db/**/<name>.yaml`` = a datacooker *engine config*
(recipe / reader / writer + source + target) plus a single ``schema:`` tag.

``structcooker build <name>`` hands the config to ``datacooker pipeline``, which
sizes the items, snake-balances them across nodes, and submits shard arrays → merge →
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

import click
from datacooker.config import load_config
from omegaconf import OmegaConf

from structcooker import preflight, schemas, workflow_ops
from structcooker.paths import biomol_root, mmcif_root, seq_db_root

REPO = Path(__file__).resolve().parents[2]          # StructCooker/
DB_ROOT = REPO / "db"
PIPELINE_CLI = [sys.executable, "-u", "-m", "datacooker.cli.lmdb", "pipeline"]
_SCHEMA_LINE = re.compile(r"^\s*schema\s*:.*$", re.MULTILINE)
# Inputs a release build may fail on, per stage, each with its reason and SHA256: the PDB
# entries the parser cannot map and the OpenFold alignment.npz files that are empty or
# corrupt in the release (db/release-exceptions.json, a superset of pdb-exceptions.json).
RELEASE_EXCEPTIONS = DB_ROOT / "release-exceptions.json"


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
    """:func:`structcooker.workflow_ops.op_of`, as a CLI error."""
    try:
        return workflow_ops.op_of(cfg)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


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
@click.option("--strict", is_flag=True,
              help="Exit nonzero for missing inputs or invalid manifest dependencies.")
@click.argument("name", required=False)
def inspect_cmd(manifest: Path | None, strict: bool, name: str | None) -> None:
    """Preflight: check every build-all node's external inputs are already in place.

    Reports the env roots, the external tools the recipes need, and per-node readiness
    -- which raw inputs (mmCIF, CCD, OpenFold, SabDab, ...) are on disk vs missing, and
    for the missing ones the ``structcooker download`` target that fetches them. Inputs
    built by another node (upstream) are noted, not flagged -- build-all builds those.
    """
    root = biomol_root()
    names = [name] if name else _manifest_names(manifest)
    configs = {n: load_config(_resolve_db(n)) for n in names}
    reports = preflight.inspect(configs)
    dependencies = {} if name else OmegaConf.to_container(
        OmegaConf.load(manifest or DB_ROOT / "MANIFEST.yaml"), resolve=True,
    )
    deps_map = {str(k): [str(d) for d in (v or [])] for k, v in dependencies.items()} if isinstance(dependencies, dict) else {}
    errors = preflight.dependency_errors(configs, deps_map)

    env = preflight.check_env()
    click.echo("[env]")
    click.echo(f"  BIOMOL_ROOT = {env['BIOMOL_ROOT'] or f'(unset -> {root})'}")
    click.echo(f"  SEQ_DB_ROOT = {env['SEQ_DB_ROOT'] or f'(unset -> {seq_db_root()})'}")
    seed = root / "materials" / "reference" / "seq_id_map_seed.tsv"
    click.echo(f"  seq_id seed = {'present -> match production' if seed.exists() else 'absent -> fresh id space'}  ({seed})")

    click.echo("[tools available for optional workflows; not all are needed by every node]")
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

    missing = preflight.missing_externals(reports, str(root))
    ready = sum(1 for r in reports if r.ready)
    click.echo(f"[summary] {ready}/{len(reports)} nodes ready; "
               f"{len(missing)} distinct external inputs missing")
    if missing:
        click.echo("[missing external inputs -- provide before build-all]")
        for path, hint in sorted(missing.items()):
            action = f"structcooker download {hint}" if hint else "provide (tool output / lab-supplied)"
            click.echo(f"  {path}\n      -> {action}")
    for error in errors:
        click.echo(f"[invalid dependency] {error}")
    if strict and (missing or errors):
        msg = "Preflight failed; resolve the inputs/dependencies listed above."
        raise click.ClickException(msg)


_MANUAL_FIXES_PATH = REPO / "db" / "pdb" / "manual_cif_fixes.txt"


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

    mmcif_dir = mmcif_root(biomol_root())
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
    """Fetch a raw external input into BIOMOL_ROOT/materials/raw.

    ``ccd`` / ``sabdab`` are small and download directly; ``mmcif`` (~90 GB+) needs
    ``--yes``; ``openfold`` (TB-scale, portal-hosted) only prints instructions. The
    ``seq_id_map`` seed is NOT here -- it is the Hugging Face dataset ``biomol/seq-id-map``
    (see the README); omit it for a fresh id space.
    """
    from structcooker import downloads

    root = biomol_root()
    try:
        if target == "ccd":
            dest = downloads.download_ccd(root)
        elif target == "sabdab":
            dest = downloads.download_sabdab(root)
        elif target == "mmcif":
            dest = downloads.download_mmcif(root, confirmed=yes)
        else:
            dest = downloads.download_openfold(root, confirmed=yes)
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


def _submit_workflow(
    name: str,
    workdir: Path | None,
    *,
    op: str,
    dry_run: bool,
    depends_on: tuple[str, ...] = (),
    cfg: dict | None = None,
) -> str | None:
    """Submit a workflow op (materialize / extract / parallel) as SLURM job(s).

    These configs produce files (a TSV/fasta projection, or a parallel op's output dir),
    not an LMDB, so they run through ``datacooker.cli.workflow`` (``run`` /
    ``extract-lmdb`` / ``parallel-run``), not the planning-first tier pipeline. Returns an
    afterok-chainable job id; submission raises on failure, so a returned id is always
    real. Pass ``cfg`` to reuse an already-loaded config.
    """
    from datacooker.executors.slurm import SlurmExecutor

    db_path = _resolve_db(name)
    if cfg is None:
        cfg = load_config(db_path)
    target = cfg.get("output_data_path")
    wd = Path(workdir) if workdir else REPO / "logs" / "workflow" / Path(name).name
    click.echo(f"[structcooker] {db_path.relative_to(REPO)}  op={op}  -> {target}"
               + (f"  afterok={','.join(depends_on)}" if depends_on else ""))
    execu = SlurmExecutor(workdir=wd, repo=REPO, submit=not dry_run)
    return workflow_ops.submit(execu, Path(name).name, op, db_path, cfg, depends_on)


@cli.command("patch")
@click.argument("name")
@click.option("--keys", "keys", type=click.Path(exists=True, path_type=Path), default=None,
              help="Keys to rebuild, one per line (a rebuild op, e.g. *_cif_attached).")
@click.option("--files", "files", type=click.Path(exists=True, path_type=Path), default=None,
              help="Input files to build, one per line (a build op, e.g. *_cif).")
@click.option("--workdir", type=click.Path(path_type=Path), default=None,
              help="Scratch for the patch scripts (default: logs/patch/<name>/<stamp>).")
@click.option("--nodes", type=int, default=1, show_default=True,
              help="Build the patch as an N-task array, one whole node each.")
@click.option("--dry-run", is_flag=True, help="Write the scripts, do not sbatch.")
def patch(name: str, keys: Path | None, files: Path | None, workdir: Path | None,
          nodes: int, dry_run: bool) -> None:
    """Rebuild only some records of an already published DB and swap them in.

    For records that were built wrong or are missing, when rebuilding the whole DB is not
    worth it (e.g. 111 of 16M OFD long entries built from partially downloaded models). The
    records are built into a separate LMDB next to the target, the target's copies are
    replaced (a shard collection gains it as one more shard), and the index is rewritten.
    A full `build` resumes only within its own run, and skips keys that already exist --
    neither fixes a record that exists but is wrong.
    """
    from datetime import datetime

    from datacooker.executors.slurm import SlurmExecutor

    if (keys is None) == (files is None):
        msg = "Give exactly one of --keys (rebuild op) or --files (build op)."
        raise click.ClickException(msg)
    db_path = _resolve_db(name)
    cfg = load_config(db_path)
    op = "rebuild" if cfg.get("old_env_path") else "build"
    if (op == "rebuild") != (keys is not None):
        msg = f"{name} is a {op} op: use {'--keys' if op == 'rebuild' else '--files'}."
        raise click.ClickException(msg)
    target = Path(str(cfg.get("new_env_path") or cfg.get("env_path")))
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # noqa: DTZ005
    wd = Path(workdir) if workdir else REPO / "logs" / "patch" / Path(name).name / stamp
    engine = _write_engine_yaml(db_path, wd)
    out = target.parent / f".{target.stem}.build" / f"patch_{stamp}" / f"{target.stem}_patch.lmdb"
    lmdb_cli = [sys.executable, "-u", "-m", "datacooker.cli.lmdb"]
    build = [*lmdb_cli, op, str(engine), "--key-list" if keys else "--file-list", str(keys or files),
             "--output", str(out), "--n-jobs", "112"]
    execu = SlurmExecutor(workdir=wd, repo=REPO, submit=not dry_run)
    if nodes > 1:
        # Each task builds its key-hash slice into <out stem>_shard<i>; all are applied at once.
        b = execu.run_tier_shards(name=f"{target.stem}_patch_build", n_shards=nodes, mem_gb=490, cores=112,
                                  argv_for_shard=lambda i: [*build, "--shard-idx", str(i), "--n-shards", str(nodes)])
        outs = [out.with_name(f"{out.stem}_shard{i}{out.suffix}") for i in range(nodes)]
    else:
        b = execu.run_once(name=f"{target.stem}_patch_build", argv=build, mem_gb=490, cores=112)
        outs = [out]
    apply = [*lmdb_cli, "patch", "--target", str(target)]
    for o in outs:
        apply += ["--patch", str(o)]
    a = execu.run_once(name=f"{target.stem}_patch_apply", argv=apply, mem_gb=64, cores=8,
                       depends_on=(b.job_id,))
    # The apply step also rewrites the index from the old rows plus the patch's sizes;
    # a full re-index of a multi-TB DB reads every value (hours) for no new information.
    click.echo(f"[patch] {name}: build {b.job_id} -> apply+index {a.job_id}  ({out})")
    click.echo(f"PIPELINE_TERMINAL_JOB={a.job_id}")


@cli.command("build")
@click.argument("name")
@click.option("--workdir", type=click.Path(path_type=Path), default=None,
              help="Scratch for shard item-lists + sbatch scripts "
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
    if op in workflow_ops.WORKFLOW_CMD:
        if show:
            argv = [*workflow_ops.WORKFLOW_CLI, workflow_ops.WORKFLOW_CMD[op], str(db_path)]
            click.echo(f"  op={op} -> {' '.join(argv)}")
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
    """Submit or resume verified release stages through finite native SLURM jobs."""
    from structcooker.release import start

    man_path = Path(manifest) if manifest else DB_ROOT / "MANIFEST.yaml"
    if dry_run:
        raw = OmegaConf.to_container(OmegaConf.load(man_path))
        if not isinstance(raw, dict):
            msg = "Manifest must be a mapping"
            raise click.ClickException(msg)
        dependencies = {str(n): [str(p) for p in parents] for n, parents in raw.items()}
        configs = {name: load_config(_resolve_db(name)) for name in dependencies}
        errors = preflight.dependency_errors(configs, dependencies)
        if errors:
            raise click.ClickException("\n".join(errors))
        click.echo(f"Validated {len(raw)} release stages; no jobs submitted")
        return
    if force:
        msg = "Use a fresh BIOMOL_ROOT and run directory for a new release; --force cannot bypass ownership"
        raise click.ClickException(msg)
    if workdir is None:
        msg = "--workdir is required for persistent, verified release state"
        raise click.ClickException(msg)
    try:
        start(REPO, man_path, workdir.resolve(), RELEASE_EXCEPTIONS)
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Release jobs submitted/resumed; inspect {workdir / 'summary.json'} for verified completion")


@cli.command("release-check")
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_all.yaml",
              show_default=True)
def release_check(manifest: Path) -> None:
    """Check a release tree is complete and self-contained (SLURM job).

    Every stage's output at its configured path under BIOMOL_ROOT, shards inside their own
    collection, entry counts agreeing with .meta.json and .index.tsv, and a sample of
    records decoding cleanly. Report: BIOMOL_ROOT/release_check.json.
    """
    from datacooker.executors.slurm import SlurmExecutor

    out_root = biomol_root()
    execu = SlurmExecutor(workdir=REPO / "logs" / "release_check", repo=REPO)
    job = execu.run_once(name="release_check", argv=[sys.executable, "-u", "-m", "structcooker.release_check",
                                                       str(manifest.resolve()), str(out_root / "release_check.json")],
                         mem_gb=256, cores=16).job_id
    click.echo(f"release-check job {job}; report -> {out_root / 'release_check.json'}")


@cli.command("release-docs")
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_all.yaml",
              show_default=True)
def release_docs(manifest: Path) -> None:
    """Write BIOMOL_ROOT/docs: README, provenance, databases, what-to-copy scenarios, check.

    The database list, the scenario sizes and the check summary are generated from the tree
    (configs, .meta.json, data.mdb sizes, release_check.json); rerun after a build or a move.
    """
    from structcooker import release_docs as docs

    for path in docs.write(biomol_root(), manifest.resolve()):
        click.echo(f"wrote {path}")


@cli.group("e2e-sample")
def e2e_sample() -> None:
    """Sampled end-to-end check: the release DAG on a few hundred inputs.

    prepare (SLURM job) -> run (build-all on the sample) -> status -> verify (SLURM job).
    See structcooker.sample_e2e and docs/e2e-build.md.
    """


def _sample_job(root: Path, name: str, argv: list[str], mem_gb: int) -> str:
    from datacooker.executors.slurm import SlurmExecutor

    execu = SlurmExecutor(workdir=root / "jobs" / name, repo=REPO)
    return execu.run_once(name=f"e2e_{name}", argv=[sys.executable, "-u", "-m", "structcooker.sample_e2e", *argv],
                          mem_gb=mem_gb, cores=8).job_id


@e2e_sample.command("prepare")
@click.argument("root", type=click.Path(path_type=Path))
@click.option("--count", "counts", multiple=True, metavar="SOURCE=N",
              help="Sample size per source, e.g. --count pdb=50 (sources: see sample_e2e.DEFAULT_COUNTS).")
def e2e_sample_prepare(root: Path, counts: tuple[str, ...]) -> None:
    """Pick the sample and build ROOT, a small BIOMOL_ROOT."""
    import json

    from structcooker.sample_e2e import DEFAULT_COUNTS

    sizes = {}
    for item in counts:
        source, _, n = item.partition("=")
        if source not in DEFAULT_COUNTS or not n.isdigit():
            msg = f"--count {item}: expected SOURCE=N with SOURCE in {sorted(DEFAULT_COUNTS)}"
            raise click.ClickException(msg)
        sizes[source] = int(n)
    root = root.resolve()
    if (root / "sample.json").exists():
        msg = f"{root} already holds a sample; use a new root"
        raise click.ClickException(msg)
    root.mkdir(parents=True, exist_ok=True)
    job = _sample_job(root, "prepare", ["prepare", str(root), json.dumps(sizes)], 128)
    click.echo(f"prepare job {job}; then: structcooker e2e-sample run {root}")


@e2e_sample.command("run")
@click.argument("root", type=click.Path(exists=True, path_type=Path))
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_all.yaml",
              show_default=True)
def e2e_sample_run(root: Path, manifest: Path) -> None:
    """build-all on the sample (ROOT/run), with every path confined to ROOT."""
    from structcooker.release import start
    from structcooker.sample_e2e import escaping_paths, sample_environment

    root = root.resolve()
    if not (root / "sample.json").exists():
        msg = f"{root} has no sample.json; run `structcooker e2e-sample prepare` first (and let it finish)"
        raise click.ClickException(msg)
    env = sample_environment(root)
    os.environ.clear()
    os.environ.update(env)
    names = _manifest_names(manifest)
    bad = escaping_paths({n: load_config(_resolve_db(n)) for n in names}, root)
    if bad:
        raise click.ClickException("paths outside the sample root:\n  " + "\n  ".join(bad))
    try:
        start(REPO, manifest, root / "run", RELEASE_EXCEPTIONS)
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"sample build submitted/resumed; structcooker e2e-sample status {root}")


_FAILED_JOB_STATES = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL", "DEADLINE"}


@e2e_sample.command("status")
@click.argument("root", type=click.Path(exists=True, path_type=Path))
def e2e_sample_status(root: Path) -> None:
    """Stage phases of the sampled build (no scanning, no jobs)."""
    import json
    from collections import Counter

    summary_path = root / "run" / "summary.json"
    if not summary_path.exists():
        msg = f"no build under {root / 'run'} yet"
        raise click.ClickException(msg)
    phases = json.loads(summary_path.read_text())["nodes"]
    states = {name: json.loads((root / "run" / "nodes" / name / "state.json").read_text())
              for name, phase in phases.items() if phase not in {"complete", "waiting"}}
    for name, state in states.items():
        # A stage whose work job failed keeps a verifier that can never start; the phase
        # stays "submitted" until `run` reconciles it.
        if phases[name] == "submitted" and state.get("attempt"):
            jobs = [p.read_text().strip() for p in Path(state["attempt"]).rglob("*.job-id")]
            shown = subprocess.run(["sacct", "-n", "-P", "-o", "State", "-j", ",".join(jobs)],  # noqa: S603, S607
                                   capture_output=True, text=True, check=False).stdout if jobs else ""
            if any(s.split()[0] in _FAILED_JOB_STATES for s in shown.split("\n") if s):
                phases[name] = "stalled"
                state["error"] = "a job failed; see the attempt logs, fix, then rerun `e2e-sample run`"
    click.echo("  ".join(f"{p}={n}" for p, n in sorted(Counter(phases.values()).items())))
    for name, state in sorted(states.items()):
        click.echo(f"  {phases[name]:<20} {name}  {(state.get('error') or '')[:160]}")


@e2e_sample.command("verify")
@click.argument("root", type=click.Path(exists=True, path_type=Path))
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_all.yaml",
              show_default=True)
def e2e_sample_verify(root: Path, manifest: Path) -> None:
    """Compare every sampled output with the release (SLURM job; report in ROOT/verify.json)."""
    root = root.resolve()
    job = _sample_job(root, "verify", ["verify", str(root), str(root / "run"), str(manifest.resolve())], 128)
    click.echo(f"verify job {job}; report -> {root / 'verify.json'}")


def main() -> None:
    """Entry point for ``python -m structcooker`` / the ``structcooker`` script."""
    cli()


@cli.command("pdb-build")
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_pdb.yaml")
@click.option("--run-dir", type=click.Path(path_type=Path), required=True)
@click.option("--policy", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "pdb-exceptions.json")
@click.option("--rebuild-untracked", is_flag=True, help="Explicitly allow rebuilding existing outputs without receipts.")
def pdb_build(manifest: Path, run_dir: Path, policy: Path, rebuild_untracked: bool) -> None:
    """Build PDB outputs with content identities, verified receipts, and safe restart."""
    from structcooker.production import execute_manifest

    try:
        results = execute_manifest(REPO, manifest, run_dir.resolve(), policy,
                                   rebuild_untracked=rebuild_untracked)
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Verified {len(results)} nodes: {results}")


@cli.command("pdb-status")
@click.option("--run-dir", type=click.Path(exists=True, path_type=Path), required=True)
def pdb_status(run_dir: Path) -> None:
    """Read durable node states without scanning datasets or submitting any jobs."""
    import json

    for path in sorted((run_dir / "nodes").rglob("state.json")):
        state = json.loads(path.read_text())
        name = path.parent.relative_to(run_dir / "nodes")
        click.echo(f"{name}: {state.get('phase', 'unknown')} job={state.get('job', '-')} "
                   f"error={state.get('error') or '-'}")


@cli.command("release-build")
@click.option("--manifest", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "MANIFEST_distillation.yaml")
@click.option("--run-dir", type=click.Path(path_type=Path), required=True)
@click.option("--policy", type=click.Path(exists=True, path_type=Path), default=DB_ROOT / "distillation-exceptions.json")
def release_build(manifest: Path, run_dir: Path, policy: Path) -> None:
    """Submit or resume a release DAG without a resident scheduler controller."""
    from structcooker.release import start

    try:
        start(REPO, manifest, run_dir.resolve(), policy)
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Submitted/resumed release; status: {run_dir / 'summary.json'}")


@cli.group("foldbench")
def foldbench_group() -> None:
    """FoldBench inference DBs: key-list subsets of the full-release LMDBs (see docs/foldbench.md)."""


@foldbench_group.command("keys")
@click.argument("package", type=click.Path(exists=True, path_type=Path))
@click.option("--out", type=click.Path(path_type=Path), default=None, help="Also write the key lists as JSON.")
def foldbench_keys(package: Path, out: Path | None) -> None:
    """Count (and save) the protein / RNA / template / CCD keys a FoldBench package looks up."""
    import json
    from dataclasses import asdict

    from structcooker import foldbench

    keys = foldbench.load_keys(package)
    click.echo(json.dumps(keys.counts()))
    if out is not None:
        out.write_text(json.dumps(asdict(keys), indent=1))


@foldbench_group.command("build")
@click.argument("package", type=click.Path(exists=True, path_type=Path))
@click.argument("out", type=click.Path(path_type=Path))
@click.option("--release", type=click.Path(exists=True, path_type=Path), default=None,
              help="Release tree to cut from (default: BIOMOL_ROOT).")
@click.option("--ccd", type=click.Path(exists=True, path_type=Path), default=None,
              help="Inference CCD LMDB (default: RELEASE/lmdb/ccd/ccd_inference.lmdb, "
                   "scripts/staging/build_inference_ccd.py).")
@click.option("--fallback", type=click.Path(exists=True, path_type=Path), multiple=True,
              help="DB holding the from-scratch (edge) MSAs (default: RELEASE/lmdb/foldbench/msa_edge.lmdb).")
@click.option("--local", is_flag=True, help="Run here instead of as a SLURM job (a few thousand key reads).")
def foldbench_build(package: Path, out: Path, release: Path | None, ccd: Path | None, fallback: tuple[Path, ...],
                    local: bool) -> None:
    """Cut the FoldBench MSA / template / CCD LMDBs under OUT (never replaces an existing DB)."""
    if release is None:
        release = biomol_root()
    # The release keeps what FoldBench needs beyond the PDB-wide DBs under foldbench/.
    if ccd is None and (release / "lmdb" / "ccd" / "ccd_inference.lmdb").is_dir():
        ccd = release / "lmdb" / "ccd" / "ccd_inference.lmdb"
    if not fallback and (release / "lmdb" / "foldbench" / "msa_edge.lmdb").is_dir():
        fallback = (release / "lmdb" / "foldbench" / "msa_edge.lmdb",)
    argv = ["build", str(package.resolve()), str(out.resolve()), "--release", str(release.resolve())]
    if ccd is not None:
        argv += ["--ccd", str(ccd.resolve())]
    for f in fallback:
        argv += ["--fallback", str(f.resolve())]
    if local:
        from structcooker import foldbench

        raise SystemExit(foldbench.main(argv))
    from datacooker.executors.slurm import SlurmExecutor

    execu = SlurmExecutor(workdir=REPO / "logs" / "foldbench", repo=REPO)
    job = execu.run_once(name="foldbench_build", argv=[sys.executable, "-u", "-m", "structcooker.foldbench", *argv],
                         mem_gb=32, cores=4).job_id
    click.echo(f"foldbench build job {job}; report -> {out / 'foldbench_build.json'}")


@foldbench_group.command("verify")
@click.argument("package", type=click.Path(exists=True, path_type=Path))
@click.argument("out", type=click.Path(exists=True, path_type=Path))
def foldbench_verify(package: Path, out: Path) -> None:
    """Check every a3m_key / template lookup in the package's data.yaml resolves in OUT, and decode a sample."""
    from structcooker import foldbench

    raise SystemExit(foldbench.main(["verify", str(package), str(out)]))


if __name__ == "__main__":
    main()
