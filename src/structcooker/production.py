"""Verified native builds with immutable code snapshots and persistent receipts."""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import ExitStack
from graphlib import TopologicalSorter
from importlib.metadata import distributions
from pathlib import Path
from threading import Semaphore
from typing import Any
from uuid import uuid4

from datacooker.config import load_config
from datacooker.executors.slurm import wait_for_jobs
from datacooker.utils.paths import scan_paths
from omegaconf import OmegaConf

from structcooker import preflight
from structcooker.build_state import (
    BuildFailedError,
    BuildState,
    Fingerprints,
    claim_output,
    digest,
    exclusive,
    write_json,
)
from structcooker.build_verify import verify_output


def encode_engine_paths(original: Any, resolved: Any) -> Any:
    """Preserve nested Path types when round-tripping frozen config through YAML."""
    if isinstance(original, Path):
        return "${p:" + json.dumps(str(resolved)) + "}"
    if isinstance(original, dict):
        return {k: encode_engine_paths(v, resolved[k]) for k, v in original.items()}
    if isinstance(original, (list, tuple)):
        return [encode_engine_paths(v, r) for v, r in zip(original, resolved, strict=True)]
    return resolved


def saved_exit(attempt: Path, job: str) -> bool:
    """Use a native script's exit receipt when SLURM has expired its job record."""
    paths = list(attempt.rglob(f"{job}.exit.json"))
    if len(paths) != 1:
        return False
    receipt = json.loads(paths[0].read_text())
    if receipt.get("job_id") != job or not isinstance(receipt.get("exit_code"), int):
        msg = f"Invalid terminal receipt: {paths[0]}"
        raise RuntimeError(msg)
    if receipt["exit_code"]:
        msg = f"Recorded job {job} exited {receipt['exit_code']}"
        raise BuildFailedError(msg)
    return True


def snapshot(repo: Path, workdir: Path) -> tuple[Path, str]:
    """Freeze source packages and dependency declarations before submitting any job."""
    fingerprints = Fingerprints()
    sources = [repo / "src", repo / "libs/datacooker/src", repo / "pyproject.toml", repo / "pixi.lock"]
    contents = {str(p.relative_to(repo)): fingerprints.path(p) for p in sources}
    runtime = {"python": sys.version, "packages": sorted((str(d.metadata.get("Name", "unknown")), d.version) for d in distributions())}
    signature = digest({"contents": contents, "runtime": runtime})
    destination = workdir / "snapshots" / signature
    if not destination.exists():
        temporary = destination.with_name(signature + "." + uuid4().hex)
        temporary.mkdir(parents=True)
        for source in sources:
            if not source.exists():
                continue
            target = temporary / source.relative_to(repo)
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.is_dir():
                shutil.copytree(source, target, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            else:
                shutil.copy2(source, target)
        temporary.rename(destination)
        write_json(destination / "environment.json", runtime)
    captured = {str(p.relative_to(repo)): fingerprints.path(destination / p.relative_to(repo)) for p in sources}
    # Optional dependency declaration files may be absent in test repositories.
    for source in sources:
        if not source.exists():
            captured[str(source.relative_to(repo))] = contents[str(source.relative_to(repo))]
    if captured != contents:
        msg = f"Code snapshot is incomplete or was modified: {destination}"
        raise RuntimeError(msg)
    return destination, signature


def freeze_config(cfg: dict[str, Any], repo: Path, frozen: Path) -> dict[str, Any]:
    """Resolve path objects and pin all recipe paths to the captured source tree."""
    def convert(value: Any) -> Any:
        if isinstance(value, dict):
            return {k: convert(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [convert(v) for v in value]
        if callable(value):
            module = getattr(value, "__module__", None)
            name = getattr(value, "__qualname__", "")
            if not module or not name or "<" in name:
                msg = f"Config callable must be importable by name: {value}"
                raise ValueError(msg)
            return f"{module}.{name}"
        if isinstance(value, Path):
            if not value.is_absolute():
                value = repo / value
            value = str(value)
        if isinstance(value, str) and value.startswith(str(repo / "src") + "/"):
            return str(frozen / Path(value).relative_to(repo))
        return value
    return convert(cfg)


def execute_manifest(
    repo: Path, manifest: Path, workdir: Path, policy: Path,
    *, rebuild_untracked: bool = False,
) -> dict[str, str]:
    """Build ready nodes, verify before reuse, and resume recorded scheduler jobs."""
    if not os.environ.get("SLURM_JOB_ID"):
        msg = "pdb-build must run in a compute allocation; hashing/auditing does not run on login nodes."
        raise RuntimeError(msg)
    raw = OmegaConf.to_container(OmegaConf.load(manifest), resolve=True)
    if not isinstance(raw, dict) or any(not isinstance(v, list) for v in raw.values()):
        msg = "Manifest must map node names to dependency lists"
        raise ValueError(msg)
    dependencies = {str(k): [str(d) for d in v] for k, v in raw.items()}
    if any(Path(name).is_absolute() or ".." in Path(name).parts for name in dependencies):
        msg = "Manifest node names must stay inside db/"
        raise ValueError(msg)
    configs = {name: load_config(repo / "db" / f"{name}.yaml") for name in dependencies}
    errors = preflight.dependency_errors(configs, dependencies)
    errors.extend(f"{report.name}: missing prerequisite {path}"
                  for report in preflight.inspect(configs) for path in report.missing)
    if errors:
        raise ValueError("\n".join(errors))
    exceptions = json.loads(policy.read_text())
    workdir.mkdir(parents=True, exist_ok=True)
    with exclusive(workdir / "run.lock"):
        frozen, code_hash = snapshot(repo, workdir)
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(frozen / "src") + os.pathsep + str(frozen / "libs/datacooker/src")
        verification_gate = Semaphore(1)

        def stage(name: str) -> str:
            cfg = freeze_config(configs[name], repo, frozen)
            directory = workdir / "nodes" / name
            directory.mkdir(parents=True, exist_ok=True)
            record = BuildState(directory / "state.json")
            outputs = preflight.output_paths(cfg)
            if not outputs:
                msg = f"{name}: no auditable outputs"
                raise ValueError(msg)
            failures = exceptions.get(name, {})
            if not isinstance(failures, dict) or any(
                not isinstance(v, dict) or not v.get("reason") or not re.fullmatch(r"[0-9a-f]{64}", v.get("sha256", ""))
                for v in failures.values()
            ):
                msg = f"{name}: every accepted failure needs an exact key, input SHA256, and reason"
                raise ValueError(msg)
            hashes = Fingerprints()

            def input_identity() -> str:
                inputs = {}
                for path in preflight.input_paths(cfg):
                    if str(path) == str(cfg.get("data_dir")):
                        inputs[str(path)] = hashes.files(sorted(scan_paths(path, pattern=cfg.get("file_pattern", "*"))))
                    elif str(path) == str(cfg.get("file_list")) and not cfg.get("keyed"):
                        inputs[str(path)] = hashes.files(list(map(Path, path.read_text().splitlines())))
                    else:
                        inputs[str(path)] = hashes.path(path)
                return digest({"config": cfg, "code": code_hash, "failures": failures,
                               "inputs": inputs})

            def output_signature() -> object:
                paths = list(outputs)
                for path in outputs:
                    if (path / "data.mdb").exists() or (path / "shards.json").exists():
                        paths.extend([Path(str(path) + ".index.tsv"), Path(str(path) + ".meta.json")])
                return {str(p): hashes.path(p) for p in paths}

            identity = input_identity()
            with ExitStack() as stack:
                for output in sorted(outputs):
                    stack.enter_context(exclusive(Path(str(output) + ".build.lock")))
                if not record.data and any(p.exists() for p in outputs) and not rebuild_untracked:
                    msg = f"{name}: existing output has no receipt; use a fresh OUTPUT_ROOT or explicitly rebuild untracked outputs."
                    raise RuntimeError(msg)
                for output in sorted(outputs):
                    claim_output(output, record.path)

                def submit() -> str:
                    attempt = directory / "attempts" / uuid4().hex
                    attempt.mkdir(parents=True)
                    engine = attempt / "engine.yaml"
                    typed_cfg = encode_engine_paths(configs[name], cfg)
                    engine_cfg = {k: v for k, v in typed_cfg.items() if k != "schema"}
                    OmegaConf.save(OmegaConf.create(engine_cfg), engine)
                    record.save(attempt=str(attempt))
                    if cfg.get("output_data_path"):
                        from datacooker.executors.slurm import SlurmExecutor

                        op = "extract-lmdb" if cfg.get("db_path") or cfg.get("extract_recipe") else "run"
                        executor = SlurmExecutor(workdir=attempt, repo=repo)
                        # Pin PYTHONPATH in the command, not in process-global env shared by threads.
                        command = ["env", "PYTHONPATH=" + environment["PYTHONPATH"], sys.executable,
                                   "-m", "datacooker.cli.workflow", op, str(engine)]
                        handle = executor.run_once(name=Path(name).name, argv=command, mem_gb=128, cores=4)
                        if handle.job_id is None:
                            msg = "Projection submission returned no job ID"
                            raise RuntimeError(msg)
                        return handle.job_id
                    command = [sys.executable, "-u", "-m", "datacooker.cli.lmdb", "pipeline", str(engine),
                               "--workdir", str(attempt), "--schema", str(cfg["schema"]), "--repo", str(repo)]
                    result = subprocess.run(command, env=environment, capture_output=True, text=True, check=False)  # noqa: S603 - argv, no shell
                    (attempt / "submit.out").write_text(result.stdout)
                    (attempt / "submit.err").write_text(result.stderr)
                    result.check_returncode()
                    matched = re.search(r"^PIPELINE_TERMINAL_JOB=(\d+)$", result.stdout, re.MULTILINE)
                    if not matched:
                        msg = "Missing terminal job ID; reconcile the saved submission log"
                        raise RuntimeError(msg)
                    return matched.group(1)

                def wait_job(job: str) -> None:
                    try:
                        wait_for_jobs((job,))
                    except RuntimeError as exc:
                        state = subprocess.run(["scontrol", "show", "job", "-o", job], capture_output=True, text=True, check=False)  # noqa: S603, S607 - scheduler argv
                        if state.returncode and saved_exit(Path(record.data["attempt"]), job):
                            return
                        if (re.search(r"JobState=(FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|PREEMPTED)\b", state.stdout)
                                or "Reason=DependencyNeverSatisfied" in state.stdout):
                            raise BuildFailedError(str(exc)) from exc
                        raise
                    result = subprocess.run(["scontrol", "show", "job", job], capture_output=True, text=True, check=False)  # noqa: S603, S607 - scheduler argv
                    if result.returncode and not saved_exit(Path(record.data["attempt"]), job):
                        result.check_returncode()
                    (directory / "terminal.txt").write_text(result.stdout)

                def verify() -> dict[str, Any]:
                    with verification_gate:
                        return verify_output(cfg, Path(record.data["attempt"]), failures)

                return record.run(identity=identity, submit=submit, wait=wait_job, verify=verify,
                                  unchanged=lambda: input_identity() == identity, output_signature=output_signature)

        graph = TopologicalSorter(dependencies)
        graph.prepare()
        results: dict[str, str] = {}
        try:
            with ThreadPoolExecutor(max_workers=2) as pool:
                pending = {}
                while graph.is_active():
                    for name in graph.get_ready():
                        pending[pool.submit(stage, name)] = name
                    completed, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in completed:
                        name = pending.pop(future)
                        results[name] = future.result()
                        graph.done(name)
                        write_json(workdir / "summary.json", {"phase": "running", "nodes": results})
        except BaseException as exc:
            write_json(workdir / "summary.json", {"phase": "failed", "nodes": results, "error": str(exc)})
            raise
        write_json(workdir / "summary.json", {"phase": "complete", "nodes": results})
        return results
