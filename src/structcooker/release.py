"""Verified releases as finite native SLURM jobs, with no resident polling process."""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any
from uuid import uuid4

from datacooker.config import load_config
from datacooker.config.runtime import activate_environment
from datacooker.executors.slurm import SlurmExecutor
from datacooker.utils.paths import scan_paths
from omegaconf import OmegaConf

from structcooker import preflight
from structcooker.build_state import (
    Fingerprints,
    claim_output,
    digest,
    exclusive,
    write_json,
)
from structcooker.build_verify import verify_output
from structcooker.production import encode_engine_paths, freeze_config, snapshot

_ACTIVE = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED"}


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text()) if path.exists() else {}


def _state(run: Path, name: str) -> Path:
    return run / "nodes" / name / "state.json"


def input_identity(cfg: dict[str, Any]) -> str:
    """Include raw payloads, auxiliary structures and declared reader references."""
    hashes = Fingerprints()
    inputs = {}
    for path in preflight.input_paths(cfg):
        if str(path) == str(cfg.get("data_dir")):
            patterns = cfg.get("input_patterns", [cfg.get("file_pattern", "*")])
            paths = sorted({p for pattern in patterns for p in scan_paths(path, pattern=pattern)})
            inputs[str(path)] = hashes.files(paths)
        elif str(path) == str(cfg.get("file_list")) and not cfg.get("keyed"):
            inputs[str(path)] = {"list": hashes.file(path),
                                 "payloads": hashes.files([Path(s) for s in path.read_text().splitlines() if s])}
        else:
            inputs[str(path)] = hashes.path(path)
    return digest({"config": freeze_config(cfg, Path.cwd(), Path.cwd()), "inputs": inputs})


def output_identity(cfg: dict[str, Any]) -> dict[str, Any]:
    """Include physical shard payloads and required publication sidecars."""
    paths = preflight.output_paths(cfg)
    for p in list(paths):
        if (p / "data.mdb").exists() or (p / "shards.json").exists():
            paths.extend([Path(str(p) + ".index.tsv"), Path(str(p) + ".meta.json")])
    return {str(p): Fingerprints().path(p) for p in paths}


def _executor(run: Path, folder: Path) -> SlurmExecutor:
    definition = _read(run / "release.json")
    os.environ.setdefault("DATACOOKER_EXCLUDE_NODES", "node20")
    return SlurmExecutor(workdir=folder, repo=Path(definition["repo"]),
                         pixi_manifest=(Path(definition["repo"]) / "pyproject.toml").resolve())


def _submit(run: Path, name: str, action: str, depends: tuple[str, ...] = ()) -> str:
    definition = _read(run / "release.json")
    frozen = Path(definition["frozen"])
    folder = run / "jobs" / uuid4().hex
    command = ["env", f"PYTHONPATH={frozen}/src:{frozen}/libs/datacooker/src",
               sys.executable, "-m", "structcooker.release", action, str(run), name]
    handle = _executor(run, folder).run_once(
        name=f"release_{action}_{name.replace('/', '_')}", argv=command,
        mem_gb=128, cores=8, depends_on=depends,
    )
    if not handle.job_id:
        msg = "Scheduler did not return a job ID"
        raise RuntimeError(msg)
    return handle.job_id


def dispatch(run: Path) -> None:
    """Submit ready stage planners once, then exit; verification jobs dispatch children."""
    definition = _read(run / "release.json")
    with exclusive(run / "dispatch.lock"):
        states = {name: _read(_state(run, name)) for name in definition["graph"]}
        for name, parents in definition["graph"].items():
            state = states[name]
            if state.get("phase") != "waiting" or any(states[p].get("phase") != "complete" for p in parents):
                continue
            state.update(phase="submitting")
            write_json(_state(run, name), state)
            job = _submit(run, name, "plan")
            state.update(phase="queued", job=job, action="plan")
            write_json(_state(run, name), state)
        write_json(run / "summary.json", {"complete": all(s.get("phase") == "complete" for s in states.values()),
                                         "nodes": {n: s.get("phase") for n, s in states.items()}})


def start(repo: Path, manifest: Path, run: Path, policy: Path) -> None:
    """Freeze a portable definition or reconcile a previous run without waiting for jobs."""
    graph = OmegaConf.to_container(OmegaConf.load(manifest), resolve=True)
    if not isinstance(graph, dict) or any(not isinstance(v, list) for v in graph.values()):
        msg = "Manifest must map node names to dependency lists"
        raise ValueError(msg)
    graph = {str(n): [str(p) for p in parents] for n, parents in graph.items()}
    if any(Path(str(n)).is_absolute() or ".." in Path(str(n)).parts for n in graph):
        msg = "Manifest names must stay inside db/"
        raise ValueError(msg)
    configs = {str(n): load_config(repo / "db" / f"{n}.yaml") for n in graph}
    errors = preflight.dependency_errors(configs, graph)
    if errors:
        raise ValueError("\n".join(errors))
    run.mkdir(parents=True, exist_ok=True)
    with exclusive(run / "definition.lock"):
        frozen, code = snapshot(repo, run)
        encoded = {n: encode_engine_paths(cfg, freeze_config(cfg, repo, frozen)) for n, cfg in configs.items()}
        exceptions = _read(policy)
        for name, entries in exceptions.items():
            if not isinstance(entries, dict) or any(
                not isinstance(v, dict) or not v.get("reason")
                or not re.fullmatch(r"[0-9a-f]{64}", v.get("sha256", "")) for v in entries.values()
            ):
                msg = f"{name}: exceptions require exact input keys, SHA256 and reasons"
                raise ValueError(msg)
        definition = {"repo": str(repo), "frozen": str(frozen), "code": code,
                      "graph": graph, "configs": encoded, "policy": exceptions}
        signature = digest(definition)
        previous = _read(run / "release.json")
        if previous and previous.get("signature") != signature:
            msg = "Release definition changed; use a new run directory and output root"
            raise RuntimeError(msg)
        if not previous:
            write_json(run / "release.json", {**definition, "signature": signature})
            for name, cfg in encoded.items():
                folder = _state(run, name).parent
                folder.mkdir(parents=True, exist_ok=True)
                OmegaConf.save(OmegaConf.create(cfg), folder / "config.yaml")
                write_json(_state(run, name), {"phase": "waiting"})
        else:
            with exclusive(run / "dispatch.lock"):
                reconcile(run)
    dispatch(run)


def reconcile(run: Path) -> None:
    """One scheduler observation per outstanding stage; unknown submissions fail closed."""
    definition = _read(run / "release.json")
    states = {n: _read(_state(run, n)) for n in definition["graph"]}
    active = False
    for name, state in states.items():
        if state.get("phase") == "submitting":
            msg = f"{name}: ambiguous submission; reconcile its job manifests before retrying"
            raise RuntimeError(msg)
        if state.get("phase") in {"waiting", "complete"}:
            continue
        jobs = {str(state["job"])} if state.get("job") else set()
        attempt = Path(state["attempt"]) if state.get("attempt") else None
        if attempt:
            jobs.update(p.read_text().strip() for p in attempt.rglob("*.job-id"))
        if attempt and (attempt / "submit.out").exists():
            jobs.update(re.findall(r"\[slurm\].*?: job (\d+)", (attempt / "submit.out").read_text()))
        stage_active = False
        for job in jobs:
            result = subprocess.run(["scontrol", "show", "job", "-o", job], capture_output=True, text=True, check=False)  # noqa: S603, S607 - scheduler argv
            found = re.findall(r"JobState=(\S+)", result.stdout)
            if result.returncode or not found:
                receipts = list(run.rglob(f"{job}.exit.json"))
                if len(receipts) != 1:
                    msg = f"{name}: cannot establish terminal state for {job}"
                    raise RuntimeError(msg)
                found = ["COMPLETED" if _read(receipts[0])["exit_code"] == 0 else "FAILED"]
            if any(s in _ACTIVE for s in found) and "Reason=DependencyNeverSatisfied" not in result.stdout:
                stage_active = True
        if stage_active:
            active = True
            continue
        state.update(phase="waiting", error="Replay a terminal stage using its recorded attempt")
        write_json(_state(run, name), state)
    if not active:
        # A repeat invocation must check completed input/output identities on compute.
        for name, state in states.items():
            if state.get("phase") in {"complete", "failed", "verification_failed"}:
                state["phase"] = "waiting"
                write_json(_state(run, name), state)


def _pin_references(cfg: dict[str, Any], attempt: Path) -> dict[str, Any]:
    """Freeze small reader maps/recovery references while retaining original input identity."""
    import shutil

    cfg = dict(cfg)
    runtime = dict(cfg.get("runtime_environment", {}))
    for key, value in runtime.items():
        if not value:
            continue
        source = Path(value)
        if not source.is_file():
            msg = f"Missing declared reader input: {source}"
            raise FileNotFoundError(msg)
        target = attempt / "references" / (key + source.suffix)
        target.parent.mkdir(exist_ok=True)
        if key == "OPENFOLD_STRUCTURE_OVERRIDES":
            mapping = _read(source)
            pinned = {}
            for original, replacement in mapping.items():
                payload = Path(replacement)
                dest = target.parent / (Fingerprints().file(payload) + payload.suffix)
                shutil.copy2(payload, dest)
                pinned[original] = str(dest)
            write_json(target, pinned)
        else:
            shutil.copy2(source, target)
            if Fingerprints().file(source) != Fingerprints().file(target):
                msg = f"Reader reference changed while freezing: {source}"
                raise RuntimeError(msg)
        runtime[key] = str(target)
    runtime["STRUCTCOOKER_IMMUTABLE_REFERENCES"] = "1"
    cfg["runtime_environment"] = runtime
    return cfg


def plan_stage(run: Path, name: str) -> None:
    """Fingerprint on a compute node, submit a native pipeline and verifier, then exit."""
    definition = _read(run / "release.json")
    state_path = _state(run, name)
    state = _read(state_path)
    cfg = load_config(state_path.parent / "config.yaml")
    activate_environment(cfg)
    state.update(phase="planning")
    write_json(state_path, state)
    identity = digest({"inputs": input_identity(cfg), "definition": definition["signature"]})
    outputs = preflight.output_paths(cfg)
    if not outputs:
        msg = f"{name}: no declared outputs"
        raise ValueError(msg)
    for output in outputs:
        with exclusive(Path(str(output) + ".build.lock")):
            if not state.get("attempt") and output.exists():
                msg = f"Untracked output {output}; use a fresh output root"
                raise RuntimeError(msg)
            claim_output(output, state_path)
    if state.get("completed_identity") == identity and state.get("output") == output_identity(cfg):
        state.update(phase="complete", reused=True, error=None)
        write_json(state_path, state)
        dispatch(run)
        return
    attempt = Path(state["attempt"]) if state.get("identity") == identity and state.get("attempt") else None
    if attempt and (attempt / "submitted.json").exists():
        terminal = repair_attempt(run, cfg, attempt)
    else:
        if attempt and (attempt / "submit.out").exists():
            msg = "Ambiguous previous pipeline submission; reconcile its recorded job IDs first"
            raise RuntimeError(msg)
        attempt = state_path.parent / "attempts" / uuid4().hex
        attempt.mkdir(parents=True)
        engine_cfg = _pin_references(cfg, attempt)
        OmegaConf.save(OmegaConf.create(encode_engine_paths(engine_cfg, freeze_config(engine_cfg, Path(definition["repo"]), Path(definition["frozen"])))), attempt / "engine.yaml")
        state.update(attempt=str(attempt), identity=identity)
        write_json(state_path, state)
        if cfg.get("output_data_path"):
            op = "extract-lmdb" if cfg.get("db_path") or cfg.get("extract_recipe") else "run"
            handle = _executor(run, attempt).run_once(
                name="projection", argv=[sys.executable, "-m", "datacooker.cli.workflow", op, str(attempt / "engine.yaml")],
                mem_gb=128, cores=8,
            )
            terminal = str(handle.job_id)
        else:
            result = subprocess.run([sys.executable, "-m", "datacooker.cli.lmdb", "pipeline",  # noqa: S603 - captured config argv
                                     str(attempt / "engine.yaml"), "--workdir", str(attempt),
                                     "--schema", str(cfg["schema"]), "--repo", definition["repo"]],
                                    capture_output=True, text=True, check=False)
            (attempt / "submit.out").write_text(result.stdout)
            (attempt / "submit.err").write_text(result.stderr)
            result.check_returncode()
            match = re.search(r"^PIPELINE_TERMINAL_JOB=(\d+)$", result.stdout, re.MULTILINE)
            if not match:
                msg = "Ambiguous pipeline submission; inspect attempt manifests"
                raise RuntimeError(msg)
            terminal = match.group(1)
        write_json(attempt / "submitted.json", {"terminal": terminal})
    with exclusive(run / "dispatch.lock"):
        job = _submit(run, name, "verify", (terminal,))
        state.update(phase="submitted", job=job, terminal=terminal, identity=identity, action="verify")
        write_json(state_path, state)


def repair_attempt(run: Path, cfg: dict[str, Any], attempt: Path) -> str:
    """Replay only incomplete worker commands, then publication/index; keep built shards."""
    executor = _executor(run, attempt / "repairs" / uuid4().hex)
    workers = []
    publish = index = None
    for path in sorted(attempt.glob("*/*.cmds")):
        lines = path.read_text().splitlines()
        if path.stem.endswith("_publish"):
            publish = shlex.split(lines[0])
        elif path.stem.endswith("_index"):
            index = shlex.split(lines[0])
        else:
            for i, line in enumerate(lines):
                command = shlex.split(line)
                target = Path(command[command.index("--output") + 1])
                report = _read(Path(str(target) + ".build-report.json"))
                if report and not report.get("failed") and not report.get("skipped_empty"):
                    continue
                handle = executor.run_once(name=f"repair_{path.stem}_{i}", argv=command, mem_gb=128, cores=8)
                workers.append(str(handle.job_id))
    if publish is None or index is None:
        if cfg.get("output_data_path"):
            # Projection retries overwrite their own declared outputs through the same recipe.
            op = "extract-lmdb" if cfg.get("db_path") or cfg.get("extract_recipe") else "run"
            handle = executor.run_once(name="projection_retry", argv=[sys.executable, "-m", "datacooker.cli.workflow", op,
                                                                       str(attempt / "engine.yaml")], mem_gb=128, cores=8)
            return str(handle.job_id)
        msg = "Incomplete submission manifests; reconcile before retrying"
        raise RuntimeError(msg)
    merged = executor.run_once(name="repair_publish", argv=publish, mem_gb=128, cores=8, depends_on=tuple(workers))
    indexed = executor.run_once(name="repair_index", argv=index, mem_gb=64, cores=8, depends_on=(str(merged.job_id),))
    return str(indexed.job_id)


def verify_stage(run: Path, name: str) -> None:
    """Verify before releasing children; failed verification retains the stage attempt."""
    definition = _read(run / "release.json")
    path = _state(run, name)
    state = _read(path)
    cfg = load_config(path.parent / "config.yaml")
    state["phase"] = "verifying"
    write_json(path, state)
    attempt = Path(state["attempt"])
    frozen_cfg = load_config(attempt / "engine.yaml")
    activate_environment(frozen_cfg)
    report = verify_output(frozen_cfg, attempt, definition["policy"].get(name, {}))
    identity = digest({"inputs": input_identity(cfg), "definition": definition["signature"]})
    if identity != state["identity"]:
        msg = "Inputs changed during stage execution"
        raise RuntimeError(msg)
    state.update(phase="complete", completed_identity=identity, verification=report,
                 output=output_identity(cfg), error=None)
    write_json(path, state)
    dispatch(run)


def main() -> None:
    """Run one finite job phase; public submission is through the StructCooker CLI."""
    action, directory, name = sys.argv[1:]
    run = Path(directory)
    if not os.environ.get("SLURM_JOB_ID"):
        msg = "Release stage work requires a compute allocation"
        raise RuntimeError(msg)
    # A short submission lock prevents a fast child from reading half a submission.
    with (run / "dispatch.lock").open("a") as handle:
        import fcntl

        fcntl.flock(handle, fcntl.LOCK_EX)
        state = _read(_state(run, name))
        if state.get("job") != os.environ["SLURM_JOB_ID"]:
            msg = "Stale stage job; a different submission owns this state"
            raise RuntimeError(msg)
    try:
        (plan_stage if action == "plan" else verify_stage)(run, name)
    except BaseException as exc:
        path = _state(run, name)
        state = _read(path)
        state.update(phase="failed" if action == "plan" else "verification_failed", error=str(exc))
        write_json(path, state)
        dispatch(run)
        raise


if __name__ == "__main__":
    main()
