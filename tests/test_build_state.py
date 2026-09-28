import json

import pytest

from structcooker.build_state import (
    BuildFailedError,
    BuildState,
    Fingerprints,
    exclusive,
)


def callbacks(tmp_path, events):
    output = tmp_path / "output"

    def submit():
        events.append("submit")
        output.write_text("built")
        return "123"

    def wait(job):
        events.append(("wait", job))

    return {"submit": submit, "wait": wait, "verify": lambda: {"checked": True},
            "unchanged": lambda: True,
            "output_signature": lambda: output.read_text() if output.exists() else None}


def test_verified_output_reused_but_changed_output_rebuilt(tmp_path):
    events = []
    state = BuildState(tmp_path / "state.json")
    options = callbacks(tmp_path, events)
    assert state.run(identity="one", **options) == "built"
    assert BuildState(state.path).run(identity="one", **options) == "reused"
    (tmp_path / "output").write_text("corrupted")
    assert state.run(identity="one", **options) == "built"
    assert events.count("submit") == 2


def test_input_change_invalidates_receipt(tmp_path):
    events = []
    state = BuildState(tmp_path / "state.json")
    options = callbacks(tmp_path, events)
    state.run(identity="old", **options)
    assert state.run(identity="new", **options) == "built"
    assert events.count("submit") == 2


def test_interrupted_wait_resumes_same_job_without_resubmission(tmp_path):
    events = []
    options = callbacks(tmp_path, events)
    normal_wait = options["wait"]

    def interrupt(_job):
        raise KeyboardInterrupt

    options["wait"] = interrupt
    state = BuildState(tmp_path / "state.json")
    with pytest.raises(KeyboardInterrupt):
        state.run(identity="same", **options)
    options["wait"] = normal_wait
    BuildState(state.path).run(identity="same", **options)
    assert events.count("submit") == 1


def test_interrupted_verification_does_not_rebuild(tmp_path):
    events = []
    options = callbacks(tmp_path, events)

    def interrupt():
        raise KeyboardInterrupt

    options["verify"] = interrupt
    state = BuildState(tmp_path / "state.json")
    with pytest.raises(KeyboardInterrupt):
        state.run(identity="same", **options)
    assert json.loads(state.path.read_text())["phase"] == "built"
    options["verify"] = lambda: {"checked": True}
    BuildState(state.path).run(identity="same", **options)
    assert events.count("submit") == 1


def test_unknown_job_does_not_trigger_duplicate_submit(tmp_path):
    events = []
    options = callbacks(tmp_path, events)

    def unknown(_job):
        msg = "Scheduler record unavailable"
        raise RuntimeError(msg)

    options["wait"] = unknown
    state = BuildState(tmp_path / "state.json")
    for _ in range(2):
        with pytest.raises(RuntimeError, match="unavailable"):
            state.run(identity="same", **options)
    assert events.count("submit") == 1
    with pytest.raises(RuntimeError, match="Inputs changed"):
        state.run(identity="different", **options)


def test_confirmed_terminal_failure_can_retry(tmp_path):
    events = []
    options = callbacks(tmp_path, events)
    normal_wait = options["wait"]

    def failed(_job):
        msg = "Terminal job failed"
        raise BuildFailedError(msg)

    options["wait"] = failed
    state = BuildState(tmp_path / "state.json")
    with pytest.raises(BuildFailedError):
        state.run(identity="same", **options)
    options["wait"] = normal_wait
    state.run(identity="same", **options)
    assert events.count("submit") == 2


def test_ambiguous_submission_stops_retry(tmp_path):
    options = callbacks(tmp_path, [])

    def lost_response():
        msg = "Connection lost after submission"
        raise RuntimeError(msg)

    options["submit"] = lost_response
    state = BuildState(tmp_path / "state.json")
    with pytest.raises(RuntimeError, match="Connection"):
        state.run(identity="same", **options)
    with pytest.raises(RuntimeError, match="ambiguous"):
        state.run(identity="same", **options)


def test_mutating_inputs_never_get_complete_receipt(tmp_path):
    options = callbacks(tmp_path, [])
    options["unchanged"] = lambda: False
    state = BuildState(tmp_path / "state.json")
    with pytest.raises(RuntimeError, match="changed during"):
        state.run(identity="same", **options)
    assert state.data["phase"] == "verification_failed"


def test_exclusive_output_lock(tmp_path):
    path = tmp_path / "output.lock"
    with exclusive(path), pytest.raises(RuntimeError, match="Another build"), exclusive(path):
        pass
    with exclusive(path):
        pass


def test_fingerprints_detect_nested_content_changes_and_ignore_lmdb_locks(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    source = root / "source"
    source.write_text("abc")
    hashes = Fingerprints()
    old = hashes.path(root)
    source.write_text("xyz")
    assert hashes.path(root) != old
    (root / "data.mdb").write_bytes(b"payload")
    old = hashes.path(root)
    (root / "lock.mdb").write_bytes(b"readers")
    assert hashes.path(root) == old


def test_code_snapshot_is_immutable_and_captures_new_version(tmp_path):
    from structcooker.production import snapshot

    repo = tmp_path / "repo"
    source = repo / "src/package/module.py"
    source.parent.mkdir(parents=True)
    source.write_text("VALUE = 1\n")
    (repo / "libs/datacooker/src").mkdir(parents=True)
    (repo / "pyproject.toml").write_text("[project]\nname='example'\n")
    work = tmp_path / "run"
    first, signature = snapshot(repo, work)
    source.write_text("VALUE = 2\n")
    second, changed = snapshot(repo, work)
    assert signature != changed
    assert first != second
    assert (first / source.relative_to(repo)).read_text() == "VALUE = 1\n"
    (second / source.relative_to(repo)).write_text("tampered")
    with pytest.raises(RuntimeError, match="modified"):
        snapshot(repo, work)


def test_resolved_hooks_and_relative_recipes_are_frozen(tmp_path):
    from pathlib import Path

    from structcooker.instructions.readers.cif import get_cif_data
    from structcooker.production import freeze_config

    repo = tmp_path / "repo"
    frozen = tmp_path / "frozen"
    cfg = freeze_config({"recipe": Path("src/structcooker/recipe.py"),
                         "reader": {"loader": get_cif_data}}, repo, frozen)
    assert cfg["recipe"] == str(frozen / "src/structcooker/recipe.py")
    assert cfg["reader"]["loader"] == "structcooker.instructions.readers.cif.get_cif_data"


@pytest.mark.parametrize("exit_code", [0, 7])
def test_native_exit_receipt_preserves_command_status(tmp_path, monkeypatch, exit_code):
    import os
    import subprocess

    from datacooker.executors.slurm import SlurmExecutor

    from structcooker.build_state import BuildFailedError
    from structcooker.production import saved_exit

    binary = tmp_path / "pixi"
    binary.write_text('#!/bin/bash\nshift 4\nexec "$@"\n')
    binary.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    executor = SlurmExecutor(workdir=tmp_path / "attempt", repo=tmp_path, submit=False)
    executor.run_once(name="terminal", argv=["bash", "-c", f"exit {exit_code}"], mem_gb=1, cores=1)
    result = subprocess.run(["bash", str(executor.workdir / "terminal.sbatch")], check=False)  # noqa: S603, S607 - generated test script
    assert result.returncode == exit_code
    if exit_code:
        with pytest.raises(BuildFailedError, match="exited 7"):
            saved_exit(executor.workdir, "12345")
    else:
        assert saved_exit(executor.workdir, "12345")


def test_missing_or_wrong_exit_receipt_is_not_success(tmp_path):
    import json

    from structcooker.production import saved_exit

    assert not saved_exit(tmp_path, "12345")
    (tmp_path / "12345.exit.json").write_text(json.dumps({"job_id": "999", "exit_code": 0}))
    with pytest.raises(RuntimeError, match="Invalid terminal receipt"):
        saved_exit(tmp_path, "12345")


def test_output_ownership_survives_driver_lock_release(tmp_path):
    from structcooker.build_state import claim_output, exclusive

    output = tmp_path / "db"
    state = tmp_path / "run1/state.json"
    lock = tmp_path / "db.build.lock"
    with exclusive(lock):
        claim_output(output, state)
    with exclusive(lock):
        claim_output(output, state)
        with pytest.raises(RuntimeError, match="another run"):
            claim_output(output, tmp_path / "run2/state.json")


def test_engine_yaml_preserves_nested_path_types(tmp_path):
    from datacooker.config import load_config
    from omegaconf import OmegaConf

    from structcooker.production import encode_engine_paths, freeze_config

    cfg = {"metadata_input": {"tsv": tmp_path / "reference.tsv",
                              "paths": [tmp_path / "other.tsv"]}}
    frozen = freeze_config(cfg, tmp_path, tmp_path)
    engine = tmp_path / "engine.yaml"
    OmegaConf.save(OmegaConf.create(encode_engine_paths(cfg, frozen)), engine)
    actual = load_config(engine)
    assert actual == cfg
