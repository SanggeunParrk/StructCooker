"""Finite scheduler DAG recovery without production data or scheduler access."""
import json
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from structcooker import release
from structcooker.build_state import write_json


def fixture_run(tmp_path):
    write_json(tmp_path / "release.json", {"graph": {"base": [], "cap": ["base"]}})
    for name in ("base", "cap"):
        write_json(tmp_path / "nodes" / name / "state.json", {"phase": "waiting"})
    return tmp_path


def test_dispatch_only_ready_nodes_once(tmp_path, monkeypatch):
    run = fixture_run(tmp_path)
    submit = Mock(return_value="123")
    monkeypatch.setattr(release, "_submit", submit)
    release.dispatch(run)
    release.dispatch(run)
    submit.assert_called_once_with(run, "base", "plan")
    assert json.loads((run / "summary.json").read_text())["complete"] is False


def test_active_repair_prevents_duplicate_submission(tmp_path, monkeypatch):
    run = fixture_run(tmp_path)
    attempt = run / "attempt"
    attempt.mkdir()
    (attempt / "repair.job-id").write_text("456\n")
    write_json(run / "nodes/base/state.json", {"phase": "failed", "job": "123", "attempt": str(attempt)})
    def observe(argv, **_kwargs):
        return SimpleNamespace(returncode=0, stdout="JobState=" + ("RUNNING" if argv[-1] == "456" else "FAILED"))
    monkeypatch.setattr(release.subprocess, "run", observe)
    release.reconcile(run)
    assert json.loads((run / "nodes/base/state.json").read_text())["phase"] == "failed"


def test_completed_run_must_revalidate(tmp_path):
    run = fixture_run(tmp_path)
    for name in ("base", "cap"):
        write_json(run / "nodes" / name / "state.json", {"phase": "complete"})
    release.reconcile(run)
    assert json.loads((run / "nodes/base/state.json").read_text())["phase"] == "waiting"
    assert json.loads((run / "nodes/cap/state.json").read_text())["phase"] == "waiting"


def test_unknown_scheduler_state_fails_closed(tmp_path, monkeypatch):
    run = fixture_run(tmp_path)
    write_json(run / "nodes/base/state.json", {"phase": "failed", "job": "123"})
    monkeypatch.setattr(release.subprocess, "run", lambda *_a, **_k: SimpleNamespace(returncode=1, stdout=""))
    with pytest.raises(RuntimeError, match="cannot establish terminal"):
        release.reconcile(run)


def test_identity_includes_payload_and_callable_stably(tmp_path):
    raw = tmp_path / "input"
    raw.mkdir()
    item = raw / "a.txt"
    item.write_text("one")
    cfg = {"data_dir": raw, "file_pattern": "*.txt", "key_builder": str}
    first = release.input_identity(cfg)
    assert release.input_identity(dict(cfg)) == first
    item.write_text("two")
    assert release.input_identity(cfg) != first
