"""How build and build-all submit workflow ops (projections and parallel searches)."""
from pathlib import Path
from types import SimpleNamespace

import pytest

from structcooker import workflow_ops
from structcooker.build_verify import verify_output


class FakeExecutor:
    def __init__(self):
        self.calls = []

    def run_once(self, **kw):
        self.calls.append(("once", kw))
        return SimpleNamespace(job_id="1")

    def run_tier_shards(self, **kw):
        self.calls.append(("array", kw))
        return SimpleNamespace(job_id="2")


@pytest.mark.parametrize(("cfg", "op"), [
    ({"old_env_path": "a", "env_path": "b"}, "rebuild"),
    ({"split_recipe_path": "s.py", "recipe_path": "r.py"}, "parallel"),
    ({"output_data_path": "o.tsv", "extract_recipe_path": "e.py"}, "extract"),
    ({"output_data_path": "o.tsv"}, "materialize"),
    ({"env_path": "x.lmdb"}, "build"),
])
def test_op_of(cfg, op):
    assert workflow_ops.op_of(cfg) == op


def test_parallel_with_node_count_is_an_array(monkeypatch):
    monkeypatch.delenv(workflow_ops.MAX_NODES_ENV, raising=False)
    execu = FakeExecutor()
    job = workflow_ops.submit(execu, "hmmsearch", "parallel", Path("c.yaml"), {"node_count": 8})
    kind, kw = execu.calls[0]
    assert (job, kind, kw["n_shards"]) == ("2", "array", 8)
    assert kw["argv_for_shard"](3)[-4:] == ["--node-rank", "3", "--node-count", "8"]


def test_node_cap_runs_a_parallel_op_on_one_node(monkeypatch):
    monkeypatch.setenv(workflow_ops.MAX_NODES_ENV, "1")
    execu = FakeExecutor()
    workflow_ops.submit(execu, "hmmsearch", "parallel", Path("c.yaml"), {"node_count": 8})
    kind, kw = execu.calls[0]
    assert kind == "once"
    assert kw["argv"][-2:] == ["parallel-run", "c.yaml"]


def test_materialize_is_sized_by_its_n_jobs(monkeypatch):
    monkeypatch.delenv(workflow_ops.MAX_NODES_ENV, raising=False)
    assert workflow_ops.resources("materialize", {"inputs": {"n_jobs": 112}}) == (490, 112, 1)
    assert workflow_ops.resources("materialize", {"inputs": {"n_jobs": 16}}) == (128, 16, 1)
    assert workflow_ops.resources("materialize", {}) == (128, 8, 1)


def test_parallel_output_directory_verifies_when_filled(tmp_path):
    out = tmp_path / "msa_wo_lower"
    cfg = {"split_recipe_path": "s.py", "inputs": {"output_dir": out}}
    out.mkdir()
    with pytest.raises(ValueError, match="missing or empty"):
        verify_output(cfg, tmp_path, {})
    (out / "000").mkdir()
    (out / "000" / "P1.a3m").write_text(">q\nA\n")
    assert verify_output(cfg, tmp_path, {})["outputs"] == [str(out)]
