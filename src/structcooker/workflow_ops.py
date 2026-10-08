"""How a db config runs: its op, and for workflow ops the SLURM job that runs it.

Shared by ``structcooker build`` (one config) and ``build-all`` (the release DAG,
:mod:`structcooker.release`) so both submit a config the same way.
"""
from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path

    from datacooker.executors.slurm import SlurmExecutor

# op -> `datacooker.cli.workflow` subcommand; ops not listed (build/rebuild) are LMDB pipelines.
WORKFLOW_CMD = {"materialize": "run", "extract": "extract-lmdb", "parallel": "parallel-run"}
# Single-node resource sizing per projection op. Config keys can't carry hints (the
# workflow CLI re-reads the yaml and rejects unknown keys), so size by op: parallel
# fans out (full node); extract runs joblib over a whole DB; materialize is mostly
# single-process (a couple stream a multi-GB fasta, hence the mem headroom).
WORKFLOW_RESOURCES = {"parallel": (490, 112), "extract": (200, 32), "materialize": (128, 8)}
WORKFLOW_CLI = [sys.executable, "-u", "-m", "datacooker.cli.workflow"]
# Caps every config's node_count (a sampled e2e run sets 1: a few hundred items do not
# need an 8-node array). Splitting is by item key, so the cap changes no result.
MAX_NODES_ENV = "STRUCTCOOKER_MAX_NODES"


def op_of(cfg: Mapping) -> str:
    """Infer the op from a db config's shape.

    ``old_env_path`` -> rebuild; ``split_recipe`` -> parallel; ``output_data_path`` -> a
    projection op (``extract`` when it reads a source DB, else ``materialize``);
    ``env_path`` -> build. Projection and parallel ops produce files, not an LMDB, so they
    run through ``datacooker.cli.workflow`` instead of the planning-first pipeline.
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
    raise ValueError(msg)


def resources(op: str, cfg: Mapping) -> tuple[int, int, int]:
    """``(mem_gb, cores, nodes)`` for a workflow op."""
    mem_gb, cores = WORKFLOW_RESOURCES[op]
    # A materialize recipe that fans out itself (joblib) declares it as inputs.n_jobs;
    # without this it ran n_jobs workers on the 8 default cores (template_candidates:
    # 112 on 8, ~14x slower).
    n_jobs = int((cfg.get("inputs") or {}).get("n_jobs") or 0)
    if op == "materialize" and n_jobs > cores:
        mem_gb, cores = WORKFLOW_RESOURCES["parallel"] if n_jobs >= 112 else (mem_gb, n_jobs)
    nodes = int(cfg.get("node_count") or 1) if op == "parallel" else 1
    cap = os.environ.get(MAX_NODES_ENV)
    if cap:
        nodes = max(1, min(nodes, int(cap)))
    return mem_gb, cores, nodes


def submit(executor: SlurmExecutor, name: str, op: str, config: Path, cfg: Mapping,
           depends_on: tuple[str, ...] = ()) -> str:
    """Submit one workflow op; returns an afterok-chainable job id.

    A parallel op that declares ``node_count`` fans out as an N-task array, one whole node
    each; ``parallel-run`` splits the items by key, so tasks may start at different times.
    """
    argv = [*WORKFLOW_CLI, WORKFLOW_CMD[op], str(config)]
    mem_gb, cores, nodes = resources(op, cfg)
    if nodes > 1:
        handle = executor.run_tier_shards(
            name=name, n_shards=nodes, mem_gb=mem_gb, cores=cores, depends_on=depends_on,
            argv_for_shard=lambda i: [*argv, "--node-rank", str(i), "--node-count", str(nodes)],
        )
    else:
        handle = executor.run_once(name=name, argv=argv, mem_gb=mem_gb, cores=cores, depends_on=depends_on)
    return handle.job_id
