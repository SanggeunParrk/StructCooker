"""Foldseek clustering of the design backbones (stage ``design/<set>_clusters``)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

CLUSTER_HEADER = ["representative", "member"]


def run_foldseek(args: list[str]) -> None:
    """Run ``foldseek <args>``; a failure raises (a stale or partial table must never be mistaken for a result)."""
    subprocess.run(["foldseek", *args], check=True)  # noqa: S603, S607 (foldseek on PATH, argv from the config)


def read_cluster_table(path: Path) -> list[list[str]]:
    """``createtsv`` output (representative, member per line) as rows; a ``.pdb`` suffix on a name is dropped."""
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            rep, member = line.rstrip("\n").split("\t")[:2]
            rows.append([rep.removesuffix(".pdb"), member.removesuffix(".pdb")])
    return rows


def foldseek_clusters(
    backbone_dir: Path,
    work_dir: Path,
    coverage: float = 0.9,
    threads: int = 8,
    split_memory_limit: str = "400G",
) -> list[list[str]]:
    """Cluster every backbone PDB in ``backbone_dir/part_*.tar`` with ``foldseek cluster`` (defaults, ``-c coverage``).

    Returns the table ``[["representative", "member"], ...]``, one row per item (a representative lists itself). The
    work directory is removed before and after; foldseek's tmp directory is several times the size of the input.
    """
    tars = sorted(Path(backbone_dir).glob("part_*.tar"))
    if not tars:
        msg = f"no part_*.tar in {backbone_dir}"
        raise FileNotFoundError(msg)
    work = Path(work_dir)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    db, clu, tmp, tsv = work / "db", work / "clu", work / "tmp", work / "clusters.tsv"
    try:
        run_foldseek(
            ["createdb", *map(str, tars), str(db), "--threads", str(threads), "-v", "1"]
        )
        run_foldseek(
            [
                "cluster",
                str(db),
                str(clu),
                str(tmp),
                "--threads",
                str(threads),
                "-c",
                str(coverage),
                "--split-memory-limit",
                split_memory_limit,
                "-v",
                "1",
            ]
        )
        run_foldseek(["createtsv", str(db), str(db), str(clu), str(tsv), "-v", "1"])
        rows = read_cluster_table(tsv)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return [CLUSTER_HEADER, *rows]
