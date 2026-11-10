import tarfile
from pathlib import Path

import pytest

from structcooker.instructions.design import cluster


def _tar(path, names):
    with tarfile.open(path, "w") as tar:
        for name in names:
            (path.parent / name).write_text("END\n")
            tar.add(path.parent / name, arcname=name)


@pytest.fixture
def fake_foldseek(monkeypatch):
    calls = []

    def run(args):
        calls.append(args)
        if args[0] == "createtsv":
            Path(args[4]).write_text("a\ta\na\tb.pdb\nc\tc\n")

    monkeypatch.setattr(cluster, "run_foldseek", run)
    return calls


def test_table_has_a_header_and_one_row_per_member(tmp_path, fake_foldseek):
    bb = tmp_path / "bb"
    bb.mkdir()
    _tar(bb / "part_00000.tar", ["a.pdb"])
    _tar(bb / "part_00001.tar", ["b.pdb", "c.pdb"])
    rows = cluster.foldseek_clusters(bb, tmp_path / "work", coverage=0.8, threads=3)
    assert rows == [["representative", "member"], ["a", "a"], ["a", "b"], ["c", "c"]]


def test_foldseek_gets_every_tar_in_order_and_the_parameters(tmp_path, fake_foldseek):
    bb = tmp_path / "bb"
    bb.mkdir()
    _tar(bb / "part_00001.tar", ["b.pdb"])
    _tar(bb / "part_00000.tar", ["a.pdb"])
    cluster.foldseek_clusters(
        bb, tmp_path / "work", coverage=0.75, threads=7, split_memory_limit="9G"
    )
    create, run, tsv = fake_foldseek
    assert create[:3] == [
        "createdb",
        str(bb / "part_00000.tar"),
        str(bb / "part_00001.tar"),
    ]
    assert create[-4:-2] == ["--threads", "7"]
    assert run[0] == "cluster"
    assert run[run.index("-c") + 1] == "0.75"
    assert run[run.index("--split-memory-limit") + 1] == "9G"
    assert tsv[0] == "createtsv"


def test_work_dir_is_removed_afterwards_and_stale_content_is_not_reused(
    tmp_path, fake_foldseek
):
    bb = tmp_path / "bb"
    bb.mkdir()
    _tar(bb / "part_00000.tar", ["a.pdb"])
    work = tmp_path / "work"
    work.mkdir()
    (work / "stale").write_text("old")
    cluster.foldseek_clusters(bb, work)
    assert not work.exists()


def test_failing_foldseek_raises_and_cleans_up(tmp_path, monkeypatch):
    bb = tmp_path / "bb"
    bb.mkdir()
    _tar(bb / "part_00000.tar", ["a.pdb"])

    def run(args):
        if args[0] == "cluster":
            raise RuntimeError(args[0])

    monkeypatch.setattr(cluster, "run_foldseek", run)
    with pytest.raises(RuntimeError, match="cluster"):
        cluster.foldseek_clusters(bb, tmp_path / "work")
    assert not (tmp_path / "work").exists()


def test_no_archives_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        cluster.foldseek_clusters(tmp_path, tmp_path / "work")


def test_read_cluster_table_drops_pdb_suffixes(tmp_path):
    path = tmp_path / "c.tsv"
    path.write_text("a.pdb\ta.pdb\na.pdb\tb.pdb\n")
    assert cluster.read_cluster_table(path) == [["a", "a"], ["a", "b"]]
