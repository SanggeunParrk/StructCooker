import tarfile

import lmdb
import pytest
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.instructions.design import stats


def _ligand_only():
    lig = chain(
        "L", [residue("XYZ", "X", [("C1", "C", (0, 0, 9))])], etype="non-polymer"
    )
    return make_mol([lig])


def _source(tmp_path, keys):
    """Make a source DB of dummy records plus the index sidecar a rebuild leaves."""
    path = tmp_path / "src.lmdb"
    env = lmdb.open(str(path), map_size=1 << 24, subdir=True)
    with env.begin(write=True) as txn:
        for key in keys:
            txn.put(key.encode(), b"x")
    env.close()
    index = tmp_path / "src.lmdb.index.tsv"
    index.write_text(
        "# db\tsrc.lmdb\n" + "".join(f"{k}\t10\t0.0\n" for k in sorted(keys))
    )
    return path


@pytest.fixture
def mols(monkeypatch):
    table = {
        "b": protein_chain("A", 5),
        "a": protein_chain("A", 3),
    }

    def fake_load(raw):
        raw = bytes(raw)
        if raw == b"bad":
            msg = "cannot decode"
            raise ValueError(msg)
        return make_mol([table[raw.decode()]]) if raw != b"x" else _ligand_only()

    monkeypatch.setattr(stats, "load_item", fake_load)


def test_item_stats_columns_and_backbone():
    columns, pdb = stats.item_stats(make_mol([protein_chain("A", 4)]))
    assert len(columns) == len(stats.STAT_COLUMNS) - 1
    assert columns[0] == "4"  # length: residues with a complete backbone
    assert float(columns[1]) == pytest.approx(15.0, abs=0.01)  # mean b factor
    assert all(0.0 <= float(c) <= 1.0 for c in columns[2:5])
    assert pdb is not None
    assert pdb.endswith("END\n")
    assert len(pdb.splitlines()) == 4 * 4 + 1


def test_item_without_a_protein_backbone_has_no_pdb():
    columns, pdb = stats.item_stats(_ligand_only())
    assert columns[0] == "0"
    assert pdb is None


def test_keys_come_from_the_index_sorted(tmp_path):
    src = _source(tmp_path, ["b", "a", "c"])
    assert stats.source_keys(src) == ["a", "b", "c"]


def test_keys_are_scanned_when_there_is_no_index(tmp_path):
    src = _source(tmp_path, ["b", "a"])
    (tmp_path / "src.lmdb.index.tsv").unlink()
    assert stats.source_keys(src) == ["a", "b"]


def test_parts_cover_the_keys_in_order(tmp_path):
    src = _source(tmp_path, [f"k{i}" for i in range(5)])
    parts = stats.pending_parts(src, tmp_path / "stats", 2)
    assert [p["part"] for p in parts] == ["part_00000", "part_00001", "part_00002"]
    assert [p["keys"] for p in parts] == [["k0", "k1"], ["k2", "k3"], ["k4"]]


def test_a_finished_part_is_not_listed_again(tmp_path):
    src = _source(tmp_path, [f"k{i}" for i in range(5)])
    out = tmp_path / "stats"
    out.mkdir()
    (out / "part_00001.tsv").write_text("")
    parts = stats.pending_parts(src, out, 2)
    assert [p["part"] for p in parts] == ["part_00000", "part_00002"]


def test_empty_source_is_an_error(tmp_path):
    env = lmdb.open(str(tmp_path / "e.lmdb"), subdir=True)
    env.close()
    with pytest.raises(ValueError, match="no records"):
        stats.pending_parts(tmp_path / "e.lmdb", tmp_path / "o", 10)


def test_write_part_writes_rows_tar_and_failures(tmp_path, mols):
    src = tmp_path / "src.lmdb"
    env = lmdb.open(str(src), map_size=1 << 24, subdir=True)
    with env.begin(write=True) as txn:
        for key, value in {"a": b"a", "b": b"b", "bad": b"bad", "lig": b"x"}.items():
            txn.put(key.encode(), value)
    env.close()
    out, bb = tmp_path / "stats", tmp_path / "backbones"
    status = stats.write_part(
        src, out, bb, "part_00000", ["a", "b", "bad", "lig", "gone"]
    )
    assert status == "part_00000: 3 items, 2 failed"
    rows = [ln.split("\t") for ln in (out / "part_00000.tsv").read_text().splitlines()]
    assert [r[0] for r in rows] == ["a", "b", "lig"]
    assert [r[1] for r in rows] == ["3", "5", "0"]
    failed = (out / "part_00000.failed.tsv").read_text().splitlines()
    assert [ln.split("\t")[0] for ln in failed] == ["bad", "gone"]
    assert "cannot decode" in failed[0]
    with tarfile.open(bb / "part_00000.tar") as tar:
        assert tar.getnames() == ["a.pdb", "b.pdb"]  # no backbone, no member
    assert not list(out.glob("*.tmp")) + list(bb.glob("*.tmp"))


def test_part_without_any_backbone_leaves_no_tar_and_no_failure_file(tmp_path, mols):
    src = _source(tmp_path, ["lig"])
    out, bb = tmp_path / "stats", tmp_path / "backbones"
    stats.write_part(src, out, bb, "part_00000", ["lig"])
    assert not (bb / "part_00000.tar").exists()
    assert not (out / "part_00000.failed.tsv").exists()
    assert (out / "part_00000.tsv").read_text().startswith("lig\t0\t")
