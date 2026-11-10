import pickle
from types import SimpleNamespace

import pytest

from structcooker.instructions.design import items

HEADER = "\t".join(items.ITEM_COLUMNS)


def _stats(length="100", plddt="90.0", coil="0.3", run="5"):
    return [length, plddt, "0.4", "0.3", coil, run]


@pytest.mark.parametrize(
    ("stats", "expected"),
    [
        (_stats(), True),
        (_stats(plddt="80.0"), True),  # at the threshold
        (_stats(plddt="79.99"), False),
        (_stats(length="32"), True),
        (_stats(length="31"), False),
        (_stats(length="512"), True),
        (_stats(length="513"), False),
        (_stats(coil="0.4999"), True),
        (_stats(coil="0.5"), False),  # strictly below
        (_stats(run="20"), True),
        (_stats(run="21"), False),
    ],
)
def test_la_proteina_filters(stats, expected):
    assert items.passes(stats, 80.0, 32, 512, 0.5, 20) is expected


def _write_stats(tmp_path, rows, failed=()):
    out = tmp_path / "stats"
    out.mkdir(exist_ok=True)
    (out / "part_00000.tsv").write_text(
        "".join("\t".join([k, *_stats(**kw)]) + "\n" for k, kw in rows)
    )
    if failed:
        (out / "part_00000.failed.tsv").write_text(
            "".join(f"{k}\terr\n" for k in failed)
        )
    return out


def _write_clusters(tmp_path, pairs):
    path = tmp_path / "clusters.tsv"
    path.write_text(
        "representative\tmember\n" + "".join(f"{r}\t{m}\n" for r, m in pairs)
    )
    return path


def _source(tmp_path, keys):
    path = tmp_path / "src.lmdb"
    path.mkdir(exist_ok=True)
    (tmp_path / "src.lmdb.index.tsv").write_text(
        "# db\tsrc\n" + "".join(f"{k}\t10\t0\n" for k in keys)
    )
    return path


def test_items_table_joins_statistics_clusters_and_flag(tmp_path):
    stats_dir = _write_stats(
        tmp_path,
        [("a", {}), ("b", {"plddt": "50"}), ("c", {"length": "0"}), ("d", {})],
    )
    clusters = _write_clusters(tmp_path, [("a", "a"), ("a", "b"), ("d", "d")])
    src = _source(tmp_path, ["a", "b", "c", "d"])
    rows = items.build_items(stats_dir, clusters, src)
    assert rows[0] == list(items.ITEM_COLUMNS)
    by_id = {r[0]: r for r in rows[1:]}
    assert [r[0] for r in rows[1:]] == ["a", "b", "c", "d"]
    assert by_id["a"][7:] == ["a", "2", "1"]
    assert by_id["b"][7:] == ["a", "2", "0"]  # low pLDDT: in the cluster, does not pass
    assert by_id["c"][7:] == ["c", "1", "0"]  # no backbone: its own cluster
    assert by_id["d"][7:] == ["d", "1", "1"]


def test_apply_filters_false_passes_everything(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {"plddt": "1", "length": "3"})])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    rows = items.build_items(
        stats_dir, clusters, _source(tmp_path, ["a"]), apply_filters=False
    )
    assert rows[1][-1] == "1"


def test_thresholds_are_parameters(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {"plddt": "70"})])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    src = _source(tmp_path, ["a"])
    assert items.build_items(stats_dir, clusters, src)[1][-1] == "0"
    assert items.build_items(stats_dir, clusters, src, min_plddt=60.0)[1][-1] == "1"


def test_a_source_key_without_statistics_refuses_the_build(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {})])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    with pytest.raises(ValueError, match="1 keys missing"):
        items.build_items(stats_dir, clusters, _source(tmp_path, ["a", "z"]))


def test_statistics_of_an_item_not_in_the_source_refuse_the_build(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {}), ("x", {})])
    clusters = _write_clusters(tmp_path, [("a", "a"), ("x", "x")])
    with pytest.raises(ValueError, match="not in the source"):
        items.build_items(stats_dir, clusters, _source(tmp_path, ["a"]))


def test_failed_records_are_refused_unless_allowed(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {})], failed=["z"])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    src = _source(tmp_path, ["a", "z"])
    with pytest.raises(ValueError, match="1 records failed"):
        items.build_items(stats_dir, clusters, src)
    rows = items.build_items(stats_dir, clusters, src, max_failed=1)
    assert [r[0] for r in rows[1:]] == ["a"]


def test_an_item_with_a_backbone_in_no_cluster_refuses_the_build(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {}), ("b", {})])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    with pytest.raises(ValueError, match="b has a backbone but is in no cluster"):
        items.build_items(stats_dir, clusters, _source(tmp_path, ["a", "b"]))


def test_an_item_in_two_parts_refuses_the_build(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {})])
    (stats_dir / "part_00001.tsv").write_text("\t".join(["a", *_stats()]) + "\n")
    clusters = _write_clusters(tmp_path, [("a", "a")])
    with pytest.raises(ValueError, match="two parts"):
        items.build_items(stats_dir, clusters, _source(tmp_path, ["a"]))


def test_source_without_index_is_an_error(tmp_path):
    stats_dir = _write_stats(tmp_path, [("a", {})])
    clusters = _write_clusters(tmp_path, [("a", "a")])
    with pytest.raises(ValueError, match="no index sidecar"):
        items.build_items(stats_dir, clusters, tmp_path / "nothing.lmdb")


def _table(tmp_path, rows):
    path = tmp_path / "items.tsv"
    path.write_text(HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return items.load_item_table(path)


def _row(key, keep="1"):
    return [key, "120", "91.5", "0.5", "0.25", "0.25", "7", "rep", "4", keep]


def _mol(key):
    return SimpleNamespace(id=[key], to_dict=lambda: {"metadata": {"id": [key]}})


def test_item_table_returns_rows_by_column(tmp_path):
    table = _table(tmp_path, [_row("a")])
    assert table.get("a")["cluster_size"] == "4"
    assert table.get("zz") is None
    assert pickle.loads(pickle.dumps(table)).get("a")["length"] == "120"  # noqa: S301


def test_item_table_refuses_other_columns(tmp_path):
    path = tmp_path / "items.tsv"
    path.write_text("id\tx\n")
    with pytest.raises(ValueError, match="unexpected columns"):
        items.load_item_table(path)


def test_annotate_adds_the_metadata_and_keeps_the_rest(tmp_path):
    table = _table(tmp_path, [_row("a")])
    record = items.annotate_item(_mol("a"), table, "distillation_long")
    assert record["metadata"] == {
        "id": ["a"],
        "source": "distillation_long",
        "length": 120,
        "mean_plddt": 91.5,
        "ss_helix": 0.5,
        "ss_strand": 0.25,
        "ss_coil": 0.25,
        "ss_max_coil_run": 7,
        "fs_cluster": "rep",
        "fs_cluster_size": 4,
    }


def test_annotate_drops_failing_items_and_rejects_unknown_ones(tmp_path):
    table = _table(tmp_path, [_row("a", keep="0")])
    assert items.annotate_item(_mol("a"), table, "s") is None
    with pytest.raises(KeyError, match="not in the items table"):
        items.annotate_item(_mol("b"), table, "s")
