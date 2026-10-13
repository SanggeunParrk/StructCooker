"""release_docs: sizes come from the tree; scenarios name only paths that exist."""
from pathlib import Path

import lmdb
import yaml
from datacooker.lmdb.sharded import publish_shards

from structcooker import release_docs

REPO = Path(__file__).resolve().parents[1]


def _db(path, n, prefix="k"):
    with lmdb.open(str(path), map_size=1 << 22) as env, env.begin(write=True) as txn:
        for i in range(n):
            txn.put(f"{prefix}{i}".encode(), b"x" * 100)


def test_disk_bytes_counts_lmdb_files_and_shards(tmp_path):
    _db(tmp_path / "one.lmdb", 3)
    assert release_docs.disk_bytes(tmp_path / "one.lmdb") == (tmp_path / "one.lmdb" / "data.mdb").stat().st_size
    (tmp_path / "coll.lmdb" / "shards").mkdir(parents=True)
    for i in range(2):
        _db(tmp_path / "coll.lmdb" / "shards" / f"{i:03d}.lmdb", 2, prefix=f"s{i}_")
    publish_shards([tmp_path / "coll.lmdb" / "shards" / f"{i:03d}.lmdb" for i in range(2)], tmp_path / "coll.lmdb",
                   overwrite=True)
    assert release_docs.disk_bytes(tmp_path / "coll.lmdb") == sum(
        (tmp_path / "coll.lmdb" / "shards" / f"{i:03d}.lmdb" / "data.mdb").stat().st_size for i in range(2))
    assert release_docs.disk_bytes(tmp_path) is None


def test_every_scenario_path_is_a_release_path(monkeypatch, tmp_path):
    """Scenario paths name configured outputs, the out-of-DAG DBs, or materials/ directories."""
    from datacooker.config import load_config

    from structcooker import preflight

    monkeypatch.setenv("BIOMOL_ROOT", str(tmp_path))
    outputs = set()
    for name in yaml.safe_load((REPO / "db/MANIFEST_all.yaml").read_text()):
        outputs |= {p.relative_to(tmp_path).as_posix()
                    for p in preflight.output_paths(load_config(REPO / "db" / f"{name}.yaml"))}
    extra = {"lmdb/ccd/ccd_inference.lmdb", "lmdb/foldbench/msa_edge.lmdb"}
    for item in yaml.safe_load((REPO / "release_docs/scenarios.yaml").read_text()):
        for entry in [*(item.get("take") or []), *(item.get("optional") or [])]:
            path = entry["path"]
            assert path in outputs or path in extra or path.startswith("materials/"), (item["id"], path)


def test_html_and_markdown_render_the_same_rows(tmp_path):
    rows = [{"path": "lmdb/pdb/cif/cif.lmdb", "group": "lmdb/pdb", "stage": "pdb/cif", "schema": "A",
             "title": "structure", "entries": 1234, "size": 5_000_000, "op": "build", "source": "materials/raw/cif"},
            {"path": "lmdb/ccd/x<y>.lmdb", "group": "lmdb/ccd", "stage": "", "schema": "", "title": "",
             "entries": None, "size": None, "op": "staging", "source": ""}]
    items = [{"id": "s", "title": "PDB", "when": "training", "take": [{"path": "lmdb/pdb/cif/cif.lmdb",
              "exists": True, "size": 5_000_000, "note": "structures"}], "take_total": 5_000_000}]
    report = {"stages": 94, "with_problems": 0, "when": "2026-10-13 10:00", "reports": []}
    page = release_docs.render_html(tmp_path, rows, items, report)
    assert "94/94 OK" in page and "5.0 MB" in page and "1,234" in page
    assert "x&lt;y&gt;" in page and "x<y>" not in page            # paths are escaped
    md = release_docs.databases_md(rows)
    assert "## `lmdb/ccd`" in md and "| 1,234 |" in md
    assert "**needed** — 5.0 MB" in release_docs.scenarios_md(items)
    assert "not run" in release_docs.render_html(tmp_path, rows, items, None)
