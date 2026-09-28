import importlib
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from structcooker.build_state import Fingerprints
from structcooker.build_verify import verify_output


def test_native_failure_ledger_and_verified_coverage(tmp_path):
    native = importlib.import_module("datacooker.cli.lmdb")
    source = tmp_path / "input"
    source.mkdir()
    (source / "good.a3m").write_text(">query\nACD\n>hit\nAaCD\n")
    (source / "bad.a3m").write_text(">query\nACD\n>hit\nA\n")
    target = tmp_path / "msa.lmdb"
    recipe = Path(__file__).resolve().parents[1] / "src/structcooker/workflows/ingest/a3m.py"
    cfg = {"env_path": str(target), "data_dir": str(source), "file_pattern": "*.a3m",
           "recipe": str(recipe), "n_jobs": 1, "test_run": False,
           "reader": {"loader": "structcooker.instructions.readers.a3m.get_a3m_data"},
           "writer": {"serializer": "structcooker.instructions.transforms.codecs.to_bytes"}}
    engine = tmp_path / "engine.json"
    engine.write_text(json.dumps(cfg))
    runner = CliRunner()
    result = runner.invoke(native.cli, ["build", str(engine)])
    assert result.exit_code == 0, result.output
    report = json.loads(Path(str(target) + ".build-report.json").read_text())
    assert report["failed"] == 1
    assert report["written"] == 1
    result = runner.invoke(native.cli, ["index", str(target), "--schema", "E"])
    assert result.exit_code == 0, result.output
    items = tmp_path / "items"
    items.mkdir()
    (items / "merge_shards.txt").write_text(str(target) + "\n")
    (items / "items_all_shard0.txt").write_text("\n".join(map(str, sorted(source.glob("*.a3m")))) + "\n")
    cfg["schema"] = "E"
    with pytest.raises(ValueError, match="Unapproved"):
        verify_output(cfg, tmp_path, {})
    policy = {"bad": {"reason": "Intentionally malformed regression fixture", "sha256": Fingerprints().file(source / "bad.a3m")}}
    audit = verify_output(cfg, tmp_path, policy)
    assert audit["records"] == 1
    with pytest.raises(ValueError, match="depth exceeds"):
        verify_output({**cfg, "parameters": {"max_depth": 1}}, tmp_path, policy)
    from datacooker.lmdb.sharded import publish_shards

    logical = tmp_path / "logical.lmdb"
    publish_shards([target], logical)
    cap = tmp_path / "cap.json"
    cap.write_text(json.dumps({
        "old_env_path": str(logical), "new_env_path": str(tmp_path / "cap.lmdb"),
        "recipe": str(recipe.parents[1] / "exports/cap_msa_depth.py"),
        "parameters": {"max_depth": 1},
        "reader": {"deserializer": "structcooker.instructions.transforms.codecs.from_bytes",
                   "adapter": "structcooker.instructions.transforms.openfold.adapt_msa_for_cap"},
        "writer": {"serializer": "structcooker.instructions.transforms.codecs.to_bytes"},
    }))
    result = runner.invoke(native.cli, ["rebuild", str(cap), "--n-jobs", "1"])
    assert result.exit_code == 0, result.output
    import lmdb

    from structcooker.instructions.transforms.codecs import from_bytes

    with lmdb.open(str(tmp_path / "cap.lmdb"), readonly=True) as env, env.begin() as txn:
        capped = from_bytes(txn.get(b"good"))
        assert capped["msa_dict"]["sequences"]["aligned_sequences"].shape[0] == 1
    report_path = Path(str(target) + ".build-report.json")
    report["written"] += 1
    report_path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="accounting"):
        verify_output(cfg, tmp_path, policy)
    report["written"] -= 1
    report_path.write_text(json.dumps(report))
    (items / "items_all_shard0.txt").write_text(str(source / "good.a3m") + "\n")
    with pytest.raises(ValueError, match="coverage"):
        verify_output(cfg, tmp_path, policy)
