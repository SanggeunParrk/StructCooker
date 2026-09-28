from pathlib import Path

# Test transaction rollback and injected loader failures at the actual engine boundary.
# ruff: noqa: SIM117, SLF001
from unittest.mock import patch

import lmdb
import pytest
from datacooker._ray import Admission
from datacooker.conditions import output_absent
from datacooker.config.runtime import canonical_inputs
from datacooker.lmdb import core

from structcooker import preflight, schemas
from structcooker.instructions.readers import openfold
from structcooker.instructions.transforms import template


def test_stale_sidecar_cannot_prove_completion(tmp_path):
    target = tmp_path / "missing.lmdb"
    Path(str(target) + ".index.tsv").write_text("stale")
    with pytest.raises(RuntimeError, match="receipt"):
        output_absent({"env_path": target})


def test_duplicate_writes_fail_in_both_completion_orders(tmp_path):
    for i, values in enumerate([(b"a", b"b"), (b"b", b"a")]):
        with lmdb.open(str(tmp_path / str(i)), map_size=1048576) as env:
            with patch.object(core, "iter_parallel_chunks", return_value=iter([[(b"key", v, None) for v in values]])):
                with pytest.raises(ValueError, match="Duplicate output key"):
                    core._parallel_write(env=env, items=[1, 2], chunk_size=2, n_jobs=1, process_item=None)
            assert env.stat()["entries"] == 0


def test_canonical_source_does_not_depend_on_enumeration(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for paths in ([a, b], [b, a]):
        selected, excluded = canonical_inputs(paths, lambda _p: "same", "first_path")
        assert selected == [a]
        assert excluded == [{"key": "same", "selected": str(a), "excluded": str(b)}]
        with pytest.raises(ValueError, match="Duplicate input key"):
            canonical_inputs(paths, lambda _p: "same")


def test_reader_cache_changes_with_path_and_content(tmp_path, monkeypatch):
    monkeypatch.delenv("STRUCTCOOKER_IMMUTABLE_REFERENCES", raising=False)
    monkeypatch.setattr(openfold, "_MONOMER_SEQID", None)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text("entry\tP1\n")
    b.write_text("entry\tP2\n")
    monkeypatch.setenv("MONOMER_SEQID_MAP", str(a))
    assert openfold.openfold_seqid_key(tmp_path / "entry/alignment.npz") == "P1"
    monkeypatch.setenv("MONOMER_SEQID_MAP", str(b))
    assert openfold.openfold_seqid_key(tmp_path / "entry/alignment.npz") == "P2"
    b.write_text("entry\tP333\n")
    assert openfold.openfold_seqid_key(tmp_path / "entry/alignment.npz") == "P333"


def test_template_corruption_fails_and_absence_is_accounted(tmp_path):
    hits = {"1abc_A": ("A", "A")}
    with patch.object(template, "load_raw_data", return_value=b"broken"):
        with patch.object(template, "load_bytes", side_effect=ValueError("corrupt")):
            with pytest.raises(ValueError, match="1abc_A"):
                template.load_templates_with_report(tmp_path, hits, 20)
    with patch.object(template, "load_raw_data", return_value=None):
        mols, report = template.load_templates_with_report(tmp_path, hits, 20)
        assert mols == {}
        assert report["missing_chain_hits"] == ["1abc_A"]
        assert not schemas.get("D").validate({"template_mols": mols, "template_report": report})
        report["candidate_count"] = 2
        assert schemas.get("D").validate({"template_mols": mols, "template_report": report})
        with pytest.raises(ValueError, match="Missing reference"):
            template.load_templates_from_chain_db(tmp_path, hits, 20)


def test_declared_reader_references_include_recovery_payload(tmp_path):
    import json

    manifest, recovered = tmp_path / "overrides.json", tmp_path / "recovered.npz"
    manifest.write_text(json.dumps({"original": str(recovered)}))
    cfg = {"runtime_environment": {"OPENFOLD_STRUCTURE_OVERRIDES": str(manifest),
                                   "MONOMER_SEQID_MAP": str(tmp_path / "ids.tsv")}}
    assert set(preflight.input_paths(cfg)) == {manifest, recovered, tmp_path / "ids.tsv"}


def test_admission_recovers_on_shared_node():
    gate = Admission(32)
    for _ in range(1000):
        gate.observe(0.62, 1)
    assert gate.target == 32
    gate.observe(0.8, 0)
    assert gate.target == 16
    assert not gate.admits(1, 0.8)
