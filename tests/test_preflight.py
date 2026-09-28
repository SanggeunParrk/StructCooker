from pathlib import Path

import pytest
from click.testing import CliRunner
from datacooker.config import load_config
from omegaconf import OmegaConf

from structcooker import preflight
from structcooker.cli import DB_ROOT, cli
from structcooker.paths import distillation_root, mmcif_root


def test_nested_and_top_level_inputs(tmp_path):
    source = tmp_path / "input"
    cfg = {
        "env_path": tmp_path / "output",
        "recipe": tmp_path / "recipe.py",
        "ccd_db_path": source / "ccd.lmdb",
        "inputs": {"fasta_paths": [source / "a.fasta", source / "b.fasta"],
                   "out_path": tmp_path / "output.tsv",
                   "old_seq_id_map_path": source / "optional-seed.tsv"},
        "parameters": {"reference_selection_path": source / "choices.lmdb"},
    }
    assert set(preflight.input_paths(cfg)) == {
        source / "ccd.lmdb", source / "a.fasta", source / "b.fasta", source / "choices.lmdb",
    }


def test_dependency_validation(tmp_path):
    configs = {
        "producer": {"env_path": tmp_path / "db"},
        "middle": {"inputs": {"db_path": tmp_path / "db"}},
        "consumer": {"inputs": {"db_path": tmp_path / "db"}},
    }
    assert preflight.dependency_errors(configs, {
        "producer": [], "middle": ["producer"], "consumer": ["middle"],
    }) == []
    assert any("without depending" in e for e in preflight.dependency_errors(configs, {}))
    assert any("cycle" in e for e in preflight.dependency_errors(configs, {
        "producer": ["middle"], "middle": ["producer"],
    }))
    assert any("unknown" in e for e in preflight.dependency_errors(configs, {"producer": ["absent"]}))


def test_output_collision(tmp_path):
    configs = {"a": {"env_path": tmp_path / "same"}, "b": {"env_path": tmp_path / "same"}}
    assert any("both write" in e for e in preflight.dependency_errors(configs, {}))


@pytest.mark.parametrize("overrides", [False, True])
def test_portable_config_paths(monkeypatch, tmp_path, overrides):
    for key in preflight.check_env():
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("DATA_ROOT", str(tmp_path))
    if overrides:
        monkeypatch.setenv("MMCIF_ROOT", str(tmp_path / "custom-cif"))
        monkeypatch.setenv("DISTILLATION_ROOT", str(tmp_path / "custom-openfold"))
        monkeypatch.setenv("OUTPUT_ROOT", str(tmp_path / "custom-output"))
        monkeypatch.setenv("SEQ_CLUSTER30_PATH", str(tmp_path / "provided-clusters.tsv"))
    configs = {str(p.relative_to(DB_ROOT)): load_config(p)
               for p in DB_ROOT.rglob("*.yaml") if not p.name.startswith("MANIFEST")}
    for cfg in configs.values():
        for path in preflight.input_paths(cfg) + preflight.output_paths(cfg):
            assert path.is_relative_to(tmp_path), path
    assert configs["pdb/cif.yaml"]["data_dir"] == mmcif_root(tmp_path)
    assert configs["distillation/short_cif.yaml"]["data_dir"].is_relative_to(distillation_root(tmp_path))
    # Every PDB consumer reads the one PDB clustering (cPDB_*; docs/seq-id-and-cluster-scheme.md).
    cluster_paths = [configs[p]["metadata_input"]["seqcluster_path"]
                     for p in ("pdb/cif_attached.yaml", "valid/valid1_attach.yaml", "valid/valid2.yaml")]
    cluster_paths.append(configs["metadata/interacting_seq_clusters.yaml"]["inputs"]["seqcluster_path"])
    assert len(set(cluster_paths)) == 1
    if not overrides:
        assert cluster_paths[0] == configs["metadata/pdb_seq_cluster30.yaml"]["output_data_path"]
    assert preflight.download_hint(mmcif_root(tmp_path), str(tmp_path)) == "mmcif"
    assert preflight.download_hint(Path("/elsewhere/materials/raw/cif"), str(tmp_path)) is None


@pytest.mark.parametrize("manifest", sorted(DB_ROOT.glob("MANIFEST*.yaml")))
def test_shipped_manifest_dependencies(manifest):
    deps = OmegaConf.to_container(OmegaConf.load(manifest))
    configs = {name: load_config(DB_ROOT / f"{name}.yaml") for name in deps}
    assert preflight.dependency_errors(configs, deps) == []


def test_strict_inspect_missing_input(monkeypatch, tmp_path):
    monkeypatch.setenv("DATA_ROOT", str(tmp_path))
    monkeypatch.setenv("OUTPUT_ROOT", str(tmp_path / "out"))
    result = CliRunner().invoke(cli, ["inspect", "pdb/cif", "--strict"])
    assert result.exit_code == 1
    assert "ccd/biomol_CCD.lmdb" in result.output
    assert "Preflight failed" in result.output
