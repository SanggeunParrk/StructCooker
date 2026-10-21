from click.testing import CliRunner

from structcooker.cli import cli


def test_interface_commands_print_their_sbatch(monkeypatch, tmp_path):
    monkeypatch.setenv("BIOMOL_ROOT", str(tmp_path))
    runner = CliRunner()
    out = runner.invoke(cli, ["interface", "build", "--dry-run"]).output
    assert "build_full.sbatch" in out and str(tmp_path / "lmdb/interface/pdb") in out
    out = runner.invoke(cli, ["interface", "validate", "--dry-run", "--n", "50", "--out", "/x/db"]).output
    assert "validate_db.sbatch /x/db 50 0" in out and "-c 4" in out
    for sub in ("annotate", "finalize", "split"):
        assert runner.invoke(cli, ["interface", sub, "--dry-run"]).exit_code == 0
