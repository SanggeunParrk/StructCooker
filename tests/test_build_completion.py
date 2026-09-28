import importlib

import pytest
from click.testing import CliRunner

cli_module = importlib.import_module("structcooker.cli")


@pytest.mark.parametrize("failed", [False, True])
def test_build_all_uses_verified_native_release(monkeypatch, tmp_path, failed):
    from structcooker import release

    manifest = tmp_path / "manifest.yaml"
    manifest.write_text("ccd/ccd: []\n")
    called = []

    def start(repo, manifest, run, policy):
        called.append((repo, manifest, run, policy))
        if failed:
            msg = "Untracked output"
            raise RuntimeError(msg)

    monkeypatch.setattr(release, "start", start)
    result = CliRunner().invoke(cli_module.cli, ["build-all", "--manifest", str(manifest),
                                                "--workdir", str(tmp_path / "run")])
    assert len(called) == 1
    assert result.exit_code == int(failed)
    assert "completed " not in result.output
    if not failed:
        assert "submitted/resumed" in result.output


def test_dry_run_does_not_submit(monkeypatch, tmp_path):
    from structcooker import release

    manifest = tmp_path / "manifest.yaml"
    manifest.write_text("ccd/ccd: []\n")
    calls = []
    monkeypatch.setattr(release, "start", lambda *args: calls.append(args))
    result = CliRunner().invoke(cli_module.cli, ["build-all", "--manifest", str(manifest), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert not calls
    assert "no jobs submitted" in result.output
