import gzip
from pathlib import Path

import pytest

from structcooker import downloads
from structcooker.paths import mmcif_root


def test_download_uses_build_input_without_deleting(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(downloads.subprocess, "run", lambda argv, **_kw: calls.append(argv))
    assert downloads.download_mmcif(tmp_path, confirmed=True) == mmcif_root(tmp_path)
    assert calls[0][-1] == str(mmcif_root(tmp_path))
    assert "--delete" not in calls[0]


def test_unconfirmed_download_has_no_side_effects(monkeypatch, tmp_path):
    with pytest.raises(RuntimeError, match="--yes"):
        downloads.download_mmcif(tmp_path, confirmed=False)
    assert not mmcif_root(tmp_path).exists()


