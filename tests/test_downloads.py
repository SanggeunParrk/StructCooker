import gzip
from pathlib import Path

import pytest

from structcooker import cif_fixes, downloads
from structcooker.paths import mmcif_root


def test_download_uses_build_input_without_deleting(monkeypatch, tmp_path):
    monkeypatch.setenv("MMCIF_ROOT", str(tmp_path / "mirror"))
    calls = []
    monkeypatch.setattr(downloads.subprocess, "run", lambda argv, **_kw: calls.append(argv))
    assert downloads.download_mmcif(tmp_path, confirmed=True) == mmcif_root(tmp_path)
    assert calls[0][-1] == str(mmcif_root(tmp_path))
    assert "--delete" not in calls[0]


def test_unconfirmed_download_has_no_side_effects(monkeypatch, tmp_path):
    monkeypatch.setenv("MMCIF_ROOT", str(tmp_path / "mirror"))
    with pytest.raises(RuntimeError, match="--yes"):
        downloads.download_mmcif(tmp_path, confirmed=False)
    assert not mmcif_root(tmp_path).exists()


@pytest.mark.parametrize("compressed", [False, True])
def test_manual_fix_replaces_the_ingested_file(tmp_path, compressed):
    source = tmp_path / "source"
    source.mkdir()
    source_file = source / ("1abc.cif.gz" if compressed else "1abc.cif")
    opener = gzip.open if compressed else Path.open
    with opener(source_file, "wb") as stream:
        stream.write(b"data_1abc\n")
    mirror = tmp_path / "mirror"
    target = mirror / "ab" / "1abc.cif.gz"
    target.parent.mkdir(parents=True)
    with gzip.open(target, "wb") as stream:
        stream.write(b"old")
    (mirror / cif_fixes.APPLIED_MARKER).write_text("2def\n")
    assert cif_fixes.apply_fixes(["1abc"], source, mirror) == (["1abc"], [])
    assert list(mirror.rglob("*.cif.gz")) == [target]
    with gzip.open(target, "rb") as stream:
        assert stream.read() == b"data_1abc\n"
    assert cif_fixes.applied_ids(mirror) == {"1abc", "2def"}


def test_corrupt_replacement_preserves_original(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "1abc.cif.gz").write_bytes(b"broken gzip")
    mirror = tmp_path / "mirror"
    mirror.mkdir()
    target = mirror / "1abc.cif.gz"
    target.write_bytes(b"original")
    with pytest.raises(gzip.BadGzipFile):
        cif_fixes.apply_fixes(["1abc"], source, mirror)
    assert target.read_bytes() == b"original"
    assert list(mirror.iterdir()) == [target]
