import json

import pytest

from structcooker.build_state import (
    Fingerprints,
    exclusive,
)


def callbacks(tmp_path, events):
    output = tmp_path / "output"

    def submit():
        events.append("submit")
        output.write_text("built")
        return "123"

    def wait(job):
        events.append(("wait", job))

    return {"submit": submit, "wait": wait, "verify": lambda: {"checked": True},
            "unchanged": lambda: True,
            "output_signature": lambda: output.read_text() if output.exists() else None}


def test_exclusive_output_lock(tmp_path):
    path = tmp_path / "output.lock"
    with exclusive(path), pytest.raises(RuntimeError, match="Another build"), exclusive(path):
        pass
    with exclusive(path):
        pass


def test_fingerprints_detect_nested_content_changes_and_ignore_lmdb_locks(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    source = root / "source"
    source.write_text("abc")
    hashes = Fingerprints()
    old = hashes.path(root)
    source.write_text("xyz")
    assert hashes.path(root) != old
    (root / "data.mdb").write_bytes(b"payload")
    old = hashes.path(root)
    (root / "lock.mdb").write_bytes(b"readers")
    assert hashes.path(root) == old


def test_code_snapshot_is_immutable_and_captures_new_version(tmp_path):
    from structcooker.release import snapshot

    repo = tmp_path / "repo"
    source = repo / "src/package/module.py"
    source.parent.mkdir(parents=True)
    source.write_text("VALUE = 1\n")
    (repo / "libs/datacooker/src").mkdir(parents=True)
    (repo / "pyproject.toml").write_text("[project]\nname='example'\n")
    work = tmp_path / "run"
    first, signature = snapshot(repo, work)
    source.write_text("VALUE = 2\n")
    second, changed = snapshot(repo, work)
    assert signature != changed
    assert first != second
    assert (first / source.relative_to(repo)).read_text() == "VALUE = 1\n"
    (second / source.relative_to(repo)).write_text("tampered")
    with pytest.raises(RuntimeError, match="modified"):
        snapshot(repo, work)


def test_resolved_hooks_and_relative_recipes_are_frozen(tmp_path):
    from pathlib import Path

    from structcooker.instructions.readers.cif import get_cif_data
    from structcooker.release import freeze_config

    repo = tmp_path / "repo"
    frozen = tmp_path / "frozen"
    cfg = freeze_config({"recipe": Path("src/structcooker/recipe.py"),
                         "reader": {"loader": get_cif_data}}, repo, frozen)
    assert cfg["recipe"] == str(frozen / "src/structcooker/recipe.py")
    assert cfg["reader"]["loader"] == "structcooker.instructions.readers.cif.get_cif_data"


def test_output_ownership_survives_driver_lock_release(tmp_path):
    from structcooker.build_state import claim_output, exclusive

    output = tmp_path / "db"
    state = tmp_path / "run1/state.json"
    lock = tmp_path / "db.build.lock"
    with exclusive(lock):
        claim_output(output, state)
    with exclusive(lock):
        claim_output(output, state)
        with pytest.raises(RuntimeError, match="another run"):
            claim_output(output, tmp_path / "run2/state.json")


def test_engine_yaml_preserves_nested_path_types(tmp_path):
    from datacooker.config import load_config
    from omegaconf import OmegaConf

    from structcooker.release import encode_engine_paths, freeze_config

    cfg = {"metadata_input": {"tsv": tmp_path / "reference.tsv",
                              "paths": [tmp_path / "other.tsv"]}}
    frozen = freeze_config(cfg, tmp_path, tmp_path)
    engine = tmp_path / "engine.yaml"
    OmegaConf.save(OmegaConf.create(encode_engine_paths(cfg, frozen)), engine)
    actual = load_config(engine)
    assert actual == cfg


def test_exclusive_waits_when_asked(tmp_path):
    import threading
    import time

    from structcooker.build_state import exclusive

    lock = tmp_path / "dispatch.lock"
    entered = []

    def peer():
        with exclusive(lock, wait=10):
            entered.append(True)

    with exclusive(lock):
        thread = threading.Thread(target=peer)
        thread.start()
        time.sleep(1)
        assert not entered
        with pytest.raises(RuntimeError, match="Another build owns"):
            with exclusive(lock):
                pass
    thread.join()
    assert entered


def test_directory_fingerprint_is_content_keyed_and_rejects_dir_symlinks(tmp_path):
    from structcooker.build_state import Fingerprints

    root = tmp_path / "msa"
    (root / "P" / "001").mkdir(parents=True)
    (root / "P" / "001" / "a.a3m").write_text(">q\nA\n")
    (root / "b.a3m").write_text(">q\nC\n")
    first = Fingerprints().path(root)
    assert sorted(first) == ["P/001/a.a3m", "b.a3m"]
    (root / "b.a3m").write_text(">q\nG\n")
    assert Fingerprints().path(root)["b.a3m"] != first["b.a3m"]
    other = tmp_path / "other"
    other.mkdir()
    (root / "linked").symlink_to(other)
    with pytest.raises(ValueError, match="Directory symlink"):
        Fingerprints().path(root)
