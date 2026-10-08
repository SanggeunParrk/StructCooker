"""Every db/*.yaml must point at code that exists: recipe files and dotted hooks.

The end-to-end build (db/MANIFEST_all.yaml) runs every one of these configs; a recipe path or
reader/writer hook that no longer resolves would otherwise surface only when its stage runs,
possibly days into a build.
"""
import importlib
import re
from pathlib import Path

import pytest
from datacooker.config import load_config

REPO = Path(__file__).resolve().parents[1]
DB = REPO / "db"
CONFIGS = sorted(p for p in DB.rglob("*.yaml") if not p.name.startswith("MANIFEST"))
_DOTTED = re.compile(r"^structcooker(\.[A-Za-z_][A-Za-z0-9_]*)+$")


def _walk(value):
    if isinstance(value, dict):
        for v in value.values():
            yield from _walk(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from _walk(v)
    else:
        yield value


@pytest.mark.parametrize("path", CONFIGS, ids=lambda p: str(p.relative_to(DB)))
def test_config_references_resolve(path):
    cfg = load_config(path)
    for value in _walk(cfg):
        text = str(value)
        if text.endswith(".py") and "/src/structcooker/" in text:
            assert Path(text).is_file(), f"{path.name}: recipe {text} is missing"
        elif _DOTTED.match(text):
            module, _, attr = text.rpartition(".")
            assert hasattr(importlib.import_module(module), attr), f"{path.name}: {text} does not resolve"


@pytest.mark.parametrize("path", CONFIGS, ids=lambda p: str(p.relative_to(DB)))
def test_lmdb_configs_declare_a_schema(path):
    cfg = load_config(path)
    if cfg.get("env_path") or cfg.get("new_env_path"):
        assert cfg.get("schema"), f"{path.name} builds an LMDB but has no schema:"
