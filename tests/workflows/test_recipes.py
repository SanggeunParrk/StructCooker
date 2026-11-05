"""Every recipe module loads, declares what the engine reads, and is used by a config."""

import importlib
import pkgutil
from pathlib import Path

import pytest

from structcooker import workflows

REPO = Path(__file__).resolve().parents[2]
RECIPES = sorted(
    m.name
    for m in pkgutil.walk_packages(workflows.__path__, "structcooker.workflows.")
    if not m.ispkg
)
CONFIG_TEXT = "\n".join(p.read_text() for p in (REPO / "db_configs").rglob("*.yaml"))


@pytest.mark.parametrize("name", RECIPES)
def test_recipe_module_declares_recipe_and_targets(name):
    module = importlib.import_module(name)
    assert module.RECIPE is not None
    assert module.TARGETS


@pytest.mark.parametrize("name", RECIPES)
def test_recipe_is_used_by_a_config(name):
    path = "src/" + name.replace(".", "/") + ".py"
    assert path in CONFIG_TEXT, (
        f"{path} is not the recipe of any db_configs/ config: dead code"
    )
