"""Extract per-chain records, preserving provided reference model choices."""
from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import extract_selected_chain

recipe = RecipeBook()
recipe.step(
    outputs=(("biomoldict", dict),), instruction=extract_selected_chain,
    kwargs={"record": ("record", dict), "chain_id": ("chain_id", str),
            "cif_key": ("cif_key", str),
            "model_cache": ("model_cache", dict | None),
            "reference_selection_path": ("reference_selection_path", Path | None)},
)
RECIPE = recipe
TARGETS = ["biomoldict"]
