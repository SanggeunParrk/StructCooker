from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.release_patch import load_release_map

"""Load the pdbid -> release_date table as a shared metadata resource for the
release-date patch rebuilds."""

metadata_recipe = RecipeBook()

metadata_recipe.step(
    outputs=(("release_map", dict),),
    instruction=load_release_map,
    kwargs={
        "release_table_path": ("release_table_path", Path),
    },
)

RECIPE = metadata_recipe
TARGETS = ["release_map"]
