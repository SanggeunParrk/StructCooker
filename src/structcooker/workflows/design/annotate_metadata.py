from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.design.items import load_item_table

"""Load the items table of a design set once, as a node-local LMDB (the rebuild's metadata recipe)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("design_items", object),),
    instruction=load_item_table,
    kwargs={"items_path": ("items_path", Path)},
)

RECIPE = recipe
TARGETS = ["design_items"]
