"""Recipe: the manifest row (type, sizes, node hashes, ...) of one interface-unit record."""

from datacooker import RecipeBook

from structcooker.instructions.interface.record import manifest_row

recipe = RecipeBook()

recipe.step(
    outputs=(("manifest_row", dict),),
    instruction=manifest_row,
    kwargs={"rec": ("db_data", dict)},
)

RECIPE = recipe
TARGETS = ["manifest_row"]
