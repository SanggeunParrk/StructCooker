"""Recipe (build): one entry's Re-Glyco worklist (read by ``readers.reglyco.read_worklist``) -> its stored value."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.reglyco_worklist import keep_records

recipe = RecipeBook()
recipe.step(
    outputs=(("records", dict),),
    instruction=keep_records,
    kwargs={"records": ("records", dict)},
)
RECIPE = recipe
TARGETS = ["records"]
