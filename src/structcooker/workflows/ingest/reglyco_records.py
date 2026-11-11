"""Recipe (explode rebuild): a record's stored value, made by ``readers.reglyco.explode_reglyco_entry``, passes through."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.reglyco_worklist import keep_value

recipe = RecipeBook()
recipe.step(
    outputs=(("reglyco_record", bytes),),
    instruction=keep_value,
    kwargs={"reglyco_record": ("reglyco_record", bytes)},
)
RECIPE = recipe
TARGETS = ["reglyco_record"]
