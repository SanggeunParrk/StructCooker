"""Recipe: one unit (a sub-entry made by the reader adapter) -> its stored record.

The adapter (``structcooker.instructions.unit.units.explode_units``) turns a ``cif_pdb`` entry
into one sub-entry per unit; this step builds each unit's record.
"""

from datacooker import RecipeBook

from structcooker.instructions.unit.units import UnitInput, build_unit_record

recipe = RecipeBook()

recipe.step(
    outputs=(("unit_record", dict),),
    instruction=build_unit_record,
    kwargs={"unit_input": ("unit_input", UnitInput)},
)

RECIPE = recipe
TARGETS = ["unit_record"]
