"""cif LMDB validation analysis recipe.

A Transform-only recipe that takes one decoded cif record and produces its
list of validity issues (empty == valid). Run it per entry over a built cif
LMDB to QC the database item-by-item against the reference BioMolDB schema.

    from datacooker import execute
    from structcooker.workflows.analysis.cif_validate import RECIPE, TARGETS

    issues = execute(RECIPE, {"record": record}, targets=TARGETS)["issues"]
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.validate import validate_cif_record

validate_recipe = RecipeBook()

validate_recipe.step(
    outputs=("issues", list),
    instruction=validate_cif_record,
    kwargs={"record": ("record", dict)},
)

RECIPE = validate_recipe
TARGETS = ["issues"]
