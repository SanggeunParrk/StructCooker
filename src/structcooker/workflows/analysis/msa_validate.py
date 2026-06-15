"""msa LMDB validation analysis recipe.

A Transform-only recipe that takes one decoded msa record and produces its
list of validity issues (empty == valid). Run it per entry over a built msa
LMDB to QC the database item-by-item against the reference BioMolDB schema.

    from datacooker import execute
    from structcooker.workflows.analysis.msa_validate import RECIPE, TARGETS

    issues = execute(RECIPE, {"record": record}, targets=TARGETS)["issues"]
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.validate import validate_msa_record

validate_recipe = RecipeBook()

validate_recipe.step(
    outputs=("issues", list),
    instruction=validate_msa_record,
    kwargs={"record": ("record", dict)},
)

RECIPE = validate_recipe
TARGETS = ["issues"]
