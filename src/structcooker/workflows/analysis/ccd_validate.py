"""CCD validation analysis recipe.

A Transform-only recipe that takes one parsed ``chem_comp`` entry and produces
its list of validity issues (empty == valid). Run it per entry over a built CCD
LMDB to QC the database item-by-item — this is post-processing/analysis, not a
regression diff against a previous build.

    from datacooker import execute
    from structcooker.workflows.analysis.ccd_validate import RECIPE, TARGETS

    issues = execute(RECIPE, {"chem_comp": entry}, targets=TARGETS)["issues"]
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.ccd import validate_chem_comp

validate_recipe = RecipeBook()

validate_recipe.step(
    outputs=("issues", list),
    instruction=validate_chem_comp,
    kwargs={"chem_comp": ("chem_comp", dict | None)},
)

RECIPE = validate_recipe
TARGETS = ["issues"]
