"""Recipe (materialize): the Re-Glyco index and example tables -> one worklist file per entry, a list and the sidecar."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.reglyco_worklist import (
    write_reglyco_worklists,
)

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=write_reglyco_worklists,
    kwargs={
        "index_path": ("index_path", str),
        "examples_path": ("examples_path", str),
        "examples_glycan_path": ("examples_glycan_path", str),
        "out_dir": ("out_dir", str),
        "out_list": ("out_list", str),
        "out_sidecar": ("out_sidecar", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
