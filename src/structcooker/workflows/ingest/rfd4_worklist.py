"""Recipe (materialize): an RFD4-Proteina component's index and example tables -> its worklist LMDB and sidecar."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.rfd4_worklist import write_worklist

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=write_worklist,
    kwargs={
        "name": ("name", str),
        "root": ("root", str),
        "out_lmdb": ("out_lmdb", str),
        "out_sidecar": ("out_sidecar", str),
        "out_list": ("out_list", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
