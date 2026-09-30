"""Recipe: the sorted files matching several ``<dir>::<glob>`` sources, one per line."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import glob_file_list

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=glob_file_list,
    kwargs={
        "sources": ("sources", list),
        "out_path": ("out_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
