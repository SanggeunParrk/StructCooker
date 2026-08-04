"""Recipe: column-0 keys of a TSV as a one-per-line list (e.g. seq_ids)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import tsv_key_list

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=tsv_key_list,
    kwargs={
        "tsv_path": ("tsv_path", str),
        "out_path": ("out_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
