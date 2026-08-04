from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import tsv_key_list

"""Small list/fasta projection for the template pipeline."""

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
