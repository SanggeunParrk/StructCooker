"""Recipe: date split of the interface units with leakage flags (``split.tsv``, ``split_summary.md``)."""

from datacooker import RecipeBook

from structcooker.instructions.interface.tables import write_split

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=write_split,
    kwargs={
        "manifest_full_path": ("manifest_full_path", str),
        "out_path": ("out_path", str),
        "out_summary": ("out_summary", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
