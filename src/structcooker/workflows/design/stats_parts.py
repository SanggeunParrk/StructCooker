from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.design.stats import pending_parts

"""List the parts of the source DB whose statistics are not written yet (the split recipe of ``design/<set>_stats``)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("data_list", list),),
    instruction=pending_parts,
    kwargs={
        "db_path": ("db_path", Path),
        "output_dir": ("output_dir", Path),
        "chunk_keys": ("chunk_keys", int),
    },
)

RECIPE = recipe
TARGETS = ["data_list"]
