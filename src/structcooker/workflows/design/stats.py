from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.design.stats import write_part

"""Statistics and backbone archive of one part of the source DB (the work recipe of ``design/<set>_stats``)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("status", str),),
    instruction=write_part,
    kwargs={
        "db_path": ("db_path", Path),
        "output_dir": ("output_dir", Path),
        "backbone_dir": ("out_backbone_dir", Path),
        "part": ("part", str),
        "keys": ("keys", list),
    },
)

RECIPE = recipe
TARGETS = ["status"]
