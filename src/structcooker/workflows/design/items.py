from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.design.items import build_items

"""Join statistics and clusters into the items table with the ``pass`` flag (``design/<set>_items``)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("rows", list),),
    instruction=build_items,
    kwargs={
        "stats_dir": ("stats_dir", Path),
        "clusters_path": ("clusters_path", Path),
        "db_path": ("db_path", Path),
        "min_plddt": ("min_plddt", float),
        "min_length": ("min_length", int),
        "max_length": ("max_length", int),
        "max_coil": ("max_coil", float),
        "max_coil_run": ("max_coil_run", int),
        "apply_filters": ("apply_filters", bool),
        "max_failed": ("max_failed", int),
    },
)

RECIPE = recipe
TARGETS = ["rows"]
