from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.design.cluster import foldseek_clusters

"""Cluster the backbone archives of a design set with foldseek (``design/<set>_clusters``)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("rows", list),),
    instruction=foldseek_clusters,
    kwargs={
        "backbone_dir": ("backbone_dir", Path),
        "work_dir": ("work_dir", Path),
        "coverage": ("coverage", float),
        "threads": ("n_jobs", int),
        "split_memory_limit": ("split_memory_limit", str),
    },
)

RECIPE = recipe
TARGETS = ["rows"]
