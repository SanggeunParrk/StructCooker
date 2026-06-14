from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.graph import (
    build_interacting_seq_clusters,
)
from structcooker.instructions.transforms.metadata import (
    load_tsv,
)

"""Extract metadata from cifmol"""

recipe = RecipeBook()

recipe.step(
    outputs=(("interacting_seq_ids", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("interacting_seq_ids_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)

recipe.step(
    outputs=(("seqclusters2seqids", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seqcluster_path", Path),
    },
    params={
        "split_by_comma": True,
    },
)


recipe.step(
    outputs=(("interacting_seq_clusters", dict),),
    instruction=build_interacting_seq_clusters,
    kwargs={
        "interacting_seq_ids": ("interacting_seq_ids", dict),
        "seqclusters2seqids": ("seqclusters2seqids", dict),
    },
)

RECIPE = recipe
TARGETS = ["interacting_seq_clusters"]
