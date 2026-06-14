from pathlib import Path

from datacooker import RecipeBook

from structcooker.ops.graph_instructions import (
    build_interacting_seq_clusters,
)
from structcooker.ops.metadata_instructions import (
    load_tsv,
)

"""Extract metadata from cifmol"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("interacting_seq_ids", dict),),
        (("seqclusters2seqids", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
                "tsv_file_path": ("interacting_seq_ids_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
        {
            "kwargs": {
                "tsv_file_path": ("seqcluster_path", Path),
            },
            "params": {
                "split_by_comma": True,
            },
        },
    ],
)


recipe.add(
    targets=[
        (("interacting_seq_clusters", dict),),  # seq+moltype -> seqid
    ],
    instruction=build_interacting_seq_clusters,
    inputs=[
        {
            "kwargs": {
                "interacting_seq_ids": ("interacting_seq_ids", dict),
                "seqclusters2seqids": ("seqclusters2seqids", dict),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["interacting_seq_clusters"]
