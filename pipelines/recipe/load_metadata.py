from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.metadata_instructions import (
    load_signalp,
    load_tsv,
    reverse_dict,
)

"""Rebuild a CIF lmdb to train AF3"""

metadata_recipe = RecipeBook()


metadata_recipe.add(
    targets=[
        (("seqid2seq", dict),),
        (("clusterid2seqid", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
                "tsv_file_path": ("seqid2seq_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
        {
            "kwargs": {
                "tsv_file_path": ("clusterid2seqid_path", Path),
            },
            "params": {
                "split_by_comma": True,
            },
        },
    ],
)

metadata_recipe.add(
    targets=[
        (("seq2seqid", dict),),
        (("seqid2clusterid", dict),),
    ],
    instruction=reverse_dict,
    inputs=[
        {
            "kwargs": {
                "input_dict": ("seqid2seq", dict),
            },
        },
        {
            "kwargs": {
                "input_dict": ("clusterid2seqid", dict),
            },
        },
    ],
)


metadata_recipe.add(
    targets=[(("signalp_dict", dict),)],
    instruction=load_signalp,
    inputs=[
        {
            "kwargs": {
                "signalp_dir": ("signalp_dir", Path | None),
            },
        },
    ],
)

RECIPE = metadata_recipe
TARGETS = ["seq2seqid", "seqid2clusterid", "signalp_dict"]

