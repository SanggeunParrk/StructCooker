from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.metadata_instructions import (
    build_seqid_map,
    load_signalp,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

metadata_recipe = RecipeBook()

metadata_recipe.add(
    targets=[
        (("seqid2seq", dict),),
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
    ],
)


metadata_recipe.add(
    targets=[
        (("seqid_map", dict),),  # seq+moltype -> seqid
    ],
    instruction=build_seqid_map,
    inputs=[
        {
            "kwargs": {
                "seqid2seq": ("seqid2seq", dict),
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
TARGETS = ["seqid_map", "signalp_dict"]
