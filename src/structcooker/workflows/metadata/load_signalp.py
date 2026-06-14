from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    build_seqid_map,
    load_signalp,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

metadata_recipe = RecipeBook()

metadata_recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seqid2seq_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)


metadata_recipe.step(
    outputs=(("seqid_map", dict),),
    instruction=build_seqid_map,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
    },
)


metadata_recipe.step(
    outputs=(("signalp_dict", dict),),
    instruction=load_signalp,
    kwargs={
        "signalp_dir": ("signalp_dir", Path | None),
    },
)

RECIPE = metadata_recipe
TARGETS = ["seqid_map", "signalp_dict"]
