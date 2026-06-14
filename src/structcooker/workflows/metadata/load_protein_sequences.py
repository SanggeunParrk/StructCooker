from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    extract_protein_seqs,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seq_id_map_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)


recipe.step(
    outputs=(("data_list", list),),
    instruction=extract_protein_seqs,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
    },
)


RECIPE = recipe
TARGETS = ["data_list"]
