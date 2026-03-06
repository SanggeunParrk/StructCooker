from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.metadata_instructions import (
    extract_protein_seqs,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("seqid2seq", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
				"tsv_file_path": ("seq_id_map_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
    ],
)


recipe.add(
    targets=[
        (("data_list", list),),
    ],
    instruction=extract_protein_seqs,
    inputs=[
        {
            "kwargs": {
				"seqid2seq": ("seqid2seq", dict),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["data_list"]
