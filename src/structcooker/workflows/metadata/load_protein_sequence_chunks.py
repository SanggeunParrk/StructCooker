from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    chunk_protein_seqs,
    load_tsv,
)

"""Split a seq_id map into chunks of queries for a batched MMseqs2 search.

The HHblits split (``load_protein_sequences``) makes one work item per sequence, because
HHblits searches one sequence at a time. MMseqs2 pays its cost per search, not per query,
so its work item is a chunk that shares one search.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={"tsv_file_path": ("seq_id_map_path", Path)},
    params={"split_by_comma": False},
)

recipe.step(
    outputs=(("data_list", list),),
    instruction=chunk_protein_seqs,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
        "chunk_size": ("chunk_size", int),
    },
)

RECIPE = recipe
TARGETS = ["data_list"]
