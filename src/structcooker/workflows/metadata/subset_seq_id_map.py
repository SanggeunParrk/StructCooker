from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.sequence import subset_seq_id_map

"""Materialize one DB's rows of the shared seq_id_map (its per-DB work list)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("seq_id_map", dict),),
    instruction=subset_seq_id_map,
    kwargs={
        "fasta_path": ("fasta_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
    },
)

RECIPE = recipe
TARGETS = ["seq_id_map"]
