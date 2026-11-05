from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_final import load_seq_id_maps
from structcooker.instructions.readers.sequence import load_fastas
from structcooker.instructions.transforms.sequence import build_seq_id_map

"""``seq_id_map`` with several maps as the seed: the release's map and the ids already handed out for ASD (one id space, no id ever moves)."""

abag_seq_id_map_recipe = RecipeBook()

abag_seq_id_map_recipe.step(
    outputs=(("fasta_dict", dict),),
    instruction=load_fastas,
    kwargs={"fasta_paths": ("fasta_paths", list)},
)

abag_seq_id_map_recipe.step(
    outputs=(("old_seq_id_map", dict),),
    instruction=load_seq_id_maps,
    kwargs={"seq_id_map_paths": ("old_seq_id_map_paths", list)},
)

abag_seq_id_map_recipe.step(
    outputs=(("seq_id_map", dict),),
    instruction=build_seq_id_map,
    kwargs={
        "fasta_dict": ("fasta_dict", dict | None),
        "old_seq_id_map": ("old_seq_id_map", dict | None),
    },
)

RECIPE = abag_seq_id_map_recipe
TARGETS = ["seq_id_map"]
