from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_chains import load_pair_chains

"""The pair table with the chain status of each pair: paired / scfv / vhh / a light chain recovered from AntiBinder / ambiguous /
unresolved (``instructions/readers/abag_chains.py``). Filtering on the status is the consumer's choice."""

abag_pairs_chains_recipe = RecipeBook()

abag_pairs_chains_recipe.step(
    outputs=(("rows", list),),
    instruction=load_pair_chains,
    kwargs={
        "pairs_path": ("pairs_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "antibinder_dir": ("antibinder_dir", Path),
    },
)

RECIPE = abag_pairs_chains_recipe
TARGETS = ["rows"]
