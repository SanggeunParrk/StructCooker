"""Load the compact signal-peptide sequence index once per workflow run."""
from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    index_signalp_sequences,
    load_signalp,
)

metadata_recipe = RecipeBook()
metadata_recipe.step(
    outputs=(("signalp_dict", dict),), instruction=load_signalp,
    kwargs={"signalp_dir": ("signalp_dir", Path | None)},
)
metadata_recipe.step(
    outputs=(("signalp_by_sequence", dict),), instruction=index_signalp_sequences,
    kwargs={"seqid2seq_path": ("seqid2seq_path", Path), "signalp_dict": ("signalp_dict", dict)},
)
RECIPE = metadata_recipe
TARGETS = ["signalp_by_sequence"]
