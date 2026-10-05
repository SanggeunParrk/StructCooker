from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.sequence import fasta_new_sequences

"""Materialize the records of one fasta whose sequence another fasta lacks."""

recipe = RecipeBook()

recipe.step(
    outputs=(("fasta_dict", dict),),
    instruction=fasta_new_sequences,
    kwargs={
        "fasta_path": ("fasta_path", Path),
        "base_fasta_path": ("base_fasta_path", Path),
    },
)

RECIPE = recipe
TARGETS = ["fasta_dict"]
