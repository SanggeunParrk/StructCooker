from datacooker import RecipeBook

from structcooker.instructions.readers.sequence import load_fastas

"""Concatenate several FASTAs into one -- a DB whose structures come from more than one
source still gets exactly one FASTA (docs/seq-id-and-cluster-scheme.md)."""

concat_fasta_recipe = RecipeBook()

concat_fasta_recipe.step(
    outputs=(("fasta_dict", dict),),
    instruction=load_fastas,
    kwargs={"fasta_paths": ("fasta_paths", list)},
)

RECIPE = concat_fasta_recipe
TARGETS = ["fasta_dict"]
