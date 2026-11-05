from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag import load_abag_fasta

"""The distinct chain sequences of ASD, AbRank and POISK as one FASTA -- ``chain_set: antigen`` or ``antibody`` (heavy and light). One FASTA per kind of chain: a
seq_id is assigned from them and the clusters of the final chains come from them."""

abag_fasta_recipe = RecipeBook()

abag_fasta_recipe.step(
    outputs=(("fasta_dict", dict),),
    instruction=load_abag_fasta,
    kwargs={
        "asd_dir": ("asd_dir", str),
        "abrank_path": ("abrank_path", Path),
        "poisk_path": ("poisk_path", Path),
        "chain_set": ("chain_set", str),
    },
)

RECIPE = abag_fasta_recipe
TARGETS = ["fasta_dict"]
