"""Recipe: the polypeptide(L) subset of a fasta (the hmmsearch template DB)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import filter_fasta_polypeptide_l

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=filter_fasta_polypeptide_l,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "out_path": ("out_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
