"""Recipe: the protein-chain work list (pdbid_chain of polypeptide(L) chains)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import fasta_chain_list

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=fasta_chain_list,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "out_path": ("out_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
