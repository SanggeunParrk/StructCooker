from datacooker import RecipeBook

from structcooker.models import CIFMol
from structcooker.ops.seq_instructions import build_fasta
from structcooker.readers.cifmol import convert_to_cifmol_dict

"""Build a CIFMol->fasta Cooker."""

fasta_recipe = RecipeBook()

fasta_recipe.add(
    targets=[
        (("cifmol_dict", dict[str, dict[str, CIFMol]]),),
    ],
    instruction=convert_to_cifmol_dict,
    inputs=[
        {
            "kwargs": {
                "value": ("db_data", dict),
            },
        },
    ],
)

fasta_recipe.add(
    targets=[
        (("fasta", str),),
    ],
    instruction=build_fasta,
    inputs=[
        {
            "kwargs": {
                "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]] | None),
            },
        },
    ],
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
