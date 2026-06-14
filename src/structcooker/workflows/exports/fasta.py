from datacooker import RecipeBook

from structcooker.mols import CIFMol
from structcooker.instructions.transforms.sequence import build_fasta
from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict

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
