from datacooker import RecipeBook

from pipelines.cifmol import CIFMol
from pipelines.instructions.seq_instructions import build_fasta
from pipelines.transforms.cifmol_transforms import convert_to_cifmol_dict

"""Build a CIFMol->fasta Cooker."""

fasta_recipe = RecipeBook()

fasta_recipe.add(
    targets=[
        (("cifmol_dict", dict[str, CIFMol]),),
    ],
    instruction=convert_to_cifmol_dict,
    inputs=[
        {
            "kwargs": {
                "value": ("inputs", dict),
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
                "cifmol_dict": ("cifmol_dict", dict[str, CIFMol] | None),
            },
        },
    ],
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
