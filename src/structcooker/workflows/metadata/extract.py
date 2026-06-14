from datacooker import RecipeBook

from structcooker.mols import CIFMol
from structcooker.instructions.transforms.metadata import extract_metadata
from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict

"""Extract metadata from cifmol"""

metadata_recipe = RecipeBook()

metadata_recipe.add(
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

metadata_recipe.add(
    targets=[
        (("metadata_dict", dict[str, dict[str, str]]),),
    ],
    instruction=extract_metadata,
    inputs=[
        {
            "kwargs": {
                "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]] | None),
            },
        },
    ],
)

RECIPE = metadata_recipe
TARGETS = ["metadata_dict"]
