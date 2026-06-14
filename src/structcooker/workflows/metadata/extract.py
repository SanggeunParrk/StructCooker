from datacooker import RecipeBook

from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict
from structcooker.instructions.transforms.metadata import extract_metadata
from structcooker.mols import CIFMol

"""Extract metadata from cifmol"""

metadata_recipe = RecipeBook()

metadata_recipe.step(
    outputs=(("cifmol_dict", dict[str, dict[str, CIFMol]]),),
    instruction=convert_to_cifmol_dict,
    kwargs={
        "value": ("db_data", dict),
    },
)

metadata_recipe.step(
    outputs=(("metadata_dict", dict[str, dict[str, str]]),),
    instruction=extract_metadata,
    kwargs={
        "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]] | None),
    },
)

RECIPE = metadata_recipe
TARGETS = ["metadata_dict"]
