from datacooker import RecipeBook

from structcooker.instructions.transforms.filtering import (
    filter_cifmol_by_clusters,
)
from structcooker.mols import CIFMolAttached

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=filter_cifmol_by_clusters,
    kwargs={
        "cifmol": ("cifmol", CIFMolAttached),
        "filtered_clusters": ("filtered_valid_2_clusters", set),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
