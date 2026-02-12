from datacooker import RecipeBook

from pipelines.cifmol import CIFMolAttached
from pipelines.instructions.filter_instructions import (
    filter_cifmol_by_clusters,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.add(
    targets=[
        (("cifmol_attached_dict", dict),),
    ],
    instruction=filter_cifmol_by_clusters,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol", CIFMolAttached),
                "filtered_clusters": ("filtered_valid_2_clusters", set),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
