from datacooker import RecipeBook

from pipelines.cifmol import CIFMolAttached
from pipelines.instructions.filter_instructions import (
    filter_cifmol_by_clusters,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.add(
    targets=[
        (("filtered_cifmol_dict", dict),),
    ],
    instruction=filter_cifmol_by_clusters,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol", CIFMolAttached),
                "filtered_clusters": ("filtered_valid_2_clusters", set),
                "raw_fasta_dict": ("raw_fasta_dict", dict),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["filtered_cifmol_dict"]
