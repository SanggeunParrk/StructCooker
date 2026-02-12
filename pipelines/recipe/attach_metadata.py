from datacooker import RecipeBook

from pipelines.cifmol import CIFMol
from pipelines.instructions.metadata_instructions import attach_metadata

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.add(
    targets=[(("cifmol_attached_dict", dict),)],
    instruction=attach_metadata,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol", CIFMol),
                "seq_metadata_map": ("seq_metadata_map", dict),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
