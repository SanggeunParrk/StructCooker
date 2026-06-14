from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import attach_metadata
from structcooker.mols import CIFMol

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=attach_metadata,
    kwargs={
        "cifmol": ("cifmol", CIFMol),
        "seq_metadata_map": ("seq_metadata_map", dict),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
