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
        # object, not dict: seq_metadata_map is an LmdbDict (mmap-backed) so workers
        # share one node-local copy instead of each holding a multi-GB dict.
        "seq_metadata_map": ("seq_metadata_map", object),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
