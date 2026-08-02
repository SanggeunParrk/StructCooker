from datacooker import RecipeBook

from structcooker.instructions.transforms.release_patch import (
    keep_cad_if_release_in_range,
)

"""Cheap date-scoped subset of cif_pdb_attached_release.lmdb (split_entries):
keep only assemblies whose stored release_date is in [start, end). Uses
adapter=adapt_cif_attached (renames cifmol_attached_dict -> _cad per assembly).
No CIFMol reconstruction -> out-of-range giants are dropped after a cheap
from_bytes, never reconstructed. The resulting subset feeds the train filter so
the expensive reconstruction only ever runs on in-range structures."""

recipe = RecipeBook()

recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=keep_cad_if_release_in_range,
    kwargs={
        "cad": ("_cad", dict),
        "start_date": ("start_date", str | None),
        "end_date": ("end_date", str | None),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
