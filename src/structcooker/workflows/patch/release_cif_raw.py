from datacooker import RecipeBook

from structcooker.instructions.transforms.release_patch import (
    inject_release_metadata_dict,
    passthrough,
)

"""Patch cif_pdb.lmdb: inject release_date into metadata_dict, pass assembly_dict
through unchanged. Used with adapter=adapt_cif_raw (renames raw keys to _*)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("assembly_dict", dict),),
    instruction=passthrough,
    kwargs={"value": ("_assembly_dict", dict)},
)

recipe.step(
    outputs=(("metadata_dict", dict),),
    instruction=inject_release_metadata_dict,
    kwargs={
        "metadata_dict": ("_metadata_dict", dict),
        "release_map": ("release_map", dict),
    },
)

RECIPE = recipe
TARGETS = ["assembly_dict", "metadata_dict"]
