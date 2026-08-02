from datacooker import RecipeBook

from structcooker.instructions.transforms.release_patch import inject_release_cad

"""Patch cif_pdb_attached.lmdb (split_entries): inject release_date into each
assembly's metadata. Used with adapter=adapt_cif_attached (renames the stored
cifmol_attached_dict to _cad per assembly) and split_entries=True."""

recipe = RecipeBook()

recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=inject_release_cad,
    kwargs={
        "cad": ("_cad", dict),
        "release_map": ("release_map", dict),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
