from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold_structure import build_structure_dict

"""Build an OpenFold3 distillation structure (CIFMol) Cooker.

Reconstructs a unified ``CIFMol`` BioMol payload from each per-entry
``structure.npz`` flat atom table (protein and RNA share one schema), filling
features absent from the distillation output with empty fields.
"""

structure_recipe = RecipeBook()

structure_recipe.add(
    targets=(("cifmol_dict", dict),),
    instruction=build_structure_dict(),
    inputs={
        "kwargs": {
            "atom_arrays": ("atom_arrays", dict),
        },
    },
)

RECIPE = structure_recipe
TARGETS = ["cifmol_dict"]
