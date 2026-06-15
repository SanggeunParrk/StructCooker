from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import build_msa_dict

"""Build an OpenFold3 distillation MSA Cooker.

Parses per-entry ``alignment.npz`` files (protein monomers and RNA monomers
share one schema) into an ``msa_dict`` keyed by MSA source name.
"""

msa_recipe = RecipeBook()

msa_recipe.add(
    targets=(("msa_dict", dict),),
    instruction=build_msa_dict(),
    inputs={
        "kwargs": {
            "msa_sources": ("msa_sources", dict),
        },
    },
)

RECIPE = msa_recipe
TARGETS = ["msa_dict"]
