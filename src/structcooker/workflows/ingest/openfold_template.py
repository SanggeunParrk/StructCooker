from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import build_template_dict

"""Build an OpenFold3 distillation template Cooker.

Parses per-entry ``template.npz`` files (protein monomers only) into a
``template_dict`` keyed by ``<pdb_id>_<chain_id>`` template hit.
"""

template_recipe = RecipeBook()

template_recipe.add(
    targets=(("template_dict", dict),),
    instruction=build_template_dict(),
    inputs={
        "kwargs": {
            "template_hits": ("template_hits", dict),
        },
    },
)

RECIPE = template_recipe
TARGETS = ["template_dict"]
