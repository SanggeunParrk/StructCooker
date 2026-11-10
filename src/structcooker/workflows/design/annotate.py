from datacooker import RecipeBook

from structcooker.instructions.design.items import annotate_item
from structcooker.instructions.mols import CIFMolAttached

"""Keep the items that pass the filters and add their statistics and cluster to the metadata (``design/<set>``)."""

recipe = RecipeBook()

recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=annotate_item,
    kwargs={
        "cifmol": ("cifmol", CIFMolAttached),
        "items": ("design_items", object),
        "source": ("source", str),
    },
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
