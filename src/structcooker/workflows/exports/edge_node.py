from datacooker import RecipeBook

from structcooker.instructions.transforms.analysis import extract_edge_node
from structcooker.mols import CIFMolAttached

recipe = RecipeBook()


recipe.step(
    outputs=(("monomer_map", dict), ("interface_map", dict)),
    instruction=extract_edge_node,
    kwargs={
        "cifmol_dict": ("db_data", dict[str, dict[str, CIFMolAttached]]),
    },
)

RECIPE = recipe
TARGETS = ["monomer_map", "interface_map"]
