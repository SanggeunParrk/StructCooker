from datacooker import RecipeBook

from structcooker.mols import CIFMolAttached
from structcooker.instructions.transforms.analysis import extract_edge_node

recipe = RecipeBook()


recipe.add(
    targets=[
        (("monomer_map", dict), ("interface_map", dict)),
    ],
    instruction=extract_edge_node,
    inputs=[
        {
            "kwargs": {
                "cifmol_dict": ("db_data", dict[str, dict[str, CIFMolAttached]]),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["monomer_map", "interface_map"]
