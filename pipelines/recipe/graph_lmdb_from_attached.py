from datacooker import RecipeBook

from pipelines.cifmol import CIFMol
from pipelines.instructions.graph_cluster_instructions import (
    extract_graph_per_cifmol_attached,
)

"""Build a CIFMol->fasta Cooker."""

gc_recipe = RecipeBook()

gc_recipe.add(
    targets=[
        (("graph_bytes", dict[str, bytes]),),
    ],
    instruction=extract_graph_per_cifmol_attached(),
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol", list[CIFMol] | None),
            },
        },
    ],
)


RECIPE = gc_recipe
TARGETS = ["graph_bytes"]
