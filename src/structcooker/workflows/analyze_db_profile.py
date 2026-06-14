from datacooker import RecipeBook

from structcooker.models import CIFMolAttached
from structcooker.ops.analyze_instructions import analyze_db_profile

"""Build a CIFMol->fasta Cooker."""

recipe = RecipeBook()


recipe.add(
    targets=[
        (("monomer_clusters", set), ("interface_clusters", set)),
    ],
    instruction=analyze_db_profile,
    inputs=[
        {
            "kwargs": {
                "cifmol_dict": ("db_data", dict[str, dict[str, CIFMolAttached]]),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["monomer_clusters", "interface_clusters"]
