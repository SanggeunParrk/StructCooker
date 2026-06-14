from datacooker import RecipeBook

from structcooker.instructions.transforms.analysis import analyze_db_profile
from structcooker.mols import CIFMolAttached

"""Build a CIFMol->fasta Cooker."""

recipe = RecipeBook()


recipe.step(
    outputs=(("monomer_clusters", set), ("interface_clusters", set)),
    instruction=analyze_db_profile,
    kwargs={
        "cifmol_dict": ("db_data", dict[str, dict[str, CIFMolAttached]]),
    },
)

RECIPE = recipe
TARGETS = ["monomer_clusters", "interface_clusters"]
