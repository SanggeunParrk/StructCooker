from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.ccd import split_each_cif_files

"""Build a CIFMol->fasta Cooker."""

recipe = RecipeBook()

recipe.step(
    outputs=(("each_cif_lines", dict[str, str]),),
    instruction=split_each_cif_files,
    kwargs={
        "cif_path": ("cif_path", Path),
    },
)


RECIPE = recipe
TARGETS = ["each_cif_lines"]
