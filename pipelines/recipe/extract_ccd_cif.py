from datacooker import RecipeBook

from pathlib import Path
from pipelines.instructions.ccd_instructions import split_each_cif_files

"""Build a CIFMol->fasta Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("each_cif_lines", dict[str, str]),),
    ],
    instruction=split_each_cif_files,
    inputs=[
        {
            "kwargs": {
                "cif_path": ("cif_path", Path),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["each_cif_lines"]
