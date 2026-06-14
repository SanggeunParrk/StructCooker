from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import run_hhmake

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.step(
    outputs=(("done_result", str),),
    instruction=run_hhmake,
    kwargs={
        "input_a3m_path": ("input_a3m_path", Path),
        "output_path": ("output_path", Path),
    },
)

RECIPE = recipe
TARGETS = ["done_result"]
