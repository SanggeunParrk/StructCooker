from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import run_hhmake

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.add(
	targets=[
		(
            ("done_result", str),
		),
    ],
	instruction=run_hhmake,
    inputs=[
        {
            "kwargs": {
                "input_a3m_path": ("input_fasta", Path),
                "output_path": ("output_path", Path),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["done_result"]
