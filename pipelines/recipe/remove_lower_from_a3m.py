from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import remove_lower_from_a3m

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("results", str),),
    ],
    instruction=remove_lower_from_a3m,
    inputs=[
        {
            "kwargs": {
                "input_a3m_path": ("input_a3m_path", Path),
                "output_path": ("output_path", Path),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["results"]
