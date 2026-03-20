from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import run_hmmbuild, run_hmmsearch

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("hmmbuild_results", dict),),
    ],
    instruction=run_hmmbuild,
    inputs=[
        {
            "kwargs": {
                "input_a3m_path": ("input_a3m_path", Path),
                "hmm_path": ("output_path", Path),
            },
        },
    ],
)


recipe.add(
    targets=[
        (("hmmsearch_results", dict),),
    ],
    instruction=run_hmmsearch,
    inputs=[
        {
            "kwargs": {
				"output_dir" : ("hmm_output_dir", Path),
                "hmm_path": ("output_path", Path),
                "fasta_path": ("fasta_path", Path),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["hmmsearch_results", "hmmbuild_results"]
