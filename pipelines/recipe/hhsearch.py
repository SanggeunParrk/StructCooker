from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import run_template_search

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("hhsearch_results", dict),),
    ],
    instruction=run_template_search,
    inputs=[
        {
            "kwargs": {
                "msa_path": ("msa_path", Path),
                "hhr_path": ("hhr_path", Path),
                "cpu": ("cpu_per_job", int),
                "mem": ("mem_per_job", int),
                "db_template": ("db_template", str | Path),
                "hhsuite_bin_dir": ("hhsuite_bin_dir", str | Path),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["hhsearch_results"]
