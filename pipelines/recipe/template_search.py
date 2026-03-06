from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_search_instructions import run_template_search_batch

"""Build an MSA->TemplateSearch (HHsearch) Cooker."""

template_search_recipe = RecipeBook()

template_search_recipe.add(
    targets=[
        (("template_results", dict),),
    ],
    instruction=run_template_search_batch,
    inputs=[
        {
            "kwargs": {
                "msa_root": ("msa_root", str | Path),
                "cpu_per_job": ("cpu_per_job", int),
                "mem": ("mem", int),
                "num_jobs": ("num_jobs", int),
                "db_template": ("db_template", str | Path),
                "hhsuite_bin_dir": ("hhsuite_bin_dir", str | Path),
                "msa_file_name": ("msa_file_name", str),
            },
        },
    ],
)

RECIPE = template_search_recipe
TARGETS = ["template_results"]
