from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import load_a3m_list

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("data_list", list),),
    ],
    instruction=load_a3m_list,
    inputs=[
        {
            "kwargs": {
				"data_dir": ("data_dir", Path),
				"output_dir": ("output_dir", Path),
				"output_pattern": ("output_pattern", str),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["data_list"]
