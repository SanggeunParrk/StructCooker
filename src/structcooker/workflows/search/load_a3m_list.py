from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import load_a3m_list

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.step(
    outputs=(("data_list", list),),
    instruction=load_a3m_list,
    kwargs={
        "data_dir": ("data_dir", Path),
        "output_dir": ("output_dir", Path),
        "output_pattern": ("output_pattern", str),
    },
)


RECIPE = recipe
TARGETS = ["data_list"]
