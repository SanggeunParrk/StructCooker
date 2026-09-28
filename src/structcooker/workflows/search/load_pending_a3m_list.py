from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import load_a3m_list

"""List the a3m files whose final product (done_dir/<shard>/<stem><done_suffix>) is missing."""

recipe = RecipeBook()

recipe.step(
    outputs=(("data_list", list),),
    instruction=load_a3m_list,
    kwargs={
        "data_dir": ("data_dir", Path),
        "output_dir": ("output_dir", Path),
        "output_pattern": ("output_pattern", str),
        "done_dir": ("done_dir", Path),
        "done_suffix": ("done_suffix", str),
    },
)


RECIPE = recipe
TARGETS = ["data_list"]
