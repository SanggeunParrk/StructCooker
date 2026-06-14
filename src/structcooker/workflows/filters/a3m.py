from datacooker import RecipeBook

from structcooker.instructions.transforms.filtering import (
    filter_a3m,
)

"""Rebuild a a3m lmdb to train AF3"""

recipe = RecipeBook()

recipe.step(
    outputs=(
        ("sequences", dict),
        ("headers", dict),
    ),
    instruction=filter_a3m(max_msa_depth=4096),
    kwargs={
        "sequences": ("_sequences", dict),
        "headers": ("_headers", dict),
    },
)


RECIPE = recipe
TARGETS = ["sequences", "headers"]
