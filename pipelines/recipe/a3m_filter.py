from biomol.core import FeatureContainer
from datacooker import RecipeBook

from pipelines.instructions.filter_instructions import (
    filter_a3m,
)

"""Rebuild a a3m lmdb to train AF3"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (
            ("residue_container", FeatureContainer),
            ("chain_container", FeatureContainer),
        ),
    ],
    instruction=filter_a3m(max_msa_depth=16_384),
    inputs=[
        {
            "kwargs": {
                "residue_container": ("_residue_container", FeatureContainer),
                "chain_container": ("_chain_container", FeatureContainer),
            },
        },
    ],
)


RECIPE = recipe
TARGETS = ["residue_container", "chain_container"]
