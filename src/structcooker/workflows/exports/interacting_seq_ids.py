"""Extract interacting sequence IDs using the already-attached chain contact graph."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.graph import (
    filter_seq_ids,
    interacting_seq_ids_from_attached,
)

recipe = RecipeBook()
recipe.step(
    outputs=(("interacting_seq_ids", set[tuple[str, str]]),),
    instruction=interacting_seq_ids_from_attached,
    kwargs={"record": ("db_data", dict)},
)
recipe.step(
    outputs=(("filtered_seq_ids", set[tuple[str, str]]),),
    instruction=filter_seq_ids,
    kwargs={"interacting_seq_ids": ("interacting_seq_ids", set[tuple[str, str]])},
    params={"valid_entity_types": {"P", "Q", "D", "R", "N"}},
)
RECIPE = recipe
TARGETS = ["filtered_seq_ids"]
