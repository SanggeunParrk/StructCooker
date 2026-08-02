from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import cap_msa_dict

"""Lightweight an existing MSA DB by capping each record's depth.

Pure re-cap of a built ``{msa_dict}`` record (no source re-read): keeps the
query + first ``max_depth - 1`` hits and recomputes deletion_mean / profile.
"""


def cap_msa_record(db_data: dict, max_depth: int) -> dict:
    """Cap one MSA record to ``max_depth`` rows, preserving the wrapper."""
    return {"msa_dict": cap_msa_dict(db_data["msa_dict"], max_depth)}


lightweight_msa_recipe = RecipeBook()
lightweight_msa_recipe.step(
    outputs=(("msa_dict_record", dict),),
    instruction=cap_msa_record,
    kwargs={"db_data": ("db_data", dict), "max_depth": ("max_depth", int)},
)
RECIPE = lightweight_msa_recipe
TARGETS = ["msa_dict_record"]
