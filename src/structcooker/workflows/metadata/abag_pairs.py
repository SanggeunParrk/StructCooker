from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_pairs import load_abag_pairs

"""The antibody-antigen pairs of ASD, AbRank and POISK that bind with Kd <= 1 uM, as seq_ids (S2, S3, antigen length, the binding filter)."""

abag_pairs_recipe = RecipeBook()

abag_pairs_recipe.step(
    outputs=(("rows", list),),
    instruction=load_abag_pairs,
    kwargs={
        "asd_dir": ("asd_dir", str),
        "abrank_path": ("abrank_path", Path),
        "poisk_path": ("poisk_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "min_antigen": ("min_antigen", int),
        "max_antigen": ("max_antigen", int),
    },
)

RECIPE = abag_pairs_recipe
TARGETS = ["rows"]
