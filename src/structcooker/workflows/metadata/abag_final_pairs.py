from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_final import load_final_pairs

"""The final ASD pair table: the kept pairs with a seq_id and a cluster for every chain."""

abag_final_pairs_recipe = RecipeBook()

abag_final_pairs_recipe.step(
    outputs=(("rows", list),),
    instruction=load_final_pairs,
    kwargs={
        "pairs_path": ("pairs_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "antigen_cluster_path": ("antigen_cluster_path", Path),
        "antibody_cluster_paths": ("antibody_cluster_paths", dict),
        "antigen_keep_path": ("antigen_keep_path", Path | None),
    },
)

RECIPE = abag_final_pairs_recipe
TARGETS = ["rows"]
