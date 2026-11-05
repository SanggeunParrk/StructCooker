from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_final import load_antigen_worklist

"""The antigens of the kept pairs, as seq_id <TAB> sequence: the work list of the MSA and template searches."""

abag_ag_worklist_recipe = RecipeBook()

abag_ag_worklist_recipe.step(
    outputs=(("seq_id_map", dict),),
    instruction=load_antigen_worklist,
    kwargs={
        "pairs_path": ("pairs_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
    },
)

RECIPE = abag_ag_worklist_recipe
TARGETS = ["seq_id_map"]
