from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_basic import load_basic_pairs

"""The ASD pairs with the basic filters as columns (duplicates merged, polyprotein, experimental structure of this very antibody, usable chain status) and
the ``keep`` flag: one auditable table (docs/abag_distillation/pair-filters.md)."""

abag_pairs_basic_recipe = RecipeBook()

abag_pairs_basic_recipe.step(
    outputs=(("rows", list),),
    instruction=load_basic_pairs,
    kwargs={
        "pairs_chains_path": ("pairs_chains_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "sabdab_summary_path": ("sabdab_summary_path", Path),
        "pdb_fasta_path": ("pdb_fasta_path", Path),
        "structure_match": ("structure_match", str),
    },
)

RECIPE = abag_pairs_basic_recipe
TARGETS = ["rows"]
