from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.readers.abag_final import load_final_fasta

"""The FASTA of the chains of the kept ASD pairs -- ``chain_set: antigen`` or ``antibody`` (heavy and light, recovered light chains included)."""

abag_final_fasta_recipe = RecipeBook()

abag_final_fasta_recipe.step(
    outputs=(("fasta_dict", dict),),
    instruction=load_final_fasta,
    kwargs={
        "pairs_path": ("pairs_path", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "chain_set": ("chain_set", str),
        "antigen_keep_path": ("antigen_keep_path", Path | None),
    },
)

RECIPE = abag_final_fasta_recipe
TARGETS = ["fasta_dict"]
