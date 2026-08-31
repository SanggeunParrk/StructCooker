from datacooker import RecipeBook

from structcooker.instructions.transforms.monomer_seqid import build_monomer_seqid_map

"""Materialize the distillation monomer entry -> seq_id map (feeds openfold_seqid_key)."""

monomer_seqid_recipe = RecipeBook()

monomer_seqid_recipe.step(
    outputs=(("monomer_seqid", dict),),
    instruction=build_monomer_seqid_map,
    kwargs={
        "fasta_paths": ("fasta_paths", list),
        "seq_id_map_path": ("seq_id_map_path", str),
    },
)

RECIPE = monomer_seqid_recipe
TARGETS = ["monomer_seqid"]
