from datacooker import RecipeBook

from structcooker.instructions.transforms.afm_msa_seqid import build_afm_msa_seqid_map

"""Materialize the AFM MSA entity -> seq_id map (feeds readers.afdb.afm_msa_seqid_key)."""

afm_msa_seqid_recipe = RecipeBook()

afm_msa_seqid_recipe.step(
    outputs=(("afm_msa_seqid", dict),),
    instruction=build_afm_msa_seqid_map,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "seq_id_map_path": ("seq_id_map_path", str),
        "chain_msa_path": ("chain_msa_path", str),
    },
)

RECIPE = afm_msa_seqid_recipe
TARGETS = ["afm_msa_seqid"]
