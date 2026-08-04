"""Template Phase 3 seq maps — seqid_to_seq + chain_to_seq (the sequences kalign needs)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.template_candidates import precompute_seqs

recipe = RecipeBook()

recipe.step(
    outputs=(("status", str),),
    instruction=precompute_seqs,
    kwargs={
        "seqid_chains": ("seqid_chains", str),
        "seq_id_map": ("seq_id_map", str),
        "cif_fasta": ("cif_fasta", str),
        "out_seqid_seq": ("out_seqid_seq", str),
        "out_chain_seq": ("out_chain_seq", str),
    },
)

RECIPE = recipe
TARGETS = ["status"]
