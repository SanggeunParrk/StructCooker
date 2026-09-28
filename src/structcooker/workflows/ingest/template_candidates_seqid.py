"""Template Phase 1/2 for sets without query release dates -- per seq_id candidates.

Materialize recipe: writes seqid_to_templates / the Phase 3 work list / chain_to_seq /
reduced hmm as side effects (output_data_path = seqid_to_templates.tsv is the marker).
"""
from datacooker import RecipeBook

from structcooker.instructions.transforms.template_candidates import (
    precompute_seqid_candidates,
)

recipe = RecipeBook()

recipe.step(
    outputs=(("status", str),),
    instruction=precompute_seqid_candidates,
    kwargs={
        "seq_ids_path": ("seq_ids_path", str),
        "pdb_dates": ("pdb_dates", str),
        "hmm_dir": ("hmm_dir", str),
        "cif_fasta": ("cif_fasta", str),
        "out_seqid_templates": ("out_seqid_templates", str),
        "out_seqid_list": ("out_seqid_list", str),
        "out_chain_seq": ("out_chain_seq", str),
        "out_reduced_hmm_dir": ("out_reduced_hmm_dir", str),
        "date_cutoff": ("date_cutoff", str),
        "max_candidates": ("max_candidates", int),
        "n_jobs": ("n_jobs", int),
    },
)

RECIPE = recipe
TARGETS = ["status"]
