from datacooker import RecipeBook

from structcooker.instructions.transforms.template_candidates import (
    precompute_candidates,
)

"""Template Phase 1/2 — seq_id -> chains and chain -> <=topk templates (+ reduced hmm).

Materialize recipe: the instruction reads the hmm outputs + fasta + maps and writes
seqid_to_chains / chain_to_templates / reduced_hmm as side effects, returning a status
(output_data_path = chain_to_templates.tsv is the build-done marker).
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("status", str),),
    instruction=precompute_candidates,
    kwargs={
        "cif_fasta": ("cif_fasta", str),
        "chain_list": ("chain_list", str),
        "seq_id_map": ("seq_id_map", str),
        "pdb_dates": ("pdb_dates", str),
        "hmm_dir": ("hmm_dir", str),
        "out_seqid_chains": ("out_seqid_chains", str),
        "out_chain_templates": ("out_chain_templates", str),
        "out_reduced_hmm_dir": ("out_reduced_hmm_dir", str),
        "date_cutoff": ("date_cutoff", str),
        "day_diff": ("day_diff", int),
        "topk": ("topk", int),
        "n_jobs": ("n_jobs", int),
    },
)

RECIPE = recipe
TARGETS = ["status"]
