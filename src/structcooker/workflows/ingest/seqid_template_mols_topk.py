from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import build_seqid_template_mols

"""Phase 3 with a final top-k (max_keep): seq_id -> template mols LMDB (+ template_ids).

For a set with no per-chain Phase 4 (teddymer): keeps the first max_keep hits, in
e-value order, that pass the coverage filter. Otherwise as seqid_template_mols.py.


Driven per seq_id. The reduced hmm (Phase 2) already holds the union of the
date-filtered hits selected by any chain of this seq_id, so the expensive
kalign + template-mol build runs ONCE per seq_id (deduped ~5.7x vs per chain).
Templates are read from the per-chain cif LMDB (memory-safe). Uses lean seq maps
(query_seqs, template_seqs) so metadata is small and n_jobs can be high.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("template_mols", dict), ("template_ids", list)),
    instruction=build_seqid_template_mols,
    kwargs={
        "file_path": ("file_path", Path),
        "query_seqs": ("query_seqs", dict),
        "template_seqs": ("template_seqs", dict),
        "reduced_hmm_dir": ("reduced_hmm_dir", Path),
        "cif_chain_db_path": ("cif_chain_db_path", Path),
        "max_keep": ("max_keep", int),
        "allow_missing_chains": ("allow_missing_chains", bool),
    },
)

RECIPE = recipe
TARGETS = ["template_mols", "template_ids"]
