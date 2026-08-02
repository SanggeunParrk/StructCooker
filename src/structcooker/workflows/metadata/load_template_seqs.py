from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import load_seq_tsv

"""Lean metadata for Phase 3: just the sequences kalign needs.

query_seqs    : seq_id -> sequence   (only the query seq_ids we build for)
template_seqs : chain  -> sequence   (polypeptide(L) template chains)

Replaces load_template_metadata (~8 GB: full 16.9 M seq_id map + 3.6 M
chain->{date,seq}) which Phase 3 mostly did not use -- it needs no dates and
only the sequences of the seq_ids/templates actually involved.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("query_seqs", dict),),
    instruction=load_seq_tsv,
    kwargs={"path": ("seqid_seq_path", Path)},
)

recipe.step(
    outputs=(("template_seqs", dict),),
    instruction=load_seq_tsv,
    kwargs={"path": ("chain_seq_path", Path)},
)

RECIPE = recipe
TARGETS = ["query_seqs", "template_seqs"]
