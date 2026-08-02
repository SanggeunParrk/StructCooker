from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.metadata import (
    invert_seqid_to_chains,
    load_chain_templates,
)

"""Metadata for Phase 4: the two small chain-level maps from Phase 1/2 outputs.

chain2seqid    : chain -> seq_id  (inverted from seqid_to_chains.tsv)
chain2templates: chain -> [template ids]  (from chain_to_templates.tsv, Phase 2)
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("chain2seqid", dict),),
    instruction=invert_seqid_to_chains,
    kwargs={"path": ("seqid_chains_path", Path)},
)

recipe.step(
    outputs=(("chain2templates", dict),),
    instruction=load_chain_templates,
    kwargs={"path": ("chain_templates_path", Path)},
)

RECIPE = recipe
TARGETS = ["chain2seqid", "chain2templates"]
