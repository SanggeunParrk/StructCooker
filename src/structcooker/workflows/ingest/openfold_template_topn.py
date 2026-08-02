from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import (
    reconstruct_template_alignments,
)
from structcooker.instructions.transforms.template import (
    load_templates_from_chain_db,
    rank_template_hits_by_coverage,
)

"""Build a top-N OpenFold3 distillation template Cooker.

Same light path as :mod:`openfold_template`, but first orders each target's
template hits by query coverage (cheap -- ``idx_map`` only, no decode), then
decodes them best-first and stops at ``max_keep`` successes. Cuts both stored
templates per target (~167 -> max_keep) and build cost.
"""

template_topn_recipe = RecipeBook()

template_topn_recipe.add(
    targets=(("ranked_hits", dict),),
    instruction=rank_template_hits_by_coverage,
    inputs={
        "kwargs": {
            "template_hits": ("template_hits", dict),
            "query_len": ("query_len", int),
            "min_coverage": ("min_coverage", float),
        },
    },
)

template_topn_recipe.add(
    targets=(("align_results", dict),),
    instruction=reconstruct_template_alignments,
    inputs={
        "kwargs": {
            "template_hits": ("ranked_hits", dict),
            "query_len": ("query_len", int),
        },
    },
)

template_topn_recipe.add(
    targets=(("template_mols", dict),),
    instruction=load_templates_from_chain_db,
    inputs={
        "kwargs": {
            "cif_chain_db_path": ("cif_chain_db_path", Path),
            "align_results": ("align_results", dict),
            "max_keep": ("max_keep", int),
        },
    },
)

RECIPE = template_topn_recipe
TARGETS = ["template_mols"]
