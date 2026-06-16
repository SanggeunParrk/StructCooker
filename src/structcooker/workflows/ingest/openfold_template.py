from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import (
    reconstruct_template_alignments,
)
from structcooker.instructions.transforms.template import load_templates_from_chain_db

"""Build an OpenFold3 distillation template Cooker.

The distillation ``template.npz`` already lists the selected template hits and
their query/template residue mapping (``idx_map``). Those mappings are expressed
as placeholder alignments, then resolved against a prebuilt per-chain cif LMDB
(``cif_chain``) so each hit is a light keyed read instead of re-decoding +
rebuilding every assembly per hit.
"""

template_recipe = RecipeBook()

template_recipe.add(
    targets=(("align_results", dict),),
    instruction=reconstruct_template_alignments,
    inputs={
        "kwargs": {
            "template_hits": ("template_hits", dict),
            "query_len": ("query_len", int),
        },
    },
)

template_recipe.add(
    targets=(("template_mols", dict),),
    instruction=load_templates_from_chain_db,
    inputs={
        "kwargs": {
            "cif_chain_db_path": ("cif_chain_db_path", Path),
            "align_results": ("align_results", dict),
        },
    },
)

RECIPE = template_recipe
TARGETS = ["template_mols"]
