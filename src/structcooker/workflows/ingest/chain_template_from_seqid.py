from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import build_chain_from_seqid_db

"""Phase 4: chain-keyed template DB by cheap lookup into the Phase 3 seq_id DB.

Each chain's per-chain template ids (Phase 2) are looked up in the seq_id union
mol DB (Phase 3) and copied out -- no kalign, no CIF decode. Produces the final
``{pdbid}_{chain} -> {template_id: mol}`` template DB.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("template_mols", dict),),
    instruction=build_chain_from_seqid_db,
    kwargs={
        "file_path": ("file_path", Path),
        "chain2seqid": ("chain2seqid", dict),
        "chain2templates": ("chain2templates", dict),
        "seqid_template_db_path": ("seqid_template_db_path", Path),
    },
)

RECIPE = recipe
TARGETS = ["template_mols"]
