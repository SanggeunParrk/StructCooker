from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.template import build_chain_template

"""Build the PDB template DB keyed by chain_id ({pdbid}_{chain}).

Driven per chain (``file_path`` is a ``{pdbid}_{chain}`` work item): the chain's
seq_id hmmsearch hits are filtered against THIS chain's own deposition date, so
templates keep their per-structure date instead of collapsing to the sequence's
earliest occurrence (which seq_id keying forced). The standard ``lmdb build``
handles the chunked LMDB writes and sharding; ``load_template_metadata`` supplies
the maps via ``metadata_recipe``.
"""

recipe = RecipeBook()

recipe.step(
    outputs=(("template_mols", dict),),
    instruction=build_chain_template,
    kwargs={
        "file_path": ("file_path", Path),
        "template_metadata_map": ("template_metadata_map", dict),
        "chain2seqid": ("chain2seqid", dict),
        "seqid2earliest_date": ("seqid2earliest_date", dict),
        "filtered_seqid2seq": ("filtered_seqid2seq", dict),
        "hmm_dir": ("hmm_dir", Path),
        "cif_db_path": ("cif_db_path", Path),
    },
)

RECIPE = recipe
TARGETS = ["template_mols"]
