"""Rebuild recipe: ``chain/cif_chain`` -- one per-chain BioMol dict per key.

Planning-first replacement for ``scripts/maintenance/build_cif_chain.py``. A 1:N
explode of the base cif DB: the reader adapter
(``template.adapt_cif_record_to_chains``) splits each PDB record into
``{base_chain: {biomoldict}}`` sub-entries (best-occupancy assembly per chain), and
this recipe passes each chain's ``biomoldict`` through so explode writes it under
``<pdbid>_<base_chain>`` (Schema C). Template ingest then resolves a hit with a light
keyed read instead of re-decoding and rebuilding every assembly.

The chain extraction lives in the adapter, not here, because ``split_entries``
requires the deserialized entry to already be a mapping of sub-entries; this step is
therefore a passthrough of the single per-chain value.
"""
from datacooker import RecipeBook

recipe = RecipeBook()


def _passthrough(biomoldict: dict) -> dict:
    return biomoldict


recipe.step(
    outputs=(("biomoldict", dict),),
    instruction=_passthrough,
    kwargs={"biomoldict": ("biomoldict", dict)},
)

RECIPE = recipe
TARGETS = ["biomoldict"]
