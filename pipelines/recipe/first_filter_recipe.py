from datetime import date

from datacooker import RecipeBook

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.instructions.filter_instructions import (
    filter_by_resolution_and_date,
    filter_cifmol_by_chain_count,
    filter_cifmol_by_token_count,
    filter_signalp,
    filter_water,
)
from pipelines.instructions.metadata_instructions import (
    attach_metadata,
)

"""Rebuild a CIF lmdb to train AF3"""

filter_recipe = RecipeBook()

filter_recipe.add(
    targets=[(("cifmol_filtered_by_resolution_date", CIFMol),)],
    instruction=filter_by_resolution_and_date,
    inputs=[
        {
            "kwargs": {
                "resolution_cutoff": ("resolution_cutoff", float),
                "start_date": ("start_date", date|str),
                "end_date": ("end_date", date|str),
                "cifmol": ("cifmol", CIFMol),
            },
        },
    ],
)

filter_recipe.add(
    targets=[(("cifmol_filtered_by_chain_count", CIFMol),)],
    instruction=filter_cifmol_by_chain_count,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_filtered_by_resolution_date", CIFMol),
                "max_chain_count": ("max_chain_count", int),
            },
        },
    ],
)

filter_recipe.add(
    targets=[(("cifmol_filtered_by_token_count", CIFMol),)],
    instruction=filter_cifmol_by_token_count,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_filtered_by_chain_count", CIFMol),
                "max_token_count": ("max_token_count", int),
            },
        },
    ],
)


filter_recipe.add(
    targets=[(("cifmol_wo_water", CIFMol),)],
    instruction=filter_water,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_filtered_by_token_count", CIFMol),
            },
        },
    ],
)


filter_recipe.add(
    targets=[(("cifmol_attached", CIFMolAttached),)],
    instruction=attach_metadata,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_wo_water", CIFMol),
                "seq2seqid": ("seq2seqid", dict),
                "seqid2clusterid": ("seqid2clusterid", dict),
            },
        },
    ],
)

filter_recipe.add(
    targets=[(("cifmol_dict", dict),)],
    instruction=filter_signalp,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_attached", CIFMolAttached),
                "seq2seqid_dict": ("seq2seqid", dict),
                "signalp_dict": ("signalp_dict", dict),
            },
        },
    ],
)



RECIPE = filter_recipe
TARGETS = ["cifmol_dict"]
