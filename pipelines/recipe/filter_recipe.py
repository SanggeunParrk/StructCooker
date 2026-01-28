from datetime import date

from datacooker import RecipeBook

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.instructions.filter_instructions import (
    filter_by_resolution_and_date,
    filter_signalp,
    filter_water,
)
from pipelines.instructions.metadata_instructions import (
    attach_metadata,
)

"""Rebuild a CIF lmdb to train AF3"""

af3_recipe = RecipeBook()

af3_recipe.add(
    targets=[(("cifmol_filtered_by_resolution_date", CIFMol),)],
    instruction=filter_by_resolution_and_date(resolution_cutoff=9.0, date_cutoff=date(2021,9,30)),
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol", CIFMol),
            },
        },
    ],
)


af3_recipe.add(
    targets=[(("cifmol_no_water", CIFMol),)],
    instruction=filter_water,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_filtered_by_resolution_date", CIFMol),
            },
        },
    ],
)


af3_recipe.add(
    targets=[(("cifmol_attached", CIFMolAttached),)],
    instruction=attach_metadata,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_no_water", CIFMol),
                "seq2seqID": ("seq2seqID", dict),
                "seqID2clusterID": ("seqID2clusterID", dict),
            },
        },
    ],
)

af3_recipe.add(
    targets=[(("cifmol_dict", dict),)],
    instruction=filter_signalp,
    inputs=[
        {
            "kwargs": {
                "cifmol": ("cifmol_attached", CIFMolAttached),
                "seq2seqID_dict": ("seq2seqID", dict),
                "signalp_dict": ("signalp_dict", dict),
            },
        },
    ],
)



RECIPE = af3_recipe
TARGETS = ["cifmol_dict"]
