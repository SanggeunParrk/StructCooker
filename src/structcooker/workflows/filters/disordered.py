from datetime import date

from datacooker import RecipeBook

from structcooker.instructions.transforms.filtering import (
    filter_by_resolution_and_date,
    filter_cifmol_by_polymer_chain_count,
    filter_cifmol_by_token_count,
    filter_min_unresolved_residues,
    filter_signalp,
)
from structcooker.instructions.transforms.sequence import filter_water
from structcooker.mols import CIFMol

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()


recipe.step(
    outputs=(("cifmol_filtered_by_min_unresolved_residues", CIFMol),),
    instruction=filter_min_unresolved_residues,
    kwargs={
        "cifmol": ("cifmol", CIFMol),
        "min_unresolved_residues": ("min_unresolved_residues", int),
    },
)


recipe.step(
    outputs=(("cifmol_filtered_by_token_count", CIFMol),),
    instruction=filter_cifmol_by_token_count,
    kwargs={
        "cifmol": ("cifmol_filtered_by_min_unresolved_residues", CIFMol),
        "min_token_count": ("min_token_count", int),
        "max_token_count": ("max_token_count", int),
    },
)


recipe.step(
    outputs=(("cifmol_filtered_by_chain_count", CIFMol),),
    instruction=filter_cifmol_by_polymer_chain_count,
    kwargs={
        "cifmol": ("cifmol_filtered_by_token_count", CIFMol),
        "max_polymer_chain_count": ("max_polymer_chain_count", int),
    },
)


recipe.step(
    outputs=(("cifmol_filtered_by_resolution_date", CIFMol),),
    instruction=filter_by_resolution_and_date,
    kwargs={
        "resolution_cutoff": ("resolution_cutoff", float),
        "start_date": ("start_date", date | str),
        "end_date": ("end_date", date | str),
        "cifmol": ("cifmol_filtered_by_chain_count", CIFMol),
    },
)


recipe.step(
    outputs=(("cifmol_wo_water", CIFMol),),
    instruction=filter_water,
    kwargs={
        "cifmol": ("cifmol_filtered_by_resolution_date", CIFMol),
    },
)


recipe.step(
    outputs=(("cifmol_dict", dict),),
    instruction=filter_signalp,
    kwargs={
        "cifmol": ("cifmol_wo_water", CIFMol),
        "seqid_map": ("seqid_map", dict),
        "signalp_dict": ("signalp_dict", dict),
    },
)


RECIPE = recipe
TARGETS = ["cifmol_dict"]
