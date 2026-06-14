from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict
from structcooker.instructions.transforms.graph import (
    filter_seq_ids,
    interacting_seq_ids,
)
from structcooker.instructions.transforms.metadata import (
    build_seqid_map,
    load_tsv,
)
from structcooker.mols import CIFMol

"""Extract metadata from cifmol"""

recipe = RecipeBook()

recipe.step(
    outputs=(("seqid2seq", dict),),
    instruction=load_tsv,
    kwargs={
        "tsv_file_path": ("seqid2seq_path", Path),
    },
    params={
        "split_by_comma": False,
    },
)


recipe.step(
    outputs=(("seqid_map", dict),),
    instruction=build_seqid_map,
    kwargs={
        "seqid2seq": ("seqid2seq", dict),
    },
)

recipe.step(
    outputs=(("cifmol_dict", dict[str, dict[str, CIFMol]]),),
    instruction=convert_to_cifmol_dict,
    kwargs={
        "value": ("db_data", dict),
    },
)

recipe.step(
    outputs=(("interacting_seq_ids", set[tuple[str, str]]),),
    instruction=interacting_seq_ids,
    kwargs={
        "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]]),
        "seqid_map": ("seqid_map", dict),
    },
)

recipe.step(
    outputs=(("filtered_seq_ids", set[tuple[str, str]]),),
    instruction=filter_seq_ids,
    kwargs={
        "interacting_seq_ids": ("interacting_seq_ids", set[tuple[str, str]]),
    },
    params={
        "valid_entity_types": ["P", "Q", "D", "R", "N"],
    },
)

RECIPE = recipe
TARGETS = ["filtered_seq_ids"]
