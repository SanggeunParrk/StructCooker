from pathlib import Path

from datacooker import RecipeBook

from pipelines.cifmol import CIFMol
from pipelines.instructions.graph_instructions import (
    filter_seq_ids,
    interacting_seq_ids,
)
from pipelines.instructions.metadata_instructions import (
    build_seqid_map,
    load_tsv,
)
from pipelines.transforms.cifmol_transforms import convert_to_cifmol_dict

"""Extract metadata from cifmol"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("seqid2seq", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
                "tsv_file_path": ("seqid2seq_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
    ],
)


recipe.add(
    targets=[
        (("seqid_map", dict),),  # seq+moltype -> seqid
    ],
    instruction=build_seqid_map,
    inputs=[
        {
            "kwargs": {
                "seqid2seq": ("seqid2seq", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("cifmol_dict", dict[str, dict[str, CIFMol]]),),
    ],
    instruction=convert_to_cifmol_dict,
    inputs=[
        {
            "kwargs": {
                "value": ("db_data", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("interacting_seq_ids", set[tuple[str, str]]),),
    ],
    instruction=interacting_seq_ids,
    inputs=[
        {
            "kwargs": {
                "cifmol_dict": ("cifmol_dict", dict[str, dict[str, CIFMol]]),
                "seqid_map": ("seqid_map", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("filtered_seq_ids", set[tuple[str, str]]),),
    ],
    instruction=filter_seq_ids,
    inputs=[
        {
            "kwargs": {
                "interacting_seq_ids": ("interacting_seq_ids", set[tuple[str, str]]),
            },
            "params": {
                "valid_entity_types": ["P", "Q", "D", "R", "N"],
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["filtered_seq_ids"]
