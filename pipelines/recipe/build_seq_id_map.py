from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.convert import load_fasta, load_seq_id_map
from pipelines.instructions.seq_instructions import build_seq_id_map

"""Build a CIFMol->fasta Cooker."""

hash_map_recipe = RecipeBook()

hash_map_recipe.add(
    targets=[
        (("fasta_dict", dict),),
    ],
    instruction=load_fasta,
    inputs=[
        {
            "kwargs": {
                "fasta_path": ("fasta_path", str | Path),
            },
        },
    ],
)

hash_map_recipe.add(
    targets=[
        (("old_seq_id_map", dict),),
    ],
    instruction=load_seq_id_map,
    inputs=[
        {
            "kwargs": {
                "seq_id_map_path": ("old_seq_id_map_path", str | Path),
            },
        },
    ],
)

hash_map_recipe.add(
    targets=[
        (("seq_id_map", dict),),
    ],
    instruction=build_seq_id_map,
    inputs=[
        {
            "kwargs": {
                "fasta_dict": ("fasta_dict", dict | None),
                "old_seq_id_map": ("old_seq_id_map", dict | None),
            },
        },
    ],
)

RECIPE = hash_map_recipe
TARGETS = ["seq_id_map"]
