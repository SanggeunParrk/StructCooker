from datacooker import RecipeBook

from pipelines.instructions.seq_instructions import build_seq_id_map

"""Build a CIFMol->fasta Cooker."""

hash_map_recipe = RecipeBook()

hash_map_recipe.add(
    targets=[
        (("seq_id_map", dict),),
    ],
    instruction=build_seq_id_map(),
    inputs=[
        {
            "kwargs": {
                "fasta_dict": ("fasta_dict", dict | None),
            },
        },
    ],
)

RECIPE = hash_map_recipe
TARGETS = ["seq_id_map"]
