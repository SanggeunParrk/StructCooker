import numpy as np
from datacooker import RecipeBook

from structcooker.instructions.transforms.msa import build_dict, parse_sequence
from structcooker.instructions.transforms.openfold import (
    merge_msa_sources,
    parse_openfold_msa_headers,
    reconstruct_a3m_sequences,
)

"""Build an OpenFold3 distillation MSA Cooker.

Adapts the distillation ``alignment.npz`` (aligned matrix + deletion counts,
one source for proteins / several for RNA) back into raw a3m strings so the
canonical ``parse_sequence`` / ``build_dict`` instructions can be reused. The
resulting ``msa_dict`` matches the schema of the existing a3m LMDB
(``query_sequence`` / ``aligned_sequences`` / ``deletions`` / ``deletion_mean``
/ ``profile`` + parsed headers).
"""

msa_recipe = RecipeBook()

msa_recipe.add(
    targets=(
        ("msa", np.ndarray),
        ("deletion_matrix", np.ndarray),
        ("metadata", np.ndarray),
    ),
    instruction=merge_msa_sources,
    inputs={"kwargs": {"msa_sources": ("msa_sources", dict)}},
)

msa_recipe.add(
    targets=(("raw_sequences", list),),
    instruction=reconstruct_a3m_sequences,
    inputs={
        "kwargs": {
            "msa": ("msa", np.ndarray),
            "deletion_matrix": ("deletion_matrix", np.ndarray),
        },
    },
)

msa_recipe.add(
    targets=(("parsed_sequences", dict),),
    instruction=parse_sequence,
    inputs={
        "kwargs": {
            "raw_sequences": ("raw_sequences", list),
            "a3m_type": ("a3m_type", str | None),
        },
    },
)

msa_recipe.add(
    targets=(("parsed_headers", dict),),
    instruction=parse_openfold_msa_headers,
    inputs={"kwargs": {"metadata": ("metadata", np.ndarray)}},
)

msa_recipe.add(
    targets=(("msa_dict", dict),),
    instruction=build_dict,
    inputs={
        "kwargs": {
            "sequences": ("parsed_sequences", dict),
            "headers": ("parsed_headers", dict),
        },
    },
)

RECIPE = msa_recipe
TARGETS = ["msa_dict"]
