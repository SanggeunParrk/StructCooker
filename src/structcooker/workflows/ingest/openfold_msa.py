import numpy as np
from datacooker import RecipeBook

from structcooker.instructions.transforms.msa import build_dict
from structcooker.instructions.transforms.openfold import (
    build_msa_features,
    cap_msa_depth,
    merge_msa_sources,
    parse_openfold_msa_headers,
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
        ("msa_all", np.ndarray),
        ("deletion_all", np.ndarray),
        ("metadata_all", np.ndarray),
    ),
    instruction=merge_msa_sources,
    inputs={"kwargs": {"msa_sources": ("msa_sources", dict)}},
)

# Cap deep alignments (e.g. long monomers) to bound memory; no-op when max_depth is null.
msa_recipe.add(
    targets=(
        ("msa", np.ndarray),
        ("deletion_matrix", np.ndarray),
        ("metadata", np.ndarray),
    ),
    instruction=cap_msa_depth,
    inputs={
        "kwargs": {
            "msa": ("msa_all", np.ndarray),
            "deletion_matrix": ("deletion_all", np.ndarray),
            "metadata": ("metadata_all", np.ndarray),
            "max_depth": ("max_depth", int | None),
        },
    },
)

msa_recipe.add(
    targets=(("parsed_sequences", dict),),
    instruction=build_msa_features,
    inputs={
        "kwargs": {
            "msa": ("msa", np.ndarray),
            "deletion_matrix": ("deletion_matrix", np.ndarray),
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
