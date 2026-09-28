"""Rebuild recipe: cap an existing ``{msa_dict}`` record to ``max_depth`` rows.

Planning-first replacement for the legacy ``scripts/maintenance/lightweight_msa.py``:
run via ``structcooker build`` → ``datacooker pipeline`` (rebuild op) so it inherits
bounded Ray execution and input-size balancing.

The reader adapter ``adapt_msa_for_cap`` renames the deserialized ``msa_dict`` to
``msa_src`` so the step can write its result back as ``msa_dict`` (matching the source
layout, ``TARGETS = ["msa_dict"]`` → ``to_bytes({"msa_dict": ...})``) without an
input/output name collision — a same-name target silently returns the *uncapped* input.
``max_depth`` comes from the config's ``parameters``. Unlike the older ``lightweight_msa``
recipe, this neither double-wraps (``msa_dict_record``) nor no-ops the deep records.
"""
from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold import cap_msa_dict

cap_msa_recipe = RecipeBook()
cap_msa_recipe.step(
    outputs=(("msa_dict", dict),),
    instruction=cap_msa_dict,
    kwargs={
        "msa_dict": ("msa_src", dict),
        "max_depth": ("max_depth", int),
    },
)
RECIPE = cap_msa_recipe
TARGETS = ["msa_dict"]
