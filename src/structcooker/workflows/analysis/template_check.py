from datacooker import RecipeBook

from structcooker.instructions.transforms.cifmol import convert_to_templatemol_dict
from structcooker.instructions.transforms.template import length_check
from structcooker.mols import TemplateMol

"""Build a TemplateMol unittest Cooker."""

recipe = RecipeBook()

recipe.step(
    outputs=(("templatemol_dict", dict[str, TemplateMol]),),
    instruction=convert_to_templatemol_dict,
    kwargs={
        "value": ("db_data", dict),
    },
)

recipe.step(
    outputs=(("results", str),),
    instruction=length_check,
    kwargs={
        "query_seq_id": ("key", str),
        "seqid2seq": ("seqid2seq", dict),
        "templatemol_dict": (
            "templatemol_dict",
            dict[str, TemplateMol] | None,
        ),
    },
)

RECIPE = recipe
TARGETS = ["results"]
