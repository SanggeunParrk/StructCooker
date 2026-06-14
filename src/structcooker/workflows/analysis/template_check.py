from datacooker import RecipeBook

from structcooker.mols import TemplateMol
from structcooker.instructions.transforms.template import length_check
from structcooker.instructions.transforms.cifmol import convert_to_templatemol_dict

"""Build a TemplateMol unittest Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("templatemol_dict", dict[str, TemplateMol]),),
    ],
    instruction=convert_to_templatemol_dict,
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
        (("results", str),),
    ],
    instruction=length_check,
    inputs=[
        {
            "kwargs": {
                "query_seq_id": ("key", str),
                "seqid2seq": ("seqid2seq", dict),
                "templatemol_dict": (
                    "templatemol_dict",
                    dict[str, TemplateMol] | None,
                ),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["results"]
