"""Recipe: one template link per chain id of the disordered set (see ``disordered_template_links``)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.filelists import disordered_template_links

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=disordered_template_links,
    kwargs={
        "structure_dir": ("structure_dir", str),
        "template_dir": ("template_dir", str),
        "out_dir": ("out_dir", str),
        "out_list": ("out_list", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
