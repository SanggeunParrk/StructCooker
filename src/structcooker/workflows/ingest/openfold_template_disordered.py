from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold_structure import (
    build_disordered_template_mols,
)

"""Build the disordered template DB: per query, ``{template_mols: {hit: CIFMol}}``.

Disordered templates ship as per-chain atom tables (``<id>/<id>_<chain>.npz``);
this groups a query's chains into one ``template_mols`` record so the layout
matches the long/short/PDB template DBs.
"""

disordered_template_recipe = RecipeBook()
disordered_template_recipe.step(
    outputs=(("template_mols", dict),),
    instruction=build_disordered_template_mols,
    kwargs={
        "templates_atom_arrays": ("templates_atom_arrays", dict),
        "ccd_db_path": ("ccd_db_path", Path),
    },
)
RECIPE = disordered_template_recipe
TARGETS = ["template_mols"]
