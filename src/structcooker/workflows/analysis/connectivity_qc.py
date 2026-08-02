"""Connectivity QC analysis recipe.

Transform-only recipe: one decoded cif record + the CCD linker reference path
produce a per-entry connectivity annotation (broken / missing / clash /
chain_mismatch with locations). Empty dict == clean. Run per entry over a built
cif LMDB.

    from datacooker import execute
    from structcooker.workflows.analysis.connectivity_qc import RECIPE, TARGETS

    qc = execute(RECIPE, {"record": rec, "ccd_linker_db_path": path}, targets=TARGETS)["qc"]
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.connectivity import connectivity_qc_record

connectivity_qc_recipe = RecipeBook()

connectivity_qc_recipe.step(
    outputs=("qc", dict),
    instruction=connectivity_qc_record,
    kwargs={
        "record": ("record", dict),
        "ccd_linker_db_path": ("ccd_linker_db_path", str),
    },
)

RECIPE = connectivity_qc_recipe
TARGETS = ["qc"]
