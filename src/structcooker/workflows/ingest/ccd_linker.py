"""CCD linker-reference ingest recipe.

Turns the raw CIF chemical-component tables of one component into a linker
descriptor (atoms that can form inter-residue bonds + intra ideal bond lengths)
via :func:`build_ccd_linker`. Hydrogens are kept while grouping so that
leaving hydrogens (which mark the amino-N / anomeric-C linkers) are visible.
Read/Write boundaries live in ``configs/ingest/ccd_linker.yaml``.
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.connectivity import build_ccd_linker
from structcooker.instructions.transforms.tables import get_smaller_dict

ccd_linker_recipe = RecipeBook()
_group_rows = get_smaller_dict(dtype=str)

ccd_linker_recipe.step(
    outputs=("_chem_comp_atom_dict", dict),
    instruction=_group_rows,
    kwargs={"cif_raw_dict": ("_chem_comp_atom", str | None)},
    params={
        "tied_to": "comp_id",
        "columns": [
            "atom_id",
            "type_symbol",
            "charge",
            "model_Cartn_x",
            "model_Cartn_y",
            "model_Cartn_z",
            "pdbx_leaving_atom_flag",
        ],
    },
)

ccd_linker_recipe.step(
    outputs=("_chem_comp_bond_dict", dict),
    instruction=_group_rows,
    kwargs={"cif_raw_dict": ("_chem_comp_bond", str | None)},
    params={
        "tied_to": "comp_id",
        "columns": ["atom_id_1", "atom_id_2", "value_order"],
    },
)

ccd_linker_recipe.step(
    outputs=("ccd_linker", dict),
    instruction=build_ccd_linker,
    kwargs={
        "chem_comp_atom_dict": ("_chem_comp_atom_dict", dict | None),
        "chem_comp_bond_dict": ("_chem_comp_bond_dict", dict | None),
    },
    params={"unwrap": True},
)

RECIPE = ccd_linker_recipe
TARGETS = ["ccd_linker"]
