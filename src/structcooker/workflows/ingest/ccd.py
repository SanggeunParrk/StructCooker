"""CCD ingest recipe.

Declares the *Transform* graph that turns the raw CIF chemical-component tables
(`_chem_comp`, `_chem_comp_atom`, `_chem_comp_bond`) into the parsed
`chem_comp_dict`. The Read (loader) and Write (serializer) boundaries are not
part of the recipe; they are configured in `configs/ingest/ccd_lmdb.yaml`.
"""

from datacooker import RecipeBook

from structcooker.instructions.transforms.cif import parse_chem_comp
from structcooker.instructions.transforms.tables import get_smaller_dict

ccd_recipe = RecipeBook()

# One grouping instruction is shared by the three table-grouping steps.
_group_rows = get_smaller_dict(dtype=str)

ccd_recipe.step(
    outputs=("_chem_comp_dict", dict),
    instruction=_group_rows,
    kwargs={"cif_raw_dict": ("_chem_comp", str | None)},
    params={"tied_to": "id", "columns": ["name", "formula"]},
)

ccd_recipe.step(
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
            "pdbx_aromatic_flag",
            "pdbx_stereo_config",
        ],
    },
)

ccd_recipe.step(
    outputs=("_chem_comp_bond_dict", dict),
    instruction=_group_rows,
    kwargs={"cif_raw_dict": ("_chem_comp_bond", str | None)},
    params={
        "tied_to": "comp_id",
        "columns": [
            "atom_id_1",
            "atom_id_2",
            "value_order",
            "pdbx_aromatic_flag",
            "pdbx_stereo_config",
        ],
    },
)

ccd_recipe.step(
    outputs=("chem_comp_dict", dict),
    instruction=parse_chem_comp,
    kwargs={
        "chem_comp_dict": ("_chem_comp_dict", dict | None),
        "chem_comp_atom_dict": ("_chem_comp_atom_dict", dict | None),
        "chem_comp_bond_dict": ("_chem_comp_bond_dict", dict | None),
    },
    params={"remove_hydrogen": True, "unwrap": True},
)


RECIPE = ccd_recipe
TARGETS = ["chem_comp_dict"]
