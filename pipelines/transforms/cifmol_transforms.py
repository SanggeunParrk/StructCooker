from typing import TYPE_CHECKING, cast

from pipelines.cifmol import CIFMol

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict


def convert_to_cifmol_dict(
    value: dict,
) -> dict[str, CIFMol]:
    """Convert a dictionary containing CIFMol data into a dictionary of CIFMol objects."""
    value, metadata = value["assembly_dict"], value["metadata_dict"]

    cifmol_dict: dict[str, CIFMol] = {}
    for cif_key, _item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")

        md = dict(metadata)
        md["assembly_id"] = assembly_id
        md["model_id"] = model_id
        md["alt_id"] = alt_id

        item = dict(_item)
        item["metadata"] = md
        item = cast("BioMolDict", item)

        cifmol_dict[cif_key] = CIFMol.from_dict(item)

    return cifmol_dict
