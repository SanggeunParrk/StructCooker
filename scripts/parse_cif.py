from pathlib import Path

import click
from datacooker import parse_file

from pipelines.cifmol import CIFMol, to_cif
from pipelines.transforms.cif_transforms import dot_transform, get_cif_data
from pipelines.utils.convert import to_dict


@click.command()
@click.argument(
    "ccd_db_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
)
@click.argument(
    "cif_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.option(
    "--recipe-path",
    "-r",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    default=Path("pipelines/recipe/cif_recipe_book.py"),
    show_default=True,
    help="Path to the Cooker recipe file.",
)
@click.option(
    "--target",
    "-t",
    "targets",
    multiple=True,
    help="Targets to serve from the recipe. If omitted, all targets are served.",
)
def cli(
    ccd_db_path: Path,
    cif_path: Path,
    recipe_path: Path,
    targets: tuple[str, ...],
) -> None:
    """
    Parse a CIF file using a Cooker recipe and export CIFMol objects.

    CCD_DB_PATH: Path to the CCD LMDB/database directory.
    CIF_PATH:    Path to the input .cif or .cif.gz file.
    """
    targets_list: list[str] | None = list(targets) if targets else None

    result = parse_file(
        load_func=get_cif_data,
        transform_func=dot_transform,
        recipe_path=recipe_path,
        file_path=cif_path,
        targets=targets_list,
        ccd_db_path=ccd_db_path,
    )
    result = to_dict(result)

    value, metadata = result["assembly_dict"], result["metadata_dict"]
    cifmol_dict: dict[str, CIFMol] = {}

    for cif_key, item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")
        metadata["assembly_id"] = assembly_id
        metadata["model_id"] = model_id
        metadata["alt_id"] = alt_id
        item["metadata"] = metadata

        cifmol = CIFMol.from_dict(item)
        cifmol_dict[cif_key] = cifmol
        to_cif(cifmol, Path(f"test_{cif_key}.cif"))


if __name__ == "__main__":
    cli()
    # python scripts/parse_cif.py /public_data/CCD/biomol_CCD_202602.lmdb /home/psk6950/data/BioMolDB_20260224/cif/raw/hc/1hcu.cif.gz
    # python scripts/parse_cif.py /public_data/CCD/biomol_CCD_202602.lmdb /home/psk6950/data/BioMolDB_20260224/cif/raw/ud/6udr.cif.gz
    # python scripts/parse_cif.py /public_data/CCD/biomol_CCD_202602.lmdb /home/psk6950/data/BioMolDB_20260224/cif/raw/nm/4nmg.cif.gz
