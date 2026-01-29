import logging
import sys
from pathlib import Path

import click
import lmdb
from datacooker.core import ConvertFunc, ProjectFunc, parse_dict
from joblib import Parallel, delayed
from omegaconf import OmegaConf

from pipelines.utils import load_config
from pipelines.utils.lmdb import extract_key_list, read_lmdb

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    stream=sys.stdout,
    force=True,
)
OmegaConf.register_new_resolver("p", lambda x: Path(x))

_ENV_CACHE: dict[str, lmdb.Environment] = {}


@click.group()
def cli() -> None:
    """Build and merge LMDB databases from CIF files."""


@cli.command("data_transform")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
def data_transform(
    config: Path,
) -> None:
    """Load data and project down."""
    config_dict = load_config(config)

    # necessary fields : input_data_path, recipe_path, project_func, output_data_path
    necessary_fields = [
        "input_data_path",
        "recipe_path",
        "project_func",
        "output_data_path",
    ]
    for field in necessary_fields:
        if field not in config_dict:
            msg = f"Missing necessary field '{field}' in config."
            raise KeyError(msg)
    if not isinstance(config_dict["project_func"], ProjectFunc):
        msg = "'project_func' must be a ProjectFunc callable."
        raise TypeError(msg)

    results = parse_dict(
        datadict={"data_path": config_dict["input_data_path"]},
        recipe_path=config_dict["recipe_path"],
        transform_func=config_dict.get("transform_func", None),
    )
    config_dict["project_func"](
        data=results,
        output_path=config_dict["output_data_path"],
    )


@cli.command("db_extract")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--n_chunk", type=int, default=100)
def db_extract(
    config: Path,
    n_chunk: int,
) -> None:
    """Load data from db and project down collectively."""
    config_dict = load_config(config)

    # necessary fields : input_db_path, extract_recipe_path, project_func, output_data_path
    necessary_fields = [
        "input_db_path",
        "extract_recipe_path",
        "project_func",
        "output_data_path",
    ]
    for field in necessary_fields:
        if field not in config_dict:
            msg = f"Missing necessary field '{field}' in config."
            raise KeyError(msg)
    if not isinstance(config_dict["project_func"], ProjectFunc):
        msg = "'project_func' must be a ProjectFunc callable."
        raise TypeError(msg)

    input_db_path: Path = config_dict["input_db_path"]
    convert_func: ConvertFunc | None = config_dict.get("convert_func", None)
    transform_func = config_dict.get("transform_func", None)

    key_list = extract_key_list(config_dict["input_db_path"])

    def _process_chunk(keys: list[str]) -> dict[tuple[str, str], dict]:
        output_dict = {}
        for key in keys:
            data = read_lmdb(input_db_path, key)
            data = convert_func(data) if convert_func is not None else data
            datadict = {"inputs": data}
            datadict.update(**config_dict)

            output_dict[key] = parse_dict(
                datadict=datadict,
                recipe_path=config_dict["extract_recipe_path"],
                transform_func=transform_func,
            )
        return output_dict

    # make chunks
    key_chunks = [key_list[i : i + n_chunk] for i in range(0, len(key_list), n_chunk)]
    click.echo(f"Processing {len(key_chunks)} chunks of size {n_chunk}...")

    # parallel batch
    chunk_results = Parallel(n_jobs=-1, verbose=10)(
        delayed(_process_chunk)(chunk) for chunk in key_chunks
    )

    # merge
    merging_input_dict = {}
    for d in chunk_results:
        if d is None:
            continue
        merging_input_dict.update(d)

    if "merge_recipe_path" in config_dict:
        results = parse_dict(
            datadict={"data_dict": merging_input_dict},
            recipe_path=config_dict["merge_recipe_path"],
        )
    else:
        results = merging_input_dict
    config_dict["project_func"](
        data=results,
        output_path=config_dict["output_data_path"],
    )


if __name__ == "__main__":
    cli()
