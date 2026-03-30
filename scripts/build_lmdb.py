import logging
import sys
from pathlib import Path

import click
import lmdb
from omegaconf import OmegaConf

from pipelines.utils import load_config, load_data_list
from pipelines.utils.lmdb import build_lmdb, extract_key_list, read_lmdb, rebuild_lmdb

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    stream=sys.stdout,
    force=True,
)
OmegaConf.register_new_resolver("p", lambda x: Path(x))


# ==============================================================
# Command Group
# ==============================================================
@click.group()
def cli() -> None:
    """Build and merge LMDB databases from CIF files."""


# ==============================================================
# 1. Build Command
# ==============================================================
@cli.command("build")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", "-m", type=float, default=1e10, show_default=True)
@click.option(
    "--shard-idx",
    "-i",
    type=int,
    default=None,
    help="Index of the shard to process (0-based).",
)
@click.option("--n-shards", "-n", type=int, default=1, show_default=True)
def build(
    config: Path,
    map_size: float,
    shard_idx: int | None,
    n_shards: int,
) -> None:
    """
    Build an LMDB database from CIF_DIR into ENV_PATH.

    Example:
        python build_lmdb.py build --config configs/ccd_lmdb.yaml --map-size 1e12 --shard-idx 0 --n-shards 4
    """
    map_size = int(map_size)
    config_dict = load_config(config)
    click.echo("Start load data")
    data_list = load_data_list(
        config_dict["data_dir"],
        pattern=config_dict["file_pattern"],
    )
    click.echo("End load data")

    if shard_idx is not None:
        if shard_idx < 0 or shard_idx >= n_shards:
            msg = f"Invalid shard index {shard_idx} for {n_shards} shards."
            raise click.BadParameter(msg)
        # Split CIF list into n_shards and take only the shard_idx part
        data_list = [
            data for i, data in enumerate(data_list) if i % n_shards == shard_idx
        ]
        click.echo(
            f"Processing shard {shard_idx}/{n_shards} with {len(data_list)} files.",
        )
        config_dict["env_path"] = Path(config_dict["env_path"]).with_name(
            f"{Path(config_dict['env_path']).stem}_shard{shard_idx}{Path(config_dict['env_path']).suffix}",
        )
    else:
        click.echo(f"Processing all {len(data_list)} files as a single shard.")

    build_lmdb(*data_list, **config_dict, map_size=map_size)

    # check the key count
    env = lmdb.open(str(config_dict["env_path"]), readonly=True, lock=False)
    with env.begin() as txn:
        cursor = txn.cursor()
        key_count = sum(1 for _ in cursor)
    env.close()
    click.echo(f"[Done] Built LMDB at {config_dict['env_path']} with {key_count} keys.")


# ==============================================================
# 2. Merge Command (auto-detect *.lmdb by stem pattern)
# ==============================================================
@cli.command("merge")
@click.argument("shard_pattern", type=str)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output merged LMDB path.",
)
@click.option("--map-size", "-m", type=float, default=1e12, show_default=True)
@click.option("--overwrite", is_flag=True, help="Overwrite existing LMDB if it exists.")
def merge(shard_pattern: str, output: Path, map_size: float, overwrite: bool) -> None:  # noqa: FBT001
    """
    Merge multiple LMDB shard databases into a single LMDB file.

    Example:
        python build_lmdb.py merge "/data/BioMolDBv2_2024Oct21/cif_shard*.lmdb" -o /data/BioMolDBv2_2024Oct21/cif_merged.lmdb
    """
    map_size = int(map_size)

    # Expand wildcard pattern
    pattern_path = Path(shard_pattern)
    shard_paths = sorted(pattern_path.parent.glob(pattern_path.name))
    if not shard_paths:
        msg = f"No LMDB files found for pattern: {shard_pattern}"
        raise click.ClickException(msg)

    if output.exists() and not overwrite:
        msg = f"{output} already exists. Use --overwrite to replace it."
        raise click.ClickException(
            msg,
        )

    click.echo(f"Found {len(shard_paths)} shards:")
    for s in shard_paths:
        click.echo(f"  - {s}")

    merged_env = lmdb.open(str(output), map_size=map_size)
    total_keys = 0

    for shard_path in shard_paths:
        click.echo(f"Merging {shard_path}")
        shard_env = lmdb.open(str(shard_path), readonly=True, lock=False)
        with shard_env.begin() as shard_txn, merged_env.begin(write=True) as merged_txn:
            cursor = shard_txn.cursor()
            for key, value in cursor:
                merged_txn.put(key, value)
                total_keys += 1
        shard_env.close()

    merged_env.sync()
    merged_env.close()

    click.echo(f"[Done] Merged {len(shard_paths)} shards into {output}")
    click.echo(f"Total keys merged: {total_keys}")


@cli.command("rebuild")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", "-m", type=float, default=1e12, show_default=True)
def rebuild(
    config: Path,
    map_size: float,
) -> None:
    """
    Rebuild an LMDB database from an existing LMDB database with transformations.

    Example:
        python build_lmdb.py rebuild --config configs/AF3_training.yaml
    """
    map_size = int(map_size)
    config_dict = load_config(config)

    rebuild_lmdb(
        **config_dict,
        map_size=map_size,
    )

    # check the key count
    key_list = extract_key_list(config_dict["new_env_path"])
    key_count = len(key_list)
    click.echo(
        f"[Done] Built LMDB at {config_dict['new_env_path']} with {key_count} keys.",
    )


@cli.command("check_db")
@click.argument("db_path", type=click.Path(exists=True, path_type=Path))
def check_db(db_path: Path) -> None:
    """Check the number of keys in the LMDB database."""
    env = lmdb.open(str(db_path), readonly=True, lock=False)
    with env.begin() as txn:
        cursor = txn.cursor()
        key_count = sum(1 for _ in cursor)
    env.close()
    click.echo(f"{db_path}: {key_count} keys")


@cli.command("test_db")
@click.argument("db_path", type=click.Path(exists=True, path_type=Path))
def test_db(db_path: Path) -> None:
    """Test reading a specific key from the LMDB database."""
    key = "P0003600"
    data = read_lmdb(db_path, key)


# ==============================================================
# Entrypoint
# ==============================================================
if __name__ == "__main__":
    cli()
