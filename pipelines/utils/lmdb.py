import logging
from pathlib import Path
from typing import Any

import lmdb
from datacooker import ConvertFunc, LoadFunc, TransformFunc, parse_dict, parse_file
from joblib import Parallel, delayed

from pipelines.utils.convert import from_bytes, to_bytes

logger = logging.getLogger(__name__)

_ENV_CACHE: dict[tuple[str, bool], lmdb.Environment] = {}


def _get_env(
    path: Path,
    *,
    readonly: bool = True,
    lock: bool = False,
) -> lmdb.Environment:
    """Get (or lazily create) an LMDB env for the current process."""
    key = (str(path), readonly)
    env = _ENV_CACHE.get(key)
    if env is None:
        env = lmdb.open(str(path), readonly=readonly, lock=lock)
        _ENV_CACHE[key] = env
    return env


def extract_key_list(env_path: Path) -> list[str]:
    """
    Retrieve all keys from the LMDB database.

    Args:
        env_path: Path to the LMDB environment.

    Returns
    -------
        set
            A set of all keys in the LMDB database.
    """
    env = lmdb.open(str(env_path), readonly=True, lock=False)
    with env.begin() as txn:
        key_set = {
            key.decode() for key in txn.cursor().iternext(keys=True, values=False)
        }
    env.close()
    return list(key_set)


def build_lmdb(  # noqa: PLR0913
    *data_list: Path,
    env_path: Path,
    recipe: Path,
    load_func: LoadFunc,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),  # ~1TB
    **extra_kwargs: Any,  # noqa: ANN401
) -> None:
    """
    Build an LMDB database from parsed data.

    Args:
        env_path: Path to the LMDB environment.
        data_list: List of paths to data files to parse.
        parser: Function to parse individual data files.
        n_jobs: Number of parallel jobs for parsing.
        map_size: Maximum size of the LMDB database in bytes.
    """
    env = lmdb.open(str(env_path), map_size=int(map_size))

    def _process_file(data_file: Path) -> tuple[bytes, bytes, Exception | None]:
        """Parse a single file and return (key, compressed_data, error)."""
        key = data_file.name.split(".")[0]
        try:
            data_dict = parse_file(
                recipe_path=recipe,
                file_path=data_file,
                load_func=load_func,
                transform_func=transform_func,
                **extra_kwargs,
            )
            zcompressed_data = to_bytes(data_dict)
            return key.encode(), zcompressed_data, None
        except Exception as error:  # noqa: BLE001
            return key.encode(), to_bytes({}), error

    # remove UNL
    filtered_data_list: list[Path]
    filtered_data_list = [p for p in data_list if p.stem != "UNL"]
    _already_parsed_keys = extract_key_list(env_path)
    logger.info("Already parsed %d entries. (%s)", len(_already_parsed_keys), env_path)
    filtered_data_list = [
        data
        for data in filtered_data_list
        if data.name.split(".")[0] not in _already_parsed_keys
    ]
    logger.info("To be parsed %d entries.", len(filtered_data_list))

    # --- Parallel processing ---
    for i in range(0, len(filtered_data_list), chunk_size):
        logger.info(
            "Processing files %d to %d / %d",
            i,
            min(i + chunk_size, len(filtered_data_list)),
            len(filtered_data_list),
        )
        data_chunk = filtered_data_list[i : i + chunk_size]

        results = Parallel(n_jobs=n_jobs, verbose=10, prefer="processes")(
            delayed(_process_file)(data_file) for data_file in data_chunk
        )
        # --- Write results to LMDB ---
        with env.begin(write=True) as txn:
            for result in results:
                if result is None:
                    continue
                key, zcompressed_data, error = result
                if error is not None:
                    # Log error message but continue
                    logger.error("Error processing %s: %s", key.decode(), error)
                    continue
                txn.put(key, zcompressed_data)

    env.close()


def rebuild_lmdb(  # noqa: PLR0913, PLR0915
    old_env_path: Path,
    new_env_path: Path,
    recipe: Path,
    parameters: dict[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: dict[str, Path] | None = None,
    convert_func: ConvertFunc | None = None,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),  # ~1TB
    **extra_kwargs: Any,  # including metadata path  # noqa: ANN401
) -> None:
    """
    Build an LMDB database from parsed data.

    Args:
        env_path: Path to the LMDB environment.
        data_list: List of paths to data files to parse.
        parser: Function to parse individual data files.
        n_jobs: Number of parallel jobs for parsing.
        map_size: Maximum size of the LMDB database in bytes.
    """
    new_env = lmdb.open(str(new_env_path), map_size=int(map_size))

    metadata_dict: dict = {}
    if metadata_recipe is not None:
        if metadata_input is None:
            msg = "metadata_input must be provided if metadata_recipe is specified."
            raise ValueError(msg)
        metadata_dict = parse_dict(
            recipe_path=metadata_recipe,
            datadict=metadata_input,
            **extra_kwargs,
        )

    def _process_file(
        inputs: tuple[bytes, bytes],
    ) -> tuple[bytes, bytes, Exception | None] | None:
        """Parse a single file and return (key, compressed_data, error)."""
        key, values = inputs
        data_bytes: bytes = bytes(values)
        data: dict = from_bytes(data_bytes)
        data = convert_func(data) if convert_func is not None else data
        output_dict = {}
        for data_key, inner_dict in data.items():
            datadict: dict = inner_dict.copy()
            if metadata_recipe is not None:
                datadict.update(metadata_dict)
            if parameters is not None:
                datadict.update(parameters)
            rebuild_data_dict = parse_dict(
                recipe_path=recipe,
                datadict=datadict,
                transform_func=transform_func,
                **extra_kwargs,
            )
            if all(value is None for value in rebuild_data_dict.values()):
                continue

            output_dict[data_key] = rebuild_data_dict
        if len(output_dict) == 0:
            return None
        return key, to_bytes(output_dict), None
        try:
            for data_key, inner_dict in data.items():
                datadict: dict = inner_dict.copy()
                if metadata_recipe is not None:
                    datadict.update(metadata_dict)
                if parameters is not None:
                    datadict.update(parameters)
                rebuild_data_dict = parse_dict(
                    recipe_path=recipe,
                    datadict=datadict,
                    transform_func=transform_func,
                    **extra_kwargs,
                )
                if all(value is None for value in rebuild_data_dict.values()):
                    continue

                output_dict[data_key] = rebuild_data_dict
            if len(output_dict) == 0:
                return None
            return key, to_bytes(output_dict), None
        except Exception as error:  # noqa: BLE001
            return key, to_bytes({}), error

    # remove UNL
    old_key_list = extract_key_list(old_env_path)
    _already_parsed_keys = extract_key_list(new_env_path)
    _already_parsed_keys = set(_already_parsed_keys)
    logger.info(
        "Already parsed %d entries. (%s)",
        len(_already_parsed_keys),
        new_env_path,
    )
    old_set = set(old_key_list)
    key_list = list(old_set - _already_parsed_keys)
    logger.info("To be parsed %d entries.", len(key_list))
    old_env = lmdb.open(str(old_env_path), readonly=True, lock=False)
    new_env = lmdb.open(str(new_env_path), map_size=int(map_size))

    # test run
    test_keys = key_list[:10]
    test_data_chunk = []
    with old_env.begin() as txn:
        for key in test_keys:
            key_bytes = key.encode()
            val_bytes = txn.get(key_bytes)
            if val_bytes is not None:
                test_data_chunk.append((key_bytes, bytes(val_bytes)))
    test_results = []
    for data in test_data_chunk:
        result = _process_file(data)
        if result is not None:
            test_results.append(result)
    logger.info("Test run processed %d / 10 entries.", len(test_results))

    # --- Parallel processing ---
    for i in range(0, len(key_list), chunk_size):
        logger.info(
            "Processing files %d to %d / %d",
            i,
            min(i + chunk_size, len(key_list)),
            len(key_list),
        )
        key_chunk = key_list[i : i + chunk_size]
        data_chunk = []
        with old_env.begin() as txn:
            for key in key_chunk:
                key_bytes = key.encode()
                val_bytes = txn.get(key_bytes)
                if val_bytes is not None:
                    data_chunk.append((key_bytes, bytes(val_bytes)))
        results = Parallel(
            n_jobs=n_jobs,
            verbose=10,
            prefer="threads",
            pre_dispatch="2*n_jobs",
            inner_max_num_threads=1,
        )(delayed(_process_file)(data) for data in data_chunk)
        # --- Write results to LMDB ---
        with new_env.begin(write=True) as txn:
            for result in results:
                if result is None:
                    continue
                key, zcompressed_data, error = result
                if error is not None:
                    # Log error message but continue
                    logger.error("Error processing %s: %s", key.decode(), error)
                    continue
                txn.put(key, zcompressed_data)

    new_env.close()


def merge_lmdb_shards(
    shard_paths: list[Path],
    merged_env_path: Path,
    map_size: int = int(1e12),
    overwrite: bool = False,  # noqa: FBT001, FBT002
) -> None:
    """
    Merge multiple LMDB shard databases into a single LMDB file.

    Args:
        shard_paths: List of LMDB shard directories to merge.
        merged_env_path: Output LMDB path for the merged database.
        map_size: Maximum size of the merged LMDB in bytes.
        overwrite: Whether to overwrite existing merged file if it exists.
    """
    if merged_env_path.exists() and not overwrite:
        msg = f"{merged_env_path} already exists. Use overwrite=True to replace it."
        raise FileExistsError(
            msg,
        )

    # --- Create a new LMDB environment for the merged database ---
    merged_env = lmdb.open(str(merged_env_path), map_size=map_size)

    total_keys = 0

    # --- Iterate through each shard and copy all entries ---
    for shard_path in shard_paths:
        logger.info("Merging shard: %s", shard_path)
        shard_env = lmdb.open(str(shard_path), readonly=True, lock=False)
        with shard_env.begin() as shard_txn, merged_env.begin(write=True) as merged_txn:
            cursor = shard_txn.cursor()
            for key, value in cursor:
                merged_txn.put(key, value)
                total_keys += 1
        shard_env.close()

    merged_env.sync()
    merged_env.close()

    logger.info("[Done] Merged %d shards into %s", len(shard_paths), merged_env_path)
    logger.info("Total keys merged: %d", total_keys)


def read_lmdb(env_path: Path, key: str) -> dict:
    """
    Read a value from the LMDB database by key.

    Args:
        env_path: Path to the LMDB environment.
        key: Key of the data to retrieve.

    Returns
    -------
        dict
            The data dictionary retrieved from the LMDB database.
    """
    env = _get_env(env_path, readonly=True, lock=False)
    with env.begin() as txn:
        value = txn.get(key.encode())
    if value is None:
        msg = f"Key '{key}' not found in LMDB database."
        raise KeyError(msg)
    value_bytes: bytes = bytes(value)
    return from_bytes(value_bytes)
