import logging
from pathlib import Path
from typing import Any

import lmdb
from datacooker import TransformFunc, parse_dict
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)



def parallel_process(  # noqa: PLR0913
    data_list: list,
    inputs: dict[str, Any],
    recipe: Path,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    **extra_kwargs: Any,  # noqa: ANN401
) -> list:
    """
    Build an LMDB database from parsed data.

    Args:
        env_path: Path to the LMDB environment.
        data_list: List of paths to data files to parse.
        parser: Function to parse individual data files.
        n_jobs: Number of parallel jobs for parsing.
        map_size: Maximum size of the LMDB database in bytes.
    """
    def _process_item(data_dict : dict) -> tuple[dict, Exception | None]:
        """Parse a single file and return (key, compressed_data, error)."""
        try:
            data_dict.update(**inputs)
            results = parse_dict(
                recipe_path=recipe,
                datadict=data_dict,
                transform_func=transform_func,
                **extra_kwargs,
            )
            return results, None
        except Exception as error:  # noqa: BLE001
            return {}, error

    logger.info("To be parsed %d entries.", len(data_list))

    if test_run:
        test_data = data_list[:10]
        test_results = []
        for data in test_data:
            result = _process_item(data)
            if result is not None:
                test_results.append(result)
            if result[1] is not None:
                raise result[1]
    logger.info("Test run successful. Proceeding with full processing...")

    # --- Parallel processing ---
    results = []
    for i in range(0, len(data_list), chunk_size):
        logger.info(
            "Processing files %d to %d / %d",
            i,
            min(i + chunk_size, len(data_list)),
            len(data_list),
        )
        data_chunk = data_list[i : i + chunk_size]

        _results = Parallel(n_jobs=n_jobs, verbose=10, prefer="processes")(
            delayed(_process_item)(data) for data in data_chunk
        )
        results.extend(_results)
    return results
