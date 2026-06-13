import logging
from pathlib import Path
from typing import Any

from datacooker import TransformFunc, parse_dict
from datacooker.utils.sharding import resolve_node_config, shard_items
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


def parallel_process(  # noqa: PLR0913
    data_list: list,
    inputs: dict[str, Any],
    recipe: Path,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,  # noqa: FBT001, FBT002
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,  # noqa: ANN401
) -> list:
    """
    Build an LMDB database from parsed data.

    Args:
        data_list: List of paths to data files to parse.
        inputs: Shared input dictionary merged into each data entry.
        recipe: Recipe path used by ``parse_dict``.
        transform_func: Optional transform function for parsing.
        chunk_size: Number of entries to process per parallel batch.
        n_jobs: Number of parallel jobs for parsing.
        test_run: Whether to run a small serial smoke test before full run.
        node_rank: 0-based node index for multi-node processing.
        node_count: Total number of nodes for multi-node processing.
    """
    if chunk_size < 1:
        msg = f"chunk_size must be >= 1, got {chunk_size}."
        raise ValueError(msg)

    resolved_node_rank, resolved_node_count = resolve_node_config(
        node_rank=node_rank,
        node_count=node_count,
    )
    sharded_data_list = shard_items(
        data_list,
        node_rank=resolved_node_rank,
        node_count=resolved_node_count,
    )

    logger.info(
        "Node %d/%d assigned %d/%d entries.",
        resolved_node_rank,
        resolved_node_count,
        len(sharded_data_list),
        len(data_list),
    )

    if not sharded_data_list:
        logger.warning(
            "No data assigned to node %d/%d. Exiting early.",
            resolved_node_rank,
            resolved_node_count,
        )
        return []

    def _process_item(data_dict: dict[str, Any]) -> tuple[dict, Exception | None]:
        """Parse a single file and return (key, compressed_data, error)."""
        try:
            process_input = dict(data_dict)
            process_input.update(**inputs)
            results = parse_dict(
                recipe_path=recipe,
                datadict=process_input,
                transform_func=transform_func,
                **extra_kwargs,
            )
        except Exception as error:  # noqa: BLE001
            return {}, error
        return results, None

    logger.info("To be parsed %d entries.", len(sharded_data_list))

    if test_run:
        test_data = sharded_data_list[:10]
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
    for i in range(0, len(sharded_data_list), chunk_size):
        logger.info(
            "Processing files %d to %d / %d",
            i,
            min(i + chunk_size, len(sharded_data_list)),
            len(sharded_data_list),
        )
        data_chunk = sharded_data_list[i : i + chunk_size]

        _results = Parallel(n_jobs=n_jobs, verbose=10, prefer="processes")(
            delayed(_process_item)(data) for data in data_chunk
        )
        results.extend(_results)

    error_count = sum(1 for _, error in results if error is not None)
    logger.info(
        "Processing completed with %d errors out of %d entries.",
        error_count,
        len(results),
    )
    if error_count > 0:
        logger.info("Sample errors:")
        max_logged_errors = 10
        for i, (_, error) in enumerate(results):
            if error is not None:
                logger.info("Error %d: %s", i + 1, str(error))
                if i >= max_logged_errors - 1:
                    break

    return results
