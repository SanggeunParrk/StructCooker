from pathlib import Path
from typing import Any

from datacooker import ConvertFunc, LoadFunc, TransformFunc
from datacooker.utils.db import (
    build_lmdb as _build_lmdb,
)
from datacooker.utils.db import (
    extract_lmdb_keys,
)
from datacooker.utils.db import (
    read_lmdb as _read_lmdb,
)
from datacooker.utils.db import (
    rebuild_lmdb as _rebuild_lmdb,
)

from pipelines.utils.convert import from_bytes, to_bytes


def extract_key_list(env_path: Path) -> list[str]:
    """Backward-compatible alias for ``extract_lmdb_keys``."""
    return extract_lmdb_keys(env_path)


def build_lmdb(  # noqa: PLR0913
    *data_list: Path,
    env_path: Path,
    recipe: Path,
    inputs: dict[str, float | Path] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: dict[str, Path] | None = None,
    load_func: LoadFunc | None = None,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    test_run: bool = True,
    **extra_kwargs: Any,  # noqa: ANN401
) -> None:
    """Build an LMDB database using the shared DataCooker utility."""
    _build_lmdb(
        *data_list,
        env_path=env_path,
        recipe=recipe,
        serialize=to_bytes,
        inputs=inputs,
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        load_func=load_func,
        transform_func=transform_func,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        map_size=map_size,
        test_run=test_run,
        **extra_kwargs,
    )


def rebuild_lmdb(  # noqa: PLR0913
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
    map_size: int = int(1e12),
    test_run: bool = True,  # noqa: FBT001, FBT002
    **extra_kwargs: Any,  # noqa: ANN401
) -> None:
    """Rebuild an LMDB database using the shared DataCooker utility."""
    _rebuild_lmdb(
        old_env_path=old_env_path,
        new_env_path=new_env_path,
        recipe=recipe,
        serialize=to_bytes,
        deserialize=from_bytes,
        parameters=parameters,
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        convert_func=convert_func,
        transform_func=transform_func,
        split_entries=True,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        map_size=map_size,
        test_run=test_run,
        **extra_kwargs,
    )


def read_lmdb(env_path: Path, key: str) -> dict[str, Any]:
    """Read a StructCooker LMDB entry through the shared DataCooker utility."""
    return _read_lmdb(env_path, key, deserialize=from_bytes)
