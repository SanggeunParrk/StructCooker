"""Table-grouping transform helper for StructCooker ingest recipes.

The generic, dependency-free instruction helpers that used to live here
(``single_value_instruction``, ``extract_float_single``, ``key_stack``,
``merge_dict``) now live in :mod:`datacooker.transforms`. ``get_smaller_dict``
stays here because it depends on numpy.
"""

from collections.abc import Callable
from typing import TypeVar

import numpy as np
from numpy.typing import NDArray

InputType = TypeVar("InputType", str, int, float)


def get_smaller_dict(
    *,
    dtype: type[InputType],
) -> Callable[..., dict[str, dict[str, NDArray]]]:
    """
    Group rows of cif_raw_dict by tied_to column(s).

    Parameters
    ----------
    cif_raw_dict : dict[str, list[str]]
        Column -> values (all lists must be same length)
    tied_to : str | tuple[str, str]
        Column(s) used for grouping.
        If tuple, the two values are joined with "|" to form the key.

    Returns
    -------
    dict[str, dict[str, list[str]]]
        Outer dict: group key -> inner dict (col -> list of values).
    """

    def _group_rows(
        cif_raw_dict: dict[str, list[str]],
        tied_to: str | tuple[str, str],
        columns: list[str] | None = None,
    ) -> dict[str, dict[str, NDArray]]:
        if not cif_raw_dict:
            return {}

        n = len(next(iter(cif_raw_dict.values())))
        cols = set(cif_raw_dict.keys())
        if columns is not None:
            cols = cols & set(columns)
        cols = list(
            cols - {tied_to}
            if isinstance(tied_to, str)
            else cols - {tied_to[0], tied_to[1]},
        )

        result: dict[str, dict[str, NDArray]] = {}

        for i in range(n):
            row = {col: cif_raw_dict[col][i] for col in cols}

            if isinstance(tied_to, str):
                key = cif_raw_dict[tied_to][i]
            else:
                key = (cif_raw_dict[tied_to[0]][i], cif_raw_dict[tied_to[1]][i])

            if key not in result:
                result[key] = {col: [] for col in cols}

            for col in cols:
                result[key][col].append(row[col])

        # convert lists to arrays
        for outer_dict in result.values():
            for col in outer_dict:
                outer_dict[col] = np.array(outer_dict[col], dtype=dtype)
        return result

    return _group_rows
