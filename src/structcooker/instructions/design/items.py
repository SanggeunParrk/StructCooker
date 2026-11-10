"""The items table of a design set and the DB rebuilt from it (stages ``design/<set>_items`` and ``design/<set>``).

The table has one row per source item: its statistics (:mod:`structcooker.instructions.design.stats`), the foldseek
cluster it belongs to, and the ``pass`` flag of the La-Proteina filters. The DB keeps the passing items and adds
the row to the metadata of each.
"""

from __future__ import annotations

import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import lmdb
from datacooker.lmdb.index import read_index

from structcooker.instructions.design.stats import STAT_COLUMNS

if TYPE_CHECKING:
    from structcooker.instructions.mols import CIFMolAttached

ITEM_COLUMNS = (*STAT_COLUMNS, "cluster", "cluster_size", "pass")
_SEP = "\t"


def passes(
    stats: list[str],
    min_plddt: float,
    min_length: int,
    max_length: int,
    max_coil: float,
    max_coil_run: int,
) -> bool:
    """La-Proteina filters on ``[length, mean_plddt, helix, strand, coil, max_coil_run]`` (strings, as written).

    pLDDT at least ``min_plddt``, ``min_length <= length <= max_length``, coil fraction below ``max_coil`` and the longest
    coil run at most ``max_coil_run``.
    """
    length, plddt, _, _, coil, run = stats
    return (
        float(plddt) >= min_plddt
        and min_length <= int(length) <= max_length
        and float(coil) < max_coil
        and int(run) <= max_coil_run
    )


def read_stats(stats_dir: Path) -> tuple[dict[str, list[str]], list[str]]:
    """``({id: columns after id}, failed ids)`` of every ``part_*.tsv`` / ``part_*.failed.tsv`` in ``stats_dir``."""
    stats: dict[str, list[str]] = {}
    for path in sorted(Path(stats_dir).glob("part_*.tsv")):
        if path.name.endswith(".failed.tsv"):
            continue
        for line in path.read_text().splitlines():
            key, *columns = line.split(_SEP)
            if key in stats:
                msg = f"{key} is in two parts"
                raise ValueError(msg)
            stats[key] = columns
    failed = [
        line.split(_SEP)[0]
        for path in sorted(Path(stats_dir).glob("part_*.failed.tsv"))
        for line in path.read_text().splitlines()
    ]
    return stats, failed


def read_representatives(clusters_path: Path) -> dict[str, str]:
    """``{member: representative}`` of the foldseek table written by the clusters stage."""
    rep: dict[str, str] = {}
    with Path(clusters_path).open() as handle:
        next(handle)  # header
        for line in handle:
            representative, member = line.rstrip("\n").split(_SEP)
            rep[member] = representative
    return rep


def build_items(
    stats_dir: Path,
    clusters_path: Path,
    db_path: Path,
    min_plddt: float = 80.0,
    min_length: int = 32,
    max_length: int = 512,
    max_coil: float = 0.5,
    max_coil_run: int = 20,
    apply_filters: bool = True,
    max_failed: int = 0,
) -> list[list[str]]:
    """Build the items table as rows (header first), sorted by id.

    Refuses the build when a source key has neither statistics nor a failure record, when more than ``max_failed``
    records failed, or when an item with a backbone is in no cluster. An item without any residue of complete backbone
    is not clustered (it has no structure to cluster): it is its own representative and does not pass. With
    ``apply_filters`` false every item passes (a set that was filtered when it was made). ``cluster_size`` counts the
    members of the cluster among all items, passing or not.
    """
    stats, failed = read_stats(stats_dir)
    expected = {item.key for item in read_index(db_path)}
    if not expected:
        msg = f"{db_path} has no index sidecar to check the statistics against"
        raise ValueError(msg)
    missing = expected - stats.keys() - set(failed)
    unknown = stats.keys() - expected
    if missing or unknown:
        msg = (
            f"statistics do not cover the source: {len(missing)} keys missing (e.g. {sorted(missing)[:3]}), "
            f"{len(unknown)} not in the source (e.g. {sorted(unknown)[:3]})"
        )
        raise ValueError(msg)
    if len(failed) > max_failed:
        msg = f"{len(failed)} records failed (allowed {max_failed}), e.g. {sorted(failed)[:3]}; see part_*.failed.tsv"
        raise ValueError(msg)
    rep = read_representatives(clusters_path)
    cluster: dict[str, str] = {}
    for key, columns in stats.items():
        if key in rep:
            cluster[key] = rep[key]
        elif int(columns[0]) == 0:
            cluster[key] = key
        else:
            msg = f"{key} has a backbone but is in no cluster"
            raise ValueError(msg)
    size = Counter(cluster.values())
    rows: list[list[str]] = [list(ITEM_COLUMNS)]
    for key in sorted(stats):
        keep = not apply_filters or passes(
            stats[key], min_plddt, min_length, max_length, max_coil, max_coil_run
        )
        rows.append(
            [key, *stats[key], cluster[key], str(size[cluster[key]]), str(int(keep))]
        )
    return rows


class ItemTable:
    """Read-only ``id -> row`` map of the items table, backed by an on-disk LMDB.

    The table has millions of rows; a Python dict in every worker would not fit. The driver writes the LMDB once and
    the workers on the node mmap that file; the object pickles to its path.
    """

    def __init__(self, path: str | Path) -> None:
        self._path = str(path)
        self._env: Any = None

    def __getstate__(self) -> dict[str, str]:
        """Pickle the path only: the environment handle is per process."""
        return {"_path": self._path}

    def __setstate__(self, state: dict[str, str]) -> None:
        """Restore from the path; the environment opens on first use."""
        self._path = state["_path"]
        self._env = None

    def get(self, key: str) -> dict[str, str] | None:
        """Return the row of ``key`` by column name (without ``id``), or None."""
        if self._env is None:
            self._env = lmdb.open(
                self._path, readonly=True, lock=False, subdir=True, readahead=False
            )
        with self._env.begin() as txn:
            raw = txn.get(key.encode())
        if raw is None:
            return None
        return dict(zip(ITEM_COLUMNS[1:], bytes(raw).decode().split(_SEP), strict=True))


def load_item_table(items_path: Path) -> ItemTable:
    """Write the items table to a node-local LMDB (under ``TMPDIR``) and return the map over it."""
    scratch = tempfile.mkdtemp(prefix="design_items_", dir=os.environ.get("TMPDIR"))
    path = Path(scratch) / "items.lmdb"
    env = lmdb.open(str(path), map_size=64 * (1 << 30), subdir=True)
    with Path(items_path).open() as handle, env.begin(write=True) as txn:
        header = handle.readline().rstrip("\n").split(_SEP)
        if header != list(ITEM_COLUMNS):
            msg = f"{items_path}: unexpected columns {header}"
            raise ValueError(msg)
        for line in handle:
            key, rest = line.rstrip("\n").split(_SEP, 1)
            txn.put(key.encode(), rest.encode())
    env.close()
    return ItemTable(path)


def annotate_item(cifmol: CIFMolAttached, items: Any, source: str) -> dict | None:
    """Return the design-DB record of one source item, or None when it does not pass the filters.

    Adds to the CIFMolAttached ``metadata``: ``source``, ``length``, ``mean_plddt``, ``ss_helix`` / ``ss_strand`` /
    ``ss_coil`` (fractions of the assessed residues), ``ss_max_coil_run``, ``fs_cluster`` (the foldseek representative)
    and ``fs_cluster_size``. An item that is not in the table is an error, not a rejection.
    """
    key = cifmol.id[0]
    row = items.get(key)
    if row is None:
        msg = f"{key} is not in the items table"
        raise KeyError(msg)
    if row[ITEM_COLUMNS[-1]] != "1":
        return None
    record = cast("dict", cifmol.to_dict())
    record["metadata"].update(
        {
            "source": source,
            "length": int(row["length"]),
            "mean_plddt": float(row["mean_plddt"]),
            "ss_helix": float(row["helix"]),
            "ss_strand": float(row["strand"]),
            "ss_coil": float(row["coil"]),
            "ss_max_coil_run": int(row["max_coil_run"]),
            "fs_cluster": row["cluster"],
            "fs_cluster_size": int(row["cluster_size"]),
        }
    )
    return record
