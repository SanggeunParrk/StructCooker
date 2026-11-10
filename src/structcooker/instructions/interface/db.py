"""Read access to the interface-unit DB (``interface/interface_pdb``)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import lmdb
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes


def _is_db(path: Path) -> bool:
    return (path / "data.mdb").is_file() or (path / "shards.json").is_file()


class InterfaceDB:
    """``InterfaceDB(root).get("1tsr_1_1_0")``; keys are ``<pdb>_<assembly>_<model>_<unit index>``.

    ``root`` is the DB directory (``$BIOMOL_ROOT/lmdb/interface/pdb``) or the ``interface_pdb.lmdb``
    inside it, a merged LMDB or a published shard collection (``shards.json``), both read through
    ``datacooker.lmdb.sharded.open_env``. A directory of bare ``shards/shard_NNN.lmdb`` (the layout
    of the earlier driver, before the DB became a stage) is read too; there ``shard`` is the number
    in the file name, which a manifest of that layout lists, and ``get`` can use it to skip the search.
    """

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        coll = self.root / "interface_pdb.lmdb"
        self._path = next((p for p in (self.root, coll) if _is_db(p)), None)
        self.shards: dict[int, Path] = {}
        if self._path is None:
            shard_dir = (
                coll / "shards" if (coll / "shards").is_dir() else self.root / "shards"
            )
            self.shards = {
                int(p.stem.split("_")[1]): p
                for p in sorted(shard_dir.glob("shard_*.lmdb"))
            }
        self._env: dict[int, lmdb.Environment] = {}
        self._whole: Any = None

    def _open(self, i: int) -> lmdb.Environment:
        if i not in self._env:
            self._env[i] = lmdb.open(
                str(self.shards[i]), readonly=True, lock=False, readahead=False
            )
        return self._env[i]

    def get(self, unit_id: str, shard: int | None = None) -> dict:
        """Decode one record; ``shard`` only matters for the bare-shards layout."""
        if self._path is not None:
            if self._whole is None:
                self._whole = open_env(
                    str(self._path), readonly=True, lock=False, readahead=False
                )
            with self._whole.begin() as txn:
                raw = txn.get(unit_id.encode())
            if raw is None:
                raise KeyError(unit_id)
            return from_bytes(bytes(raw))
        for i in [shard] if shard is not None else self.shards:
            with self._open(i).begin() as txn:
                raw = txn.get(unit_id.encode())
            if raw is not None:
                return from_bytes(bytes(raw))
        raise KeyError(unit_id)
