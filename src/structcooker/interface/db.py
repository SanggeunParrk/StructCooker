"""Read access to the interface-unit DB (a directory of LMDB shards + manifest.tsv)."""

from __future__ import annotations

from pathlib import Path

import lmdb

from structcooker.instructions.transforms.codecs import from_bytes


class InterfaceDB:
    """``InterfaceDB(root).get("1tsr_1_1_0")``; keys are ``<pdb>_<assembly>_<model>_<unit index>``.

    ``root`` is the DB directory (``$BIOMOL_ROOT/lmdb/interface/pdb``): ``interface_pdb.lmdb/`` (the
    published shard collection; before ``finalize`` a plain ``shards/``) next to ``manifest.tsv``.
    """

    def __init__(self, root: str | Path):
        self.root = Path(root)
        coll = self.root / "interface_pdb.lmdb" / "shards"
        shard_dir = coll if coll.is_dir() else self.root / "shards"
        # the manifest's ``shard`` column is the number in the file name
        self.shards = {int(p.stem.split("_")[1]): p for p in sorted(shard_dir.glob("shard_*.lmdb"))}
        self._env: dict[int, lmdb.Environment] = {}

    def _open(self, i: int) -> lmdb.Environment:
        if i not in self._env:
            self._env[i] = lmdb.open(str(self.shards[i]), readonly=True, lock=False, readahead=False)
        return self._env[i]

    def get(self, unit_id: str, shard: int | None = None) -> dict:
        """Decode one record; pass the manifest's ``shard`` to avoid searching all shards."""
        for i in ([shard] if shard is not None else self.shards):
            with self._open(i).begin() as txn:
                raw = txn.get(unit_id.encode())
            if raw is not None:
                return from_bytes(bytes(raw))
        raise KeyError(unit_id)
