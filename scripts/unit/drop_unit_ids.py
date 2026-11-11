"""Delete unit records from a single-LMDB unit DB and bring its sidecars (index.tsv, meta.json) in line.

usage: drop_unit_ids.py <unit_pdb.lmdb> <ids.txt>

For the tail of ``structcooker patch unit/unit_pdb --keys``: a patch replaces the units it rebuilt, but an entry
that now has fewer (or renumbered) units leaves old ids behind that no rebuilt unit overwrites. Those ids are
listed in ``ids.txt`` (old ids minus new ids of the patched entries) and removed here; each must exist.
Run through sbatch after the patch's apply job; the manifest / manifest_full / split stages are rerun afterwards.
"""

import sys
from dataclasses import replace
from pathlib import Path

import lmdb
from datacooker.lmdb.index import read_index, read_meta, write_index, write_meta


def main() -> None:
    db = Path(sys.argv[1])
    ids = [x.strip() for x in Path(sys.argv[2]).read_text().splitlines() if x.strip()]
    gone = set(ids)
    with (
        lmdb.open(str(db), map_size=1 << 40, lock=True) as env,
        env.begin(write=True) as txn,
    ):
        missing = [k for k in ids if txn.get(k.encode()) is None]
        if missing:
            sys.exit(
                f"{len(missing)} ids are not in the DB, e.g. {missing[:3]}; nothing deleted"
            )
        for k in ids:
            txn.delete(k.encode())
        entries = txn.stat()["entries"]
    old = read_index(db)
    keep = [s for s in old if s.key not in gone]
    if len(old) - len(keep) != len(gone):
        sys.exit(f"index holds {len(old) - len(keep)} of the {len(gone)} ids")
    stats = write_index(db, keep)
    meta = read_meta(db)
    if meta is not None:
        write_meta(
            db,
            replace(
                meta,
                entries=stats.entries,
                total_bytes=stats.total_bytes,
                max_bytes=stats.max_bytes,
            ),
        )
    assert stats.entries == entries, (stats.entries, entries)
    print(f"dropped {len(gone)} units; {entries} left (index and meta agree)")


if __name__ == "__main__":
    main()
