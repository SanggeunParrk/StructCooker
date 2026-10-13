"""Gather a published shard collection's shards inside its own directory.

A shard_lmdb collection (``<db>.lmdb/shards.json``) may point at shards left in the build,
release and patch directories it was assembled from (``.cif_long.build/<attempt>/``,
``.release_20260915/``, ``.msa_long.build/patch_*/``, ...). For a deployable tree each
collection holds its own shards: ``<db>.lmdb/shards/<NNN>.lmdb``.

Each shard's files are hard-linked (same filesystem: no copy, instant), the manifest is
replaced atomically, and the old paths are left in place, so readers never see a broken
collection and the old directories can be deleted afterwards. Default is a dry run.

    python scripts/staging/consolidate_shards.py <db>.lmdb [...] [--apply]
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

import lmdb

MANIFEST = "shards.json"


def plan(target: Path) -> list[tuple[Path, Path]]:
    manifest = json.loads((target / MANIFEST).read_text())
    moves = []
    for i, shard in enumerate(manifest["shards"]):
        src = Path(os.path.realpath(target / shard["path"]))
        moves.append((src, target / "shards" / f"{i:03d}.lmdb"))
    return moves


def entries(path: Path) -> int:
    with lmdb.open(str(path), readonly=True, lock=False) as env:
        return env.stat()["entries"]


def consolidate(target: Path, *, apply: bool) -> None:
    manifest = json.loads((target / MANIFEST).read_text())
    moves = plan(target)
    inside = sum(src == dst.resolve() for src, dst in moves)
    print(f"{target}: {len(moves)} shards, {inside} already inside, {manifest['entries']} entries")
    if not apply or inside == len(moves):
        for src, dst in moves[:3]:
            print(f"  {src} -> {dst}")
        return
    for (src, dst), shard in zip(moves, manifest["shards"], strict=True):
        if dst.exists():
            if os.path.samefile(dst / "data.mdb", src / "data.mdb"):
                continue
            msg = f"{dst} exists and is not {src}"
            raise FileExistsError(msg)
        dst.mkdir(parents=True)
        os.link(src / "data.mdb", dst / "data.mdb")
        if entries(dst) != shard["entries"]:
            msg = f"{dst}: entry count differs from the manifest"
            raise ValueError(msg)
    new = {**manifest, "shards": [{**s, "path": os.path.relpath(dst, target)}
                                  for (src, dst), s in zip(moves, manifest["shards"], strict=True)]}
    with tempfile.NamedTemporaryFile("w", dir=target, delete=False) as handle:
        json.dump(new, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    Path(handle.name).replace(target / MANIFEST)
    total = sum(entries(target / s["path"]) for s in new["shards"])
    if total != manifest["entries"]:
        msg = f"{target}: {total} entries after consolidation, manifest says {manifest['entries']}"
        raise ValueError(msg)
    print(f"  consolidated: {total} entries in {target / 'shards'}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="+", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    for target in args.targets:
        consolidate(target.resolve(), apply=args.apply)


if __name__ == "__main__":
    main()
