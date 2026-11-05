"""Durable build receipts and fail-closed restart handling."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
import time
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path


def digest(value: object) -> str:
    """Hash a canonical JSON description."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Publish a complete, flushed receipt atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, sort_keys=True, default=str)
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def exclusive(path: Path, wait: float = 0) -> Iterator[None]:
    """Reject concurrent writers; kernel locks are released on process death.

    ``wait`` > 0 queues for up to that many seconds instead: for short critical sections
    that peers enter routinely (two release stages finishing at once both dispatch).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        deadline = time.monotonic() + wait
        while True:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError as exc:
                if time.monotonic() >= deadline:
                    msg = f"Another build owns {path}"
                    raise RuntimeError(msg) from exc
                time.sleep(0.5)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def claim_output(output: Path, state_path: Path) -> None:
    """Keep ownership across driver death while native jobs can still publish.

    Caller must hold the output's filesystem lock. Reuse the original run directory
    to resume; changing run directories cannot bypass its outstanding job state.
    """
    owner = Path(str(output) + ".build-owner.json")
    identity = str(state_path.resolve())
    if owner.exists():
        if json.loads(owner.read_text()).get("state") != identity:
            msg = f"Output belongs to another run; resume its recorded state: {owner}"
            raise RuntimeError(msg)
    else:
        write_json(owner, {"state": identity})


class Fingerprints:
    """Content identity independent of filesystem timestamp resolution."""

    def file(self, path: Path) -> str:
        """Hash content on every check; coarse filesystem timestamps cannot prove identity."""
        before = path.stat()
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024**2), b""):
                hasher.update(block)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino,
        ):
            msg = f"Input changed while hashing: {path}"
            raise RuntimeError(msg)
        return hasher.hexdigest()

    def files(self, paths: list[Path]) -> dict[str, str]:
        """Overlap bounded file I/O without trusting size/mtime as a content hash."""
        # Network-filesystem reads overlap well: 64 in flight hash ~3.5x faster than 16
        # (OpenFold long raw: 3,200 vs 900 files/s, 2026-10-08).
        with ThreadPoolExecutor(max_workers=64) as pool:
            values = {}
            for offset in range(0, len(paths), 1024):
                batch = paths[offset:offset + 1024]
                values.update(zip(map(str, batch), pool.map(self.file, batch), strict=False))
            return values

    def path(self, path: Path) -> object:
        """Hash files/directories, excluding mutable LMDB reader lock files."""
        if not path.exists():
            return {"missing": str(path)}
        if path.is_file():
            return self.file(path)
        if (path / "shards.json").is_file():
            from datacooker.lmdb.sharded import shard_paths

            return {"manifest": self.file(path / "shards.json"),
                    "shards": {str(p): self.path(p) for p in shard_paths(path)}}
        if (path / "data.mdb").is_file():
            return {"data.mdb": self.file(path / "data.mdb")}
        from datacooker.utils.paths import scan_paths

        children = []
        for child in scan_paths(path, pattern="*"):        # parallel listing, sorted
            if "__pycache__" in child.parts or child.suffix == ".pyc":
                continue
            if child.is_dir():                              # listed as a non-directory: a symlink
                msg = f"Directory symlink requires an explicit input root: {child}"
                raise ValueError(msg)
            children.append(child)
        hashes = self.files(children)
        return {str(child.relative_to(path)): hashes[str(child)] for child in children}


