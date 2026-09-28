"""Durable build receipts and fail-closed restart handling."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
from collections.abc import Callable, Iterator, Mapping
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from typing import Any


class BuildFailedError(RuntimeError):
    """The recorded terminal job can no longer publish an output."""


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
def exclusive(path: Path) -> Iterator[None]:
    """Reject concurrent writers; kernel locks are released on process death."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            msg = f"Another build owns {path}"
            raise RuntimeError(msg) from exc
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
        with ThreadPoolExecutor(max_workers=16) as pool:
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
        entries = {}
        for child in sorted(path.rglob("*")):
            if "__pycache__" in child.parts or child.suffix == ".pyc":
                continue
            if child.is_symlink() and child.is_dir():
                msg = f"Directory symlink requires an explicit input root: {child}"
                raise ValueError(msg)
            if child.is_file():
                entries[str(child.relative_to(path))] = self.file(child)
        return entries


class BuildState:
    """One node's persistent submit/wait/verify state; never trust output existence."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.data: dict[str, Any] = json.loads(path.read_text()) if path.exists() else {}

    def save(self, **fields: Any) -> None:
        """Persist a state transition."""
        self.data.update(fields)
        write_json(self.path, self.data)

    def run(
        self, *, identity: str, submit: Callable[[], str], wait: Callable[[str], None],
        verify: Callable[[], Mapping[str, object]], unchanged: Callable[[], bool],
        output_signature: Callable[[], object],
    ) -> str:
        """Reuse verified outputs, resume a recorded job, or submit a new attempt."""
        previous = self.data.get("identity")
        phase = self.data.get("phase")
        if phase == "submitting":
            msg = "Submission outcome is ambiguous; reconcile scheduler jobs before retrying."
            raise RuntimeError(msg)
        if phase == "submitted" and previous != identity:
            msg = "Inputs changed while a recorded job is outstanding; finish or cancel it first."
            raise RuntimeError(msg)
        if phase == "complete" and previous == identity and self.data.get("output") == output_signature():
            return "reused"
        if phase not in {"submitted", "built"} or previous != identity:
            self.save(phase="submitting", identity=identity, error=None)
            try:
                job = submit()
                if not job:
                    msg = "Submission returned no job ID"
                    raise RuntimeError(msg)  # noqa: TRY301 - retain ambiguous submission state
            except Exception as exc:
                # A submitter may have created some jobs before losing its response.
                self.save(error=str(exc))
                raise
            self.save(phase="submitted", job=job)
        job = str(self.data["job"])
        # Keep phase=submitted on an interrupted/unknown wait, so retry cannot
        # submit a duplicate writer. Explicit scheduler reconciliation is required.
        if self.data["phase"] != "built":
            try:
                wait(job)
            except BuildFailedError as exc:
                self.save(phase="failed", error=str(exc))
                raise
            except Exception as exc:
                self.save(error=str(exc))
                raise
            self.save(phase="built")
        try:
            report = verify()
            if not unchanged():
                msg = "Build inputs changed during execution"
                raise RuntimeError(msg)  # noqa: TRY301 - persist verification failure
            signature = output_signature()
        except Exception as exc:
            self.save(phase="verification_failed", error=str(exc))
            raise
        self.save(phase="complete", output=signature, verification=dict(report), error=None)
        return "built"
