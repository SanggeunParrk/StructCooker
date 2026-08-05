"""Preflight readiness check for ``structcooker build-all``.

Answers one question: *is everything a build needs already in place?* For each db
config it resolves the paths the recipe **reads** and sorts them into

* **upstream** -- produced by another node in ``db/MANIFEST.yaml`` (build-all will
  build it first; not something the user provides), and
* **external** -- a raw input under ``DATA_ROOT`` that must already exist (mmCIF, CCD,
  OpenFold sets, SabDab, provided MSAs/SignalP).

plus the external tools the recipes shell out to and the two env roots. The report
tells you exactly which external inputs are missing and, when known, the
``structcooker download`` target that fetches them -- so a newcomer knows what to get
before running build-all.
"""
from __future__ import annotations

import importlib.util
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

# Top-level keys that hold this node's primary OUTPUT.
_PRIMARY_OUTPUT_KEYS = ("new_env_path", "env_path", "output_data_path", "output_path")
# Keys inside `inputs` that are this node's *side-effect* outputs (a multi-output
# materialize writes out_* files; a parallel op fans out to a directory), so a
# downstream node reading one of them depends on this node, not on a raw input.
_SIDE_OUTPUT_KEYS = frozenset({"output_dir", "hmm_output_dir"})   # + any `out_*` key
# Keys inside `inputs`/`metadata_input` that are NOT prerequisites to check: this
# node's own scratch/outputs, or the optional seq_id seed (absent == fresh id space).
_NON_INPUT_KEYS = frozenset(
    {"tmp_dir", "output_dir", "hmm_output_dir", "out_path", "old_seq_id_map_path"},
)

# External binaries the recipes shell out to (name -> optional/licensed).
_BIN_TOOLS: tuple[tuple[str, bool], ...] = (
    ("hmmbuild", False), ("hmmsearch", False), ("mmseqs", False),
    ("cd-hit", False), ("signalp6", True),
)
# Tools imported as Python modules (kalign-python, anarci), not PATH binaries.
_PY_TOOLS: tuple[tuple[str, bool], ...] = (("kalign", False), ("anarci", False))

# DATA_ROOT-relative prefixes a downloader can fetch -> `structcooker download` target.
_DOWNLOADABLE: tuple[tuple[str, str], ...] = (
    ("mmcif_files_latest", "mmcif"),
    ("openfold_distillation", "openfold"),   # huge -- code-only, warns before fetching
    ("external/SabDab", "sabdab"),
)


@dataclass
class NodeReport:
    """Readiness of one db config."""

    name: str
    upstream: list[str] = field(default_factory=list)     # inputs built by other nodes
    present: list[str] = field(default_factory=list)      # external inputs that exist
    missing: list[str] = field(default_factory=list)      # external inputs absent

    @property
    def ready(self) -> bool:
        """True when no external input is missing (upstream deps are build-all's job)."""
        return not self.missing


def _is_pathish(value: object) -> bool:
    return isinstance(value, Path) or (isinstance(value, str) and "/" in value)


def input_paths(cfg: Mapping[str, object]) -> list[Path]:
    """Return the filesystem paths a config *reads* (its prerequisites)."""
    paths: list[Path] = []
    for key in ("data_dir", "file_list", "old_env_path", "db_path"):
        value = cfg.get(key)
        if value is not None and _is_pathish(value):
            paths.append(Path(str(value)))
    for block in ("inputs", "metadata_input", "additional_inputs"):
        section = cfg.get(block)
        if not isinstance(section, dict):
            continue
        for key, value in section.items():
            if key in _NON_INPUT_KEYS or key.startswith("out_"):
                continue
            if _is_pathish(value):
                paths.append(Path(str(value)))
    return paths


def output_paths(cfg: Mapping[str, object]) -> list[Path]:
    """Return every path this config produces -- primary target + side-effect files.

    A multi-output materialize also writes its ``out_*`` files, and a parallel op fans
    out to ``output_dir`` / ``hmm_output_dir``; downstream nodes read those, so they
    count as this node's outputs (not raw inputs) for the upstream/external split.
    """
    outs: list[Path] = []
    for key in _PRIMARY_OUTPUT_KEYS:
        value = cfg.get(key)
        if value is not None and _is_pathish(value):
            outs.append(Path(str(value)))
    section = cfg.get("inputs")
    if isinstance(section, dict):
        for key, value in section.items():
            if (key.startswith("out_") or key in _SIDE_OUTPUT_KEYS) and _is_pathish(value):
                outs.append(Path(str(value)))
    return outs


def download_hint(path: Path, data_root: str) -> str | None:
    """Return the ``structcooker download`` target for a missing external path, if any."""
    rel = str(path).removeprefix(data_root.rstrip("/") + "/")
    for marker, target in _DOWNLOADABLE:
        if marker in rel:
            return target
    return None


def inspect(configs: Mapping[str, Mapping[str, object]]) -> list[NodeReport]:
    """Classify every config's inputs as upstream (a node output) or external.

    ``configs`` maps node name -> loaded config. An input path that equals another
    node's output (primary or side-effect) is *upstream*; anything else is an *external*
    raw input whose on-disk presence is checked here.
    """
    produced: dict[str, str] = {}
    for name, cfg in configs.items():
        for out in output_paths(cfg):
            produced[str(out)] = name

    reports: list[NodeReport] = []
    for name, cfg in configs.items():
        report = NodeReport(name)
        for path in input_paths(cfg):
            key = str(path)
            if key in produced and produced[key] != name:
                report.upstream.append(produced[key])
            elif path.exists():
                report.present.append(key)
            else:
                report.missing.append(key)
        reports.append(report)
    return reports


def check_tools() -> list[tuple[str, bool, bool]]:
    """Return ``(tool, found, optional)`` for each external tool the recipes need.

    Binaries are looked up on ``PATH``; kalign / anarci are Python modules.
    """
    found: list[tuple[str, bool, bool]] = [
        (tool, shutil.which(tool) is not None, optional) for tool, optional in _BIN_TOOLS
    ]
    found += [
        (f"{mod} (py)", importlib.util.find_spec(mod) is not None, optional)
        for mod, optional in _PY_TOOLS
    ]
    return found


def check_env() -> dict[str, str | None]:
    """Return the two roots + the optional seq_id seed from the environment."""
    return {name: os.environ.get(name) for name in ("DATA_ROOT", "OUTPUT_ROOT", "SEQID_SEED")}


def missing_externals(reports: Iterable[NodeReport], data_root: str) -> dict[str, str | None]:
    """Return every missing external path -> its download hint (deduped)."""
    out: dict[str, str | None] = {}
    for report in reports:
        for path in report.missing:
            out.setdefault(path, download_hint(Path(path), data_root))
    return out
