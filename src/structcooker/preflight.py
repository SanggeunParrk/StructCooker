"""Preflight readiness check for ``structcooker build-all``.

Answers one question: *is everything a build needs already in place?* For each db
config it resolves the paths the recipe **reads** and sorts them into

* **upstream** -- produced by another node in ``db/MANIFEST.yaml`` (build-all will
  build it first; not something the user provides), and
* **external** -- a raw input under ``DATA_ROOT`` that must already exist (mmCIF, CCD,
  OpenFold sets, SabDab, provided MSAs/SignalP).

plus the external tools the recipes shell out to and the configured input and output roots. The report
tells you exactly which external inputs are missing and, when known, the
``structcooker download`` target that fetches them -- so a newcomer knows what to get
before running build-all.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from structcooker.paths import distillation_root, mmcif_root

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
    # done_dir: the resume marker tree, which is this node's own output (hmm_output_dir).
    {"tmp_dir", "output_dir", "hmm_output_dir", "done_dir", "out_path", "old_seq_id_map_path"},
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
    ("materials/raw/cif", "mmcif"),
    ("openfold_distillation", "openfold"),   # huge -- code-only, warns before fetching
    ("external/SabDab", "sabdab"),
    ("materials/raw/ccd", "ccd"),
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

    def collect(value: object) -> None:
        if isinstance(value, (list, tuple)):
            for item in value:
                collect(item)
        elif isinstance(value, dict):
            for key, item in value.items():
                if key not in _NON_INPUT_KEYS and not key.startswith("out_"):
                    collect(item)
        elif _is_pathish(value):
            paths.append(Path(str(value)))

    for key, value in cfg.items():
        if key in (*_PRIMARY_OUTPUT_KEYS, "recipe", "recipe_path", "metadata_recipe") or key.endswith("_recipe_path"):
            continue
        if key in ("data_dir", "file_list") or key.endswith(("_path", "_paths")):
            collect(value)
    for block in ("inputs", "metadata_input", "additional_inputs", "parameters", "runtime_environment"):
        section = cfg.get(block)
        if isinstance(section, dict):
            collect(section)
    runtime = cfg.get("runtime_environment", {})
    if isinstance(runtime, dict) and runtime.get("OPENFOLD_STRUCTURE_OVERRIDES"):
        manifest = Path(str(runtime["OPENFOLD_STRUCTURE_OVERRIDES"]))
        if manifest.is_file():
            replacements = json.loads(manifest.read_text())
            if not isinstance(replacements, dict) or any(not isinstance(v, str) for v in replacements.values()):
                msg = "Invalid structure override manifest"
                raise ValueError(msg)
            paths.extend(Path(v) for v in replacements.values())
    return list(dict.fromkeys(paths))


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
    for root, target in ((mmcif_root(data_root), "mmcif"),
                         (distillation_root(data_root), "openfold")):
        if path == root or root in path.parents:
            return target
    try:
        rel = path.relative_to(data_root).as_posix()
    except ValueError:
        return None
    for marker, target in _DOWNLOADABLE:
        if rel == marker or rel.startswith(marker + "/"):
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


def dependency_errors(
    configs: Mapping[str, Mapping[str, object]],
    dependencies: Mapping[str, list[str]],
) -> list[str]:
    """Find invalid DAGs and reads whose producers are not declared ancestors."""
    errors: list[str] = []
    ancestors: dict[str, set[str]] = {}

    def visit(name: str, active: frozenset[str]) -> set[str]:
        if name in active:
            errors.append(f"dependency cycle at {name}")
            return set()
        if name not in configs:
            errors.append(f"unknown dependency {name}")
            return set()
        if name not in ancestors:
            parents: set[str] = set()
            for parent in dependencies.get(name, []):
                parents.add(parent)
                parents.update(visit(parent, active | {name}))
            ancestors[name] = parents
        return ancestors[name]

    produced: dict[Path, str] = {}
    for name, cfg in configs.items():
        visit(name, frozenset())
        for path in output_paths(cfg):
            if path in produced and produced[path] != name:
                errors.append(f"{name} and {produced[path]} both write {path}")
            produced[path] = name
    for name, cfg in configs.items():
        for path in input_paths(cfg):
            producer = produced.get(path)
            if producer and producer != name and producer not in ancestors[name]:
                errors.append(f"{name} reads {path} without depending on {producer}")
    return sorted(set(errors))


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
    """Return deployment roots and supplied-input overrides."""
    return {name: os.environ.get(name) for name in (
        "DATA_ROOT", "OUTPUT_ROOT", "MMCIF_ROOT", "DISTILLATION_ROOT", "MSA_ROOT", "SEQ_CLUSTER30_PATH", "SEQ_CLUSTER40_PATH",
    )}


def missing_externals(reports: Iterable[NodeReport], data_root: str) -> dict[str, str | None]:
    """Return every missing external path -> its download hint (deduped)."""
    out: dict[str, str | None] = {}
    for report in reports:
        for path in report.missing:
            out.setdefault(path, download_hint(Path(path), data_root))
    return out
