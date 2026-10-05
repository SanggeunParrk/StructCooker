"""The release manifest: every per-set manifest under db/ merged into one build DAG.

Each db/MANIFEST_<set>.yaml stays the unit a person builds or reviews (PDB, distillation,
teddymer, AFDB multimer, ...). A node can appear in several of them -- a set lists the shared
nodes it reads (seq_id_map, cif_fasta, cif_chain) as roots so it can be built alone -- so the
union merges each node's dependency lists. ``structcooker build-all --manifest
db/MANIFEST_all.yaml`` is then the end-to-end build; ``write`` regenerates that file and a
test keeps it equal to the union of its sources.
"""
from __future__ import annotations

from pathlib import Path

import yaml

DB_ROOT = Path(__file__).resolve().parents[2] / "db"
RELEASE = DB_ROOT / "MANIFEST_all.yaml"


def source_manifests(db_root: Path = DB_ROOT) -> list[Path]:
    """Every per-set manifest (MANIFEST.yaml and MANIFEST_<set>.yaml), not the release one."""
    return sorted(p for p in db_root.glob("MANIFEST*.yaml") if p.name != RELEASE.name)


def merged(db_root: Path = DB_ROOT) -> dict[str, list[str]]:
    """Union of the per-set manifests: node -> sorted union of its dependencies."""
    deps: dict[str, set[str]] = {}
    for path in source_manifests(db_root):
        for node, needs in (yaml.safe_load(path.read_text()) or {}).items():
            deps.setdefault(node, set()).update(needs or [])
    return {node: sorted(deps[node]) for node in sorted(deps)}


def write(db_root: Path = DB_ROOT) -> Path:
    """Regenerate db/MANIFEST_all.yaml from the per-set manifests."""
    nodes = merged(db_root)
    header = (
        "# Release DAG: the union of every db/MANIFEST*.yaml (structcooker.manifests.write).\n"
        "# Do not edit by hand -- edit the per-set manifest and regenerate:\n"
        "#   python -c 'from structcooker.manifests import write; write()'\n"
        "# End to end: structcooker build-all --manifest db/MANIFEST_all.yaml (docs/e2e-build.md).\n"
    )
    body = "".join(f"{node}: [{', '.join(needs)}]\n" for node, needs in nodes.items())
    out = db_root / RELEASE.name
    out.write_text(header + body, encoding="utf-8")
    return out
