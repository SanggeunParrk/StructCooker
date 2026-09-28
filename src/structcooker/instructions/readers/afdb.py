"""Key builders for the AFDB multimer (AFM) DBs.

Every AFDB complex-release file names its entity: ``AF-<n>-model_v1.pdb``,
``AF-<n>-model_v1.cif.gz``, ``AF-<n>-msa_v1.a3m.zst``. Structures are keyed by that entity
(``AF-<n>``), the id the release's metadata uses. MSAs are keyed by seq_id instead, like
every other MSA DB in the set: an MSA depends only on its sequence, and seq_id is the one
id space shared across DBs (docs/seq-id-and-cluster-scheme.md).
"""

from __future__ import annotations

import os
import re
from pathlib import Path

_ENTITY = re.compile(r"^(AF-\d+)-(?:model|msa)_v1\.")


def afdb_entity_key(path: Path) -> str:
    """Key a release file by its entity (``AF-<n>``), dropping the ``-model_v1`` suffix."""
    m = _ENTITY.match(Path(path).name)
    if m is None:
        msg = f"not an AFDB release file name: {Path(path).name}"
        raise ValueError(msg)
    return m.group(1)


_AFM_MSA_SEQID: dict[str, str] | None = None


def _afm_msa_seqid_map() -> dict[str, str]:
    """Load (and cache) the MSA entity -> seq_id map built by metadata/afm_msa_seqid.

    The path comes from ``AFM_MSA_SEQID_MAP`` (defaulting under ``$OUTPUT_ROOT/metadata``),
    the same pattern ``openfold_seqid_key`` uses for the distillation MSAs.
    """
    global _AFM_MSA_SEQID  # noqa: PLW0603 - process-local cache
    if _AFM_MSA_SEQID is None:
        root = os.environ.get("OUTPUT_ROOT", "/data/shared/cssb_data/BioMol_clean")
        map_path = Path(os.environ.get(
            "AFM_MSA_SEQID_MAP", str(Path(root) / "metadata" / "afm_msa_seqid.tsv"),
        ))
        mapping: dict[str, str] = {}
        with map_path.open(encoding="utf-8") as handle:
            for line in handle:
                entity, _, seq_id = line.rstrip("\n").partition("\t")
                if seq_id:
                    mapping[entity] = seq_id
        _AFM_MSA_SEQID = mapping
    return _AFM_MSA_SEQID


def afm_msa_seqid_key(path: Path) -> str:
    """Key an AFM MSA file by its query sequence's seq_id.

    Two entities with the same sequence collapse onto one seq_id; the build keeps the
    first (``duplicate_input_policy: first_path``), as the distillation MSA DBs do.
    """
    return _afm_msa_seqid_map()[afdb_entity_key(path)]
