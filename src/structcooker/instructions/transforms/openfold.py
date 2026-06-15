"""Transforms for OpenFold3 distillation ingest recipes."""

from collections.abc import Callable

import numpy as np

_MSA_FIELDS = ("msa", "deletion_matrix", "metadata")
_TEMPLATE_FIELDS = ("index", "release_date", "idx_map")


def build_msa_dict() -> Callable[..., dict[str, dict[str, np.ndarray]]]:
    """Normalise per-source MSA payloads into a single ``msa_dict``.

    Protein monomers carry one source, RNA monomers carry several; both are
    stored under ``msa_dict`` keyed by source name so a single reader/writer
    pair serves every distillation dataset.
    """

    def _worker(
        msa_sources: dict[str, dict[str, np.ndarray]],
    ) -> dict[str, dict[str, np.ndarray]]:
        return {
            source: {field: payload[field] for field in _MSA_FIELDS}
            for source, payload in msa_sources.items()
        }

    return _worker


def build_template_dict() -> Callable[..., dict[str, dict[str, object]]]:
    """Normalise per-hit template payloads into a single ``template_dict``.

    Each hit (``<pdb_id>_<chain_id>``) keeps its ``index`` / ``release_date`` /
    ``idx_map`` payload. Only protein monomers carry templates.
    """

    def _worker(
        template_hits: dict[str, dict[str, object]],
    ) -> dict[str, dict[str, object]]:
        return {
            hit: {field: payload[field] for field in _TEMPLATE_FIELDS}
            for hit, payload in template_hits.items()
        }

    return _worker
