"""Transforms for OpenFold3 distillation ingest recipes.

OpenFold distillation alignments are stored split into an already-aligned
character matrix plus a separate deletion-count matrix, whereas the canonical
MSA pipeline parses raw a3m strings (lowercase = insertions). These adapters
express the distillation arrays back as a3m strings / standard headers so the
**existing** ``parse_sequence`` / ``parse_headers`` / ``build_dict``
instructions can be reused verbatim, guaranteeing identical record content.
"""

import re
from collections.abc import Callable

import numpy as np

_TEMPLATE_FIELDS = ("index", "release_date", "idx_map")


def merge_msa_sources() -> Callable[..., dict[str, np.ndarray]]:
    """Collapse the per-source MSA payloads into one alignment.

    Protein monomers expose a single source; RNA monomers expose several
    (``rnacentral_hits`` / ``nt_hits`` / ``rfam_hits``) that share the query
    columns, so their hit rows are concatenated below a single query row.
    """

    def _worker(
        msa_sources: dict[str, dict[str, np.ndarray]],
    ) -> dict[str, np.ndarray]:
        sources = list(msa_sources.values())
        first = sources[0]
        if len(sources) == 1:
            return {k: first[k] for k in ("msa", "deletion_matrix", "metadata")}
        # keep the query (row 0) once, then every source's hit rows
        msa = [first["msa"]]
        deletion = [first["deletion_matrix"]]
        metadata = [first["metadata"]]
        for src in sources[1:]:
            msa.append(src["msa"][1:])
            deletion.append(src["deletion_matrix"][1:])
            metadata.append(src["metadata"][1:])
        return {
            "msa": np.concatenate(msa, axis=0),
            "deletion_matrix": np.concatenate(deletion, axis=0),
            "metadata": np.concatenate(metadata, axis=0),
        }

    return _worker


def reconstruct_a3m_sequences() -> Callable[..., list[str]]:
    """Express the aligned matrix + deletion counts as raw a3m strings.

    Each query column keeps its (upper-cased) residue; ``deletion_matrix[i, j]``
    lower-case placeholders are inserted before column ``j`` so that the reused
    ``parse_sequence`` recovers exactly the same aligned residues and deletions.
    """

    def _worker(alignment: dict[str, np.ndarray]) -> list[str]:
        upper = np.char.upper(alignment["msa"].astype("<U1"))
        deletions = alignment["deletion_matrix"].astype(np.int64)
        raw_sequences = []
        for row, dels in zip(upper, deletions, strict=True):
            raw_sequences.append(
                "".join("a" * int(d) + c for c, d in zip(row, dels, strict=True)),
            )
        return raw_sequences

    return _worker


def parse_openfold_msa_headers() -> Callable[..., dict[str, np.ndarray]]:
    r"""Parse distillation MSA metadata into the standard header fields.

    Protein rows look like ``UniRef100_<id>\t<mmseqs stats>``; RNA rows look
    like ``<accession>/<range>``. The distillation metadata carries no species /
    representative id, so those are filled with ``N/A`` (as the a3m BFD path
    already does).
    """
    uniref = re.compile(r"^(?P<db>UniRef\d+)_(?P<id>\S+)")

    def _worker(alignment: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        database, database_id, species, rep_id = [], [], [], []
        for raw in alignment["metadata"].tolist():
            head = str(raw).split("\t", 1)[0].strip()
            m = uniref.match(head)
            if m:
                db, did = m.group("db"), m.group("id")
            elif "/" in head:
                db, did = "rnacentral", head.split("/", 1)[0]
            else:
                db, did = "query", head
            database.append(db)
            database_id.append(did)
            species.append("N/A")
            rep_id.append(did)
        return {
            "database": np.array(database),
            "database_id": np.array(database_id),
            "species": np.array(species),
            "rep_id": np.array(rep_id),
        }

    return _worker


def build_template_dict() -> Callable[..., dict[str, dict[str, object]]]:
    """Normalise per-hit template payloads into a single ``template_dict``."""

    def _worker(
        template_hits: dict[str, dict[str, object]],
    ) -> dict[str, dict[str, object]]:
        return {
            hit: {field: payload[field] for field in _TEMPLATE_FIELDS}
            for hit, payload in template_hits.items()
        }

    return _worker
