"""Transforms for OpenFold3 distillation ingest recipes.

OpenFold distillation alignments are stored split into an already-aligned
character matrix plus a separate deletion-count matrix, whereas the canonical
MSA pipeline parses raw a3m strings (lowercase = insertions). These adapters
express the distillation arrays back as a3m strings / standard headers so the
**existing** ``parse_sequence`` / ``parse_headers`` / ``build_dict``
instructions can be reused verbatim, guaranteeing identical record content.
"""

import re

import numpy as np

_UNIREF = re.compile(r"^(?P<db>UniRef\d+)_(?P<id>\S+)")


def merge_msa_sources(
    msa_sources: dict[str, dict[str, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Collapse the per-source MSA payloads into the ``(msa, deletions, metadata)`` arrays.

    Protein monomers expose a single source; RNA monomers expose several
    (``rnacentral_hits`` / ``nt_hits`` / ``rfam_hits``) that share the query
    columns, so their hit rows are concatenated below a single query row. The
    three arrays are returned as distinct outputs so each downstream instruction
    depends only on the columns it actually consumes.
    """
    sources = list(msa_sources.values())
    first = sources[0]
    if len(sources) == 1:
        return first["msa"], first["deletion_matrix"], first["metadata"]
    # keep the query (row 0) once, then every source's hit rows
    msa = [first["msa"]]
    deletion = [first["deletion_matrix"]]
    metadata = [first["metadata"]]
    for src in sources[1:]:
        msa.append(src["msa"][1:])
        deletion.append(src["deletion_matrix"][1:])
        metadata.append(src["metadata"][1:])
    return (
        np.concatenate(msa, axis=0),
        np.concatenate(deletion, axis=0),
        np.concatenate(metadata, axis=0),
    )


def reconstruct_a3m_sequences(msa: np.ndarray, deletion_matrix: np.ndarray) -> list[str]:
    """Express the aligned matrix + deletion counts as raw a3m strings.

    Each query column keeps its (upper-cased) residue; ``deletion_matrix[i, j]``
    lower-case placeholders are inserted before column ``j`` so that the reused
    ``parse_sequence`` recovers exactly the same aligned residues and deletions.
    """
    upper = np.char.upper(msa.astype("<U1"))
    deletions = deletion_matrix.astype(np.int64)
    return [
        "".join("a" * int(d) + c for c, d in zip(row, dels, strict=True))
        for row, dels in zip(upper, deletions, strict=True)
    ]


def parse_openfold_msa_headers(metadata: np.ndarray) -> dict[str, np.ndarray]:
    r"""Parse distillation MSA metadata into the standard header fields.

    Protein rows look like ``UniRef100_<id>\t<mmseqs stats>``; RNA rows look
    like ``<accession>/<range>``. The distillation metadata carries no species /
    representative id, so those are filled with ``N/A`` (as the a3m BFD path
    already does). Stored as bytes (``|S``) to match the existing a3m headers.
    """
    database, database_id, species, rep_id = [], [], [], []
    for raw in metadata.tolist():
        head = str(raw).split("\t", 1)[0].strip()
        m = _UNIREF.match(head)
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
        "database": np.array([s.encode() for s in database]),
        "database_id": np.array([s.encode() for s in database_id]),
        "species": np.array([s.encode() for s in species]),
        "rep_id": np.array([s.encode() for s in rep_id]),
    }


def reconstruct_template_alignments(
    template_hits: dict[str, dict[str, object]],
    query_len: int,
) -> dict[str, tuple[str, str]]:
    """Express each template hit's ``idx_map`` as a placeholder query/template alignment.

    ``idx_map`` lists matched ``(query_residue, template_residue)`` index pairs.
    A gapped alignment of placeholder residues is built so the reused
    ``load_templates`` / ``to_template_mol`` recover exactly those pairs over the
    full query length (residue identity is irrelevant there — only gap masks
    are used).
    """
    align_results: dict[str, tuple[str, str]] = {}
    for hit, payload in template_hits.items():
        idx_map = np.asarray(payload["idx_map"], dtype=np.int64)
        order = np.argsort(idx_map[:, 0], kind="stable")
        pairs = idx_map[order]
        query_chars: list[str] = []
        target_chars: list[str] = []
        q_cur = t_cur = 0
        for qi, ti in pairs.tolist():
            while q_cur < qi:
                query_chars.append("A")
                target_chars.append("-")
                q_cur += 1
            while t_cur < ti:
                query_chars.append("-")
                target_chars.append("A")
                t_cur += 1
            query_chars.append("A")
            target_chars.append("A")
            q_cur += 1
            t_cur += 1
        while q_cur < query_len:
            query_chars.append("A")
            target_chars.append("-")
            q_cur += 1
        align_results[hit] = ("".join(query_chars), "".join(target_chars))
    return align_results
