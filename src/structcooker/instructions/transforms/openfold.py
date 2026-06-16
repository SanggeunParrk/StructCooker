"""Transforms for OpenFold3 distillation ingest recipes.

OpenFold distillation alignments are stored split into an already-aligned
character matrix plus a separate deletion-count matrix. ``build_msa_features``
computes the canonical MSA feature arrays (query / aligned / deletions /
deletion_mean / profile) directly from those arrays with vectorised numpy,
producing the exact same output as round-tripping through a3m strings +
``parse_sequence`` but ~100x faster (no per-row Python string work).
"""

import re

import numpy as np

from structcooker.utils.mapping import ResidueMapping

_UNIREF = re.compile(r"^(?P<db>UniRef\d+)_(?P<id>\S+)")
_DELETION_CLIP = 255  # parse_sequence stores deletion counts as uint8 (0-255)


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


def cap_msa_depth(
    msa: np.ndarray,
    deletion_matrix: np.ndarray,
    metadata: np.ndarray,
    max_depth: int | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cap the MSA to at most ``max_depth`` rows.

    Deep distillation alignments (long monomers reach millions of rows) blow up
    memory; this keeps the query (row 0) plus the first ``max_depth - 1`` hits.
    A ``None`` cap or an already-shallow alignment is returned unchanged.
    """
    if max_depth is None or len(msa) <= max_depth:
        return msa, deletion_matrix, metadata
    return msa[:max_depth], deletion_matrix[:max_depth], metadata[:max_depth]


def build_msa_features(
    msa: np.ndarray,
    deletion_matrix: np.ndarray,
    a3m_type: str | None = "protein",
) -> dict[str, np.ndarray]:
    """Compute the canonical MSA feature arrays directly from the aligned arrays.

    Equivalent to round-tripping the matrix through a3m strings and
    ``parse_sequence`` (upper-cased columns become ``aligned_sequences``,
    ``deletion_matrix`` clipped to uint8 becomes ``deletions``) but vectorised:
    a 256-entry char->index LUT (built from the same ``ResidueMapping`` so the
    encoding is identical) replaces the per-row string work.
    """
    mapping = ResidueMapping()
    n_classes = mapping.MAX_INDEX + 1
    view = mapping.rna if (a3m_type or "protein").lower() == "rna" else mapping.protein

    # char -> residue index LUT over the ASCII range, using the same mapping;
    # fold lowercase onto uppercase so the upper-casing is free (no np.char.upper).
    lut = view.map(np.array([chr(i) for i in range(256)], dtype="<U1")).astype(np.uint8)
    lut[ord("a") : ord("z") + 1] = lut[ord("A") : ord("Z") + 1]
    codes = np.ascontiguousarray(msa.astype("<U1")).view(np.uint32)
    aligned_sequences = lut[codes].astype(np.uint8)

    n_rows, length = aligned_sequences.shape
    deletions = np.clip(deletion_matrix, 0, _DELETION_CLIP).astype(np.int32)
    deletion_mean = (2 * np.arctan(deletions.astype(np.float32) / 3) / np.pi).mean(
        axis=0,
    ).astype(np.float32)
    # profile[l, k] = fraction of rows with residue k at column l (== one-hot mean),
    # via a single bincount instead of an (N, L, K) one-hot intermediate.
    flat = aligned_sequences.astype(np.int64) + np.arange(length) * n_classes
    counts = np.bincount(flat.ravel(), minlength=length * n_classes)
    profile = (counts.reshape(length, n_classes) / n_rows).astype(np.float32)

    # Query keeps the row-0 a3m string (upper columns + lowercase insertions).
    q_upper, q_del = np.char.upper(msa[0].astype("<U1")), deletion_matrix[0].astype(np.int64)
    query_string = "".join(
        "a" * int(d) + c for c, d in zip(q_upper, q_del, strict=True)
    )
    query_sequence = np.array(list(query_string))

    return {
        "query_sequence": query_sequence,
        "aligned_sequences": aligned_sequences,
        "deletions": deletions,
        "deletion_mean": deletion_mean,
        "profile": profile,
    }


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
