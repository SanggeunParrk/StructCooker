"""Template Phase 1/2 precompute (seq_id-deduplicated candidate selection).

Absorbs ``scripts/maintenance/precompute_template_candidates.py`` and
``precompute_template_seqs.py`` — the kalign-free, CIF-decode-free precompute that
turns raw hmmsearch outputs into the small per-seq_id maps Phase 3 needs. The two
public entry points (:func:`precompute_candidates`, :func:`precompute_seqs`) each
read their inputs and **write their outputs as side effects**, returning a short
status string; they are driven by a materialize recipe whose ``output_data_path``
(the primary TSV) doubles as the build-done marker.

Logic preserved verbatim from the scripts (date filters, top-k, reduced-hmm framing)
so the outputs reproduce exactly.
"""
from __future__ import annotations

import datetime as dt
import re
from collections import defaultdict
from pathlib import Path
from typing import cast

from joblib import Parallel, delayed

# hmmsearch full-sequence summary row: capture the target name.
_HIT_RE = re.compile(
    r"""
    ^\s*
    [0-9.eE+-]+ \s+[0-9.]+ \s+[0-9.]+     # full-seq E-value, score, bias
    \s+[0-9.eE+-]+ \s+[0-9.]+ \s+[0-9.]+  # best-domain E-value, score, bias
    \s+[0-9.]+ \s+\d+                      # exp, N
    \s+(\S+)                               # <-- target name (capture)
    \s+\|                                  # description starts with '|'
    """,
    re.MULTILINE | re.VERBOSE,
)


def _reduce_id(raw: str) -> str:
    """Reduce ``1B4U_A_.`` to ``1B4U_A`` (pdbid_chain)."""
    parts = raw.split("_")
    return "_".join(parts[:2])


def _parse_date(s: str) -> dt.date | None:
    try:
        return dt.date.fromisoformat(s.strip())
    except (ValueError, AttributeError):
        return None


def _load_pdb_dates(path: str | Path) -> dict[str, dt.date]:
    r"""Map pdbid(lower) -> cutoff date, preferring ``release_date``.

    Header-aware: falls back to ``deposition_date``, then to legacy positional
    (``cif_id\tresolution\tdeposition_date``).
    """
    dates: dict[str, dt.date] = {}
    with Path(path).open(encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(head)}
        id_i = col.get("pdbid", col.get("cif_id", 0))
        if "release_date" in col:
            date_i = col["release_date"]
        elif "deposition_date" in col:
            date_i = col["deposition_date"]
        else:
            date_i = 2
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(id_i, date_i):
                continue
            d = _parse_date(parts[date_i])
            if d is not None:
                dates[parts[id_i].split("_")[0].lower()] = d
    return dates


def _read_chain_list(path: str | Path) -> list[str]:
    with Path(path).open(encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def _parse_cif_fasta(path: str | Path, wanted: set[str]) -> dict[str, str]:
    """Map chain ``{pdbid}_{chain}`` -> sequence, only for chains in ``wanted``."""
    chain_seq: dict[str, str] = {}
    cur: str | None = None
    buf: list[str] = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None and cur in wanted:
                    chain_seq[cur] = "".join(buf)
                header = line[1:].split("|", 1)[0].strip()
                parts = header.split("_")
                cur = "_".join(parts[:2]) if len(parts) >= 2 else header
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None and cur in wanted:
        chain_seq[cur] = "".join(buf)
    return chain_seq


def _stream_seq_to_id(path: str | Path, needed: set[str]) -> dict[str, str]:
    """Map sequence -> seq_id, streamed from seq_id_map.tsv, filtered to ``needed``."""
    seq_to_id: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in needed:
                seq_to_id[parts[1]] = parts[0]
    return seq_to_id


def _hmm_path(hmm_dir: str | Path, seq_id: str) -> Path:
    return Path(hmm_dir) / seq_id[0] / seq_id[-3:] / f"{seq_id}.out"


def _ordered_hits(text: str) -> list[str]:
    """Return ordered (best-score-first) unique template chain ids from .out text."""
    hits: list[str] = []
    seen: set[str] = set()
    for raw in _HIT_RE.findall(text):
        tid = _reduce_id(raw)
        if tid not in seen:
            seen.add(tid)
            hits.append(tid)
    return hits


def _reduce_hmm_text(text: str, keep: set[str]) -> str:
    """Return the hmm .out reduced to only the target ids in ``keep``.

    Keeps every header / table-header / footer line untouched; drops summary rows
    and ``>>`` domain blocks for non-kept targets, so the result parses exactly like
    the full file, only smaller.
    """
    out: list[str] = []
    in_domains = False
    keep_block = True
    for line in text.splitlines(keepends=True):
        if line.startswith(">>"):
            in_domains = True
            target = _reduce_id(line[2:].split("|")[0].strip())
            keep_block = target in keep
            if keep_block:
                out.append(line)
            continue
        if in_domains:
            if line.startswith(("Internal pipeline", "//")):
                in_domains = False
                out.append(line)
            elif keep_block:
                out.append(line)
            continue
        m = _HIT_RE.match(line)
        if m is not None:
            if _reduce_id(m.group(1)) in keep:
                out.append(line)
            continue
        out.append(line)
    return "".join(out)


def _process_seqid(
    seq_id: str,
    chains: list[str],
    hmm_dir: str | Path,
    reduced_dir: str | Path,
    pdb_dates: dict[str, dt.date],
    date_cutoff: dt.date,
    day_diff: int,
    topk: int,
) -> list[tuple[str, list[str]]]:
    """Select each chain's <=topk templates + write the seq_id's reduced hmm .out."""
    hmm_path = _hmm_path(hmm_dir, seq_id)
    if not hmm_path.exists():
        return [(c, []) for c in chains]
    text = hmm_path.read_text(encoding="utf-8")
    hits = _ordered_hits(text)

    per_chain: list[tuple[str, list[str]]] = []
    union: set[str] = set()
    for chain in chains:
        q_date = pdb_dates.get(chain.split("_")[0].lower())
        if q_date is None or not hits:
            per_chain.append((chain, []))
            continue
        kept: list[str] = []
        for tid in hits:
            t_date = pdb_dates.get(tid.split("_")[0].lower())
            # AF3: template released <= max_template_date AND >= day_diff days
            # BEFORE (older than) the query structure (anti-leakage).
            if t_date is None or t_date > date_cutoff:
                continue
            if (q_date - t_date).days < day_diff:
                continue
            kept.append(tid)
            if len(kept) >= topk:
                break
        per_chain.append((chain, kept))
        union.update(kept)

    if union:
        reduced = _reduce_hmm_text(text, union)
        out_path = Path(reduced_dir) / seq_id[0] / seq_id[-3:] / f"{seq_id}.reduced.out"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(reduced, encoding="utf-8")
    return per_chain


def _phase1(
    cif_fasta: str | Path,
    chain_list: str | Path,
    seq_id_map: str | Path,
    out_seqid_chains: str | Path,
) -> dict[str, list[str]]:
    """Group protein chains by seq_id (write + reuse ``out_seqid_chains``)."""
    out_path = Path(out_seqid_chains)
    if out_path.exists():
        seqid_chains: dict[str, list[str]] = {}
        with out_path.open(encoding="utf-8") as f:
            for line in f:
                sid, cs = line.rstrip("\n").split("\t")
                seqid_chains[sid] = cs.split(",") if cs else []
        return seqid_chains

    chains = _read_chain_list(chain_list)
    chain_seq = _parse_cif_fasta(cif_fasta, set(chains))
    seq_to_id = _stream_seq_to_id(seq_id_map, set(chain_seq.values()))

    grouped: dict[str, list[str]] = defaultdict(list)
    for chain, seq in chain_seq.items():
        sid = seq_to_id.get(seq)
        if sid is not None:
            grouped[sid].append(chain)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for sid, cs in grouped.items():
            out.write(f"{sid}\t{','.join(cs)}\n")
    return dict(grouped)


def precompute_candidates(
    cif_fasta: str | Path,
    chain_list: str | Path,
    seq_id_map: str | Path,
    pdb_dates: str | Path,
    hmm_dir: str | Path,
    out_seqid_chains: str | Path,
    out_chain_templates: str | Path,
    out_reduced_hmm_dir: str | Path,
    date_cutoff: str = "2021-09-30",
    day_diff: int = 60,
    topk: int = 20,
    n_jobs: int = 112,
) -> str:
    """Phase 1 (seq_id -> chains) + Phase 2 (chain -> <=topk templates + reduced hmm).

    Writes ``out_seqid_chains``, ``out_chain_templates``, and the per-seq_id reduced
    hmm .out files under ``out_reduced_hmm_dir`` as side effects; returns a status
    line. Faithful port of ``precompute_template_candidates.py``.
    """
    cutoff = dt.date.fromisoformat(date_cutoff)
    seqid_chains = _phase1(cif_fasta, chain_list, seq_id_map, out_seqid_chains)
    dates = _load_pdb_dates(pdb_dates)

    items = sorted(seqid_chains.items())
    results = cast(
        "list[list[tuple[str, list[str]]]]",
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(_process_seqid)(
                sid, cs, hmm_dir, out_reduced_hmm_dir, dates, cutoff, day_diff, topk,
            )
            for sid, cs in items
        ),
    )

    n_chain = n_with = 0
    out_path = Path(out_chain_templates)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for per_seqid in results:
            for chain, tids in per_seqid:
                n_chain += 1
                if tids:
                    n_with += 1
                out.write(f"{chain}\t{','.join(tids)}\n")
    return f"candidates: {n_chain} chains, {n_with} with >=1 template -> {out_chain_templates}"


def precompute_seqs(
    seqid_chains: str | Path,
    seq_id_map: str | Path,
    cif_fasta: str | Path,
    out_seqid_seq: str | Path,
    out_chain_seq: str | Path,
) -> str:
    """Write the lean sequence maps Phase 3 needs (seqid->seq, chain->seq).

    ``seqid_to_seq`` covers only the query seq_ids in ``seqid_chains``;
    ``chain_to_seq`` covers polypeptide(L) template chains. Faithful port of
    ``precompute_template_seqs.py``.
    """
    needed: set[str] = set()
    with Path(seqid_chains).open(encoding="utf-8") as f:
        for line in f:
            sid = line.split("\t", 1)[0]
            if sid:
                needed.add(sid)

    n = 0
    Path(out_seqid_seq).parent.mkdir(parents=True, exist_ok=True)
    with Path(seq_id_map).open(encoding="utf-8") as f, Path(out_seqid_seq).open("w", encoding="utf-8") as out:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] in needed:
                out.write(f"{p[0]}\t{p[1]}\n")
                n += 1

    m = 0
    cur: str | None = None
    is_prot = False
    buf: list[str] = []
    Path(out_chain_seq).parent.mkdir(parents=True, exist_ok=True)
    with Path(cif_fasta).open(encoding="utf-8") as f, Path(out_chain_seq).open("w", encoding="utf-8") as out:
        for line in f:
            if line.startswith(">"):
                if cur is not None and is_prot:
                    out.write(f"{cur}\t{''.join(buf)}\n")
                    m += 1
                header = line[1:]
                cur = "_".join(header.split("|", 1)[0].strip().split("_")[:2])
                is_prot = "polypeptide(L)" in header
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None and is_prot:
            out.write(f"{cur}\t{''.join(buf)}\n")
            m += 1
    return f"seqs: seqid_to_seq={n}, chain_to_seq={m}"


def _select_seqid(
    seq_id: str,
    hmm_dir: str | Path,
    reduced_dir: str | Path,
    pdb_dates: dict[str, dt.date],
    date_cutoff: dt.date,
    max_candidates: int,
) -> tuple[str, list[str]]:
    """Date-filter one seq_id's hits (template release <= cutoff) and write its reduced hmm."""
    hmm_path = _hmm_path(hmm_dir, seq_id)
    if not hmm_path.exists():
        return seq_id, []
    text = hmm_path.read_text(encoding="utf-8")
    kept: list[str] = []
    for tid in _ordered_hits(text):
        t_date = pdb_dates.get(tid.split("_")[0].lower())
        if t_date is None or t_date > date_cutoff:
            continue
        kept.append(tid)
        if len(kept) >= max_candidates:
            break
    if kept:
        out_path = Path(reduced_dir) / seq_id[0] / seq_id[-3:] / f"{seq_id}.reduced.out"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(_reduce_hmm_text(text, set(kept)), encoding="utf-8")
    return seq_id, kept


def _select_chunk(
    seq_ids: list[str],
    hmm_dir: str | Path,
    reduced_dir: str | Path,
    pdb_dates: dict[str, dt.date],
    date_cutoff: dt.date,
    max_candidates: int,
) -> list[tuple[str, list[str]]]:
    return [_select_seqid(s, hmm_dir, reduced_dir, pdb_dates, date_cutoff, max_candidates) for s in seq_ids]


def precompute_seqid_candidates(
    seq_ids_path: str | Path,
    pdb_dates: str | Path,
    hmm_dir: str | Path,
    cif_fasta: str | Path,
    out_seqid_templates: str | Path,
    out_seqid_list: str | Path,
    out_chain_seq: str | Path,
    out_reduced_hmm_dir: str | Path,
    date_cutoff: str = "2021-09-30",
    max_candidates: int = 60,
    n_jobs: int = 112,
) -> str:
    """Phase 1/2 for a set whose queries carry no release date (predicted structures).

    The PDB path keys candidates by chain because each query chain has its own release
    date (templates must be >= 60 days older). A predicted query has none, so the only
    rule is the global template cutoff (release <= ``date_cutoff``), every chain of a
    sequence gets the same candidates, and selection is per seq_id. Up to
    ``max_candidates`` date-passing hits are kept in e-value order -- more than the 20
    finally used, so Phase 3 can backfill hits its coverage filter drops.

    ``seq_ids_path`` is a TSV whose first column is the query seq_id (e.g. the set's
    seq_id_map subset). Writes, as side effects: ``out_seqid_templates`` (seq_id ->
    candidates), ``out_seqid_list`` (seq_ids with >= 1 candidate: Phase 3's work list),
    ``out_chain_seq`` (template chain -> sequence, polypeptide(L)), and the reduced hmm
    per seq_id.
    """
    cutoff = dt.date.fromisoformat(date_cutoff)
    with Path(seq_ids_path).open(encoding="utf-8") as f:
        seq_ids = sorted({line.split("\t", 1)[0].strip() for line in f if line.strip()})
    dates = _load_pdb_dates(pdb_dates)
    # One task per chunk, not per seq_id: every task pickles ``dates`` (~250k entries),
    # and at one task per seq_id that dispatch serialised the whole run at ~200 seq_ids
    # per minute whatever the core count.
    chunk = max(1, -(-len(seq_ids) // (n_jobs * 8)))
    chunks = [seq_ids[i:i + chunk] for i in range(0, len(seq_ids), chunk)]
    per_chunk = cast(
        "list[list[tuple[str, list[str]]]]",
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(_select_chunk)(c, hmm_dir, out_reduced_hmm_dir, dates, cutoff, max_candidates)
            for c in chunks
        ),
    )
    results = [r for rs in per_chunk for r in rs]
    Path(out_seqid_templates).parent.mkdir(parents=True, exist_ok=True)
    n_with = 0
    with Path(out_seqid_templates).open("w", encoding="utf-8") as out, \
            Path(out_seqid_list).open("w", encoding="utf-8") as lst:
        for sid, tids in results:
            out.write(f"{sid}\t{','.join(tids)}\n")
            if tids:
                n_with += 1
                lst.write(f"{sid}\n")
    m = _write_chain_seqs(cif_fasta, out_chain_seq)
    return (f"candidates: {len(results)} seq_ids, {n_with} with >=1 template; "
            f"chain_to_seq={m} -> {out_seqid_templates}")


def _write_chain_seqs(cif_fasta: str | Path, out_chain_seq: str | Path) -> int:
    """Write ``{pdbid}_{chain} -> sequence`` for polypeptide(L) chains (as precompute_seqs)."""
    m = 0
    cur: str | None = None
    is_prot = False
    buf: list[str] = []
    Path(out_chain_seq).parent.mkdir(parents=True, exist_ok=True)
    with Path(cif_fasta).open(encoding="utf-8") as f, Path(out_chain_seq).open("w", encoding="utf-8") as out:
        for line in f:
            if line.startswith(">"):
                if cur is not None and is_prot:
                    out.write(f"{cur}\t{''.join(buf)}\n")
                    m += 1
                header = line[1:]
                cur = "_".join(header.split("|", 1)[0].strip().split("_")[:2])
                is_prot = "polypeptide(L)" in header
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None and is_prot:
            out.write(f"{cur}\t{''.join(buf)}\n")
            m += 1
    return m

