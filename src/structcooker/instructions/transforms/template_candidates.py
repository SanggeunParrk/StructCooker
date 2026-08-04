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
