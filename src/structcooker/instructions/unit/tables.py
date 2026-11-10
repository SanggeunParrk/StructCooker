"""Tables of the unit DB: manifest, manifest with clusters + report, date split.

Each function is the instruction of one stage of the ``unit/`` set:

=============================  =======================================================  =====================
stage                          function                                                  writes
=============================  =======================================================  =====================
``unit/manifest``         :func:`write_manifest` (over ``record.manifest_row``)     ``manifest.tsv``
``unit/manifest_full``    :func:`write_manifest_full`                               ``manifest_full.tsv``, ``report.md``
``unit/split``            :func:`write_split`                                       ``split.tsv``, ``split_summary.md``
=============================  =======================================================  =====================

Protein node hashes (md5-16 of ``"P:<sequence>"``) are mapped back to sequences through the PDB
protein FASTA, then to seq_ids through ``seq_id_map.tsv`` (one streaming pass), then to the 40%
sequence clusters of ``seq_cluster40_PDB.tsv``.
"""

from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable

Row = dict[str, Any]

TRAIN_END, VALID_END = "2021-09-30", "2023-01-12"
NO_VALUE = ("", "None", "nan")


def _text(value: Any) -> str:
    return "" if value is None else str(value)


def _write_rows(path: Path, columns: list[str], rows: list[Row]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(_text(row[c]) for c in columns) + "\n")


def read_tsv(path: Path | str) -> list[Row]:
    """Rows of a tab-separated file with a header line."""
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ----------------------------------------------------------------------------- manifest
def write_manifest(data: dict[str, dict[str, Row]], output_path: Path) -> None:
    """Materializer of ``unit/manifest``: one row per unit, sorted by unit id.

    ``data`` maps a record key to ``{"manifest_row": row}`` (``record.manifest_row`` of that record).
    """
    rows = sorted(
        (v["manifest_row"] for v in data.values()), key=lambda r: r["unit_id"]
    )
    columns = list(rows[0]) if rows else ["unit_id"]
    _write_rows(Path(output_path), columns, rows)


# ----------------------------------------------------------------------------- clusters
def protein_hashes(rows: list[Row]) -> set[str]:
    """Hashes of every protein node of the manifest rows."""
    return {
        h
        for r in rows
        for h, k in zip(
            r["node_hashes"].split("|"), r["node_kinds"].split("|"), strict=True
        )
        if k == "protein"
    }


def hash_to_sequence(fasta_path: Path | str, wanted: set[str]) -> dict[str, str]:
    """Sequences of ``wanted`` node hashes, from a one-line-per-sequence FASTA."""
    found: dict[str, str] = {}
    with Path(fasta_path).open() as handle:
        for line in handle:
            if line[0] != ">":
                s = line.strip()
                h = hashlib.md5(("P:" + s).encode()).hexdigest()[:16]  # noqa: S324  (a key, not security)
                if h in wanted:
                    found[h] = s
    return found


def sequence_ids(seq_id_map: Path | str, sequences: set[str]) -> dict[str, str]:
    """``sequence -> protein seq_id`` for ``sequences``, in one pass over ``seq_id_map.tsv``."""
    out: dict[str, str] = {}
    with Path(seq_id_map).open() as handle:
        for line in handle:
            if line[0] != "P":
                continue
            sid, _, s = line.rstrip("\n").partition("\t")
            if s in sequences:
                out[s] = sid
    return out


def cluster_of(cluster_tsv: Path | str, ids: set[str]) -> dict[str, str]:
    """``seq_id -> cluster id`` for ``ids`` from a ``cluster<TAB>member,member`` TSV."""
    out: dict[str, str] = {}
    with Path(cluster_tsv).open() as handle:
        for line in handle:
            c, _, members = line.rstrip("\n").partition("\t")
            for m in members.split(","):
                if m in ids:
                    out[m] = c
    return out


def annotate_rows(
    rows: list[Row],
    hash2id: dict[str, str | None],
    id2cluster: dict[str, str],
    including_na: dict[str, str],
) -> tuple[list[Row], int]:
    """Add node_seq_ids, node_c40, cluster_key40, n_unmapped_protein, including_NA_entry.

    Returns the new rows and the number of protein nodes without a seq_id.
    """
    unmapped_nodes = 0
    out = []
    for r in rows:
        sid, c40, unm = [], [], 0
        for h, k in zip(
            r["node_hashes"].split("|"), r["node_kinds"].split("|"), strict=True
        ):
            if k != "protein":
                sid.append("")
                c40.append("")
                continue
            i = hash2id.get(h)
            if i is None:
                unm += 1
            sid.append(i or "")
            c40.append(id2cluster.get(i, "") if i else "")
        unmapped_nodes += unm
        out.append(
            {
                **r,
                "node_seq_ids": "|".join(sid),
                "node_c40": "|".join(c40),
                "cluster_key40": "+".join(sorted({x for x in c40 if x})),
                "n_unmapped_protein": unm,
                "including_NA_entry": including_na.get(r["pdb_id"], ""),
            }
        )
    return out, unmapped_nodes


def entry_flags(cif_metadata: Path | str) -> tuple[dict[str, str], int]:
    """``pdb id -> including_NA`` (first row of the entry) and the number of entries."""
    flags: dict[str, str] = {}
    for r in read_tsv(cif_metadata):
        flags.setdefault(r["cif_id"].split("_")[0], r["including_NA"])
    return flags, len(flags)


# ----------------------------------------------------------------------------- report
def _dist(vals: list[float]) -> str:
    vals = sorted(vals)
    if not vals:
        return "-"

    def q(p: float) -> float:
        return vals[min(len(vals) - 1, int(p * len(vals)))]

    return (
        f"min {vals[0]:g} · p25 {q(0.25):g} · median {q(0.5):g} · "
        f"p75 {q(0.75):g} · p95 {q(0.95):g} · max {vals[-1]:g}"
    )


CLASS_TOTALS: list[tuple[str, Callable[[Row], bool]]] = [
    (
        "protein-protein only (no NA / ligand / ion)",
        lambda r: r["n_na"] == "0" and r["n_ligand"] == "0" and r["n_ion"] == "0",
    ),
    ("with NA", lambda r: r["n_na"] != "0"),
    ("with ligand", lambda r: r["n_ligand"] != "0"),
    ("with ion", lambda r: r["n_ion"] != "0"),
    ("with antibody / TCR", lambda r: r["has_ab"] == "1"),
    ("with a base-paired duplex", lambda r: r["n_duplex_pairs"] != "0"),
    ("truncated (member cap hit)", lambda r: r["truncated"] == "True"),
]
SIZE_COLUMNS = [
    ("n_nodes", "nodes per unit"),
    ("n_atoms", "atoms per unit"),
    ("n_residues", "residues per unit"),
    ("n_interface_res", "interface residues"),
    ("n_atom_pairs", "contacting atom pairs"),
    ("n_water", "waters per unit"),
    ("n_bridging_water", "bridging waters"),
]


def build_report(rows: list[Row], n_entries_total: int, unmapped_nodes: int) -> str:
    """Return the statistics page (``report.md``) of the annotated manifest rows."""
    n_with_units = len({r["pdb_id"] for r in rows})
    out = [
        "# Unit DB report",
        "",
        f"Units: **{len(rows):,}** from **{n_with_units:,}** PDB entries "
        f"({n_entries_total - n_with_units:,} of the {n_entries_total:,} entries of the PDB DB gave no unit).",
        "",
        "## Units by type",
        "",
        "| type | units | share |",
        "|---|---:|---:|",
    ]
    out += [
        f"| {t} | {n:,} | {100 * n / len(rows):.1f}% |"
        for t, n in Counter(r["type"] for r in rows).most_common(25)
    ]
    out += ["", "## Class totals", ""]
    out += [
        f"- {lab}: {sum(1 for r in rows if pred(r)):,}" for lab, pred in CLASS_TOTALS
    ]
    out += ["", "## Size distributions", ""]
    out += [f"- {lab}: {_dist([float(r[c]) for r in rows])}" for c, lab in SIZE_COLUMNS]
    res = [float(r["resolution"]) for r in rows if r["resolution"] not in NO_VALUE]
    out += [
        f"- resolution (A): {_dist(res)}  ({len(rows) - len(res):,} units without resolution: NMR / EM without value)",
        "",
        "## Water",
        "",
        f"- units with no water: {sum(1 for r in rows if r['n_water'] == '0'):,}",
        f"- units with >= 1 bridging water: {sum(1 for r in rows if r['n_bridging_water'] != '0'):,}",
        "",
        "## Cluster coverage",
        "",
        f"- protein nodes without a seq_id: {unmapped_nodes:,}",
        f"- units with every protein node mapped: {sum(1 for r in rows if str(r['n_unmapped_protein']) == '0'):,}",
        f"- distinct c40 cluster keys: {len({r['cluster_key40'] for r in rows}):,}",
        "",
    ]
    ccd: Counter[str] = Counter()
    ions: Counter[str] = Counter()
    for r in rows:
        for k, c in zip(
            r["node_kinds"].split("|"), r["node_ccd"].split("|"), strict=True
        ):
            if k == "ligand":
                ccd[c] += 1
            elif k == "ion":
                ions[c] += 1
    out += [
        "## Most frequent ligands / ions (nodes in units)",
        "",
        "ligands: " + ", ".join(f"{c} {n:,}" for c, n in ccd.most_common(30)),
        "",
        "ions: " + ", ".join(f"{c} {n:,}" for c, n in ions.most_common(15)),
        "",
    ]
    return "\n".join(out)


def write_manifest_full(
    manifest_path: Path | str,
    fasta_path: Path | str,
    seq_id_map_path: Path | str,
    cluster_path: Path | str,
    cif_metadata_path: Path | str,
    out_path: Path | str,
    out_report: Path | str,
) -> str:
    """Stage ``unit/manifest_full``: the manifest with seq_ids and 40% clusters, plus ``report.md``."""
    rows = read_tsv(manifest_path)
    h2s = hash_to_sequence(fasta_path, protein_hashes(rows))
    s2id = sequence_ids(seq_id_map_path, set(h2s.values()))
    id2c = cluster_of(cluster_path, set(s2id.values()))
    h2id = {h: s2id.get(s) for h, s in h2s.items()}
    including_na, n_entries = entry_flags(cif_metadata_path)
    full, unmapped = annotate_rows(rows, h2id, id2c, including_na)
    columns = [
        *(list(rows[0]) if rows else []),
        "node_seq_ids",
        "node_c40",
        "cluster_key40",
        "n_unmapped_protein",
        "including_NA_entry",
    ]
    _write_rows(Path(out_path), columns, full)
    Path(out_report).parent.mkdir(parents=True, exist_ok=True)
    Path(out_report).write_text(build_report(full, n_entries, unmapped))
    return (
        f"{len(full)} units, {len(h2id)} protein sequences, {len(id2c)} clustered ids"
    )


# ----------------------------------------------------------------------------- split
def classify_unit(release_date: str, resolution: str) -> tuple[str, str]:
    """``(split, reason)`` of one unit, AF3 style.

    train: release < 2021-09-30, resolution <= 9.0 (the rule of train_20210930); valid:
    2021-09-30 .. 2023-01-12 and test: later, both with resolution <= 4.5; ``none``: no release date,
    no resolution or above the cut.
    """
    if release_date in NO_VALUE:
        return "none", "no_release_date"
    if resolution in NO_VALUE:
        return "none", "no_resolution"
    res = float(resolution)
    if release_date < TRAIN_END:
        return ("train", "") if res <= 9.0 else ("none", "resolution>9")
    if res > 4.5:
        return "none", "resolution>4.5"
    return ("valid", "") if release_date <= VALID_END else ("test", "")


def _ligands(r: Row) -> list[str]:
    return [
        c
        for k, c in zip(
            r["node_kinds"].split("|"), r["node_ccd"].split("|"), strict=True
        )
        if k == "ligand"
    ]


def split_rows(rows: list[Row]) -> list[Row]:
    """Rows with ``split`` / ``reason`` and, for valid and test units, the leakage flags.

    Flags against the training units (nothing is dropped): ``n_prot_c40_in_train`` (protein nodes
    whose 40% cluster has a node in train), ``novel_protein`` (none of them), ``novel_pair`` (the
    unit's set of clusters does not occur in train), ``ligands_seen_in_train`` (every ligand CCD id
    occurs in train; empty without ligands).
    """
    rows = [dict(r) for r in rows]
    for r in rows:
        r["split"], r["reason"] = classify_unit(r["release_date"], r["resolution"])
    tr_c: set[str] = set()
    tr_pair: set[str] = set()
    tr_lig: set[str] = set()
    for r in rows:
        if r["split"] == "train":
            tr_c.update(c for c in r["node_c40"].split("|") if c)
            tr_pair.add(r["cluster_key40"])
            tr_lig.update(_ligands(r))
    for r in rows:
        n_in = novel_p = novel_pair = lig_seen = ""
        if r["split"] in ("valid", "test"):
            n_in = sum(c in tr_c for c in r["node_c40"].split("|") if c)
            novel_p = int(n_in == 0)
            novel_pair = int(r["cluster_key40"] not in tr_pair)
            ligs = _ligands(r)
            lig_seen = int(all(c in tr_lig for c in ligs)) if ligs else ""
        r["n_prot_c40_in_train"] = n_in
        r["novel_protein"] = novel_p
        r["novel_pair"] = novel_pair
        r["ligands_seen_in_train"] = lig_seen
    return rows


SPLIT_COLUMNS = [
    "unit_id",
    "pdb_id",
    "split",
    "reason",
    "release_date",
    "resolution",
    "type",
    "n_prot_c40_in_train",
    "n_protein",
    "novel_protein",
    "novel_pair",
    "ligands_seen_in_train",
    "cluster_key40",
]


def split_summary(rows: list[Row]) -> str:
    """Return the ``split_summary.md`` page of rows from :func:`split_rows`."""
    cnt = Counter(r["split"] for r in rows)
    excluded = Counter(r["reason"] for r in rows if r["split"] == "none")
    flags: Counter[tuple[str, str]] = Counter()
    for r in rows:
        if r["split"] in ("valid", "test"):
            flags[(r["split"], "novel_protein")] += r["novel_protein"]
            flags[(r["split"], "novel_pair")] += r["novel_pair"]
    entries = {
        s: len({r["pdb_id"] for r in rows if r["split"] == s})
        for s in ("train", "valid", "test")
    }
    lines = [
        "# Split summary",
        "",
        f"train: release < {TRAIN_END}, resolution <= 9.0; valid: {TRAIN_END} .. {VALID_END}, resolution <= 4.5; "
        f"test: after {VALID_END}, resolution <= 4.5",
        "",
        "| split | units | entries |",
        "|---|---:|---:|",
    ]
    lines += [
        f"| {s} | {cnt[s]:,} | {entries.get(s, '') if s != 'none' else ''} |"
        for s in ("train", "valid", "test", "none")
    ]
    lines += [
        "",
        "excluded: " + ", ".join(f"{k} {v:,}" for k, v in excluded.items()),
        "",
        "| split | units with no train-seen protein cluster | units with a novel cluster pair |",
        "|---|---:|---:|",
    ]
    lines += [
        f"| {s} | {flags[(s, 'novel_protein')]:,} | {flags[(s, 'novel_pair')]:,} |"
        for s in ("valid", "test")
    ]
    return "\n".join(lines) + "\n"


def write_split(
    manifest_full_path: Path | str, out_path: Path | str, out_summary: Path | str
) -> str:
    """Stage ``unit/split``: ``split.tsv`` (every unit, its split and leakage flags) + a summary."""
    rows = split_rows(read_tsv(manifest_full_path))
    _write_rows(Path(out_path), SPLIT_COLUMNS, rows)
    Path(out_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(out_summary).write_text(split_summary(rows))
    counts = Counter(r["split"] for r in rows)
    return ", ".join(f"{s} {counts[s]}" for s in ("train", "valid", "test", "none"))
