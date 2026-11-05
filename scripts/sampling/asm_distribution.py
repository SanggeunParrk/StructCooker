"""How the items of each interface class are spread over assembly x model x alt.

An item row is (entry, assembly, model, alt, chain pair). The DB keeps one cifmol per
"assembly_model_alt" key, so an entry with 20 NMR models, or with alternate locations, or with several
bioassemblies, contributes many rows for the same chain-cluster pair. This script measures, per
interface class (MiniWorld's buckets, protein_rna among them):

* per entry: how many (assembly, model, alt) keys carry rows of the class;
* per (entry, assembly): how many (model, alt) combinations;
* what share of the rows, and of the sampling mass, sits on a model != 1, an alt != '.' or an
  assembly != '1' row. Both samplers draw uniformly inside a cluster, so the mass of such rows in a
  class is the mean over the class's item clusters of the share of their rows that are such rows.

Reads the interface tables written by ``extract_items.py`` (chain-pair rows with their cluster ids).

    python asm_distribution.py --items <dir> --out report.md
"""

from __future__ import annotations

import argparse
import glob
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_sampling as cs

BINS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 10**9)]


def bname(lo: int, hi: int) -> str:
    return f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10**9 else f">={lo}")


def split_asm(key: str) -> tuple[str, str, str]:
    a, m, alt = key.rsplit("_", 2)
    return a, m, alt


def letter(cluster: str) -> str:
    return cluster.split("_", 1)[1][0] if "_" in cluster else cluster[1]


def hist(values: np.ndarray) -> list[float]:
    return [
        float(((values >= lo) & (values <= hi)).mean()) if len(values) else 0.0
        for lo, hi in BINS
    ]


class Pool:
    """Per class: (entry -> {asm key -> rows}) and per item cluster the row / duplicate counts."""

    def __init__(self) -> None:
        self.by_entry: dict[str, dict[str, dict[str, int]]] = defaultdict(
            lambda: defaultdict(lambda: defaultdict(int))
        )
        self.edge: dict[tuple[str, str, str], list[int]] = defaultdict(
            lambda: [0, 0, 0, 0, 0]
        )  # rows, model!=1, alt!=., asm!=1, any
        self.edge_entries: dict[tuple[str, str, str], set[str]] = defaultdict(set)

    def add(self, cls: str, entry: str, asm: str, c1: str, c2: str) -> None:
        self.by_entry[cls][entry][asm] += 1
        a, m, alt = split_asm(asm)
        e = self.edge[(cls, *sorted((c1, c2)))]
        e[0] += 1
        nm, na, nA = m != "1", alt != ".", a != "1"
        e[1] += nm
        e[2] += na
        e[3] += nA
        e[4] += nm or na or nA
        self.edge_entries[(cls, *sorted((c1, c2)))].add(entry)


def summarise(pool: Pool, classes: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for cls in classes:
        ents = pool.by_entry.get(cls)
        if not ents:
            continue
        rows = np.array([sum(v.values()) for v in ents.values()])
        nkeys = np.array([len(v) for v in ents.values()])
        nasm, multi_model, with_alt = [], 0, 0
        combos, models_per_asm, alts_per_model = [], [], []
        rows_nm = rows_na = rows_nA = 0
        for keys in ents.values():
            parsed = [(split_asm(k), r) for k, r in keys.items()]
            asms = {p[0][0] for p in parsed}
            nasm.append(len(asms))
            multi_model += any(p[0][1] != "1" for p in parsed)
            with_alt += any(p[0][2] != "." for p in parsed)
            per_asm: dict[str, set[tuple[str, str]]] = defaultdict(set)
            per_am: dict[tuple[str, str], set[str]] = defaultdict(set)
            for (a, m, alt), r in parsed:
                per_asm[a].add((m, alt))
                per_am[(a, m)].add(alt)
                rows_nm += r if m != "1" else 0
                rows_na += r if alt != "." else 0
                rows_nA += r if a != "1" else 0
            combos += [len(s) for s in per_asm.values()]
            models_per_asm += [len({m for m, _ in s}) for s in per_asm.values()]
            alts_per_model += [len(s) for s in per_am.values()]
        edges = [v for k, v in pool.edge.items() if k[0] == cls]
        er = np.array([e[0] for e in edges], dtype=float)
        frac = {
            name: float(np.mean([e[i] / e[0] for e in edges]))
            for name, i in (("model", 1), ("alt", 2), ("assembly", 3), ("any", 4))
        }
        n_entries_per_edge = np.array(
            [len(pool.edge_entries[k]) for k in pool.edge if k[0] == cls], dtype=float
        )
        tot = int(rows.sum())
        out[cls] = {
            "entries": len(ents),
            "rows": tot,
            "rows_per_entry_mean": float(rows.mean()),
            "nkeys_median": float(np.median(nkeys)),
            "nkeys_p90": float(np.percentile(nkeys, 90)),
            "nkeys_p99": float(np.percentile(nkeys, 99)),
            "nkeys_max": int(nkeys.max()),
            "nkeys_hist": hist(nkeys),
            "nasm_mean": float(np.mean(nasm)),
            "multi_model_entries": multi_model / len(ents),
            "alt_entries": with_alt / len(ents),
            "multi_asm_entries": float(np.mean(np.array(nasm) > 1)),
            "rows_model": rows_nm / tot,
            "rows_alt": rows_na / tot,
            "rows_assembly": rows_nA / tot,
            "combos_hist": hist(np.array(combos)),
            "combos_mean": float(np.mean(combos)),
            "models_per_asm_max": int(max(models_per_asm)),
            "alts_per_model_max": int(max(alts_per_model)),
            "edges": len(edges),
            "rows_per_edge": float(er.mean()),
            "entries_per_edge": float(n_entries_per_edge.mean()),
            "mass_on_dup": frac,
        }
    return out


def top_entries(
    pool: Pool, cls: str, n: int = 10
) -> list[tuple[str, int, int, int, int, int]]:
    rows = []
    for entry, keys in pool.by_entry[cls].items():
        parsed = [split_asm(k) for k in keys]
        rows.append(
            (
                entry,
                sum(keys.values()),
                len(keys),
                len({p[0] for p in parsed}),
                len({p[1] for p in parsed}),
                len({p[2] for p in parsed}),
            )
        )
    return sorted(rows, key=lambda r: -r[1])[:n]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--items", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    # whole-DB baseline from the chain tables
    ent_keys: dict[str, set[str]] = defaultdict(set)
    for path in sorted(glob.glob(str(args.items / "chains.*.tsv"))):
        with open(path) as handle:
            for line in handle:
                p = line.split("\t", 2)
                ent_keys[p[0]].add(p[1])
    base_nkeys = np.array([len(v) for v in ent_keys.values()])
    base_mm = float(
        np.mean([any(split_asm(k)[1] != "1" for k in v) for v in ent_keys.values()])
    )
    base_alt = float(
        np.mean([any(split_asm(k)[2] != "." for k in v) for v in ent_keys.values()])
    )
    base_asm = float(
        np.mean([len({split_asm(k)[0] for k in v}) > 1 for v in ent_keys.values()])
    )

    mw, af3 = Pool(), Pool()
    for path in sorted(glob.glob(str(args.items / "interfaces.*.tsv"))):
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                entry, asm, _c1, _c2, cl1, cl2, in5 = parts[:7]
                cls = cs.mw_edge_type_from_letters({letter(cl1), letter(cl2)})
                mw.add(cls, entry, asm, cl1, cl2)
                if in5 == "1":
                    af3.add(cls, entry, asm, cl1, cl2)

    classes = sorted(
        mw.by_entry,
        key=lambda c: -sum(sum(v.values()) for v in mw.by_entry[c].values()),
    )
    res_mw, res_af = summarise(mw, classes), summarise(af3, classes)

    out = [
        "# Assembly x model x alt, by interface class\n",
        "Row = one chain-pair contact in one (assembly, model, alt) cifmol of an entry. `1` is the first assembly / model, "
        "`.` is 'no alternate location'.\n",
        "## Whole DB (every chain table row)\n",
        f"{len(ent_keys):,} entries, {int(base_nkeys.sum()):,} (assembly, model, alt) keys, mean {base_nkeys.mean():.2f} per entry, "
        f"max {int(base_nkeys.max())}. Entries with a model != 1: {base_mm:.1%}; with an alt != '.': {base_alt:.1%}; "
        f"with more than one assembly: {base_asm:.1%}.\n",
        "Keys per entry: "
        + ", ".join(
            f"{bname(lo, hi)}: {v:.1%}"
            for (lo, hi), v in zip(BINS, hist(base_nkeys), strict=True)
        )
        + "\n",
    ]

    def tbl(title: str, res: dict[str, dict], cols) -> None:
        out.append(f"## {title}\n")
        out.append(
            "| class | "
            + " | ".join(c[0] for c in cols)
            + " |\n|---|"
            + "---:|" * len(cols)
        )
        for cls in classes:
            if cls in res:
                out.append(
                    f"| {cls} | " + " | ".join(c[1](res[cls]) for c in cols) + " |"
                )
        out.append("")

    f = lambda v: f"{v:.1%}"
    tbl(
        "Entries and keys per entry (MiniWorld pool, 6 A)",
        res_mw,
        [
            ("entries", lambda r: f"{r['entries']:,}"),
            ("rows", lambda r: f"{r['rows']:,}"),
            ("rows / entry", lambda r: f"{r['rows_per_entry_mean']:.1f}"),
            (
                "keys / entry: median / p90 / p99 / max",
                lambda r: f"{r['nkeys_median']:.0f} / {r['nkeys_p90']:.0f} / {r['nkeys_p99']:.0f} / {r['nkeys_max']}",
            ),
            ("assemblies / entry", lambda r: f"{r['nasm_mean']:.2f}"),
        ],
    )
    tbl(
        "Share of entries per key-count bin (MiniWorld pool)",
        res_mw,
        [
            (bname(lo, hi), (lambda i: lambda r: f(r["nkeys_hist"][i]))(i))  # noqa: PLC3002 (binds i per column)
            for i, (lo, hi) in enumerate(BINS)
        ],
    )
    tbl(
        "Which entries, which rows (MiniWorld pool)",
        res_mw,
        [
            ("entries with a model != 1", lambda r: f(r["multi_model_entries"])),
            ("entries with an alt != '.'", lambda r: f(r["alt_entries"])),
            ("entries with > 1 assembly", lambda r: f(r["multi_asm_entries"])),
            ("rows with model != 1", lambda r: f(r["rows_model"])),
            ("rows with alt != '.'", lambda r: f(r["rows_alt"])),
            ("rows with assembly != 1", lambda r: f(r["rows_assembly"])),
        ],
    )
    tbl(
        "(model, alt) combinations per (entry, assembly) (MiniWorld pool)",
        res_mw,
        [("mean", lambda r: f"{r['combos_mean']:.2f}")]
        + [
            (bname(lo, hi), (lambda i: lambda r: f(r["combos_hist"][i]))(i))  # noqa: PLC3002 (binds i per column)
            for i, (lo, hi) in enumerate(BINS)
        ]
        + [("max models in one assembly", lambda r: str(r["models_per_asm_max"]))],
    )
    for label, res in (("MiniWorld pool (6 A)", res_mw), ("AF3 pool (5 A)", res_af)):
        tbl(
            f"Sampling mass on a non-primary row ({label}): mean over the class's item clusters of the share of their rows",
            res,
            [
                ("model != 1", lambda r: f(r["mass_on_dup"]["model"])),
                ("alt != '.'", lambda r: f(r["mass_on_dup"]["alt"])),
                ("assembly != 1", lambda r: f(r["mass_on_dup"]["assembly"])),
                ("any of the three", lambda r: f(r["mass_on_dup"]["any"])),
                ("item clusters", lambda r: f"{r['edges']:,}"),
                ("rows / cluster", lambda r: f"{r['rows_per_edge']:.1f}"),
                ("entries / cluster", lambda r: f"{r['entries_per_edge']:.1f}"),
            ],
        )
    if "protein_rna" in mw.by_entry:
        out.append("## protein_rna: the entries with the most rows (MiniWorld pool)\n")
        out.append(
            "| entry | rows | keys | assemblies | models | alts |\n|---|---:|---:|---:|---:|---:|"
        )
        for e in top_entries(mw, "protein_rna"):
            out.append("| " + " | ".join(str(x) for x in e) + " |")
        out.append("")
    args.out.write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
