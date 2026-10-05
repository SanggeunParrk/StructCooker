"""Repair ``_entity_poly.pdbx_seq_one_letter_code`` in the staged AFDB heterodimer mmCIFs.

The one-letter code in ``_entity_poly`` is a summary of the residues; in 6 of the 80,248
release files it disagrees with the rest of the file (measured 2026-09-26, every file checked):

* 4 files list ``_entity_poly`` rows as entity 2, 1 but the sequences as entity 1, 2 -- the
  id and the sequence are crossed inside the file. ``_entity_poly_seq``, ``_struct_asym`` and
  ``_atom_site`` agree with each other and with the geometry, so only this summary is wrong.
  (The build caught these because the two lengths differ; a same-length swap would pass
  silently, which is why every file was checked rather than only the failures.)
* 2 files carry UniProt's ambiguity code ``Z`` (Glu/Gln; P01681 is an Edman-era entry) at
  positions the model built as GLU.

For a predicted model every residue is modelled, so ``_entity_poly_seq`` is the record of what
the structure contains; the one-letter code is rebuilt from it. Only the sequence token of
each ``_entity_poly`` row is replaced (that row is re-joined with single spaces, which mmCIF
reads the same); every other line is byte-identical.

The raw release files are not touched. Repaired copies are written to
``BioMol/materials/intermediate/afdb_multimer/heterodimer/repaired/cif/<last 3>/`` with
``ENTITY_POLY_REPAIRS.tsv`` (every change) and ``overrides.json`` ({raw path: repaired path}),
which db/afdb_multimer/heterodimer_cif.yaml passes to the reader as AFM_CIF_OVERRIDES.

Usage:
    python scripts/staging/repair_afdb_entity_poly.py [--list bad.tsv] [--dry-run]
      Without --list every heterodimer mmCIF under raw/ is checked (80,248 files, ~minutes);
      a file is copied only when its _entity_poly disagrees with _entity_poly_seq.
      bad.tsv: lines "BAD<TAB><file name><TAB>..." limiting the check to those files.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
from pathlib import Path

_DATA = Path(os.environ.get("DATA_ROOT", "/data/shared/cssb_data"))
RAW = _DATA / "BioMol/materials/raw/afdb_multimer/heterodimer"
OUT = _DATA / "BioMol/materials/intermediate/afdb_multimer/heterodimer/repaired"
_AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _loop(lines: list[str], cat: str) -> tuple[list[str], int, int]:
    """Return (columns, first data line, end line) of a ``loop_`` category."""
    i = 0
    while i < len(lines) and not lines[i].startswith(cat + "."):
        i += 1
    cols = []
    while i < len(lines) and lines[i].startswith(cat + "."):
        cols.append(lines[i].split(".", 1)[1].strip())
        i += 1
    start = i
    while i < len(lines) and not lines[i].startswith(("#", "loop_", "_")):
        i += 1
    return cols, start, i


def repair(path: Path, dest: Path) -> list[tuple[str, str, str]]:
    """Write a repaired copy of ``path`` to ``dest``; return (entity, before, after) per changed row."""
    with gzip.open(path, "rt") as fh:
        text = fh.read()
    lines = text.split("\n")

    cols, s, e = _loop(lines, "_entity_poly_seq")
    # mmCIF rows are token streams, not lines: a writer may pack several rows on one line.
    tokens = [t for line in lines[s:e] for t in shlex.split(line)]
    if len(tokens) % len(cols):
        msg = f"{path.name}: _entity_poly_seq token count is not a multiple of its columns"
        raise ValueError(msg)
    seq: dict[str, list[str]] = {}
    for j in range(0, len(tokens), len(cols)):
        row = dict(zip(cols, tokens[j : j + len(cols)], strict=True))
        seq.setdefault(row["entity_id"], []).append(_AA.get(row["mon_id"], "X"))

    cols, s, e = _loop(lines, "_entity_poly")
    at_id, at_code = cols.index("entity_id"), cols.index("pdbx_seq_one_letter_code")
    changes = []
    for k in range(s, e):
        tokens = lines[k].split()
        if not tokens:   # release files end some loops with a blank line rather than "#"
            continue
        if len(tokens) != len(cols):   # a row split over lines would need a real parser
            msg = f"{path.name}: _entity_poly row is not single-line; refusing to rewrite"
            raise ValueError(msg)
        ent, before = tokens[at_id], tokens[at_code]
        after = "".join(seq[ent])
        if before != after:
            tokens[at_code] = after
            lines[k] = " ".join(tokens)
            changes.append((ent, before, after))
    if changes:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(dest, "wt") as out:
            out.write("\n".join(lines))
    return changes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.list is not None:
        names = sorted({ln.split("\t")[1] for ln in a.list.read_text().splitlines() if ln.startswith("BAD\t")})
    else:
        names = sorted(p.name for p in (RAW / "cif").glob("*/AF-*-model_v1.cif.gz"))
    log = OUT / "ENTITY_POLY_REPAIRS.tsv"
    rows, overrides = [], {}
    for name in names:
        n = name.split("-")[1]
        path = RAW / "cif" / n[-3:] / name
        if a.dry_run:
            print(f"would repair {path}")
            continue
        dest = OUT / "cif" / n[-3:] / name
        changes = repair(path, dest)
        if changes:
            overrides[str(path)] = str(dest)
        for ent, before, after in changes:
            diff = [i for i, (x, y) in enumerate(zip(before, after, strict=False)) if x != y]
            kind = "ambiguity_code" if len(before) == len(after) and {before[i] for i in diff} <= set("ZBXJU") \
                else "crossed_entity"
            rows.append(f"{name}\t{ent}\t{kind}\t{len(before)}\t{len(after)}\t{before}\t{after}\n")
            print(f"{name} entity {ent}: {kind} ({len(before)} -> {len(after)})")
    if rows:
        (OUT / "overrides.json").write_text(json.dumps(overrides, indent=2) + "\n")
        with log.open("w") as fh:
            fh.write("file\tentity\tkind\tlen_before\tlen_after\tbefore\tafter\n")
            fh.writelines(rows)


if __name__ == "__main__":
    main()
