from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .abag import iter_rows, valid
from .asd_binding import pair_binds

MIN_ANTIGEN = 10  # a peptide under 10 aa has no alignment worth searching
MAX_ANTIGEN = 3000  # longer ones are polyproteins and giant single chains (and what an MSA search costs most for)
HEADER = (
    "antigen_seq_id",
    "heavy_seq_id",
    "light_seq_id",
    "nanobody",
    "scfv",
    "datasets",
    "names",
)


def load_abag_pairs(
    asd_dir: str | Path,
    abrank_path: str | Path,
    poisk_path: str | Path,
    seq_id_map_path: str | Path,
    min_antigen: int = MIN_ANTIGEN,
    max_antigen: int = MAX_ANTIGEN,
) -> list[tuple[str, ...]]:
    """The antibody-antigen pairs of ASD, AbRank and POISK that bind with Kd <= 1 uM (the rule per source: ``asd_binding``, ``abag_binding``).

    * a row needs a heavy chain and an antigen, every chain A-Z only (S2);
    * the antigen is ``min_antigen``..``max_antigen`` aa;
    * a pair is one (heavy, light, antigen) sequence triple however many rows and sources report it (S3); it binds when some row says so and none says it does not
      (``pair_binds``), the rows of every source counted together;
    * ``datasets`` lists ``source:dataset`` of the rows that say it binds, ``names`` the target names of those rows (``;``-joined, for the polyprotein test),
      ``nanobody`` / ``scfv`` are the flags of any row.

    The first row is the header. ``seq_id_map_path`` maps every chain sequence to its seq_id (``abag/seq_id_map``).
    """
    id_of: dict[str, str] = {}
    for line in Path(seq_id_map_path).open():
        seq_id, sequence = line.rstrip("\n").split("\t")
        id_of[sequence] = seq_id
    evidence: dict[tuple[str, str, str], list] = defaultdict(list)
    for row in iter_rows(asd_dir, abrank_path, poisk_path):
        if not valid(row) or not min_antigen <= len(row.antigen) <= max_antigen:
            continue
        evidence[(row.heavy, row.light, row.antigen)].append(row)
    out: list[tuple[str, ...]] = [HEADER]
    for (heavy, light, antigen), rows in evidence.items():
        if not pair_binds([r.verdict for r in rows]):
            continue
        yes = [r for r in rows if r.verdict]
        datasets = ",".join(sorted({f"{r.source}:{r.dataset}" for r in yes}))
        names = ";".join(sorted({r.name for r in yes if r.name}))[:300]
        out.append(
            (
                id_of[antigen],
                id_of[heavy],
                id_of[light] if light else "",
                str(int(any(r.nanobody for r in rows))),
                str(int(any(r.scfv for r in rows))),
                datasets,
                names,
            )
        )
    return [out[0], *sorted(out[1:])]
