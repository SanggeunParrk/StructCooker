from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from .abag import normalize as _normalize

# an antigen that is a polyprotein: ASD's own target name says so
POLYPROTEIN = re.compile(r"polyprotein|orf1a|pp1ab?\b|gag[-/ ]?pol", re.IGNORECASE)
HEADER = (
    "antigen_seq_id",
    "heavy_seq_id",
    "light_seq_id",
    "light_sequence",
    "chain_status",
    "light_source",
    "datasets",
    "polyprotein",
    "pdb_structure",
    "keep",
)
_KEEP_STATUS = (
    "paired_original",
    "scfv",
    "vhh",
    "filled_covid",
    "filled_hiv",
    "filled_biomap",
)
_K = 25


def _contains(a: str, b: str, min_length: int) -> bool:
    short, long = sorted((a, b), key=len)
    return len(short) >= min_length and short in long


class PdbComplexIndex:
    """Antibody-antigen complexes with an experimental structure: the SAbDab entries (heavy chain, light chain, antigen chains), their sequences taken
    from the PDB fasta (``polypeptide(L)`` records, matched on PDB id and author chain id).

    ``find(heavy, light, antigen)`` returns the PDB id of a complex that is **this antibody with this antigen**, else ``None``:

    * ``match="exact"`` (the rule): the heavy chain, the light chain (none for a single-chain pair) and the antigen are each **identical** to the chains of one
      SAbDab entry -- the same seq_id for every chain. A point mutant of a crystallised antibody is not found, nor is the same antibody in another construct
      (a variable domain against a Fab chain, a full UniProt chain against a fragment);
    * ``match="contained"``: the looser reading. The antibody chains are contained in the PDB chains or contain them (a PDB chain carries a constant domain or a
      tag, ASD often the variable domain), the antigen matches when the PDB antigen chain is contained in the ASD antigen, the reverse, or 80% of its 10-mers
      (sampled every 5th) occur in it. It flags about twice as many pairs (4,616 against 2,204 on the ASD pairs).
    """

    def __init__(
        self, sabdab_summary_path: str | Path, pdb_fasta_path: str | Path
    ) -> None:
        entries = []
        wanted: set[tuple[str, str]] = set()
        for line in list(Path(sabdab_summary_path).open())[1:]:
            fields = line.rstrip("\n").split("\t")
            pdb, heavy_chain, light_chain, antigen_chain, antigen_type = (
                fields[0].upper(),
                fields[1],
                fields[2],
                fields[4],
                fields[5],
            )
            if heavy_chain in ("NA", "") or (
                "protein" not in antigen_type and "peptide" not in antigen_type
            ):
                continue
            antigens = [
                c.strip()
                for c in antigen_chain.split("|")
                if c.strip() not in ("", "NA")
            ]
            if not antigens:
                continue
            entries.append(
                (
                    pdb,
                    heavy_chain,
                    None if light_chain in ("NA", "") else light_chain,
                    antigens,
                )
            )
            wanted.update((pdb, c) for c in [heavy_chain, light_chain, *antigens])
        chain: dict[tuple[str, str], str] = {}
        key = None
        for line in Path(pdb_fasta_path).open():
            if line.startswith(">"):
                parts = line[1:].split(" | ")
                key = None
                if len(parts) > 2 and "polypeptide(L)" in parts[1]:
                    candidate = (parts[0].split("_")[0], parts[2].strip()[5:])
                    key = candidate if candidate in wanted else None
            elif key:
                chain.setdefault(key, _normalize(line))
        self.entries = []
        for pdb, heavy_chain, light_chain, antigens in entries:
            heavy = chain.get((pdb, heavy_chain))
            light = chain.get((pdb, light_chain)) if light_chain else ""
            antigen_sequences = [chain[(pdb, c)] for c in antigens if (pdb, c) in chain]
            if heavy and len(heavy) >= 80 and antigen_sequences and light is not None:
                self.entries.append((pdb, heavy, light, antigen_sequences))
        self.index: dict[str, set[int]] = defaultdict(set)
        for number, (_, heavy, _, _) in enumerate(self.entries):
            for i in range(len(heavy) - _K + 1):
                self.index[heavy[i : i + _K]].add(number)
        self._kmers: dict[str, set[str]] = {}
        self.exact: dict[tuple[str, str], list[tuple[str, list[str]]]] = defaultdict(
            list
        )
        for pdb, heavy, light, antigen_sequences in self.entries:
            self.exact[(heavy, light)].append((pdb, antigen_sequences))

    def _antigen_matches(self, antigen: str, pdb_antigen: str) -> bool:
        if len(pdb_antigen) < 8:
            return False
        if pdb_antigen in antigen or (len(antigen) >= 8 and antigen in pdb_antigen):
            return True
        if len(pdb_antigen) < 30:
            return False
        kmers = self._kmers.get(antigen)
        if kmers is None:
            kmers = self._kmers[antigen] = {
                antigen[i : i + 10] for i in range(len(antigen) - 9)
            }
        sample = [pdb_antigen[i : i + 10] for i in range(0, len(pdb_antigen) - 9, 5)]
        return sum(k in kmers for k in sample) >= 0.8 * len(sample)

    def find(
        self, heavy: str, light: str, antigen: str, match: str = "exact"
    ) -> str | None:
        if match == "exact":
            for pdb, antigens in self.exact.get((heavy, light), ()):
                if antigen in antigens:
                    return pdb
            return None
        if len(heavy) < 80:
            return None
        candidates: set[int] = set()
        for position in (len(heavy) // 3, len(heavy) // 2, 2 * len(heavy) // 3):
            candidates |= self.index.get(heavy[position : position + _K], set())
        for number in sorted(candidates):
            pdb, entry_heavy, entry_light, antigens = self.entries[number]
            if not _contains(heavy, entry_heavy, 80):
                continue
            if light:
                if not entry_light or not _contains(light, entry_light, 60):
                    continue
            elif entry_light:
                continue  # the pair is a single chain (VHH / scFv-like) and the structure has two chains
            if any(self._antigen_matches(antigen, a) for a in antigens):
                return pdb
        return None


def load_basic_pairs(
    pairs_chains_path: str | Path,
    seq_id_map_path: str | Path,
    sabdab_summary_path: str | Path,
    pdb_fasta_path: str | Path,
    structure_match: str = "exact",
) -> list[tuple[str, ...]]:
    """The pair table with the chain status (``abag/pairs/pairs_chains.tsv``) after the basic filters, each one a column so nothing is lost:

    * duplicates: a pair is one (heavy, light, antigen) sequence triple; after the light chains are filled two rows can be the same triple (a filled light
      chain that a source also holds): merged, their datasets joined;
    * ``polyprotein``: 1 when a target name of the antigen says polyprotein / ORF1a(b) / gag-pol (antigens over 3,000 aa are not in the table at all);
    * ``pdb_structure``: the PDB id of a complex whose heavy chain, light chain and antigen are all identical to this pair's (``PdbComplexIndex``,
      ``structure_match="exact"``), empty if none -- a variant of a crystallised antibody is not one;
    * ``keep``: 1 when the chain status is usable (paired / scfv / vhh / filled from AntiBinder), not a polyprotein, and no structure.

    The first row is the header (``HEADER``).
    """
    sequence, id_of = {}, {}
    for line in Path(seq_id_map_path).open():
        seq_id, text = line.rstrip("\n").split("\t")
        sequence[seq_id] = text
        id_of[text] = seq_id
    pairs = [
        line.rstrip("\n").split("\t")
        for line in list(Path(pairs_chains_path).open())[1:]
    ]
    merged: dict[tuple[str, str, str], list] = {}
    for antigen, heavy, light_id, light_seq, status, source, datasets, names in pairs:
        light = sequence[light_id] if light_id else light_seq
        if (
            light_seq and light in id_of
        ):  # a recovered light chain a source also holds: it has a seq_id after all
            light_id, light_seq = id_of[light], ""
        key = (heavy, light_id or light_seq, antigen)
        if key in merged:
            merged[key][6] = ",".join(
                sorted(set(merged[key][6].split(",")) | set(datasets.split(",")))
            )
            merged[key][7] = ";".join(
                sorted(
                    set(filter(None, merged[key][7].split(";")))
                    | set(filter(None, names.split(";")))
                )
            )
        else:
            merged[key] = [
                antigen,
                heavy,
                light_id,
                light_seq,
                status,
                source,
                datasets,
                names,
            ]
    index = PdbComplexIndex(sabdab_summary_path, pdb_fasta_path)
    out: list[tuple[str, ...]] = [HEADER]
    for (
        antigen,
        heavy,
        light_id,
        light_seq,
        status,
        source,
        datasets,
        names,
    ) in merged.values():
        light = sequence[light_id] if light_id else light_seq
        polyprotein = int(bool(POLYPROTEIN.search(names)))
        pdb = (
            index.find(sequence[heavy], light, sequence[antigen], structure_match) or ""
        )
        keep = int(status in _KEEP_STATUS and not polyprotein and not pdb)
        out.append(
            (
                antigen,
                heavy,
                light_id,
                light_seq,
                status,
                source,
                datasets,
                str(polyprotein),
                pdb,
                str(keep),
            )
        )
    return out
