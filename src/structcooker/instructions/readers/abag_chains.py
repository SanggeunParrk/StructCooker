from __future__ import annotations

import csv
import re
from pathlib import Path

_LINKER = re.compile(r"(G{3,5}S){2,}")
_MIN_LIGHT = 60  # a V-domain light chain is ~107 aa; a shorter rest is not one (SelfDistillation: filled_hiv)


def _normalize(sequence: str | None) -> str:
    return "".join((sequence or "").split()).upper()


class AntiBinderIndex:
    """The AntiBinder raw tables (Brief Bioinform 26(1) bbaf008), the original source of ASD's covid-19, hiv and biomap rows, as lookups
    from (heavy chain, antigen) to the light chains the source gives.

    * covid: ``CoV-AbDab_only_sars2_filter.csv``, columns ``Heavy``, ``Light``, ``antigen``;
    * hiv: ``dataset_hiv_cls.csv``, ``antibody`` is VH + VL as one string (no linker) and ``virus_seq`` the antigen;
    * biomap: ``biomap.csv``, ``antibody_seq_a`` / ``antibody_seq_b`` (one is the heavy chain, the other the light) and ``antigen_seq``.
    """

    def __init__(self, raw_data_dir: str | Path) -> None:
        csv.field_size_limit(10**9)
        root = Path(raw_data_dir)
        self.covid: dict[tuple[str, str], set[str]] = {}
        for row in csv.DictReader((root / "COVID-19/CoV-AbDab_only_sars2_filter.csv").open(newline="", encoding="utf-8-sig")):
            heavy, light, antigen = _normalize(row["Heavy"]), _normalize(row["Light"]), _normalize(row["antigen"])
            if heavy and light and antigen:
                self.covid.setdefault((heavy, antigen), set()).add(light)
        self.hiv: dict[str, set[str]] = {}
        for row in csv.DictReader((root / "HIV/dataset_hiv_cls.csv").open(newline="", encoding="utf-8-sig")):
            antibody, antigen = _normalize(row["antibody"]), _normalize(row["virus_seq"])
            if antibody and antigen:
                self.hiv.setdefault(antigen, set()).add(antibody)
        self.biomap: dict[tuple[str, str], set[str]] = {}
        for row in csv.DictReader((root / "BioMap/biomap.csv").open(newline="", encoding="utf-8-sig")):
            first, second, antigen = _normalize(row["antibody_seq_a"]), _normalize(row["antibody_seq_b"]), _normalize(row["antigen_seq"])
            if first and second and antigen:
                self.biomap.setdefault((first, antigen), set()).add(second)
                self.biomap.setdefault((second, antigen), set()).add(first)

    def candidates(self, heavy: str, antigen: str) -> dict[str, set[str]]:
        """The light chains each source gives for this heavy chain and antigen: ``{source: {light, ...}}``, empty sources left out."""
        found = {
            "covid": set(self.covid.get((heavy, antigen), ())),
            "hiv": {antibody[len(heavy):] for antibody in self.hiv.get(antigen, ())
                    if antibody.startswith(heavy) and len(antibody) - len(heavy) >= _MIN_LIGHT},
            "biomap": set(self.biomap.get((heavy, antigen), ())),
        }
        return {source: lights for source, lights in found.items() if lights}


def classify_chains(heavy: str, light: str, *, nanobody: bool, scfv: bool, antibinder: AntiBinderIndex, antigen: str) -> tuple[str, str, str]:
    """``(chain_status, light_source, light_sequence)`` of one pair. The rules, in order (SelfDistillation, Paired-set-build, 3):

    ``paired_original`` ASD holds a light chain; ``scfv`` the scfv flag or a ``(G{3,5}S){2,}`` linker in the heavy chain; ``vhh`` the nanobody flag;
    ``filled_covid`` / ``filled_hiv`` / ``filled_biomap`` AntiBinder gives exactly one light chain for this heavy chain and antigen;
    ``ambiguous`` two or more (not filled: a wrong light chain is worse than none); ``unresolved`` none.
    ``light_sequence`` is set only for a filled pair.
    """
    if light:
        return "paired_original", "source", ""
    if nanobody:
        return "vhh", "source", ""
    if scfv or _LINKER.search(heavy):
        return "scfv", "source", ""
    found = antibinder.candidates(heavy, antigen)
    lights = {light for group in found.values() for light in group}
    if len(lights) == 1:
        names = sorted(found)
        return f"filled_{names[0]}", "antibinder_" + "+".join(names), next(iter(lights))
    return ("ambiguous" if lights else "unresolved"), "", ""


def load_pair_chains(
    pairs_path: str | Path,
    seq_id_map_path: str | Path,
    antibinder_dir: str | Path,
) -> list[tuple[str, ...]]:
    """The pair table (``abag/pairs/pairs_binding.tsv``) with the chain status of each pair:
    ``(antigen_seq_id, heavy_seq_id, light_seq_id, light_sequence, chain_status, light_source, datasets, names)``.

    The nanobody and scfv flags are columns of the pair table. A recovered light chain is a new sequence, so it has no seq_id yet: ``light_seq_id`` is
    empty and ``light_sequence`` holds it. The first row is the header.
    """
    sequence = {}
    for line in Path(seq_id_map_path).open():
        seq_id, text = line.rstrip("\n").split("\t")
        sequence[seq_id] = text
    index = AntiBinderIndex(antibinder_dir)
    out: list[tuple[str, ...]] = [("antigen_seq_id", "heavy_seq_id", "light_seq_id", "light_sequence", "chain_status", "light_source", "datasets", "names")]
    for line in list(Path(pairs_path).open())[1:]:
        antigen, heavy, light, nanobody, scfv, datasets, names = line.rstrip("\n").split("\t")
        status, source, filled = classify_chains(
            sequence[heavy], sequence[light] if light else "", nanobody=nanobody == "1", scfv=scfv == "1", antibinder=index, antigen=sequence[antigen],
        )
        out.append((antigen, heavy, light, filled, status, source, datasets, names))
    return out
