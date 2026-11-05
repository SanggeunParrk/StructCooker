from __future__ import annotations

from pathlib import Path

from .abag import abag_header

_KEEP_COLUMN = "keep"


def _read_pairs(pairs_path: str | Path) -> list[dict[str, str]]:
    lines = Path(pairs_path).read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _sequences(seq_id_map_path: str | Path) -> dict[str, str]:
    return dict(line.rstrip("\n").split("\t") for line in Path(seq_id_map_path).open())


def _antigen_keep(path: str | Path | None) -> set[str] | None:
    return None if path is None else {line.strip() for line in Path(path).open() if line.strip()}


def load_final_fasta(
    pairs_path: str | Path, seq_id_map_path: str | Path, chain_set: str, antigen_keep_path: str | Path | None = None,
) -> dict[str, str]:
    """The FASTA of the chains of the kept pairs (``abag/pairs/pairs_basic.tsv``, ``keep == 1``): ``chain_set`` is ``antigen`` or ``antibody``.

    ``antigen_keep_path`` (one seq_id a line): only pairs of these antigens (the antigens with an MSA deep enough or a template). The header is the one
    ``abag/antigen_fasta`` / ``abag/antibody_fasta`` give, so a chain that was in them has the same header; a recovered light chain (``light_sequence``) gets the
    same kind. Order: (role, length, sequence).
    """
    sequence = _sequences(seq_id_map_path)
    keep = _antigen_keep(antigen_keep_path)
    per_role: dict[str, set[str]] = {"T": set(), "H": set(), "L": set()}
    for row in _read_pairs(pairs_path):
        if row[_KEEP_COLUMN] != "1" or (keep is not None and row["antigen_seq_id"] not in keep):
            continue
        per_role["T"].add(sequence[row["antigen_seq_id"]])
        per_role["H"].add(sequence[row["heavy_seq_id"]])
        light = sequence[row["light_seq_id"]] if row["light_seq_id"] else row["light_sequence"]
        if light:
            per_role["L"].add(light)
    roles = ("T",) if chain_set == "antigen" else ("H", "L")
    fasta: dict[str, str] = {}
    for tag in roles:
        for chain in sorted(per_role[tag], key=lambda s: (len(s), s)):
            fasta.setdefault(abag_header(tag, chain), chain)
    return fasta


def load_final_pairs(
    pairs_path: str | Path,
    seq_id_map_path: str | Path,
    antigen_cluster_path: str | Path,
    antibody_cluster_paths: dict[str, str | Path],
    antigen_keep_path: str | Path | None = None,
) -> list[tuple[str, ...]]:
    """The kept pairs with every chain a seq_id and every chain a cluster: the final pair table.

    Columns: ``antigen_seq_id, heavy_seq_id, light_seq_id, chain_status, light_source, datasets, antigen_cluster40`` and, for each antibody clustering
    name in ``antibody_cluster_paths`` (``cdr100``, ``cdr90``), ``heavy_cluster_<name>, light_cluster_<name>``. A recovered light chain gets the seq_id the final
    seq_id map gave it. A chain in two clusters (a few antibody ids ANARCI could not number) takes the first.
    The first row is the header.
    """
    id_of = {text: seq_id for seq_id, text in _sequences(seq_id_map_path).items()}

    def clusters(path: str | Path) -> dict[str, str]:
        member: dict[str, str] = {}
        for line in Path(path).open():
            cluster, members = line.rstrip("\n").split("\t")
            for seq_id in members.split(","):
                member.setdefault(seq_id, cluster)
        return member

    antigen_cluster = clusters(antigen_cluster_path)
    antibody_clusters = {name: clusters(path) for name, path in antibody_cluster_paths.items()}
    header = ["antigen_seq_id", "heavy_seq_id", "light_seq_id", "chain_status", "light_source", "datasets", "antigen_cluster40"]
    for name in antibody_clusters:
        header += [f"heavy_cluster_{name}", f"light_cluster_{name}"]
    out: list[tuple[str, ...]] = [tuple(header)]
    keep = _antigen_keep(antigen_keep_path)
    for row in _read_pairs(pairs_path):
        if row[_KEEP_COLUMN] != "1" or (keep is not None and row["antigen_seq_id"] not in keep):
            continue
        light_id = row["light_seq_id"] or (id_of[row["light_sequence"]] if row["light_sequence"] else "")
        cells = [row["antigen_seq_id"], row["heavy_seq_id"], light_id, row["chain_status"], row["light_source"], row["datasets"],
                 antigen_cluster[row["antigen_seq_id"]]]
        for member in antibody_clusters.values():
            cells += [member[row["heavy_seq_id"]], member[light_id] if light_id else ""]
        out.append(tuple(cells))
    return out


def load_antigen_worklist(pairs_path: str | Path, seq_id_map_path: str | Path) -> dict[str, str]:
    """``{seq_id: sequence}`` of the antigens of the kept pairs (``keep == 1``): what the MSA and template searches run on."""
    sequence = _sequences(seq_id_map_path)
    antigens = {row["antigen_seq_id"] for row in _read_pairs(pairs_path) if row[_KEEP_COLUMN] == "1"}
    return {seq_id: sequence[seq_id] for seq_id in sorted(antigens)}


def load_seq_id_maps(seq_id_map_paths: list[str | Path]) -> dict[str, str]:
    """The union of seq_id maps (``seq_id<TAB>sequence``), keyed the way ``build_seq_id_map`` reads its seed: ``<type letter><sequence>`` -> seq_id.

    The first map to hold a sequence decides its id; a later map that gives the same sequence another id is an error (one id space means one id per sequence).
    A missing map is skipped, as ``load_seq_id_map`` does.
    """
    merged: dict[str, str] = {}
    for path in seq_id_map_paths:
        if not Path(path).exists():
            continue
        for line in Path(path).open():
            seq_id, text = line.rstrip("\n").split("\t")
            key = seq_id[0] + text
            if merged.setdefault(key, seq_id) != seq_id:
                msg = f"one sequence has two ids: {merged[key]} and {seq_id} ({path})"
                raise ValueError(msg)
    return merged
