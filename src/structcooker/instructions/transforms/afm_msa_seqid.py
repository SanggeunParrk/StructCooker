"""Build the AFM MSA entity -> seq_id map (feeds ``readers.afdb.afm_msa_seqid_key``).

The AFDB release ships one MSA per monomer entity. A homodimer uses its own entity's MSA;
a heterodimer chain uses its monomer's, via ``heterodimer/chain_msa.tsv`` (content-verified
at staging). Each MSA entity is resolved to a seq_id through the sequence of a chain that
uses it, taken from the AFM FASTA -- so the MSA DB is keyed like every other MSA DB.
"""

from __future__ import annotations

from pathlib import Path


def _chains(fasta_path: Path) -> dict[tuple[str, str], str]:
    """``(entity, chain) -> sequence`` from FASTA headers ``AF-<n>-model_v1_<chain>_<alt>``."""
    out: dict[tuple[str, str], str] = {}
    header: str | None = None
    with fasta_path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split(" | ", 1)[0]
            elif header is not None:
                entry, chain, _alt = header.rsplit("_", 2)
                out[(entry.split("-model_v1", 1)[0], chain)] = line
                header = None
    return out


def build_afm_msa_seqid_map(
    fasta_path: str | Path,
    seq_id_map_path: str | Path,
    chain_msa_path: str | Path,
) -> dict[str, str]:
    """Return ``{msa_entity: seq_id}``; raise if two chains disagree on an entity's sequence."""
    chains = _chains(Path(fasta_path))
    seq_to_id: dict[str, str] = {}
    with Path(seq_id_map_path).open(encoding="utf-8") as handle:
        for line in handle:
            seq_id, _, seq = line.rstrip("\n").partition("\t")
            if seq:
                seq_to_id[seq] = seq_id

    uses: dict[str, str] = {}  # msa entity -> sequence
    def claim(entity: str, seq: str) -> None:
        prior = uses.setdefault(entity, seq)
        if prior != seq:
            msg = f"MSA {entity} is used by chains with different sequences"
            raise ValueError(msg)

    homo = {e for (e, _c) in chains}
    with Path(chain_msa_path).open(encoding="utf-8") as handle:
        next(handle)  # header
        hetero = [line.rstrip("\n").split("\t") for line in handle]
    het_ids = {row[0] for row in hetero}
    for entity in homo - het_ids:          # a homodimer: its own entity, chain A
        claim(entity, chains[(entity, "A")])
    for het, chain, _acc, monomer, _how, verified in hetero:
        if verified != "ok":
            continue
        claim(monomer, chains[(het, chain)])

    missing = [e for e, s in uses.items() if s not in seq_to_id]
    if missing:
        msg = f"{len(missing)} MSA sequences have no seq_id (first: {missing[:3]})"
        raise KeyError(msg)
    return {e: seq_to_id[s] for e, s in uses.items()}
