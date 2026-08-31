"""Build the distillation monomer entry -> seq_id map.

The MSA depends only on the sequence, so the distillation MSA DBs are keyed by seq_id.
This joins each monomer distillation FASTA (entry -> sequence) with the seq_id_map
(sequence -> seq_id) to produce the entry -> seq_id lookup the ``openfold_seqid_key``
key builder loads at build time. Streams the FASTAs so the (large) long set stays cheap.
"""
from __future__ import annotations

from pathlib import Path


def _entry_of(header: str) -> str:
    """Recover the openfold entry key (the structure folder name) from a FASTA header.

    Headers read ``<entry>_<chain>_<ins> | <moltype> | Auth:<n>`` where ``<entry>`` may
    itself contain underscores (rna accessions), so drop the ``| ...`` tail and strip the
    trailing two ``_<chain>_<ins>`` fields.
    """
    return header.split(" | ", 1)[0].rsplit("_", 2)[0]


def build_monomer_seqid_map(
    fasta_paths: list[str],
    seq_id_map_path: str | Path,
) -> dict[str, dict[str, str]]:
    """Return ``{"monomer_seqid": {entry: seq_id}}`` for the given monomer FASTAs."""
    seq_to_id: dict[str, str] = {}
    with Path(seq_id_map_path).open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                seq_to_id[parts[1]] = parts[0]

    entry_seqid: dict[str, str] = {}
    for fasta_path in fasta_paths:
        with Path(fasta_path).open(encoding="utf-8") as handle:
            header: str | None = None
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    header = line[1:]
                elif header is not None:
                    seq_id = seq_to_id.get(line)
                    if seq_id is not None:
                        entry_seqid[_entry_of(header)] = seq_id
                    header = None
    return {"monomer_seqid": entry_seqid}
