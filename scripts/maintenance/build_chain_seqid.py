"""Build a comprehensive `<pdbid_lower>_<chain>` -> seq_id LMDB for ALL molecule
types (protein / DNA / RNA / hybrid / branched / non-polymer).

The cif chains already carry seq_id, but decoding every heavy CIFMol record to
read it is slow. That seq_id is assigned (see ``build_seq_metadata_map`` /
``attach_metadata``) purely from two flat text artifacts:

  * cif.fasta        -- per-chain header ``pdbid_chain_alt | mol_type | Auth``
                        plus the chain's sequence *representation* (linear seq for
                        polymers, ``(RES)(RES)|(bonds)`` for branched, comp id for
                        non-polymer).
  * seq_id_map.tsv   -- ``seq_id <tab> sequence`` for every unique (mol, sequence);
                        the lookup key is ``seq_id[0] + sequence`` exactly as in
                        ``readers.sequence.load_seq_id_map``.

So chain -> seq_id is a pure text join keyed by ``mol_prefix + sequence`` -- no
CIF decode. We group the fasta by that key (small) and stream seq_id_map.tsv
(large) once, so peak memory stays near the fasta size.

    pixi run python -m scripts.maintenance.build_chain_seqid \
        --fasta        /data/shared/cssb_data/BioMol/materials/raw/fasta/cif_pdb.fasta \
        --seq-id-map   /data/shared/cssb_data/BioMol/metadata/seq_id_map.tsv \
        --out          /data/shared/cssb_data/BioMol/metadata/chain_to_seqid.lmdb
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import lmdb

from structcooker.utils.mapping import mol_type_map


def _chain_key(cif_id: str) -> str | None:
    # "5KTJ_A_."  ->  "5ktj_A"
    parts = cif_id.split("_")
    if len(parts) < 2:  # noqa: PLR2004 - need pdbid + chain
        return None
    return f"{parts[0].lower()}_{parts[1]}"


def _group_fasta(fasta_path: Path) -> dict[str, list[str]]:
    """mol_prefix+sequence -> [chain_key, ...] over every fasta chain."""
    groups: dict[str, list[str]] = defaultdict(list)
    header = None
    prefix = None
    seq_parts: list[str] = []
    n = kept = 0

    def flush() -> None:
        nonlocal kept
        if header is not None and prefix is not None and seq_parts:
            ck = _chain_key(header.split("|")[0].strip())
            if ck is not None:
                groups[prefix + "".join(seq_parts)].append(ck)
                kept += 1

    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                mol_type = header.split("|")[1].strip() if "|" in header else "unknown"
                prefix = mol_type_map.get(mol_type, "X")
                seq_parts = []
                n += 1
            else:
                seq_parts.append(line.strip())
        flush()
    print(f"fasta: {n:,} chains, {kept:,} keyed, {len(groups):,} unique mol+seq", flush=True)
    return groups


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--seq-id-map", required=True, type=Path)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    groups = _group_fasta(args.fasta)

    env = lmdb.open(args.out, map_size=8 * 1024**3, subdir=True)
    txn = env.begin(write=True)
    seen_keys = 0
    written = 0
    matched_seqs = 0
    with args.seq_id_map.open() as fh:
        for raw in fh:
            row = raw.split()
            if len(row) != 2:  # noqa: PLR2004 - seq_id + sequence
                continue
            seq_id, sequence = row
            chains = groups.get(seq_id[0] + sequence)
            if chains is None:
                continue
            matched_seqs += 1
            for ck in chains:
                txn.put(ck.encode(), seq_id.encode(), overwrite=True)
                written += 1
            seen_keys += 1
            if seen_keys % 200_000 == 0:
                txn.commit()
                txn = env.begin(write=True)
                print(f"  matched {matched_seqs:,} seqs, {written:,} chains", flush=True)
    txn.commit()

    total_chains = sum(len(v) for v in groups.values())
    with env.begin() as t:
        rows = t.stat()["entries"]
    env.close()
    print(
        f"done: {written:,}/{total_chains:,} fasta chains mapped "
        f"({matched_seqs:,} unique seqs) -> {rows:,} lmdb rows -> {args.out}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
