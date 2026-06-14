from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

MIN_SEQUENCE_LENGTH = 16


def copy_tree(src: Path, dst: Path) -> None:
    """Copy one MSA directory into its remapped destination."""
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns("*.hhr", "*.atab"),
    )


def load_seq_id_map(tsv_path: Path) -> dict[str, str]:
    """Load a query-sequence to seq-id mapping from a TSV file."""
    seq_id_map: dict[str, str] = {}
    with tsv_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            seq_id, sequence = line.split("	")
            if not seq_id.startswith("P"):
                continue
            if set(sequence) == {"X"} or len(sequence) < MIN_SEQUENCE_LENGTH:
                continue
            existing = seq_id_map.get(sequence)
            if existing is not None and existing != seq_id:
                msg = (
                    f"Duplicate sequence found for {sequence!r}: "
                    f"{existing} vs {seq_id}"
                )
                raise ValueError(msg)
            seq_id_map[sequence] = seq_id
    return seq_id_map


def load_query_seq_from_a3m(a3m_path: Path) -> str:
    """Return the first query sequence from an A3M file."""
    with a3m_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(">"):
                continue
            return line
    return ""


def remap_msa_tree(seq_id_map_path: Path, old_msa_path: Path, new_msa_path: Path) -> int:
    """Copy MSA directories into a seq-id keyed destination tree."""
    seq_id_map = load_seq_id_map(seq_id_map_path)
    copied = 0
    for subdir in sorted(path for path in old_msa_path.iterdir() if path.is_dir()):
        old_msa_file = subdir / "t000_msa0.a3m"
        if not old_msa_file.exists():
            continue
        query_seq = load_query_seq_from_a3m(old_msa_file)
        seq_id = seq_id_map.get(query_seq)
        if seq_id is None:
            continue
        dst = new_msa_path / seq_id[:4] / seq_id
        dst.parent.mkdir(parents=True, exist_ok=True)
        copy_tree(subdir, dst)
        copied += 1
    return copied


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser for the remap utility."""
    parser = argparse.ArgumentParser(description="Remap MSA directories by sequence id.")
    parser.add_argument("seq_id_map", type=Path)
    parser.add_argument("old_msa_dir", type=Path)
    parser.add_argument("new_msa_dir", type=Path)
    return parser


def main() -> None:
    """Run the remap utility from the command line."""
    args = build_parser().parse_args()
    copied = remap_msa_tree(args.seq_id_map, args.old_msa_dir, args.new_msa_dir)
    sys.stdout.write(f"Copied {copied} MSA directories.\n")


if __name__ == "__main__":
    main()
