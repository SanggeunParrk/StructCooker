"""``structcooker interface ...``: submit the interface-unit DB stages to SLURM.

Each command sbatches one script of ``scripts/interface`` from the repository root, with paths under
``BIOMOL_ROOT`` (``lmdb/interface/pdb`` for the DB, ``materials/intermediate/interface`` for the
antibody cache). Logs go to ``logs/interface/``. ``--dry-run`` prints the command.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import click

from structcooker.paths import biomol_root

REPO = Path(__file__).resolve().parents[3]
SCRIPTS = REPO / "scripts" / "interface"


def default_out() -> Path:
    """Where the PDB interface DB lives in the release tree."""
    return biomol_root() / "lmdb" / "interface" / "pdb"


def _submit(script: str, args: list[str], dry_run: bool, sbatch_opts: list[str] | None = None) -> None:
    cmd = ["sbatch", *(sbatch_opts or []), str(SCRIPTS / script), *args]
    if dry_run:
        click.echo(" ".join(cmd))
        return
    (REPO / "logs" / "interface").mkdir(parents=True, exist_ok=True)
    out = subprocess.run(cmd, cwd=REPO, check=True, capture_output=True, text=True)  # noqa: S603
    click.echo(out.stdout.strip())


out_option = click.option("--out", type=click.Path(path_type=Path), default=None,
                          help="DB directory (default: $BIOMOL_ROOT/lmdb/interface/pdb).")
dry_option = click.option("--dry-run", is_flag=True, help="Print the sbatch command, do not submit.")


@click.group("interface")
def interface_group() -> None:
    """Interface-unit DB of the PDB: units of 2..N touching molecules, waters included (docs/interface-db.md)."""


@interface_group.command("annotate")
@dry_option
def annotate(dry_run: bool) -> None:
    """Classify every unique PDB protein sequence with ANARCI (antibody / TCR cache). Run once."""
    _submit("annotate_antibodies.sbatch", [], dry_run)


@interface_group.command("build")
@out_option
@dry_option
def build(out: Path | None, dry_run: bool) -> None:
    """Build all units of the PDB on one node (3 size tiers). Rerun the same command to resume."""
    _submit("build_full.sbatch", [str(out or default_out())], dry_run)


@interface_group.command("finalize")
@out_option
@dry_option
def finalize(out: Path | None, dry_run: bool) -> None:
    """Manifest with clusters, report, and the published shard collection."""
    _submit("finalize.sbatch", [str(out or default_out())], dry_run)


@interface_group.command("split")
@out_option
@dry_option
def split(out: Path | None, dry_run: bool) -> None:
    """Date split (train / valid / test) with leakage flags -> split.tsv."""
    _submit("make_split.sbatch", [str(out or default_out())], dry_run)


@interface_group.command("validate")
@out_option
@click.option("--n", type=int, default=1000, show_default=True, help="Units to recompute.")
@click.option("--seed", type=int, default=0, show_default=True)
@dry_option
def validate(out: Path | None, n: int, seed: int, dry_run: bool) -> None:
    """Recompute contacts, flags and waters of a random sample from the stored coordinates."""
    _submit("validate_db.sbatch", [str(out or default_out()), str(n), str(seed)], dry_run, ["-c", "4", "--mem=16g"])
