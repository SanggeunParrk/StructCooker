import os
import shutil
import subprocess
from pathlib import Path

from joblib import Parallel, delayed

DEFAULT_DB_UR30 = Path(
    "/data/shared/cssb_data/db_protSeq/uniref30/2022_02/UniRef30_2022_02",
)
DEFAULT_DB_BFD = Path(
    "/data/shared/cssb_data/db_protSeq/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
)
DEFAULT_HHSUITE_BIN = Path("/software/hhsuite/build/bin")


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _count_header_lines(path: Path) -> int:
    with path.open("r") as f:
        return sum(1 for line in f if line.startswith(">"))


def _run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
) -> None:
    result = subprocess.run(command, capture_output=True, text=True, env=env, check=False)
    if result.returncode != 0:
        cmd = " ".join(command)
        msg = (
            f"Command failed ({result.returncode}): {cmd}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        raise RuntimeError(msg)


def _hhblits_command(
    db_path: Path,
    cpu: int,
    mem: int,
    in_a3m: Path,
    out_a3m: Path,
    evalue: str,
) -> list[str]:
    return [
        "hhblits",
        "-o",
        "/dev/null",
        "-mact",
        "0.35",
        "-maxfilt",
        "100000000",
        "-neffmax",
        "20",
        "-cov",
        "25",
        "-cpu",
        str(cpu),
        "-nodiff",
        "-realign_max",
        "100000000",
        "-maxseq",
        "1000000",
        "-maxmem",
        str(mem),
        "-n",
        "4",
        "-d",
        str(db_path),
        "-i",
        str(in_a3m),
        "-oa3m",
        str(out_a3m),
        "-e",
        evalue,
        "-v",
        "0",
    ]


def _hhfilter_command(in_a3m: Path, out_a3m: Path, coverage: int) -> list[str]:
    return [
        "hhfilter",
        "-maxseq",
        "100000",
        "-id",
        "90",
        "-cov",
        str(coverage),
        "-i",
        str(in_a3m),
        "-o",
        str(out_a3m),
    ]

def make_input_fasta(
    seqid: str,
    sequence: str,
    output_dir: Path,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_dir = output_dir / seqid[:4] / seqid
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / f"{seqid}.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open("w") as f:
        f.write(f">{seqid}\n{sequence}\n")
    return fasta_path, out_dir

def run_signalp(  # noqa: PLR0913
    input_fasta: Path,
    out_dir: Path,
    *,
    signalp_mode: str = "fast",
) -> Path:
    """Run SignalP + HHblits/HHfilter for one FASTA."""
    out_dir.mkdir(parents=True, exist_ok=True)
    signalp_dir = out_dir / "signalp"
    signalp_dir.mkdir(parents=True, exist_ok=True)

    _run_command(
        [
            "signalp6",
            "--fastafile",
            str(input_fasta),
            "--organism",
            "other",
            "--output_dir",
            str(signalp_dir),
            "--format",
            "none",
            "--mode",
            signalp_mode,
        ],
    )

    trim_fasta = signalp_dir / "processed_entries.fasta"
    return trim_fasta if _is_nonempty(trim_fasta) else input_fasta

def run_msa_search(  # noqa: PLR0913
    input_fasta: Path,
    out_dir: Path,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_ur30: Path,
    db_bfd: Path,
    hhsuite_bin_dir: Path,
) -> str:

    hhsuite_env = os.environ.copy()
    hhsuite_env["HHLIB"] = str(hhsuite_bin_dir)
    hhsuite_env["PATH"] = f"{hhsuite_bin_dir}:{hhsuite_env.get('PATH', '')}"

    hhblits_dir = out_dir / "hhblits"
    hhblits_dir.mkdir(parents=True, exist_ok=True)
    msa0_file = out_dir / "t000_msa0.a3m"

    if not _is_nonempty(msa0_file):
        prev_a3m = input_fasta
        for evalue in ("1e-10", "1e-6", "1e-3"):
            a3m_file = hhblits_dir / f"t000_.{evalue}.a3m"
            if not _is_nonempty(a3m_file):
                _run_command(
                    _hhblits_command(
                        db_path=db_ur30,
                        cpu=cpu,
                        mem=mem,
                        in_a3m=prev_a3m,
                        out_a3m=a3m_file,
                        evalue=evalue,
                    ),
                    env=hhsuite_env,
                )

            id90cov75_file = hhblits_dir / f"t000_.{evalue}.id90cov75.a3m"
            id90cov50_file = hhblits_dir / f"t000_.{evalue}.id90cov50.a3m"
            _run_command(
                _hhfilter_command(
                    in_a3m=a3m_file,
                    out_a3m=id90cov75_file,
                    coverage=75,
                ),
                env=hhsuite_env,
            )
            _run_command(
                _hhfilter_command(
                    in_a3m=a3m_file,
                    out_a3m=id90cov50_file,
                    coverage=50,
                ),
                env=hhsuite_env,
            )
            prev_a3m = id90cov50_file

            n75 = _count_header_lines(id90cov75_file)
            n50 = _count_header_lines(id90cov50_file)

            if n75 > 2000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(id90cov75_file, msa0_file)
                    break
            elif n50 > 4000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(id90cov50_file, msa0_file)
                    break

        if not _is_nonempty(msa0_file):
            evalue = "1e-3"
            bfd_a3m_file = hhblits_dir / f"t000_.{evalue}.bfd.a3m"
            if not _is_nonempty(bfd_a3m_file):
                _run_command(
                    _hhblits_command(
                        db_path=db_bfd,
                        cpu=cpu,
                        mem=mem,
                        in_a3m=prev_a3m,
                        out_a3m=bfd_a3m_file,
                        evalue=evalue,
                    ),
                    env=hhsuite_env,
                )

            bfd_id90cov75_file = hhblits_dir / f"t000_.{evalue}.bfd.id90cov75.a3m"
            bfd_id90cov50_file = hhblits_dir / f"t000_.{evalue}.bfd.id90cov50.a3m"
            _run_command(
                _hhfilter_command(
                    in_a3m=bfd_a3m_file,
                    out_a3m=bfd_id90cov75_file,
                    coverage=75,
                ),
                env=hhsuite_env,
            )
            _run_command(
                _hhfilter_command(
                    in_a3m=bfd_a3m_file,
                    out_a3m=bfd_id90cov50_file,
                    coverage=50,
                ),
                env=hhsuite_env,
            )
            prev_a3m = bfd_id90cov50_file

            n75 = _count_header_lines(bfd_id90cov75_file)
            n50 = _count_header_lines(bfd_id90cov50_file)
            if n75 > 2000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov75_file, msa0_file)
            elif n50 > 4000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov50_file, msa0_file)

        if not _is_nonempty(msa0_file):
            shutil.copyfile(prev_a3m, msa0_file)

    return f"Done {input_fasta}"

