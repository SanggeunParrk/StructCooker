import os
import subprocess
from pathlib import Path

from joblib import Parallel, delayed

DEFAULT_DB_TEMPL = Path("/data/shared/cssb_data/pdb100/pdb100_2021Mar03/pdb100_2021Mar03")
DEFAULT_HHSUITE_BIN = Path("/software/hhsuite/build/bin")


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


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


def _hhsearch_command(
    db_template: Path,
    cpu: int,
    mem: int,
    in_msa: Path,
    out_hhr: Path,
    out_atab: Path,
) -> list[str]:
    return [
        "hhsearch",
        "-b",
        "50",
        "-B",
        "500",
        "-z",
        "50",
        "-Z",
        "500",
        "-mact",
        "0.05",
        "-cpu",
        str(cpu),
        "-maxmem",
        str(mem),
        "-aliw",
        "100000",
        "-e",
        "100",
        "-p",
        "5.0",
        "-d",
        str(db_template),
        "-i",
        str(in_msa),
        "-o",
        str(out_hhr),
        "-atab",
        str(out_atab),
        "-v",
        "0",
    ]


def run_template_search(
    msa_dir: Path | str,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_template: Path | str = DEFAULT_DB_TEMPL,
    hhsuite_bin_dir: Path | str = DEFAULT_HHSUITE_BIN,
    msa_file_name: str = "t000_msa0.a3m",
) -> str:
    """Run HHsearch for one MSA directory."""
    msa_dir = Path(msa_dir)
    db_template = Path(db_template)
    hhsuite_bin_dir = Path(hhsuite_bin_dir)

    msa_path = msa_dir / msa_file_name
    if not _is_nonempty(msa_path):
        msg = f"MSA file does not exist or is empty: {msa_path}"
        raise FileNotFoundError(msg)

    hhr_path = msa_dir / "t000_.hhr"
    if _is_nonempty(hhr_path):
        return f"Skip {msa_dir.name}"

    hhsuite_env = os.environ.copy()
    hhsuite_env["HHLIB"] = str(hhsuite_bin_dir)
    hhsuite_env["PATH"] = f"{hhsuite_bin_dir}:{hhsuite_env.get('PATH', '')}"

    _run_command(
        _hhsearch_command(
            db_template=db_template,
            cpu=cpu,
            mem=mem,
            in_msa=msa_path,
            out_hhr=hhr_path,
            out_atab=msa_dir / "t000_.atab",
        ),
        env=hhsuite_env,
    )
    return f"Done {msa_dir.name}"


def run_template_search_batch(  # noqa: PLR0913
    msa_root: Path | str,
    *,
    cpu_per_job: int = 4,
    mem: int = 20,
    num_jobs: int = 4,
    db_template: Path | str = DEFAULT_DB_TEMPL,
    hhsuite_bin_dir: Path | str = DEFAULT_HHSUITE_BIN,
    msa_file_name: str = "t000_msa0.a3m",
) -> dict[str, str]:
    """Run HHsearch for all MSA directories under msa_root."""
    msa_root = Path(msa_root)
    msa_dirs = sorted(d for d in msa_root.iterdir() if d.is_dir())

    def _process(msa_dir: Path) -> tuple[str, str]:
        try:
            result = run_template_search(
                msa_dir=msa_dir,
                cpu=cpu_per_job,
                mem=mem,
                db_template=db_template,
                hhsuite_bin_dir=hhsuite_bin_dir,
                msa_file_name=msa_file_name,
            )
            return msa_dir.name, result
        except Exception as exc:  # noqa: BLE001
            return msa_dir.name, f"error: {exc}"

    results = Parallel(n_jobs=num_jobs)(delayed(_process)(msa_dir) for msa_dir in msa_dirs)
    return dict(results)
