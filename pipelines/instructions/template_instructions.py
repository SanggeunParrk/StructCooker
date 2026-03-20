import os
from pathlib import Path
import subprocess
import fnmatch

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


def load_a3m_list(data_dir: Path, output_dir: Path, pattern: str = "P*.a3m", output_pattern: str = ".hhm") -> list[dict[str, Path]]:
    result = []
    def _scan(dir_path: Path):
        with os.scandir(dir_path) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False):
                    _scan(Path(entry.path))
                elif fnmatch.fnmatch(entry.name, pattern):
                    result.append(
                        {
                            "input_a3m_path": Path(entry.path),
                            "output_path": output_dir / f"{entry.name}{output_pattern}",
                        }
                    )
    
    _scan(data_dir)
    return result


def run_hhmake(input_a3m_path: Path, output_path: Path | None) -> str:
    # run hhmake to convert a3m to hhm
    if output_path is None:
        output_path = input_a3m_path.with_suffix(".hhm")
    command = [
        "hhmake",
        "-i",
        str(input_a3m_path),
        "-o",
        str(output_path),
    ]
    try:
        _run_command(command)
        return "hhm file created at: " + str(output_path)
    except Exception as e:
        msg = f"Error running hhmake for {input_a3m_path}: {e}"
        raise RuntimeError(msg) from e


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



def run_hhsearch(
    msa_path: Path,
    hhr_path: Path,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_template: Path,
    hhsuite_bin_dir: Path,
) -> str:
    """Run HHsearch for one MSA directory."""
    db_template = Path(db_template)
    hhsuite_bin_dir = Path(hhsuite_bin_dir)

    hhsuite_env = os.environ.copy()
    hhsuite_env["HHLIB"] = str(hhsuite_bin_dir)
    hhsuite_env["PATH"] = f"{hhsuite_bin_dir}:{hhsuite_env.get('PATH', '')}"

    if not _is_nonempty(msa_path):
        msg = f"MSA file does not exist or is empty: {msa_path}"
        raise FileNotFoundError(msg)

    if _is_nonempty(hhr_path):
        print(f"Output file {hhr_path} already exists and is non-empty. Skipping HHsearch for {msa_path}.")
        return f"Skip {hhr_path.name} (already exists and is non-empty)"


    _run_command(
        _hhsearch_command(
            db_template=db_template,
            cpu=cpu,
            mem=mem,
            in_msa=msa_path,
            out_hhr=hhr_path,
            out_atab=hhr_path.with_suffix(".atab"),
        ),
        env=hhsuite_env,
    )
    print(f"HHsearch completed for {msa_path}, output saved to {hhr_path}")
    return f"Done {hhr_path.name}"



def run_hmmbuild(input_a3m_path: Path, hmm_path: Path | None) -> str:
    # run hmmbuild to convert a3m to hmm
    if hmm_path is None:
        hmm_path = input_a3m_path.with_suffix(".hmm")
    command = [
        "hmmbuild",
        str(hmm_path),
        str(input_a3m_path),
    ]
    try:
        _run_command(command)
        return "hmm file created at: " + str(hmm_path)
    except Exception as e:
        msg = f"Error running hmmbuild for {input_a3m_path}: {e}"
        raise RuntimeError(msg) from e
    
def run_hmmsearch(
    hmm_path: Path,
    fasta_path: Path,
) -> str:
    if not _is_nonempty(hmm_path):
        msg = f"HMM file does not exist or is empty: {hmm_path}"
        raise FileNotFoundError(msg)
    if not _is_nonempty(fasta_path):
        msg = f"FASTA file does not exist or is empty: {fasta_path}"
        raise FileNotFoundError(msg)
    command = [
        "hmmsearch",
        "--noali",
        "--F1",
        "0.1",
        "--F2",
        "0.1",
        "--F3",
        "0.1",
        "--E",
        "100",
        "--incE",
        "100",
        "--domE",
        "100",
        "--incdomE",
        "100",
        str(hmm_path),
        str(fasta_path),
    ]
    try:
        _run_command(command)
        return f"hmmsearch completed for {hmm_path} against {fasta_path}"
    except Exception as e:
        msg = f"Error running hmmsearch for {hmm_path}: {e}"
        raise RuntimeError(msg) from e
