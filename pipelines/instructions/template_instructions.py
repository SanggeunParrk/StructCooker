import os
from pathlib import Path
import subprocess
import fnmatch


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


def load_a3m_list(data_dir: Path, pattern: str = "P*.a3m") -> list[dict[str, Path]]:
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
                            "output_path": Path(entry.path).with_suffix(".hhm"),
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
    