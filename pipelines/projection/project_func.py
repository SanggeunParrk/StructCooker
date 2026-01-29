from pathlib import Path
from typing import Protocol


class ProjectFunc(Protocol):
    """Protocol for data loading functions."""

    def __call__(self, data: object, output_path: Path) -> None:
        """Load a file and return a dict-like data structure."""
        ...
