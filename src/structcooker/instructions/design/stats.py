"""Per-item statistics and backbone archives of the design monomers (stage ``design/<set>_stats``).

The source DB is cut into parts of ``chunk_keys`` consecutive keys (sorted, so the cut depends on the source alone).
One part is one work item of the ``parallel`` op: it writes ``part_<n>.tsv`` (one row of statistics per item),
``part_<n>.tar`` (one backbone PDB per item, the input of ``foldseek createdb``) and, if a record could not be read,
``part_<n>.failed.tsv``. The ``.tsv`` is written last, so its presence marks a finished part and a resumed run lists
only the parts without one.
"""

from __future__ import annotations

import io
import tarfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from datacooker.lmdb.index import read_index
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.design.export import backbone_pdb_text
from structcooker.instructions.design.ss import secondary_structure
from structcooker.instructions.mols import CIFMolAttached
from structcooker.instructions.transforms.codecs import from_bytes

if TYPE_CHECKING:
    from structcooker.instructions.mols import CIFMol

STAT_COLUMNS = (
    "id",
    "length",
    "mean_plddt",
    "helix",
    "strand",
    "coil",
    "max_coil_run",
)


def source_keys(db_path: Path) -> list[str]:
    """Sorted keys of ``db_path``: its ``.index.tsv`` sidecar (listing a 450 GB LMDB reads all of it), else a scan."""
    keys = [item.key for item in read_index(db_path)]
    if not keys:
        with (
            open_env(str(db_path), readonly=True, lock=False) as env,
            env.begin() as txn,
        ):
            keys = [bytes(k).decode() for k in txn.cursor().iternext(values=False)]
    return sorted(keys)


def part_name(number: int) -> str:
    """File stem of part ``number``."""
    return f"part_{number:05d}"


def pending_parts(
    db_path: Path, output_dir: Path, chunk_keys: int
) -> list[dict[str, Any]]:
    """Work items ``{"part": stem, "keys": [...]}`` of the parts whose ``.tsv`` is not written yet."""
    keys = source_keys(Path(db_path))
    if not keys:
        msg = f"{db_path} holds no records"
        raise ValueError(msg)
    parts = []
    for number, start in enumerate(range(0, len(keys), chunk_keys)):
        stem = part_name(number)
        if not (Path(output_dir) / f"{stem}.tsv").exists():
            parts.append({"part": stem, "keys": keys[start : start + chunk_keys]})
    return parts


def load_item(raw: bytes) -> CIFMol:
    """Decode a source record (``{assembly key: {"cifmol_attached_dict": ...}}``) into its first assembly."""
    record = from_bytes(raw)
    attached = CIFMolAttached.from_dict(
        record[sorted(record)[0]]["cifmol_attached_dict"]
    )
    return cast("CIFMol", attached)


def item_stats(cm: CIFMol) -> tuple[list[str], str | None]:
    """``(columns after id, backbone PDB text)`` of one item; the PDB is None when no residue has a complete backbone.

    ``length`` counts the residues with a complete backbone, ``mean_plddt`` is the mean B factor of all atoms (OpenFold
    stores pLDDT there), the secondary-structure columns are fractions of the assessed residues
    (:func:`structcooker.instructions.design.ss.secondary_structure`).
    """
    _, s = secondary_structure(cm)
    plddt = float(np.mean(cm.atoms.b_factor.value))
    length = int(s["n_assessed"])
    columns = [
        str(length),
        str(round(plddt, 2)),
        str(round(float(s["helix"]), 4)),
        str(round(float(s["strand"]), 4)),
        str(round(float(s["coil"]), 4)),
        str(int(s["max_coil_run"])),
    ]
    return columns, backbone_pdb_text(cm) if length else None


def write_part(
    db_path: Path, output_dir: Path, backbone_dir: Path, part: str, keys: list[str]
) -> str:
    """Statistics and backbone archive of the items ``keys``; returns ``"<part>: <n> items, <m> failed"``.

    A record that cannot be read is listed in ``<part>.failed.tsv`` (key, error) and the part goes on; the items stage
    refuses the build if more records failed than the config allows.
    """
    output_dir, backbone_dir = Path(output_dir), Path(backbone_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    backbone_dir.mkdir(parents=True, exist_ok=True)
    tsv, tar_path, failed_path = (
        output_dir / f"{part}.tsv",
        backbone_dir / f"{part}.tar",
        output_dir / f"{part}.failed.tsv",
    )
    tmp = {p: p.with_name(p.name + ".tmp") for p in (tsv, tar_path, failed_path)}
    rows, failed, members = [], [], 0
    with (
        open_env(str(db_path), readonly=True, lock=False, readahead=False) as env,
        env.begin() as txn,
        tarfile.open(tmp[tar_path], "w") as tar,
    ):
        for key in keys:
            try:
                raw = txn.get(key.encode())
                if raw is None:
                    msg = "key not in the source"
                    raise KeyError(msg)  # noqa: TRY301 (caught below, recorded)
                columns, pdb = item_stats(load_item(bytes(raw)))
                if pdb is not None:
                    data = pdb.encode()
                    info = tarfile.TarInfo(f"{key}.pdb")
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))
                    members += 1
                rows.append("\t".join([key, *columns]))
            except Exception as error:  # noqa: BLE001, PERF203  (a failed record is listed, not fatal for the part)
                failed.append(
                    f"{key}\t{type(error).__name__}: {error}".replace("\n", " ")
                )
    if members:
        tmp[tar_path].replace(tar_path)
    else:
        tmp[tar_path].unlink()
    if failed:
        tmp[failed_path].write_text("\n".join(failed) + "\n")
        tmp[failed_path].replace(failed_path)
    tmp[tsv].write_text("".join(row + "\n" for row in rows))
    tmp[tsv].replace(tsv)
    return f"{part}: {len(rows)} items, {len(failed)} failed"
