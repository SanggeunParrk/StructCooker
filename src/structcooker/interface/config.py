"""Load configs/interface/interface.yaml into the dataclasses ``detect`` / ``duplex`` / ``adapt`` use."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from structcooker.interface.detect import ContactCfg, UnitCfg
from structcooker.interface.ligands import load_excluded

DEFAULT = Path(__file__).resolve().parents[3] / "configs" / "interface" / "interface.yaml"


@dataclass
class Config:
    contact: ContactCfg
    unit: UnitCfg
    excluded: frozenset
    min_heavy: int
    water_near: float
    water_bridge: float
    duplex: dict
    raw: dict


def load_config(path: Path | str = DEFAULT, overrides: dict | None = None) -> Config:
    cfg = yaml.safe_load(Path(path).read_text())
    for dotted, v in (overrides or {}).items():
        d = cfg
        *head, last = dotted.split(".")
        for h in head:
            d = d[h]
        d[last] = v
    u = cfg["unit"]
    unit = UnitCfg(
        max_polymer_chains=u["max_polymer_chains"], max_ligands=u["max_ligands"], max_ions=u["max_ions"],
        join_min_residues=u["join_min_residues"], allowed_pairs=frozenset(tuple(p) for p in u["allowed_pairs"]),
        max_polymer_chains_in_assembly=cfg["assembly"]["max_polymer_chains"], max_units=cfg["assembly"]["max_units"],
        weight_bucket=cfg["dedup"]["weight_bucket"], max_forced_polymer=u.get("max_forced_polymer", 10),
    )
    return Config(
        ContactCfg(**cfg["contact"]), unit, load_excluded(cfg["ligand"]["exclude_categories"]),
        cfg["ligand"]["min_heavy_atoms"], cfg["water"]["near"], cfg["water"]["bridge"], cfg["duplex"], cfg,
    )
