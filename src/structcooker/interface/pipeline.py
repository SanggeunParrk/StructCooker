"""CIFMol -> deduplicated interface units (everything ``detect`` needs, wired together)."""

from __future__ import annotations

import numpy as np

from structcooker.interface.adapt import build_assembly
from structcooker.interface.config import Config
from structcooker.interface.detect import dedup, drop_ion_subsets, find_units
from structcooker.interface.duplex import find_duplex


def extract_units(cm, cfg: Config, ab_lookup=None):
    """Return ``(assembly, units, duplex)``; units are ion-subset- and symmetry-deduplicated."""
    asm = build_assembly(cm, cfg.excluded, cfg.min_heavy)
    if not len(asm.kinds):
        return asm, [], {}
    ab = [ab_lookup(k[2:]) if (ab_lookup and kind == "protein") else None
          for k, kind in zip(asm.keys, asm.kinds, strict=True)]
    is_na = np.array([k == "na" for k in asm.kinds])[asm.node_idx]
    duplex = {}
    if is_na.any():
        d = cfg.duplex
        duplex = find_duplex(asm.xyz, asm.node_idx, asm.res_idx, asm.atom_name, asm.res_code, is_na,
                             d["min_bp"], d["min_hbonds"], d["hb_cutoff"], tuple(d["c1_range"]))
    units = find_units(asm.xyz, asm.node_idx, asm.res_idx, asm.kinds, cfg.contact, cfg.unit,
                       donor=asm.donor, ab_class=ab, duplex=duplex)
    units = dedup(drop_ion_subsets(units, asm.kinds), asm.keys, cfg.unit.weight_bucket)
    asm.stats["ab_class"] = ab
    return asm, units, duplex
