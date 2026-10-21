"""Run interface-unit detection over a few PDB entries and print what comes out. Run via sbatch.

usage: run_sample.py <cif_pdb.lmdb> <antibody_cache.tsv|-> <key> [<key> ...]
"""
import sys
import time
from collections import Counter

import numpy as np

from structcooker.instructions.readers.lmdb import read_lmdb
from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict
from structcooker.interface.adapt import pick_alt_keys
from structcooker.interface.antibody import AntibodyClassifier
from structcooker.interface.config import load_config
from structcooker.interface.detect import select_waters
from structcooker.interface.pipeline import extract_units

db, cache, keys = sys.argv[1], sys.argv[2], sys.argv[3:]
cfg = load_config()
ab = AntibodyClassifier(None if cache == "-" else cache)
for key in keys:
    t0 = time.time()
    try:
        rec = read_lmdb(db, key)
    except Exception as e:  # noqa: BLE001
        print(f"=== {key}: unreadable ({type(e).__name__})")
        continue
    cms = convert_to_cifmol_dict(rec)
    print(f"=== {key}: {len(cms)} keys, res={rec['metadata_dict'].get('resolution')}")
    for ck in pick_alt_keys(cms):
        asm, units, duplex = extract_units(cms[ck]["cifmol"], cfg, ab)
        abc = asm.stats.get("ab_class", [])
        print(f"  [{ck}] nodes={dict(Counter(asm.kinds))} ab={dict(Counter(a for a in abc if a))} "
              f"duplex_pairs={ {(asm.labels[a], asm.labels[b]): n for (a, b), n in duplex.items()} } "
              f"units={len(units)} water={asm.stats['n_water']}")
        for u in units[:10]:
            lab = [f"{asm.kinds[n][0]}:{asm.labels[n]}" + (f"[{abc[n]}]" if abc and abc[n] else "") for n in u.nodes]
            m = np.isin(asm.node_idx, u.nodes)
            keep, br = select_waters(asm.xyz[m], asm.node_idx[m], asm.water_xyz, cfg.water_near, cfg.water_bridge)
            print(f"     {lab} w={u.weight} seeds={len(u.seeds)} trunc={u.truncated} waters={len(keep)} bridging={int(br.sum())}")
    print(f"  ({time.time()-t0:.1f}s)")
