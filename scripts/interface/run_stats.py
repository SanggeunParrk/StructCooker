"""Sample N random cif_pdb entries, extract interface units, write summary stats (sbatch only).

usage: run_stats.py <cif_pdb.lmdb> <index.tsv> <N> <seed> <out_dir> [max_value_MB]
env:   IF_OVERRIDES="contact.cutoff=4.0,ligand.min_heavy_atoms=8"  (sensitivity runs)
"""
import json
import os
import random
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

import numpy as np

from structcooker.instructions.readers.lmdb import read_lmdb
from structcooker.instructions.transforms.cifmol import convert_to_cifmol_dict
from structcooker.interface.adapt import pick_alt_keys
from structcooker.interface.antibody import AntibodyClassifier
from structcooker.interface.config import load_config
from structcooker.interface.detect import select_waters
from structcooker.interface.pipeline import extract_units

db, index, n, seed, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), Path(sys.argv[5])
max_mb = float(sys.argv[6]) if len(sys.argv) > 6 else 200
out.mkdir(parents=True, exist_ok=True)
ov = {k: yaml_v for k, yaml_v in (kv.split("=") for kv in os.environ.get("IF_OVERRIDES", "").split(",") if kv)}
ov = {k: (float(v) if "." in v else int(v)) for k, v in ov.items()}
cfg = load_config(overrides=ov)
ab = AntibodyClassifier()

rows = [l.rstrip("\n").split("\t") for l in open(index) if not l.startswith(("#", "key"))]
rows = [r[0] for r in rows if int(r[1]) < max_mb * 1e6]
random.Random(seed).shuffle(rows)
keys = rows[:n]

C = {k: Counter() for k in ("unit_type", "n_polymer", "n_nodes", "ligand_in_units", "ion_in_units", "excluded_names",
                            "node_kinds", "skip", "waters_bucket", "truncated", "dropped", "special")}
t0 = time.time()
unit_rows = []
examples = []
for i, key in enumerate(keys):
    try:
        rec = read_lmdb(db, key)
        cms = convert_to_cifmol_dict(rec)
        res = rec["metadata_dict"].get("resolution")
        seen = set()
        for ck in pick_alt_keys(cms):
            asm, units, duplex = extract_units(cms[ck]["cifmol"], cfg, ab)
            C["node_kinds"].update(asm.kinds)
            C["excluded_names"].update(asm.stats["excluded_names"])
            C["dropped"]["covalent_merged"] += asm.stats["n_merged"]
            for k2, v in asm.stats["dropped"].items():
                C["dropped"][k2] += v
            if sum(k in ("protein", "na") for k in asm.kinds) > cfg.unit.max_polymer_chains_in_assembly:
                C["skip"]["too_many_chains"] += 1
                continue
            abc = asm.stats.get("ab_class", [])
            C["special"]["assemblies_with_duplex"] += bool(duplex)
            C["special"]["assemblies_with_na"] += "na" in asm.kinds
            for un in units:
                sig = tuple(sorted(asm.keys[x] for x in un.nodes))
                if (sig, round(un.weight / 10)) in seen:
                    continue
                seen.add((sig, round(un.weight / 10)))
                kinds = [asm.kinds[x] for x in un.nodes]
                kc = Counter(kinds)
                typ = "+".join(f"{kc[k]}{k}" for k in ("protein", "na", "ligand", "ion") if kc[k])
                C["unit_type"][typ] += 1
                C["n_polymer"][kc["protein"] + kc["na"]] += 1
                C["n_nodes"][len(un.nodes)] += 1
                C["truncated"][un.truncated] += 1
                has_ab = [abc[x] for x in un.nodes if abc and abc[x]]
                if has_ab:
                    C["special"]["units_with_ab_or_tcr"] += 1
                    C["special"]["ab_unit_ab_chain_set:" + "".join(sorted(has_ab))] += 1
                if kc["na"]:
                    nas = [x for x in un.nodes if asm.kinds[x] == "na"]
                    paired = sum(1 for (a, b) in duplex if a in nas and b in nas)
                    C["special"]["na_units_paired_strands" if paired else "na_units_single_strand"] += 1
                    if len(nas) >= 2 and not paired:
                        C["special"]["na_multi_nodes_unpaired_unit"] += 1
                        examples.append((key, ck, [asm.labels[x] for x in nas], [int((asm.node_idx == x).sum()) for x in nas]))
                for x in un.nodes:
                    if asm.kinds[x] == "ligand":
                        C["ligand_in_units"][asm.keys[x][2:]] += 1
                    if asm.kinds[x] == "ion":
                        C["ion_in_units"][asm.keys[x][2:]] += 1
                m = np.isin(asm.node_idx, un.nodes)
                keep, br = select_waters(asm.xyz[m], asm.node_idx[m], asm.water_xyz, cfg.water_near, cfg.water_bridge)
                C["waters_bucket"][min(len(keep) // 25 * 25, 200)] += 1
                unit_rows.append((key, ck, typ, len(un.nodes), un.weight, len(keep), int(br.sum()), res))
    except Exception as e:  # noqa: BLE001
        if sum(v for k, v in C["skip"].items() if k.startswith("error")) < 3:
            print(f"ERROR {key}:", traceback.format_exc()[-700:], flush=True)
        C["skip"][f"error:{type(e).__name__}"] += 1
    if (i + 1) % 50 == 0:
        print(f"{i+1}/{len(keys)} entries, {len(unit_rows)} units, {time.time()-t0:.0f}s", flush=True)

json.dump({k: dict(v.most_common()) for k, v in C.items()}, open(out / "stats.json", "w"), indent=1, default=str)
json.dump(examples, open(out / "na_unpaired_examples.json", "w"))
with open(out / "units.tsv", "w") as f:
    f.write("entry\tasm_key\ttype\tn_nodes\tweight\tn_water\tn_bridging\tresolution\n")
    for r in unit_rows:
        f.write("\t".join(map(str, r)) + "\n")
print("done", len(keys), "entries,", len(unit_rows), "units,", f"{time.time()-t0:.0f}s")
