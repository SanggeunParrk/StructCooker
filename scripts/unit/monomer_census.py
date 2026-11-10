"""How many contact and monomer units a sample of cif_pdb entries gives (sbatch only).

usage: monomer_census.py <cif_pdb.lmdb> <index.tsv> <out.tsv> (--random N SEED | --ids ID [ID ...]) [--max-mb MB]
env:   UNIT_OVERRIDES="monomer.min_residues=2"   dotted config overrides (see params/unit.yaml)

One row per entry: contact units, monomer units (n_seed_edges == 0), the node kinds of the monomers, and the
residue counts of the monomer chains. The same ``entry_units`` the ``unit/unit_pdb`` rebuild runs.
"""

import os
import random
import sys
import time
from collections import Counter

from structcooker.instructions.readers.lmdb import read_lmdb
from structcooker.instructions.unit.config import load_config
from structcooker.instructions.unit.units import entry_units

args = sys.argv[1:]
db, index, out = args[:3]
rest = args[3:]
max_mb = float(rest[rest.index("--max-mb") + 1]) if "--max-mb" in rest else 200
ov = dict(kv.split("=") for kv in os.environ.get("UNIT_OVERRIDES", "").split(",") if kv)
ov = {k: (float(v) if "." in v else int(v)) for k, v in ov.items()}
cfg = load_config(overrides=ov)
if "--ids" in rest:
    keys = [
        k.lower() for k in rest[rest.index("--ids") + 1 :] if not k.startswith("--")
    ]
else:
    n, seed = (
        int(rest[rest.index("--random") + 1]),
        int(rest[rest.index("--random") + 2]),
    )
    rows = [ln.split("\t") for ln in open(index) if not ln.startswith(("#", "key"))]
    keys = [r[0] for r in rows if int(r[1]) < max_mb * 1e6]
    random.Random(seed).shuffle(keys)
    keys = keys[:n]

tot = Counter()
t0 = time.time()
with open(out, "w") as f:
    f.write(
        "pdb_id\tn_contact\tn_monomer\tmonomer_kinds\tmonomer_residues\tcontact_types\n"
    )
    for key in keys:
        try:
            subs = entry_units(read_lmdb(db, key), cfg, lambda _s: None)
        except Exception as e:  # noqa: BLE001
            f.write(f"{key}\tERROR\t{type(e).__name__}\t\t\t\n")
            tot["error"] += 1
            continue
        mono, contact = [], []
        for u in subs.values():
            ui = u["unit_input"]
            nodes = [(ui.asm.kinds[n], ui.asm.keys[n]) for n in ui.unit.nodes]
            (contact if ui.unit.seeds else mono).append(nodes)
        tot["entries"] += 1
        tot["contact"] += len(contact)
        tot["monomer"] += len(mono)
        tot["entries_with_monomer"] += bool(mono)
        tot["entries_monomer_only"] += bool(mono) and not contact
        tot["entries_no_unit"] += not subs
        for m in mono:
            tot["mono_" + "+".join(k for k, _ in m)] += 1
        f.write(
            f"{key}\t{len(contact)}\t{len(mono)}\t"
            + ",".join("+".join(k for k, _ in m) for m in mono)
            + "\t"
            + ",".join("+".join(str(len(s) - 2) for _, s in m) for m in mono)
            + "\t"
            + ",".join("+".join(k[0] for k, _ in c) for c in contact[:6])
            + "\n"
        )
print(f"{time.time() - t0:.0f}s", dict(tot))
