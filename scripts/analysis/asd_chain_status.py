"""Light-chain status of the ASD pair table, with the SelfDistillation (kkh517) rules and POISK as the recovery source.

Per pair: paired (ASD has a light chain), scfv (linker in the heavy chain or the scfv flag), vhh (nanobody flag), else
heavy_only. For the heavy_only pairs: the POISK rows whose VH equals the heavy chain exactly give 1 distinct VL (fillable),
2+ (ambiguous: kkh517 does not fill) or none. Also the agreement of POISK's unique VL with the VL ASD already has.
Usage: asd_chain_status.py ASD_DIR PAIRS SEQ_ID_MAP POISK_CSV
"""
import collections, csv, glob, re, sys

import pyarrow.parquet as pq

asd, pairs_tsv, seq_id_map, poisk_csv = sys.argv[1:5]
csv.field_size_limit(10**9)
seq = {}; id_of = {}
for l in open(seq_id_map):
    i, s = l.rstrip("\n").split("\t"); seq[i] = s; id_of[s] = i
pairs = {}
for l in list(open(pairs_tsv))[1:]:
    a, h, li, ds = l.rstrip("\n").split("\t"); pairs[(h, li, a)] = ds.split(",")
norm = lambda x: "".join((x or "").split()).upper()
LINK = re.compile(r"(G{3,5}S){2,}")
flag = collections.defaultdict(lambda: [False, False])  # triple -> [nanobody, scfv]
light_only = collections.Counter(); heavy_only_rows = collections.Counter(); both_empty = collections.Counter()
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "chains", "nanobody", "scfv"]).to_pylist():
        ab = dict(r["chains"]["antibody_chains"] or []); tg = dict(r["chains"]["target_chains"] or [])
        h, li, a = norm(ab.get("heavy")), norm(ab.get("light")), norm(tg.get("full"))
        if not h and li: light_only[r["dataset"]] += 1
        if not h and not li: both_empty[r["dataset"]] += 1
        if not h or not a: continue
        k = (id_of.get(h), id_of.get(li) if li else "", id_of.get(a))
        if k in pairs:
            flag[k][0] |= bool(r["nanobody"]); flag[k][1] |= bool(r["scfv"])
print("ASD rows without a heavy chain: light only", sum(light_only.values()), dict(light_only.most_common(6)), "| nothing", sum(both_empty.values()))
# POISK
vh_to_vl = collections.defaultdict(set)
for x in csv.DictReader(open(poisk_csv, newline="")):
    h = norm(x["vh.sequence_from_patent"]); v = norm(x["vl.sequence_from_patent"])
    if h and v: vh_to_vl[h].add(v)
print("POISK distinct VH with a VL:", len(vh_to_vl), " with exactly one VL:", sum(len(v) == 1 for v in vh_to_vl.values()))
status = collections.defaultdict(collections.Counter); agree = collections.Counter()
for (h, li, a), ds in pairs.items():
    hs = seq[h]; d = ds[0] if len(ds) == 1 else "+".join(sorted(ds))
    if li:
        s = "paired"
        cand = vh_to_vl.get(hs)
        if cand and len(cand) == 1: agree["POISK unique VL = ASD light" if next(iter(cand)) == seq[li] else "POISK unique VL != ASD light"] += 1
        elif cand: agree["POISK ambiguous"] += 1
    elif flag[(h, li, a)][0]: s = "vhh"
    elif flag[(h, li, a)][1] or LINK.search(hs): s = "scfv"
    else:
        cand = vh_to_vl.get(hs)
        s = "heavy_only/POISK 1 VL" if cand and len(cand) == 1 else ("heavy_only/POISK 2+ VL" if cand else "heavy_only/no POISK")
    status[d][s] += 1
tot = collections.Counter()
for d, c in status.items(): tot.update(c)
print("\nstatus of", len(pairs), "pairs:", dict(tot.most_common()))
print("POISK check on paired pairs whose heavy is in POISK:", dict(agree))
print(f"\n{'dataset':<26}" + "".join(f"{k:>24}" for k in ("paired", "scfv", "vhh", "heavy_only/POISK 1 VL", "heavy_only/POISK 2+ VL", "heavy_only/no POISK")))
for d, c in sorted(status.items(), key=lambda x: -sum(x[1].values())):
    print(f"{d:<26}" + "".join(f"{c[k]:>24}" for k in ("paired", "scfv", "vhh", "heavy_only/POISK 1 VL", "heavy_only/POISK 2+ VL", "heavy_only/no POISK")))
