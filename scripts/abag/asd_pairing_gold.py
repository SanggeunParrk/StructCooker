"""Error rate of heavy-light pairings against SAbDab, whose pairs come from crystal structures (true by construction).

Gold: for each SAbDab entry with both a heavy and a light chain, the sequences of those two chains from the PDB fasta.
Matching is containment either way (the PDB chain carries a constant domain or a tag, the ASD / POISK chain is often the variable domain), >= 100 aa.
Evaluated: (A) the light chains ASD already holds in the pair table, by dataset; (B) the POISK pairs, and the 'exactly one candidate' rule the
SelfDistillation recovery uses; (C) the POISK lights that would fill the heavy-only pairs of the pair table.
(D) with a sixth argument, the AntiBinder raw_data directory: its covid-19 (Heavy, Light) and biomap (antibody_seq_a / _b) pairs.
Usage: asd_pairing_gold.py SABDAB_SUMMARY PDB_FASTA PAIRS SEQ_ID_MAP POISK_CSV [ANTIBINDER_RAW_DATA]
"""

import collections
import csv
import sys

sab, pdb_fasta, pairs_tsv, seq_id_map, poisk_csv = sys.argv[1:6]
ab_dir = sys.argv[6] if len(sys.argv) > 6 else None
csv.field_size_limit(10**9)
norm = lambda x: "".join((x or "").split()).upper()
want = collections.defaultdict(list)  # PDB ID -> [(H, L)]
for l in list(open(sab))[1:]:
    f = l.rstrip("\n").split("\t")
    pdb, hc, lc = f[0].upper(), f[1], f[2]
    if hc not in ("NA", "") and lc not in ("NA", "") and hc != lc:
        want[pdb].append((hc, lc))
chain = {}
name = None
for l in open(pdb_fasta):
    if l.startswith(">"):
        p = l[1:].split(" | ")
        pid = p[0].split("_")[0]
        name = (
            (pid, p[2].strip()[5:])
            if pid in want and "polypeptide(L)" in p[1] and len(p) > 2
            else None
        )
    elif name:
        chain.setdefault(name, norm(l))
gold = collections.defaultdict(set)  # heavy -> {light}
for pid, hl in want.items():
    for hc, lc in hl:
        h, v = chain.get((pid, hc)), chain.get((pid, lc))
        if h and v and len(h) >= 100 and len(v) >= 60:
            gold[h].add(v)
print(
    "SAbDab entries",
    len(want),
    "gold distinct heavies",
    len(gold),
    "pairs",
    sum(len(v) for v in gold.values()),
)
K = 25
idx = collections.defaultdict(set)
for h in gold:
    for i in range(len(h) - K + 1):
        idx[h[i : i + K]].add(h)


def gold_of(h):
    """Gold heavies that contain h or are contained in h"""
    if len(h) < 100:
        return []
    out = set()
    for pos in (len(h) // 3, len(h) // 2, 2 * len(h) // 3):
        for g in idx.get(h[pos : pos + K], ()):
            if h in g or g in h:
                out.add(g)
    return out


def light_ok(v, g_h):
    return any(v in L or L in v for g in g_h for L in gold[g])


seq = dict(l.rstrip("\n").split("\t") for l in open(seq_id_map))
# POISK
rows = collections.defaultdict(set)
prow = []
for x in csv.DictReader(open(poisk_csv, newline="")):
    h, v = norm(x["vh.sequence_from_patent"]), norm(x["vl.sequence_from_patent"])
    if h and v:
        rows[h].add(v)
        prow.append((h, v))
# (B) POISK rows against gold
b = collections.Counter()
cache = {}
for h, v in set(prow):
    g = cache.setdefault(h, gold_of(h))
    if not g:
        continue
    b["POISK distinct (VH,VL) pairs with VH in gold"] += 1
    b[
        "  light consistent with gold" if light_ok(v, g) else "  light NOT consistent"
    ] += 1
for h, vs in rows.items():
    g = cache.setdefault(h, gold_of(h))
    if g and len(vs) == 1:
        b["POISK VH with exactly one VL, in gold"] += 1
        b[
            "  unique VL consistent"
            if light_ok(next(iter(vs)), g)
            else "  unique VL NOT consistent"
        ] += 1
    elif g:
        b["POISK VH with 2+ VL, in gold"] += 1
print("\n(B) POISK against SAbDab")
[print(f"   {k:<55}{v:>7}") for k, v in b.items()]
# (A) and (C) over the pair table
a = collections.defaultdict(collections.Counter)
c = collections.Counter()
for l in list(open(pairs_tsv))[1:]:
    ag, h, li, ds = l.rstrip("\n").split("\t")
    hs = seq[h]
    g = cache.setdefault(hs, gold_of(hs))
    if not g:
        continue
    d = ds
    if li:
        a[d]["heavy in gold"] += 1
        a[d][
            "  ASD light consistent"
            if light_ok(seq[li], g)
            else "  ASD light NOT consistent"
        ] += 1
    else:
        vs = rows.get(hs)
        if vs and len(vs) == 1:
            c["heavy-only pairs with heavy in gold and a unique POISK VL"] += 1
            c[
                "  POISK VL consistent with gold"
                if light_ok(next(iter(vs)), g)
                else "  POISK VL NOT consistent"
            ] += 1
print("\n(A) ASD pairs that already have a light chain, heavy in gold")
for d, cc in sorted(a.items(), key=lambda x: -x[1]["heavy in gold"]):
    n = cc["heavy in gold"]
    ok = cc["  ASD light consistent"]
    print(f"   {d:<40}heavy in gold {n:>6}   light consistent {ok:>6} ({ok / n:.1%})")
print("\n(C) heavy-only pairs filled from POISK")
[print(f"   {k:<62}{v:>7}") for k, v in c.items()]

if ab_dir:
    d = collections.Counter()

    def judge(tag, h, v):
        h, v = norm(h), norm(v)
        g = cache.setdefault(h, gold_of(h)) if h else None
        d[tag + ": pairs"] += 1
        if g:
            d[tag + ": heavy in gold"] += 1
            d[tag + ":   light consistent"] += light_ok(v, g)

    seen = set()
    for x in csv.DictReader(
        open(ab_dir + "/COVID-19/CoV-AbDab_only_sars2_filter.csv", newline="")
    ):
        k = (x["Heavy"], x["Light"])
        if k not in seen:
            seen.add(k)
            judge("covid", x["Heavy"], x["Light"])
    for x in csv.DictReader(open(ab_dir + "/BioMap/biomap.csv", newline="")):
        a, b = norm(x["antibody_seq_a"]), norm(x["antibody_seq_b"])
        k = (a, b)
        if k not in seen:
            seen.add(k)
            # which one is the heavy: the one with a CH1 / VH motif is not decidable here; try both orders and keep the consistent one
            ok1 = bool(gold_of(a)) and light_ok(b, gold_of(a))
            ok2 = bool(gold_of(b)) and light_ok(a, gold_of(b))
            d["biomap: pairs"] += 1
            d["biomap: either order in gold"] += bool(gold_of(a) or gold_of(b))
            d["biomap:   consistent (either order)"] += ok1 or ok2
    print("\n(D) AntiBinder raw pairs against SAbDab")
    [print(f"   {k:<44}{v:>7}") for k, v in d.items()]
