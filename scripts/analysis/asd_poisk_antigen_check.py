"""Does requiring the antigen to match make the POISK light-chain recovery safer?

For every pair of the ASD pair table whose heavy chain is an exact VH in POISK: the light chains POISK pairs with that VH, once over all POISK rows
(what the SelfDistillation recovery does) and once over the POISK rows whose target is the pair's antigen (the POISK target sequence contains the ASD
antigen or is contained in it, >= 30 aa). Pairs that already have an ASD light chain tell how often the unique candidate is the same chain;
pairs without one tell how the candidates change.
Usage: asd_poisk_antigen_check.py PAIRS SEQ_ID_MAP POISK_CSV
"""
import collections, csv, sys

pairs_tsv, seq_id_map, poisk_csv = sys.argv[1:4]
csv.field_size_limit(10**9)
seq = dict(l.rstrip("\n").split("\t") for l in open(seq_id_map))
norm = lambda x: "".join((x or "").split()).upper()
rows = collections.defaultdict(list)
for x in csv.DictReader(open(poisk_csv, newline="")):
    h, v = norm(x["vh.sequence_from_patent"]), norm(x["vl.sequence_from_patent"])
    t = norm(x["target.sequence_from_patent"])
    if h and v: rows[h].append((v, t, x["audit.verdict"]))
def same(a, t):
    return len(t) >= 30 and (t in a or a in t)
res = collections.defaultdict(collections.Counter)
for l in list(open(pairs_tsv))[1:]:
    a, h, li, ds = l.rstrip("\n").split("\t")
    c = rows.get(seq[h])
    if not c: continue
    allc = {v for v, _, _ in c}
    ant = {v for v, t, _ in c if same(seq[a], t)}
    key = "ASD has light" if li else "heavy only"
    r = res[key]; r["pairs in POISK"] += 1
    r["all: unique" if len(allc) == 1 else "all: 2+"] += 1
    r["antigen-matched rows: none" if not ant else ("antigen-matched: unique" if len(ant) == 1 else "antigen-matched: 2+")] += 1
    if li:
        if len(allc) == 1: r["all-unique agrees" if next(iter(allc)) == seq[li] else "all-unique DISAGREES"] += 1
        if len(ant) == 1: r["antigen-unique agrees" if next(iter(ant)) == seq[li] else "antigen-unique DISAGREES"] += 1
        if len(ant) >= 1: r["ASD light is among antigen-matched" if seq[li] in ant else "ASD light NOT among antigen-matched"] += 1
    else:
        if len(allc) == 1: r["all-unique, antigen-matched rows agree"] += (next(iter(allc)) in ant); r["all-unique, antigen-matched rows exist but light differs"] += (bool(ant) and next(iter(allc)) not in ant); r["all-unique, no antigen-matched row"] += (not ant)
        if len(allc) > 1 and len(ant) == 1: r["ambiguous -> resolved by antigen"] += 1
for k, c in res.items():
    print("\n", k)
    for kk, v in sorted(c.items()): print(f"   {kk:<62}{v:>7}")
