"""How many antigens survive each stage of the ASD pair pipeline, and which datasets they come from; plus how many antigens AbRank and POISK would add.

Stages (an antigen survives while at least one of its pairs does): ASD fasta -> a row with a heavy chain (S2) -> the antigen keep list (10-3000 aa, MSA/template)
-> a pair that binds (Kd <= 1 uM) -> a usable chain status -> not a polyprotein -> no structure (= final).
Usage: asd_antigen_funnel.py ASD_DIR SEQ_ID_MAP KEEP_LIST PAIRS_KD1UM PAIRS_BASIC ABRANK_TSV POISK_CSV
"""

import collections
import csv
import glob
import re
import sys

import pyarrow.parquet as pq

asd, seq_id_map, keep_list, pairs_kd1um, pairs_basic, abrank, poisk = sys.argv[1:8]
csv.field_size_limit(10**9)
id_of = {}
for l in open(seq_id_map):
    i, s = l.rstrip("\n").split("\t")
    id_of[s] = i
norm = lambda x: "".join((x or "").split()).upper()
ok = re.compile(r"[A-Z]+")
datasets = collections.defaultdict(set)
s2 = set()
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "chains"]).to_pylist():
        ab = dict(r["chains"]["antibody_chains"] or [])
        tg = dict(r["chains"]["target_chains"] or [])
        a = id_of.get(norm(tg.get("full")))
        if a is None:
            continue
        datasets[a].add(r["dataset"])
        h, li = norm(ab.get("heavy")), norm(ab.get("light"))
        if (
            h
            and ok.fullmatch(h)
            and ok.fullmatch(norm(tg.get("full")))
            and (not li or ok.fullmatch(li))
        ):
            s2.add(a)
keep = {l.strip() for l in open(keep_list) if l.strip()}
bind = {l.split("\t")[0] for l in list(open(pairs_kd1um))[1:]}
rows = [l.rstrip("\n").split("\t") for l in open(pairs_basic)]
h = rows[0]
rows = rows[1:]
ix = {k: i for i, k in enumerate(h)}
OKC = ("paired_original", "scfv", "vhh", "filled_covid", "filled_hiv", "filled_biomap")
chain = {r[0] for r in rows if r[ix["chain_status"]] in OKC}
poly = {
    r[0] for r in rows if r[ix["chain_status"]] in OKC and r[ix["polyprotein"]] == "0"
}
final = {r[0] for r in rows if r[ix["keep"]] == "1"}
stages = [
    ("ASD antigens (fasta)", set(datasets)),
    ("a row with a heavy chain (S2)", s2),
    ("antigen keep list", s2 & keep),
    ("a pair that binds", bind),
    ("usable chain status", chain),
    ("not a polyprotein", poly),
    ("no structure = final", final),
]
prev = None
print(f"{'stage':<34}{'antigens':>9}{'lost':>7}   datasets of the lost (antigens)")
for name, st in stages:
    lost = (prev - st) if prev is not None else set()
    c = collections.Counter(d for a in lost for d in datasets[a])
    print(f"{name:<34}{len(st):>9}{len(lost):>7}   {dict(c.most_common(6))}")
    prev = st
by = collections.defaultdict(lambda: [0, 0])
for a, ds in datasets.items():
    for d in ds:
        by[d][0] += 1
        by[d][1] += a in final
print("\nantigens per dataset: all -> final")
for d, (n, k) in sorted(by.items(), key=lambda x: -x[1][0])[:14]:
    print(f"   {d:<26}{n:>6} -> {k:>5}")
pairs_per = collections.Counter(r[0] for r in rows if r[ix["keep"]] == "1")
print(
    "\nfinal antigens by pairs: top5",
    pairs_per.most_common(5),
    "| antigens with 1 pair:",
    sum(v == 1 for v in pairs_per.values()),
)
# AbRank, POISK antigens
ab = {
    norm(x["Ag_seq"])
    for x in csv.DictReader(open(abrank, newline=""), delimiter="\t")
    if x["Ag_seq"]
}
ps = set()
for x in csv.DictReader(open(poisk, newline="")):
    t = norm(x["target.sequence_from_patent"])
    if len(t) >= 10:
        ps.add(t)
asd_seqs = {s for s, i in id_of.items() if i in datasets}
print(
    f"\ndistinct antigen sequences: AbRank {len(ab)} (also in ASD {len(ab & asd_seqs)}), POISK {len(ps)} (also in ASD {len(ps & asd_seqs)}), union new vs ASD {len((ab | ps) - asd_seqs)}"
)
