"""Do the antigens of the final pair table already have an MSA and templates? (They are a subset of the antigens asd/ag_* searched, so no new search is needed.)

Usage: asd_final_antigen_check.py PAIRS_BASIC WORKLIST DEPTH_TSV TEMPLATE_COUNTS_TSV KEEP_LIST
"""

import statistics
import sys

pairs, worklist, depth_tsv, tmpl_tsv, keep_list = sys.argv[1:6]
rows = [l.rstrip("\n").split("\t") for l in open(pairs)]
ix = {k: i for i, k in enumerate(rows[0])}
final = {r[0] for r in rows[1:] if r[ix["keep"]] == "1"}
work = {l.split("\t")[0] for l in open(worklist)}
keep = {l.strip() for l in open(keep_list) if l.strip()}
depth = {k: int(v) for k, v in (l.split() for l in open(depth_tsv))}
tmpl = {k: int(v) for k, v in (l.split() for l in open(tmpl_tsv))}
d = [depth[a] for a in final if a in depth]
t = [tmpl.get(a, 0) for a in final]
print(
    f"final antigens {len(final)}: on the work list {len(final & work)}, on the keep list {len(final & keep)}, with an MSA {len(d)}, in the template table {sum(a in tmpl for a in final)}"
)
print(
    f"MSA depth: median {statistics.median(d):.0f}, min {min(d)}, <= 100: {sum(x <= 100 for x in d)}, <= 50: {sum(x <= 50 for x in d)}"
)
print(
    f"templates: none {sum(x == 0 for x in t)}, 1-4 {sum(0 < x < 5 for x in t)}, 5-19 {sum(4 < x < 20 for x in t)}, 20 {sum(x == 20 for x in t)}"
)
