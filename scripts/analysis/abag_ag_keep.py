"""The antigens the final pairs may use: the work list minus those with a shallow MSA and no template.

An antigen is dropped when MSA depth <= MAX_DEPTH and it has no template at the 2021-09-30 cutoff. Writes one seq_id a line.
Usage: abag_ag_keep.py WORKLIST DEPTH_TSV TEMPLATE_COUNTS_TSV MAX_DEPTH OUT
"""

import sys

worklist, depth_tsv, tmpl_tsv, max_depth, out = sys.argv[1:6]
ids = [l.split("\t")[0] for l in open(worklist)]
depth = {k: int(v) for k, v in (l.split() for l in open(depth_tsv))}
tmpl = {k: int(v) for k, v in (l.split() for l in open(tmpl_tsv))}
keep = sorted(
    k
    for k in ids
    if k in depth and not (depth[k] <= int(max_depth) and tmpl.get(k, 0) == 0)
)
open(out, "w").write("".join(f"{k}\n" for k in keep))
print(
    f"{len(keep)} of {len(ids)} antigens kept ({sum(k not in depth for k in ids)} without an MSA)"
)
