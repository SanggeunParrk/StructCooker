"""Where each ASD dataset loses rows on the way to the pair table (asd/pairs/pairs_kd1um.tsv), with S1 off.

Per dataset: rows; without a heavy chain or antigen / invalid characters (S2); antigen not on the kept list (length 10..3000 and MSA
depth/template); then the binding verdicts of the rows that are left (yes / no / none). Row level, so it counts rows, not pairs.
Usage: asd_dataset_funnel.py ASD_DIR SEQ_ID_MAP ANTIGEN_IDS BINDING_PY
"""
import collections, glob, importlib.util, re, sys

import pyarrow.parquet as pq

asd, seq_id_map, antigen_ids, binding_py = sys.argv[1:5]
spec = importlib.util.spec_from_file_location("asd_binding", binding_py)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
keep = {l.strip() for l in open(antigen_ids) if l.strip()}
id_of = {s: k for k, s in (l.rstrip("\n").split("\t") for l in open(seq_id_map))}
ok = re.compile(r"[A-Z]+")
norm = lambda x: "".join((x or "").split()).upper()
c = collections.defaultdict(collections.Counter)
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "affinity", "affinity_type", "chains", "metadata"]).to_pylist():
        d = c[r["dataset"]]; d["rows"] += 1
        d["has_pdb"] += bool((r["metadata"] or {}).get("target_pdb"))
        ab = dict(r["chains"]["antibody_chains"] or []); tg = dict(r["chains"]["target_chains"] or [])
        h, l, a = norm(ab.get("heavy")), norm(ab.get("light")), norm(tg.get("full"))
        if not h or not a or not ok.fullmatch(h) or not ok.fullmatch(a) or (l and not ok.fullmatch(l)):
            d["no_chain"] += 1; continue
        if id_of.get(a) not in keep:
            d["antigen_out"] += 1; continue
        v = mod.binding_verdict(r["affinity_type"], r["affinity"], r["dataset"])
        d["yes" if v else ("no" if v is False else "none")] += 1
print(f"{'dataset':<26}{'rows':>9}{'has_pdb':>9}{'no_chain':>9}{'ag_out':>8}{'yes':>9}{'no':>9}{'none':>9}")
for k, d in sorted(c.items(), key=lambda x: -x[1]["rows"]):
    print(f"{k:<26}{d['rows']:>9}{d['has_pdb']:>9}{d['no_chain']:>9}{d['antigen_out']:>8}{d['yes']:>9}{d['no']:>9}{d['none']:>9}")
