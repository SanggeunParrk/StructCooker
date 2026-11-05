"""Three ways to say "this pair has an experimental structure", counted on the pair table, and ASD's own target_pdb next to them.

  current   the rule of pairs_basic: the antibody contained in a PDB chain (constant domain / tag allowed), the antigen a fragment of / containing / 80% of 10-mers
  exact-ab  heavy and light sequence IDENTICAL to a SAbDab antibody's chains (same seq_id), the antigen as in current
  exact-all heavy, light and antigen sequence all IDENTICAL to the chains of one SAbDab entry (same seq_id for every chain)
  target_pdb  ASD's metadata.target_pdb is set on a row of the pair
Run with `pixi run python`. Usage: asd_structure_rules.py ASD_DIR PAIRS_BASIC SEQ_ID_MAP SABDAB_SUMMARY PDB_FASTA
"""

import collections
import glob
import sys

import pyarrow.parquet as pq

from structcooker.instructions.readers.asd_basic import PdbComplexIndex

asd, pairs_basic, seq_id_map, sab, fasta = sys.argv[1:6]
seq = dict(l.rstrip("\n").split("\t") for l in open(seq_id_map))
norm = lambda x: "".join((x or "").split()).upper()
idx = PdbComplexIndex(sab, fasta)
exact = collections.defaultdict(list)  # (heavy, light) -> entries
for pdb, heavy, light, antigens in idx.entries:
    exact[(heavy, light)].append((pdb, antigens))
rows = [l.rstrip("\n").split("\t") for l in open(pairs_basic)]
h = rows[0]
rows = rows[1:]
ix = {k: i for i, k in enumerate(h)}
pairs = []
for r in rows:
    heavy, antigen = seq[r[ix["heavy_seq_id"]]], seq[r[ix["antigen_seq_id"]]]
    light = (
        seq[r[ix["light_seq_id"]]] if r[ix["light_seq_id"]] else r[ix["light_sequence"]]
    )
    pairs.append((heavy, light, antigen, r))
wanted = {(a, b, c) for a, b, c, _ in pairs}
tp = set()
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for row in pq.read_table(f, columns=["chains", "metadata"]).to_pylist():
        if not (row["metadata"] or {}).get("target_pdb"):
            continue
        ab = dict(row["chains"]["antibody_chains"] or [])
        tg = dict(row["chains"]["target_chains"] or [])
        k = (norm(ab.get("heavy")), norm(ab.get("light")), norm(tg.get("full")))
        if k in wanted:
            tp.add(k)
cnt = collections.Counter()
byds = collections.defaultdict(collections.Counter)
for heavy, light, antigen, r in pairs:
    cur = bool(r[ix["pdb_structure"]])
    ent = exact.get((heavy, light), [])
    eab = any(any(idx._antigen_matches(antigen, a) for a in ants) for _, ants in ent)
    eall = any(antigen in ants for _, ants in ent)
    t = (heavy, light, antigen) in tp
    cnt["pairs"] += 1
    cnt["current"] += cur
    cnt["exact-ab"] += eab
    cnt["exact-all"] += eall
    cnt["target_pdb"] += t
    cnt[f"current={int(cur)} target_pdb={int(t)}"] += 1
    cnt["exact-all but not current"] += eall and not cur
    cnt["exact-ab but not current"] += eab and not cur
    cnt["current but not exact-ab"] += cur and not eab
    for d in r[ix["datasets"]].split(","):
        byds[d]["pairs"] += 1
        byds[d]["current"] += cur
        byds[d]["target_pdb"] += t
        byds[d]["exact-all"] += eall
for k, v in cnt.items():
    print(f"{k:<40}{v:>8}")
print("\nby dataset (pairs, current, exact-all, target_pdb)")
for d, c in sorted(byds.items(), key=lambda x: -x[1]["pairs"])[:12]:
    print(
        f"   {d:<26}{c['pairs']:>8}{c['current']:>8}{c['exact-all']:>8}{c['target_pdb']:>8}"
    )
