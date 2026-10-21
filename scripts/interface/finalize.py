"""Manifest + cluster annotation + statistics report for a finished interface DB (sbatch only).

usage: finalize.py <db_root>
Reads  <db_root>/manifest.tsv, done.tsv, failed.tsv
Writes <db_root>/manifest_full.tsv   manifest + node_seq_ids, node_c30, node_c40, cluster keys, method
       <db_root>/report.md           counts and distributions
Protein node hashes (md5-16 of "P:<sequence>") are mapped back to sequences through the PDB protein
FASTA, then to seq_ids by one streaming pass over metadata/seq_id_map.tsv (7.9 GB, all sources), then to
the 30% / 40% sequence clusters of metadata/seq_cluster30_PDB.tsv and seq_cluster40.tsv.
"""
import csv
import hashlib
import os
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
db = Path(sys.argv[1])
t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)


rows = list(csv.DictReader(open(db / "manifest.tsv"), delimiter="\t"))
log(f"{len(rows)} units")

# ---- protein node hashes -> sequences -> seq_ids -> clusters
need = set()
for r in rows:
    for h, k in zip(r["node_hashes"].split("|"), r["node_kinds"].split("|"), strict=True):
        if k == "protein":
            need.add(h)
log(f"{len(need)} distinct protein node hashes")
hash2seq = {}
with open(ROOT / "materials/intermediate/fasta/pdb_polypeptide_L.fasta") as f:
    for line in f:
        if line[0] != ">":
            s = line.strip()
            h = hashlib.md5(("P:" + s).encode()).hexdigest()[:16]  # noqa: S324
            if h in need:
                hash2seq[h] = s
log(f"sequence recovered for {len(hash2seq)} of {len(need)} hashes")
seqs = set(hash2seq.values())
seq2id = {}
with open(ROOT / "metadata/seq_id_map.tsv") as f:
    for line in f:
        if line[0] != "P":
            continue
        i, _, s = line.rstrip("\n").partition("\t")
        if s in seqs:
            seq2id[s] = i
log(f"seq_id found for {len(seq2id)} of {len(seqs)} sequences")
ids = set(seq2id.values())
id2c30, id2c40 = {}, {}
for fn, dst in (("seq_cluster30_PDB.tsv", id2c30), ("seq_cluster40.tsv", id2c40)):
    with open(ROOT / "metadata" / fn) as f:
        for line in f:
            c, _, members = line.rstrip("\n").partition("\t")
            for m in members.split(","):
                if m in ids:
                    dst[m] = c
log(f"clusters: c30 {len(id2c30)}, c40 {len(id2c40)} of {len(ids)} ids")
h2id = {h: seq2id.get(s) for h, s in hash2seq.items()}

meta = {}
with open(ROOT / "metadata/cif_metadata.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta.setdefault(r["cif_id"].split("_")[0], r)

cols = list(rows[0]) + ["node_seq_ids", "node_c30", "node_c40", "cluster_key30", "cluster_key40", "n_unmapped_protein",
                        "including_NA_entry"]
miss = Counter()
with open(db / "manifest_full.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        sid, c30, c40, unm = [], [], [], 0
        for h, k in zip(r["node_hashes"].split("|"), r["node_kinds"].split("|"), strict=True):
            if k != "protein":
                sid.append(""); c30.append(""); c40.append("")
                continue
            i = h2id.get(h)
            if i is None:
                unm += 1
                miss["unmapped_protein_nodes"] += 1
            sid.append(i or ""); c30.append(id2c30.get(i, "") if i else ""); c40.append(id2c40.get(i, "") if i else "")
        r["node_seq_ids"], r["node_c30"], r["node_c40"] = "|".join(sid), "|".join(c30), "|".join(c40)
        r["cluster_key30"] = "+".join(sorted(set(x for x in c30 if x)))
        r["cluster_key40"] = "+".join(sorted(set(x for x in c40 if x)))
        r["n_unmapped_protein"] = unm
        r["including_NA_entry"] = meta.get(r["pdb_id"], {}).get("including_NA", "")
        fo.write("\t".join(str(r[c]) for c in cols) + "\n")
log("manifest_full.tsv written")

# ---- report
done = [l.rstrip("\n").split("\t") for l in open(db / "done.tsv")]
failed = [l.rstrip("\n").split("\t") for l in open(db / "failed.tsv")] if (db / "failed.tsv").exists() else []
n_entries, n_zero, n_big = len(done), sum(1 for d in done if int(d[1]) == 0), sum(int(d[2]) for d in done)


def dist(vals):
    vals = sorted(vals)
    if not vals:
        return "-"
    q = lambda p: vals[min(len(vals) - 1, int(p * len(vals)))]  # noqa: E731
    return f"min {vals[0]:g} · p25 {q(.25):g} · median {q(.5):g} · p75 {q(.75):g} · p95 {q(.95):g} · max {vals[-1]:g}"


out = ["# Interface-unit DB report", "", f"Units: **{len(rows):,}** from **{n_entries:,}** PDB entries "
       f"({n_zero:,} entries produced no unit, {len(failed):,} failed, {n_big:,} assemblies skipped for > "
       f"max_polymer_chains).", ""]
out += ["## Units by type", "", "| type | units | share |", "|---|---:|---:|"]
tc = Counter(r["type"] for r in rows)
for t, n in tc.most_common(25):
    out.append(f"| {t} | {n:,} | {100 * n / len(rows):.1f}% |")
out += ["", "## Class totals", ""]
for lab, pred in (("protein-protein only (no NA / ligand / ion)", lambda r: r["n_na"] == "0" and r["n_ligand"] == "0" and r["n_ion"] == "0"),
                  ("with NA", lambda r: r["n_na"] != "0"), ("with ligand", lambda r: r["n_ligand"] != "0"),
                  ("with ion", lambda r: r["n_ion"] != "0"), ("with antibody / TCR", lambda r: r["has_ab"] == "1"),
                  ("with a base-paired duplex", lambda r: r["n_duplex_pairs"] != "0"),
                  ("truncated (member cap hit)", lambda r: r["truncated"] == "True")):
    out.append(f"- {lab}: {sum(1 for r in rows if pred(r)):,}")
out += ["", "## Size distributions", ""]
for c, lab in (("n_nodes", "nodes per unit"), ("n_atoms", "atoms per unit"), ("n_residues", "residues per unit"),
               ("n_interface_res", "interface residues"), ("n_atom_pairs", "contacting atom pairs"),
               ("n_water", "waters per unit"), ("n_bridging_water", "bridging waters")):
    out.append(f"- {lab}: {dist([float(r[c]) for r in rows])}")
res = [float(r["resolution"]) for r in rows if r["resolution"] not in ("", "None", "nan")]
out += [f"- resolution (A): {dist(res)}  ({len(rows) - len(res):,} units without resolution: NMR / EM without value)", ""]
out += ["## Water", "", f"- units with no water: {sum(1 for r in rows if r['n_water'] == '0'):,}",
        f"- units with >= 1 bridging water: {sum(1 for r in rows if r['n_bridging_water'] != '0'):,}", ""]
out += ["## Cluster coverage", "", f"- protein nodes without a seq_id: {miss['unmapped_protein_nodes']:,}",
        f"- units with every protein node mapped: {sum(1 for r in rows if str(r['n_unmapped_protein']) == '0'):,}",
        f"- distinct c30 cluster keys: {len(set(r['cluster_key30'] for r in rows)):,}", ""]
ccd, ions = Counter(), Counter()
for r in rows:
    for k, c in zip(r["node_kinds"].split("|"), r["node_ccd"].split("|"), strict=True):
        if k == "ligand":
            ccd[c] += 1
        elif k == "ion":
            ions[c] += 1
out += ["## Most frequent ligands / ions (nodes in units)", "", "ligands: " + ", ".join(f"{c} {n:,}" for c, n in ccd.most_common(30)),
        "", "ions: " + ", ".join(f"{c} {n:,}" for c, n in ions.most_common(15)), ""]
if failed:
    out += ["## Failed entries", "", *[f"- {k}: {m}" for k, m in failed[:30]], ""]
(db / "report.md").write_text("\n".join(out))
log("report.md written")

# ---- publish the shards as a DataCooker shard collection (shards.json) + .meta.json, which is what
# release-docs and the other tools of the release tree expect of a database
import json  # noqa: E402

from datacooker.lmdb.sharded import publish_shards  # noqa: E402

coll = db / "interface_pdb.lmdb"
if (db / "shards").is_dir() and not (coll / "shards.json").exists():
    coll.mkdir(exist_ok=True)
    (db / "shards").rename(coll / "shards")
shard_paths = sorted((coll / "shards").glob("shard_*.lmdb"))
n_published = publish_shards(shard_paths, coll, overwrite=True)  # also checks that every key is unique
if n_published != len(rows):
    raise SystemExit(f"{n_published} records in the shards but {len(rows)} units in the manifest")
(db / "interface_pdb.lmdb.meta.json").write_text(json.dumps({
    "schema": "interface-unit v1", "key_convention": "<pdb>_<assembly>_<model>_<unit index>",
    "codec": "biomol_zstd", "entries": n_published, "source": "lmdb/pdb/cif/cif_pdb.lmdb",
    "recipe": "scripts/interface/build_db.py", "extra": {"pdb_entries": n_entries, "shards": len(shard_paths)},
}, indent=2))
log(f"published {n_published} records in {len(shard_paths)} shards")
