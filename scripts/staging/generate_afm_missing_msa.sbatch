#!/bin/bash
#SBATCH -J afm_gen_msa
#SBATCH -c 112
#SBATCH --mem=490g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -t 8:00:00
#SBATCH -o afm_gen_msa_%j.out
set -uo pipefail
# MSA for an AFDB multimer entity the release ships none for (2026-09-25: only
# AF-0000000208803231). MMseqs2 vs uniref30_2302 with the teddymer MSA parameters
# (transforms/msa_mmseqs.py): search 3 iterations -s 8 -e 0.1 --max-seqs 10000 -> expandaln ->
# align -> result2msa. The release used colabfold_search with a per-taxon filter (--filter 2);
# this does not -- the lineage note goes in the README. A whole node's memory keeps the 226 GB
# index in page cache (with 128 GB each mmseqs step re-read it and timed out).
# Output: BioMol/materials/intermediate/afdb_multimer/msa_generated/ (not raw: it is ours);
# stage_afdb_multimer_msa.py finalize compresses it to .a3m.zst beside it.
D=${DATA_ROOT:-/data/shared/cssb_data}
MM=$D/seq_db/tools/mmseqs/bin/mmseqs
DB=$D/seq_db/uniref30/2302/msadb/uniref30_2302
OUT=$D/BioMol/materials/intermediate/afdb_multimer/msa_generated
ID=${1:?usage: sbatch generate_afm_missing_msa.sbatch AF-<n>  (an entity the release ships no MSA for)}
N=${ID#AF-}
W=${TMPDIR:-/tmp}/afm_gen_$SLURM_JOB_ID; mkdir -p $W $OUT
PDB=$D/BioMol/materials/raw/afdb_multimer/homodimer/pdb/${N: -3}/$ID-model_v1.pdb
echo "=== START $(date -Is) on $(hostname)  from $PDB ==="
python3 - "$PDB" "$W/query.fasta" "$ID" <<'PY'
import sys
AA={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
pdb, out, sid = sys.argv[1:]
s = "".join(AA.get(l[17:20], "X") for l in open(pdb) if l.startswith("ATOM") and l[12:16] == " CA " and l[21] == "A")
open(out, "w").write(f">{sid}\n{s}\n"); print(f"query {sid} len {len(s)}")
PY
T="--threads $SLURM_CPUS_PER_TASK --db-load-mode 2"
$MM createdb $W/query.fasta $W/qdb >/dev/null
$MM search $W/qdb $DB $W/res $W/tmp --num-iterations 3 -s 8 -e 0.1 --max-seqs 10000 -a $T >/dev/null
$MM expandaln $W/qdb ${DB}_seq $W/res ${DB}_aln $W/exp --expansion-mode 0 -e inf $T >/dev/null
$MM align $W/qdb ${DB}_seq $W/exp $W/aln -e 10 --max-accept 10000 --alt-ali 10 -a $T >/dev/null
$MM result2msa $W/qdb ${DB}_seq $W/aln $W/msa --msa-format-mode 6 $T >/dev/null
$MM unpackdb $W/msa $W/unpacked --unpack-suffix .a3m >/dev/null
cp $W/unpacked/0.a3m $OUT/$ID-msa_v1.a3m
echo "depth $(grep -c '^>' $OUT/$ID-msa_v1.a3m)  -> $OUT/$ID-msa_v1.a3m"
cat > $OUT/README.md <<EOF2
# msa_generated — MSAs the AFDB release does not ship

The AFDB complex release distributes one MSA per monomer entity (msas/251124_1.2M,
251124_8.8M, 251209_13.4M; 23,425,744 entities). Checked against every structure in
BioMol/materials/raw/afdb_multimer/{homodimer,heterodimer} on 2026-09-25, by entity id and, where the metadata could not
map an entity, by exact sequence:

- homodimers: 1,750,754 / 1,750,755 covered
- heterodimer chains: 160,496 / 160,496 covered (2,229 only findable by sequence -- their
  monomers have no homodimer_metadata row, so UniProt could not map them)

The single uncovered entity is generated here: $ID (the only structure in the
208.8M id block; every distributed MSA is in 65M-88M).

Lineage differs from the release: MMseqs2 against uniref30_2302 with the same parameters as
the teddymer MSAs (transforms/msa_mmseqs.py), without the release's per-taxon filter
(colabfold_search --filter 2). One MSA in 1.75M.
EOF2
rm -rf $W
echo "=== DONE $(date -Is) ==="
