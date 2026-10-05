#!/bin/bash
# raw = exactly what was downloaded; anything we generated or modified -> intermediate.
# One-off reorganisation, 2026-09-30. Moves are renames on one filesystem; the six repaired
# heterodimer mmCIFs are moved out and the release originals copied back (sizes checked
# against SOURCE.tsv). Nothing is deleted.
set -euo pipefail
BM=/data/shared/cssb_data/BioMol/materials
BC=/data/shared/cssb_data/BioMol_clean/materials
AFM_I=$BM/intermediate/afdb_multimer
OFD_I=$BM/intermediate/openfold_distillation
mvx() {  # mvx <src> <dst>: refuse to overwrite, create the parent
  [ -e "$1" ] || { echo "skip (absent) $1"; return 0; }
  [ -e "$2" ] && { echo "REFUSE: $2 exists"; exit 1; }
  mkdir -p "$(dirname "$2")"; mv -n "$1" "$2"; echo "mv $1 -> $2"
}

# --- BioMol/materials/raw/fasta: extracted from our DBs
mvx $BM/raw/fasta $BM/intermediate/fasta

# --- rna_alignment_arrays was downloaded (s3://openfold3-data); keep its download log with it
mvx $BM/raw/logs/rna_align_dl.log $BM/raw/rna_alignment_arrays/DOWNLOAD.log
rmdir $BM/raw/logs

# --- AFDB multimer: the one MSA we generated, the chain->MSA map, the six repairs
mvx $BM/raw/afdb_multimer/msa_generated $AFM_I/msa_generated
mvx $BM/raw/afdb_multimer/msa/231/AF-0000000208803231-msa_v1.a3m.zst \
    $AFM_I/msa_generated/AF-0000000208803231-msa_v1.a3m.zst
mvx $BM/raw/afdb_multimer/heterodimer/chain_msa.tsv $AFM_I/heterodimer/chain_msa.tsv
H=$BM/raw/afdb_multimer/heterodimer
R=$AFM_I/heterodimer/repaired
mvx $H/ENTITY_POLY_REPAIRS.tsv $R/ENTITY_POLY_REPAIRS.tsv
echo "{" > $R/overrides.json.tmp
first=1
for f in $(tail -n +2 $R/ENTITY_POLY_REPAIRS.tsv | cut -f1 | sort -u); do
  n=${f#AF-}; n=${n%%-*}; sub=${n: -3}
  mvx $H/cif/$sub/$f $R/cif/$sub/$f
  orig=$(awk -F'\t' -v e="AF-$n" '$1==e{print $3}' $H/SOURCE.tsv)
  want=$(awk -F'\t' -v e="AF-$n" '$1==e{print $4}' $H/SOURCE.tsv)
  cp -n "$orig" $H/cif/$sub/$f
  got=$(stat -c %s $H/cif/$sub/$f)
  [ "$got" = "$want" ] || { echo "size mismatch restoring $f: $got vs $want"; exit 1; }
  gzip -t $H/cif/$sub/$f
  [ $first = 1 ] || echo "," >> $R/overrides.json.tmp; first=0
  printf '  "%s": "%s"' "$H/cif/$sub/$f" "$R/cif/$sub/$f" >> $R/overrides.json.tmp
  echo "restored original $f ($got bytes)"
done
printf '\n}\n' >> $R/overrides.json.tmp; mv $R/overrides.json.tmp $R/overrides.json

# --- OpenFold distillation: our fastas, filelists/seqid tables, old LMDBs, tools and logs
for d in fasta lmdb bin _archive; do mvx $BM/raw/openfold_distillation/$d $OFD_I/$d; done
for f in $BM/raw/openfold_distillation/disordered_set/*.tsv $BM/raw/openfold_distillation/disordered_set/*.txt; do
  mvx "$f" $OFD_I/disordered_set/$(basename "$f")
done

# --- BioMol_clean: the per-DB fastas are derived, not raw
mvx $BC/raw/fasta $BC/intermediate/fasta
# TEMPORARY: afm_hmmsearch (job 275298) was started with the old fasta_path and reopens it per
# query. Remove this link (and the then-empty raw/) once that job has finished.
ln -s ../intermediate/fasta $BC/raw/fasta
echo "done $(date -Is)"
