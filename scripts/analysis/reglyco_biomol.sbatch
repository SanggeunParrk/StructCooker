#!/bin/bash
# Convert one shipped Re-Glyco shard to BioMol records, then measure the BioMol delta codec. Submit from the repo root:
#   [LIMIT=8] [PHASE=convert|evaluate|both] [LEVEL=9] sbatch scripts/analysis/reglyco_biomol.sbatch [SHARD]
#SBATCH --job-name=reglyco_biomol
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_biomol_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
R=$BIOMOL_ROOT/materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco
SHARD=${1:-$R/blob/data/rs_00005.bin}
WORK=$BIOMOL_ROOT/materials/intermediate/reglyco_biomol
EX=$BIOMOL_ROOT/materials/intermediate/reglyco_compact/hub/examples.parquet
LIMIT=${LIMIT:-0}
DIR=$WORK/$(basename "$SHARD" .bin)${LIMIT:+_n$LIMIT}
[ "$LIMIT" = 0 ] && DIR=$WORK/$(basename "$SHARD" .bin)
mkdir -p "$DIR"
if [ "${PHASE:-both}" != evaluate ]; then
  pixi run python scripts/analysis/reglyco_biomol.py convert "$SHARD" "$R/blob/index.parquet" --examples "$EX" --out "$DIR" --workers 16 --limit "$LIMIT"
fi
if [ "${PHASE:-both}" != convert ]; then
  SHIPPED=$(pixi run python -c "
import pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa, sys, os
from pathlib import Path
d=Path('$DIR'); keys={f.stem for f in d.glob('*.bin')}
t=pq.read_table('$R/blob/index.parquet', columns=['path','length'])
import re
k=pc.extract_regex(t['path'].combine_chunks(), r'(?P<k>[^/]+)\.cif\$').field('k')
m=pc.is_in(k, value_set=pa.array(sorted(keys)))
print(pc.sum(t['length'].filter(m)).as_py())
")
  pixi run python scripts/analysis/reglyco_biomol.py evaluate "$DIR" --out "$DIR/evaluate_L${LEVEL:-9}.json" --level "${LEVEL:-9}" --shipped-bytes "$SHIPPED"
fi
