#!/bin/bash
# Runs groups of the antibody-antigen recipe in order, one `structcooker build abag/<stage>` at a time (each stage is an sbatch job; the script waits for it).
#   docs/abag_distillation/run.sh ids        antigen_fasta, antibody_fasta -> seq_id_map                       (ASD + AbRank + POISK)
#   docs/abag_distillation/run.sh pairs      pairs_binding -> pairs_chains -> pairs_basic                      (about 10 minutes)
#   docs/abag_distillation/run.sh antigen [from]   ag_worklist -> MSA search -> templates (2021-09-30 and 2026) -> depth / template counts -> ag_keep.txt   (many hours)
#   docs/abag_distillation/run.sh final      the kept pairs' own fastas, seq_ids, clusters -> abag/final/pairs_final.tsv   (about 40 minutes: ANARCI twice, in parallel)
#   docs/abag_distillation/run.sh stage <name>   one stage, e.g. pairs_basic
#   docs/abag_distillation/run.sh all        ids, pairs, antigen, final
# Raw data first (RECIPE.md, step 0).
set -euo pipefail
cd "$(dirname "$0")/../.."
A=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}

build() {  # build <stage>: submit abag/<stage>, print its job id
    local out
    out=$(pixi run structcooker build "abag/$1" 2>&1 | grep -v WARN)
    echo "$out" >&2
    echo "$out" | grep -o 'job [0-9]*' | tail -1 | cut -d' ' -f2
}
wait_for() {  # wait_for <job...>: block until every task of every job ended (an array has one line per task), fail unless all are COMPLETED
    local j states tries
    for j in "$@"; do
        tries=0
        while [ -z "$(sacct -j "$j" -n -X -o State)" ] && [ $tries -lt 12 ]; do sleep 10; tries=$((tries + 1)); done  # a job just submitted is not in sacct for a moment
        while sacct -j "$j" -n -X -o State | grep -q "PENDING\|RUNNING\|REQUEUED\|COMPLETING\|CONFIGURING\|SUSPENDED"; do sleep 20; done
        states=$(sacct -j "$j" -n -X -o State | tr -d ' ' | sort -u | tr '\n' ' ')
        [ "$states" = "COMPLETED " ] || { echo "job $j ended: $states" >&2; exit 1; }
    done
}
run() { local j; j=$(build "$1"); echo "abag/$1 -> job $j"; wait_for "$j"; }
both() { local a b; a=$(build "$1"); b=$(build "$2"); wait_for "$a" "$b"; }

ids() {
    both antigen_fasta antibody_fasta
    run seq_id_map
}
pairs() {
    run pairs_binding
    run pairs_chains
    run pairs_basic
    echo "done: $A/abag/pairs/pairs_basic.tsv ($(($(wc -l < "$A/abag/pairs/pairs_basic.tsv") - 1)) pairs)"
}
antigen() {  # antigen [from-stage]: the stages in order; a stage before from-stage is skipped (to resume)
    local from=${1:-ag_worklist} on=0 st
    for st in ag_worklist ag_msa_search ag_msa_and_wo_lower ag_hmmsearch ag_template_candidates ag_template ag_keep; do
        [ "$st" = "$from" ] && on=1
        [ $on = 1 ] || continue
        case $st in
            ag_worklist|ag_msa_search|ag_hmmsearch) run "$st" ;;
            ag_msa_and_wo_lower) both ag_msa ag_msa_wo_lower ;;
            ag_template_candidates) both ag_template_candidates_20210930 ag_template_candidates_2026 ;;
            ag_template) both ag_template_20210930 ag_template_2026 ;;
            ag_keep) ag_keep ;;
        esac
    done
}
ag_keep() {
    local W=$A/abag/worklists j
    j=$(sbatch --parsable -p cpu-short --qos=cpu-short-q -c 2 --mem=16g -t 01:00:00 -J abag_agkeep --wrap "
        set -e
        python3 -I scripts/analysis/asd_msa_depth.py $A/abag/msa/ag $W/ag_msa_depth.tsv
        pixi run python scripts/analysis/asd_template_counts.py $A/abag/lmdb/template/template_ag_20210930.lmdb $W/ag_template_counts.tsv
        python3 -I scripts/analysis/abag_ag_keep.py $W/seq_id_map_AG.tsv $W/ag_msa_depth.tsv $W/ag_template_counts.tsv 50 $W/ag_keep.txt")
    wait_for "$j"
}
final() {
    both final_antigen_fasta final_antibody_fasta
    run final_seq_id_map
    local c d e
    c=$(build final_antigen_seq_cluster40); d=$(build final_antibody_seq_cluster_cdr90); e=$(build final_antibody_seq_cluster_cdr100)
    wait_for "$c" "$d" "$e"
    run final_pairs
    echo "done: $A/abag/final/pairs_final.tsv ($(($(wc -l < "$A/abag/final/pairs_final.tsv") - 1)) pairs)"
}

case "${1:-}" in
    ids) ids ;;
    pairs) pairs ;;
    antigen) antigen "${2:-ag_worklist}" ;;
    final) final ;;
    stage) run "${2:?usage: run.sh stage <name>}" ;;
    all) ids; pairs; antigen; final ;;
    *) sed -n 2,9p "$0"; exit 1 ;;
esac
