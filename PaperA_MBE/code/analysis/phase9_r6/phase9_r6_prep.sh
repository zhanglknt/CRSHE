#!/bin/bash
# Phase9 R6 prep: build WSL-local work dirs for the 400-gene subset + simulated data.
# Usage: bash phase9_r6_prep.sh
set -e
BASE=/root/hyphy_work_phase9
PROJ=/mnt/d/人类正选择基因项目
SCRIPTS=$PROJ/scripts
RESULTS=$PROJ/results/phase9_hardening
PY3=/usr/bin/python3
EVOLVER=/opt/miniconda3/envs/hyphy/bin/evolver

mkdir -p "$BASE"/{fasta,trees,fasta_cpg,trees_alt,fasta_sim,evolver_tmp} \
         "$BASE"/{results_baseline,results_srv,results_cpg,results_topo,results_sim} \
         "$BASE/logs" "$RESULTS/run_logs"

echo "[1/6] copy subset fasta + trees (local)"
awk -F, 'NR>1 {print $1}' "$RESULTS/stratified_400_genes.csv" > "$BASE/genes_400_raw.txt"
awk -F, 'NR>1 {print $1}' "$RESULTS/sim_400_genes.csv" > "$BASE/genes_sim_raw.txt"
cat "$BASE/genes_400_raw.txt" "$BASE/genes_sim_raw.txt" | sort -u | while read g; do
  cp "$PROJ/results/phase4_hyphy/fasta_v3_cleaned/$g.fasta" "$BASE/fasta/" 2>/dev/null || true
  cp "$PROJ/results/phase4_hyphy/trees_v3/$g.nwk" "$BASE/trees/" 2>/dev/null || true
done
N_FASTA=$(ls "$BASE/fasta" | wc -l)
N_TREE=$(ls "$BASE/trees" | wc -l)
echo "  copied: fasta=$N_FASTA trees=$N_TREE"

while read g; do
  [ -f "$BASE/fasta/$g.fasta" ] && [ -f "$BASE/trees/$g.nwk" ] && echo "$g"
done < "$BASE/genes_400_raw.txt" > "$BASE/genes_400.txt"
echo "  stratified subset with inputs: $(wc -l < "$BASE/genes_400.txt")"

echo "[2/6] CpG masking"
$PY3 "$SCRIPTS/phase9_r6_cpg_mask.py" "$RESULTS/stratified_400_genes.csv" \
    "$BASE/fasta" "$BASE/fasta_cpg" "$RESULTS/cpg_masking_stats.csv" \
    > "$RESULTS/run_logs/cpg_mask.log" 2>&1
tail -20 "$RESULTS/run_logs/cpg_mask.log"

echo "[3/6] alternative topology trees"
$PY3 "$SCRIPTS/phase9_r6_alt_trees.py" "$RESULTS/stratified_400_genes.csv" \
    "$BASE/trees" "$BASE/trees_alt" "$RESULTS/alt_trees_report.json" \
    > "$RESULTS/run_logs/alt_trees.log" 2>&1
tail -20 "$RESULTS/run_logs/alt_trees.log"

echo "[4/6] evolver neutral simulations"
$PY3 "$SCRIPTS/phase9_r6_evolver.py" "$RESULTS/sim_400_genes.csv" \
    "$BASE/fasta" "$BASE/trees" "$RESULTS/f3x4_codon_freqs.json" \
    "$BASE/fasta_sim" "$BASE/evolver_tmp" "$EVOLVER" \
    > "$RESULTS/run_logs/evolver.log" 2>&1
tail -25 "$RESULTS/run_logs/evolver.log"

echo "[5/6] sim gene list"
while read g; do
  [ -f "$BASE/fasta_sim/$g.fasta" ] && [ -f "$BASE/trees/$g.nwk" ] && echo "$g"
done < "$BASE/genes_sim_raw.txt" > "$BASE/genes_sim.txt"
echo "  sim subset ready: $(wc -l < "$BASE/genes_sim.txt")"

echo "[6/6] quick BUSTED smoke test on smallest gene"
SMALL=$(ls -S "$BASE/fasta" | tail -1 | sed 's/.fasta//')
timeout 120 /opt/miniconda3/envs/hyphy/bin/hyphy busted \
    --alignment "$BASE/fasta/$SMALL.fasta" --tree "$BASE/trees/$SMALL.nwk" \
    --branches Test --srv No --output "$BASE/logs/smoke_$SMALL.json" > "$BASE/logs/smoke.log" 2>&1
RC=$?
echo "  smoke gene $SMALL rc=$RC size=$(stat -c%s "$BASE/logs/smoke_$SMALL.json" 2>/dev/null || echo 0)"

echo "PREP DONE $(date)"
