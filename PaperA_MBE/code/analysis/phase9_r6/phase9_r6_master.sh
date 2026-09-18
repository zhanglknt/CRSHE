#!/bin/bash
# Phase9 R6 master: run all BUSTED groups sequentially (4 workers, 300s/gene).
# Usage: bash phase9_r6_master.sh
BASE=/root/hyphy_work_phase9
PROJ=/mnt/d/人类正选择基因项目
SCRIPTS=$PROJ/scripts
RESULTS=$PROJ/results/phase9_hardening
HELPER="$SCRIPTS/phase9_r6_busted_one.sh"
N_WORKERS=4

run_group() {
  local GROUP="$1"
  local LIST="$2"
  if [ ! -s "$LIST" ]; then
    echo "[$GROUP] empty gene list, skip"
    return 0
  fi
  local LOG="$BASE/logs/${GROUP}.log"
  echo "===== GROUP $GROUP start $(date) ====="
  local START=$(date +%s)
  cat "$LIST" | xargs -P $N_WORKERS -I{} bash "$HELPER" "$GROUP" {} > "$LOG" 2>&1
  local END=$(date +%s)
  local OK=$(grep -c '^OK,' "$LOG" || true)
  local DONE=$(grep -c '^DONE,' "$LOG" || true)
  local SKIP=$(grep -c '^SKIP,' "$LOG" || true)
  local TO=$(grep -c '^TIMEOUT,' "$LOG" || true)
  local ERR=$(grep -c '^ERROR,' "$LOG" || true)
  echo "[$GROUP] done in $(( (END-START)/60 )) min: OK=$OK DONE=$DONE SKIP=$SKIP TIMEOUT=$TO ERROR=$ERR"
}

mkdir -p "$BASE/logs" "$RESULTS/run_logs"

run_group baseline "$BASE/genes_400.txt"
run_group srv      "$BASE/genes_400.txt"
run_group cpg      "$BASE/genes_400.txt"
run_group topo     "$BASE/genes_400.txt"
run_group sim      "$BASE/genes_sim.txt"

# copy results back to Windows
echo "===== syncing results to Windows $(date) ====="
for G in baseline srv cpg topo sim; do
  mkdir -p "$RESULTS/busted_rerun_$G"
  cp "$BASE/results_$G/"*.json "$RESULTS/busted_rerun_$G/" 2>/dev/null || true
  cp "$BASE/logs/$G.log" "$RESULTS/run_logs/$G.log" 2>/dev/null || true
done
echo "MASTER DONE $(date)"
