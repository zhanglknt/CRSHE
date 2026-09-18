#!/bin/bash
# copy current group results to Windows-side phase9_hardening (partial sync for dry-run)
RESULTS=/mnt/d/人类正选择基因项目/results/phase9_hardening
BASE=/root/hyphy_work_phase9
for G in baseline srv cpg topo sim; do
  if [ -d "$BASE/results_$G" ] && [ "$(ls "$BASE/results_$G" 2>/dev/null | wc -l)" -gt 0 ]; then
    mkdir -p "$RESULTS/busted_rerun_$G"
    cp "$BASE/results_$G/"*.json "$RESULTS/busted_rerun_$G/" 2>/dev/null
    echo "$G: $(ls "$RESULTS/busted_rerun_$G" | wc -l) synced"
  fi
done
cp "$BASE/logs/"*.log "$RESULTS/run_logs/" 2>/dev/null
cp "$BASE/fasta_sim/_sim_manifest.csv" "$RESULTS/" 2>/dev/null
cp "$BASE/fasta_sim/_sim_report.json" "$RESULTS/" 2>/dev/null
echo synced
