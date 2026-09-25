#!/bin/bash
# R8 optional Task 2 master: BUSTED on all realigned neutral-sim alignments.
BASE=/root/hyphy_work_phase9
SCRIPTS=/mnt/d/人类正选择基因项目/scripts
LOG="$BASE/logs/sim_realign_busted.log"

mkdir -p "$BASE/logs" "$BASE/results_sim_realign"

echo "[simrealign-busted] $(date '+%F %T') start" | tee -a "$LOG"
ls "$BASE"/fasta_sim_realign/*.fasta 2>/dev/null | sed 's/.*\///; s/\.fasta//' | \
    xargs -P 4 -I{} bash "$SCRIPTS/phase9_r8_busted_simrealign_one.sh" {} >> "$LOG" 2>&1
N=$(ls "$BASE"/results_sim_realign/*.json 2>/dev/null | wc -l)
echo "[simrealign-busted] $(date '+%F %T') done results=$N" | tee -a "$LOG"

mkdir -p /mnt/d/人类正选择基因项目/results/phase9_hardening/busted_rerun_sim_realign
cp -u "$BASE"/results_sim_realign/*.json /mnt/d/人类正选择基因项目/results/phase9_hardening/busted_rerun_sim_realign/ 2>/dev/null
cp "$LOG" /mnt/d/人类正选择基因项目/results/phase9_hardening/sim_realign_busted.log 2>/dev/null
echo "[simrealign-busted] synced" | tee -a "$LOG"
