#!/bin/bash
# R8 optional Task 2 step 1: codon-aware realignment of the 400 neutral-sim
# (omega=1) alignments, same pipeline as the stratified-400 realign.
BASE=/root/hyphy_work_phase9
SCRIPTS=/mnt/d/人类正选择基因项目/scripts
PY=/opt/miniconda3/envs/hyphy/bin/python
LOG="$BASE/logs/sim_realign.log"

mkdir -p "$BASE/fasta_sim_realign" "$BASE/sim_realign_stats" "$BASE/realign_tmp_s" "$BASE/logs"

echo "[sim-realign] $(date '+%F %T') start" | tee -a "$LOG"
# NOTE: input is fasta_sim2 (R7 Task 2 parametric-bootstrap null set, per-gene MLE
# params) - NOT fasta_sim (old R6 fixed-param set); baseline for the calibration
# delta is busted_rerun_sim2, so the same 400 sim alignments must be realigned.
ls "$BASE"/fasta_sim2/*.fasta | sed 's/.*\///; s/\.fasta//' | xargs -P 4 -I{} \
    $PY "$SCRIPTS/phase9_r7_realign_one.py" {} \
    "$BASE/fasta_sim2" "$BASE/fasta_sim_realign" "$BASE/sim_realign_stats" "$BASE/realign_tmp_s" \
    >> "$LOG" 2>&1
N=$(ls "$BASE"/fasta_sim_realign/*.fasta 2>/dev/null | wc -l)
echo "[sim-realign] $(date '+%F %T') done fasta=$N" | tee -a "$LOG"

mkdir -p /mnt/d/人类正选择基因项目/results/phase9_hardening/sim_realign_stats
cp -u "$BASE"/sim_realign_stats/*.json /mnt/d/人类正选择基因项目/results/phase9_hardening/sim_realign_stats/ 2>/dev/null
echo "[sim-realign] synced" | tee -a "$LOG"
