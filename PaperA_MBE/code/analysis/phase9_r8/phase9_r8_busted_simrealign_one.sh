#!/bin/bash
# R8 optional Task 2: BUSTED (srv No) on realigned neutral-sim alignments.
# Usage: bash phase9_r8_busted_simrealign_one.sh <GENE_ID>
GENE="$1"
BASE=/root/hyphy_work_phase9
HYPHY=/opt/miniconda3/envs/hyphy/bin/hyphy
TIMEOUT=300

if [ -z "$GENE" ]; then echo "SKIP,,empty"; exit 0; fi

mkdir -p "$BASE/results_sim_realign"

FASTA="$BASE/fasta_sim_realign/${GENE}.fasta"
TREE="$BASE/trees/${GENE}.nwk"
OUTPUT="$BASE/results_sim_realign/${GENE}.json"

# sim genes use the sim gene id scheme; tree of the same id exists per R6 setup
if [ ! -f "$FASTA" ] || [ ! -f "$TREE" ]; then
    echo "SKIP,$GENE,no_input"
    exit 0
fi

if [ -f "$OUTPUT" ] && [ "$(stat -c%s "$OUTPUT" 2>/dev/null || echo 0)" -gt 50 ]; then
    echo "DONE,$GENE"
    exit 0
fi

rm -f "$OUTPUT" 2>/dev/null

timeout -k 30 $TIMEOUT $HYPHY busted \
    --alignment "$FASTA" \
    --tree "$TREE" \
    --branches Test \
    --srv No \
    --output "$OUTPUT" \
    > /dev/null 2>&1

exit_code=$?
if [ $exit_code -eq 0 ] && [ -f "$OUTPUT" ] && [ "$(stat -c%s "$OUTPUT" 2>/dev/null || echo 0)" -gt 50 ]; then
    echo "OK,$GENE"
elif [ $exit_code -eq 124 ] || [ $exit_code -eq 137 ]; then
    echo "TIMEOUT,$GENE"
else
    echo "ERROR,$GENE,$exit_code"
fi
