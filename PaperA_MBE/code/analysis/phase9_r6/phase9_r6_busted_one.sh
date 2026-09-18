#!/bin/bash
# Run BUSTED for a single gene in a given phase9 group.
# Usage: bash phase9_r6_busted_one.sh <GROUP> <GENE_ID>
GROUP="$1"
GENE="$2"
BASE=/root/hyphy_work_phase9
HYPHY=/opt/miniconda3/envs/hyphy/bin/hyphy
TIMEOUT=300

if [ -z "$GENE" ]; then echo "SKIP,,empty"; exit 0; fi

case "$GROUP" in
  baseline|srv) FASTA_DIR="$BASE/fasta"; TREE_DIR="$BASE/trees" ;;
  cpg)          FASTA_DIR="$BASE/fasta_cpg"; TREE_DIR="$BASE/trees" ;;
  topo)         FASTA_DIR="$BASE/fasta"; TREE_DIR="$BASE/trees_alt" ;;
  sim)          FASTA_DIR="$BASE/fasta_sim"; TREE_DIR="$BASE/trees" ;;
  *) echo "ERROR,$GENE,bad_group"; exit 1 ;;
esac

SRV=No
if [ "$GROUP" = "srv" ]; then SRV=Yes; fi

FASTA="$FASTA_DIR/${GENE}.fasta"
TREE="$TREE_DIR/${GENE}.nwk"
OUTPUT="$BASE/results_${GROUP}/${GENE}.json"

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
    --srv $SRV \
    --output "$OUTPUT" \
    > /dev/null 2>&1

exit_code=$?
if [ $exit_code -eq 0 ] && [ -f "$OUTPUT" ] && [ "$(stat -c%s "$OUTPUT" 2>/dev/null || echo 0)" -gt 50 ]; then
    echo "OK,$GENE"
elif [ $exit_code -eq 124 ]; then
    echo "TIMEOUT,$GENE"
else
    echo "ERROR,$GENE,$exit_code"
fi
