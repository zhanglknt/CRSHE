#!/bin/bash
HYPHY=/home/linux/miniforge3/envs/hyphy/bin/hyphy
FASTA_DIR=/home/linux/hyphy_work_v3/fasta
TREE_DIR=/home/linux/hyphy_work_v3/trees
GENE_LIST=/home/linux/relax_work_v3/remaining_genes_v2.txt
LOG_FILE=/home/linux/relax_work_v3/relax_remaining_progress.log

TOTAL=$(wc -l < "$GENE_LIST")
echo "Starting RELAX on $TOTAL remaining genes with 4 cores, 600s timeout"
echo "Started: $(date)"

count=0
ok=0
fail=0
timeout_count=0
skip=0

cat "$GENE_LIST" | xargs -P 4 -I {} bash -c '
    gene="{}"
    fasta="'"$FASTA_DIR"'/${gene}.fasta"
    tree="'"$TREE_DIR"'/${gene}.nwk"
    json="'"$FASTA_DIR"'/${gene}.fasta.RELAX.json"

    # Remove old empty JSON if exists
    if [ -f "$json" ] && [ ! -s "$json" ]; then
        rm -f "$json"
    fi

    # Skip if already has valid result
    if [ -f "$json" ] && [ -s "$json" ]; then
        echo "SKIP $gene"
        exit 0
    fi

    if [ ! -f "$fasta" ] || [ ! -f "$tree" ]; then
        echo "MISSING $gene"
        exit 1
    fi

    timeout 600 '"$HYPHY"' relax --alignment "$fasta" --tree "$tree" --typepermodel --test Test > /dev/null 2>&1

    rc=$?
    if [ $rc -eq 0 ] && [ -s "$json" ]; then
        echo "OK $gene"
        exit 0
    elif [ $rc -eq 124 ]; then
        echo "TIMEOUT $gene"
        exit 2
    else
        echo "FAIL $gene"
        exit 1
    fi
' 2>&1 | tee "$LOG_FILE" | while read -r line; do
    echo "$line"
    case "$line" in
        OK*) ok=$((ok+1)) ;;
        FAIL*|MISSING*) fail=$((fail+1)) ;;
        TIMEOUT*) timeout_count=$((timeout_count+1)) ;;
        SKIP*) skip=$((skip+1)) ;;
    esac
    count=$((count+1))
    if [ $((count % 50)) -eq 0 ]; then
        echo "  --- Progress: $count/$TOTAL | OK=$ok Skip=$skip Fail=$fail Timeout=$timeout_count ---"
    fi
done

echo ""
echo "=== Completed: $(date) ==="
echo "Total processed: $count | OK: $ok | Skip: $skip | Fail: $fail | Timeout: $timeout_count"
