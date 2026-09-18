#!/bin/bash
# progress check for phase9 R6 pipeline
BASE=/root/hyphy_work_phase9
echo "=== $(date) ==="
for G in baseline srv cpg topo sim; do
  N=$(ls "$BASE/results_$G" 2>/dev/null | wc -l)
  echo "  $G: $N results"
done
if [ -f "$BASE/logs/baseline.log" ]; then
  echo "baseline.log tail:"
  grep -c '^OK,' "$BASE/logs/baseline.log" 2>/dev/null | xargs echo "  OK count:"
fi
tail -5 /mnt/d/人类正选择基因项目/results/phase9_hardening/run_logs/master.log 2>/dev/null
ps aux | grep -c "[h]yphy busted" | xargs echo "  running hyphy procs:"
