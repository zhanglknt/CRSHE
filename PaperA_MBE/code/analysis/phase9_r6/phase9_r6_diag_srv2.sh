#!/bin/bash
BASE=/root/hyphy_work_phase9
echo "srv results files: $(ls $BASE/results_srv | wc -l)"
echo "srv.log totals:"
awk -F, '{c[$1]++} END {for (k in c) print k, c[k]}' "$BASE/logs/srv.log"
echo "latest 15 lines:"
tail -15 "$BASE/logs/srv.log"
