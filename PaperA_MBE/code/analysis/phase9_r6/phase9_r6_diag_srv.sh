#!/bin/bash
BASE=/root/hyphy_work_phase9
echo ===SRV_LOG_TAIL===
tail -10 "$BASE/logs/srv.log" 2>&1
echo ===RUNNING===
ps aux | grep '[h]yphy busted' | awk '{print $2, $10, $NF}' | head -10
echo ===MASTER2_LOG===
tail -10 /mnt/d/人类正选择基因项目/results/phase9_hardening/run_logs/master2.log 2>&1
