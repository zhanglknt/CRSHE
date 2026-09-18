#!/bin/bash
# systemd launcher for phase9 R6 master pipeline
exec bash /mnt/d/人类正选择基因项目/scripts/phase9_r6_master.sh \
    > /mnt/d/人类正选择基因项目/results/phase9_hardening/run_logs/master.log 2>&1
