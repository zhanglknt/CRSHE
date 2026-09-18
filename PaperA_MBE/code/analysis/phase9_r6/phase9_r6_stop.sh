#!/bin/bash
# stop current master pipeline (v1) cleanly: kill xargs + master.sh + their children
pkill -f 'phase9_r6_master.sh' 2>/dev/null
pkill -f 'xargs -P 4' 2>/dev/null
sleep 2
pkill -f 'hyphy busted' 2>/dev/null
pkill -f 'phase9_r6_busted_one.sh' 2>/dev/null
sleep 2
ps aux | grep -E 'master|xargs|hyphy busted' | grep -v grep | head
echo STOPPED
