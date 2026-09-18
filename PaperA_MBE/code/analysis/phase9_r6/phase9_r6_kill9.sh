#!/bin/bash
# force-kill any remaining hyphy/xargs/master processes
pkill -9 -f 'hyphy busted' 2>/dev/null
pkill -9 -f 'phase9_r6' 2>/dev/null
pkill -9 -f 'xargs' 2>/dev/null
sleep 2
ps aux | grep -E 'hyphy busted|xargs|master' | grep -v grep | head
echo KILLED
