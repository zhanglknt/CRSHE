#!/bin/bash
# R8 optional Task 2 babysitter: exit when busted master finishes (done results=)
# or when the driver disappears without finishing (death report). 6h cap.
LOG=/root/hyphy_work_phase9/logs/sim_realign_busted.log
for i in $(seq 1 180); do
    if grep -q 'done results=' "$LOG" 2>/dev/null; then
        echo "SIMBUSTED_DONE $(grep 'done results=' "$LOG" | tail -1)"
        exit 0
    fi
    if ! pgrep -f 'simrealign_busted_master' >/dev/null 2>&1 \
       && ! pgrep -f 'busted_simrealign_one' >/dev/null 2>&1; then
        echo "SIMBUSTED_DRIVER_GONE $(date '+%F %T')"
        exit 1
    fi
    sleep 120
done
echo "SIMBUSTED_WAIT_TIMEOUT_6H"
exit 2
