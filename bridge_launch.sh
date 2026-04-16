#!/bin/bash
# Self-daemonizing bridge script
(
    # Close stdin/stdout/stderr, start new session
    exec </dev/null >/tmp/bridge.log 2>&1
    # Wait for fix_all_resume.sh (PID 33416) to finish
    while kill -0 33416 2>/dev/null; do sleep 30; done
    echo "[$(date)] fix_all_resume.sh finished, appending DONE marker"
    echo "=== DONE ===" >> /tmp/step6_out.txt
    echo "[$(date)] Bridge complete"
) &
disown
echo "bridge launched as PID $!"