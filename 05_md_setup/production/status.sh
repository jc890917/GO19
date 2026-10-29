#!/usr/bin/env bash
# Production status for one instance.  Run on the instance: bash status.sh
cd /workspace/run 2>/dev/null || { echo "  no /workspace/run"; exit 1; }

TOTAL_NS=100
done_n=0; run_n=0; warn_n=0

for d in */prod_rep*; do
    [ -d "$d" ] || continue
    log="$d/prod.log"
    [ -f "$log" ] || { printf "  %-32s queued\n" "$d"; continue; }

    lincs=$(grep -c "LINCS WARNING" "$log" 2>/dev/null); lincs=${lincs:-0}
    warn_n=$((warn_n + lincs))

    if [ -f "$d/prod.gro" ]; then
        perf=$(awk '/^Performance:/{p=$2} END{print p}' "$log")
        printf "  %-32s DONE      100/100 ns   %s ns/day   LINCS=%s\n" "$d" "${perf:-?}" "$lincs"
        done_n=$((done_n + 1))
    else
        ps=$(awk '/^ *Step *Time/{getline; t=$2} END{printf "%.0f", t+0}' "$log")
        ns=$(awk -v t="$ps" 'BEGIN{printf "%.1f", t/1000}')
        pct=$(awk -v n="$ns" -v T="$TOTAL_NS" 'BEGIN{printf "%.0f", 100*n/T}')
        eta=$(awk -v n="$ns" -v T="$TOTAL_NS" -v log="$log" 'BEGIN{print ""}')
        printf "  %-32s running   %5s/100 ns  (%s%%)   LINCS=%s\n" "$d" "$ns" "$pct" "$lincs"
        run_n=$((run_n + 1))
    fi
done

gpu=$(nvidia-smi --query-gpu=utilization.gpu,temperature.gpu --format=csv,noheader 2>/dev/null)
alive=$(pgrep -cf "gmx mdrun")
df=$(df -h /workspace | awk 'NR==2{print $4" free"}')
printf "  -- %s done, %s running | mdrun procs: %s | GPU %s | disk %s\n" \
       "$done_n" "$run_n" "$alive" "$gpu" "$df"
[ "$warn_n" != "0" ] && printf "  !! %s LINCS warnings total -- investigate\n" "$warn_n"
exit 0
