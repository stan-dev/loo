#!/usr/bin/env bash
# Measure peak total RSS (MB) of an R run and its entire process tree
# (mclapply forks or mirai daemons), by sampling /proc every 20 ms.
#
# Usage: peak-mem.sh LOO_LIB BENCH_LABEL MODE CORES   (MODE = per-call|persist)
set -u
LOO_LIB="$1"; BENCH_LABEL="$2"; MODE="$3"; CORES="$4"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LOO_LIB="$LOO_LIB" BENCH_LABEL="$BENCH_LABEL" MODE="$MODE" CORES="$CORES" \
  Rscript "$SCRIPT_DIR/peak-mem-run.R" &
RPID=$!

tree_rss_kb() {
  # Sum RSS (KB) of the root pid and all descendants.
  local root=$1
  local pids=$root
  local frontier=$root
  while [ -n "$frontier" ]; do
    local next=""
    for p in $frontier; do
      local kids
      kids=$(pgrep -P "$p" 2>/dev/null | tr '\n' ' ')
      next="$next $kids"
    done
    frontier=$(echo "$next" | xargs)
    pids="$pids $frontier"
  done
  local total=0 rss
  for p in $pids; do
    rss=$(awk '/^VmRSS:/{print $2}' "/proc/$p/status" 2>/dev/null)
    [ -n "${rss:-}" ] && total=$((total + rss))
  done
  echo "$total"
}

PEAK=0
while kill -0 "$RPID" 2>/dev/null; do
  CUR=$(tree_rss_kb "$RPID")
  [ "$CUR" -gt "$PEAK" ] && PEAK=$CUR
  sleep 0.02
done
wait "$RPID"
PEAK_MB=$(echo "$PEAK/1024" | bc -l)

# Append a parseable record so compare.R can fold the peak-RSS results into its
# report. Override the destination with PEAK_OUT; delete it between fresh runs.
PEAK_OUT="${PEAK_OUT:-/tmp/bench-peakmem.tsv}"
printf "%s\t%s\t%s\t%.0f\n" "$BENCH_LABEL" "$MODE" "$CORES" "$PEAK_MB" >> "$PEAK_OUT"

printf "PEAK_RSS %s/%s/cores=%s : %.0f MB  (appended to %s)\n" \
  "$BENCH_LABEL" "$MODE" "$CORES" "$PEAK_MB" "$PEAK_OUT"
