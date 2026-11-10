#!/usr/bin/env bash
# FINN-model-NBI-Netherlands | calibration/v02/run_pilots.sh
# Runs the v0.2 pilots ONE AT A TIME (laptop budget: one pilot process per agent),
# each single-threaded, logging to ../../FINN-zoo-data/nbi/fits/v02-pilots/<name>.log.
# Resources and sequencing only; every scientific setting is inside the pilot scripts.
# Run from the repository root:  bash calibration/v02/run_pilots.sh [pilot names...]
set -u
LOGDIR="../../FINN-zoo-data/nbi/fits/v02-pilots"; mkdir -p "$LOGDIR"
PILOTS=("$@"); if [ ${#PILOTS[@]} -eq 0 ]; then PILOTS=(p0_v1_reference p1_v2_data p2_v2_hygiene p3_v2_harvest p4_v2_linear p5_v2_countspace p6_v2_staged); fi
for p in "${PILOTS[@]}"; do
  echo "==== $(date '+%F %T') start $p" | tee -a "$LOGDIR/chain.log"
  Rscript "calibration/v02/pilot_$p.R" > "$LOGDIR/$p.log" 2>&1
  echo "==== $(date '+%F %T') end $p status $?" | tee -a "$LOGDIR/chain.log"
done
