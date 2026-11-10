#!/usr/bin/env bash
# FINN-model-NBI-Netherlands | calibration/v02/run_pilots.sh
# Runs the v0.2 pilots ONE AT A TIME (laptop budget: one pilot process per agent),
# each single-threaded, logging to ../../FINN-zoo-data/nbi/fits/v02-pilots/<name>.log.
# Resources and sequencing only; every scientific setting is inside the pilot scripts.
# Run from the repository root:  bash calibration/v02/run_pilots.sh [pilot names...]
set -u
LOGDIR="../../FINN-zoo-data/nbi/fits/v02-pilots"; mkdir -p "$LOGDIR"
PILOTS=("$@"); if [ ${#PILOTS[@]} -eq 0 ]; then PILOTS=(p1_baseline p2_cosine p3_theta_cap p4_pertree p5_harvest p7_combo p6_floor); fi
for p in "${PILOTS[@]}"; do
  echo "==== $(date '+%F %T') start $p" | tee -a "$LOGDIR/chain.log"
  Rscript "calibration/v02/pilot_$p.R" > "$LOGDIR/$p.log" 2>&1
  echo "==== $(date '+%F %T') end $p status $?" | tee -a "$LOGDIR/chain.log"
done
