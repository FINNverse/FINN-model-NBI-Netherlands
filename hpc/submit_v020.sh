#!/bin/bash
# ===========================================================================
# FINN-model-NBI-Netherlands | submit the v020-defaults production chain on the cluster
#   bash hpc/submit_v020.sh            (from the repository root, after the smoke job)
#
# Three jobs, each one self-contained R script run by hpc/run_v020.sbatch
# (resources only; every scientific choice is in the script):
#   fit        calibration/v020/fit.R           24 h segment 1 (refuses to start over existing checkpoints)
#   continue   calibration/v020/fit_continue.R  24 h segment 2 (afterany: fit; exits at once when model.pt exists)
#   report     calibration/v020/validation.R    holdout report, 4 CPUs (afterok: segment 2)
# Sizing from the smoke (fit_smoke.R: 11.5 s/epoch on a laptop core at batch
# 100; the cluster's cores ran the v0.2 recipe stage A ~3x slower): 500 epochs
# ~ 3-5 h on the cluster, so segment 1 finishes the fit and the continuation
# exits at once; it is the safety net for a slower node. Memory: 0.7 GB peak
# on the laptop smoke, 8 G requested.
# Smoke first: sbatch --time=01:00:00 --mem=8G hpc/run_v020.sbatch calibration/v020/fit_smoke.R
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
F=$(sbatch --parsable --job-name=zoo-nl-v020-fit --time=24:00:00 --mem=8G hpc/run_v020.sbatch calibration/v020/fit.R)
C1=$(sbatch --parsable --job-name=zoo-nl-v020-c1 --time=24:00:00 --mem=8G --dependency=afterany:"$F" hpc/run_v020.sbatch calibration/v020/fit_continue.R)
V=$(sbatch --parsable --job-name=zoo-nl-v020-val --time=02:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$C1" hpc/run_v020.sbatch calibration/v020/validation.R)
echo "submitted: fit $F -> continue $C1 (afterany) -> validation $V (afterok); logs in logs/v020_<jobid>.out|err"
