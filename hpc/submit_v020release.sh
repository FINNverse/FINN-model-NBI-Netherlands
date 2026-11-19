#!/bin/bash
# ===========================================================================
# FINN-model-NBI-Netherlands | submit the v020-release production chain on the cluster
#   bash hpc/submit_v020release.sh        (from the repository root)
#
# The chain of hpc/submit_v020soil.sh (same shape and limits) for
# calibration/v020-release (the v020-soil fit refitted unchanged on the released
# FINN 0.2.0, commit 6c9156f; data pairs-v3-soil). Three jobs, each one
# self-contained R script run by hpc/run_v020release.sbatch (resources only;
# every scientific choice is in the script):
#   fit        calibration/v020-release/fit.R           24 h segment 1 (refuses to start over existing checkpoints)
#   continue   calibration/v020-release/fit_continue.R  24 h segment 2 (afterany: fit; exits at once when model.pt exists)
#   report     calibration/v020-release/validation.R    holdout report, 4 CPUs (afterok: segment 2)
# The v020-soil fit took ~3 h on one cluster core, so segment 1 should finish the
# fit and the continuation exit at once. The sizing notes of the soil chain follow.
# Sizing from measurements: the v020 fit took 22.1 s/epoch on one cluster core
# (184 min for 500 epochs, fits/pairs-v2-v020/fit_summary.md); on the laptop
# (efficiency cores, nice 15, beside another running fit) the soil smoke took
# 88 s/epoch against 75 s/epoch for v020 under the same conditions (+17 % wall
# clock incl. one checkpoint per epoch in the soil smoke only; +4 % CPU time),
# so expect 23-26 s/epoch, 500 epochs ~ 3.2-3.6 h: segment 1 finishes the fit
# and the continuation exits at once (safety net for a slower node). Memory:
# 0.45 GB peak RSS (fit), 0.43 GB (resume), 0.42 GB (validation smoke with 60
# bare-ground years); 8 G requested for the fit, 16 G for the validation as v020.
# Prerequisite: prepared/pairs-v3-soil synced to $HOME/finn-zoo/FINN-zoo-data/nbi/prepared/.
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
test -d ../../FINN-zoo-data/nbi/prepared/pairs-v3-soil/train || { echo "prepared/pairs-v3-soil missing: sync it first"; exit 1; }
F=$(sbatch --parsable --job-name=zoo-nl-v020rel-fit --time=24:00:00 --mem=8G hpc/run_v020release.sbatch calibration/v020-release/fit.R)
C1=$(sbatch --parsable --job-name=zoo-nl-v020rel-c1 --time=24:00:00 --mem=8G --dependency=afterany:"$F" hpc/run_v020release.sbatch calibration/v020-release/fit_continue.R)
V=$(sbatch --parsable --job-name=zoo-nl-v020rel-val --time=02:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$C1" hpc/run_v020release.sbatch calibration/v020-release/validation.R)
echo "submitted: fit $F -> continue $C1 (afterany) -> validation $V (afterok); logs in logs/v020rel_<jobid>.out|err"
