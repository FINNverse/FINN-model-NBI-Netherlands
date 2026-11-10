#!/bin/bash
# ===========================================================================
# FINN-model-NBI-Netherlands | submit the zoo v0.2 recipe chain on the cluster
#   bash hpc/submit_recipe.sh            (from the repository root)
#
# Three jobs, each one self-contained R script run by hpc/run_pilot.sbatch
# (resources only; every scientific choice is in the script; the wrapper puts
# the libtorch of lib-finn-0848b60 on LD_LIBRARY_PATH):
#   stage A  calibration/v02/recipe/fit_stageA.R   demography, observed recruitment
#   stage B  calibration/v02/recipe/fit_stageB.R   regeneration (afterok: stage A)
#   report   calibration/v02/recipe/validation.R   holdout report (afterok: stage B)
# Time limits from the measured cost (laptop, one thread, 2026-11-10): stage A
# 1.5 s/epoch on the 200 smoke pairs beyond the fixed cost (12 env rows, expanded
# trees) -> ~9 s/epoch on the 1,179 training pairs, 300 epochs ~45 min (A) and
# 150 epochs ~23 min (B); the cluster's cores ran the v0.2 pilots ~3x slower
# than the laptop (~2.3 h / ~1.1 h), hence 6 h and 3 h below. Validation
# simulates the 583 holdout and 1,179 training pairs plus six 500-year
# bare-ground runs.
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
A=$(sbatch --parsable --job-name=zoo-nl-A --time=06:00:00 --mem=8G hpc/run_pilot.sbatch calibration/v02/recipe/fit_stageA.R)
B=$(sbatch --parsable --job-name=zoo-nl-B --time=03:00:00 --mem=8G --dependency=afterok:"$A" hpc/run_pilot.sbatch calibration/v02/recipe/fit_stageB.R)
V=$(sbatch --parsable --job-name=zoo-nl-V --time=02:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$B" hpc/run_pilot.sbatch calibration/v02/recipe/validation.R)
echo "submitted: stage A $A -> stage B $B (afterok) -> validation $V (afterok); logs in logs/pilot_<jobid>.out|err"
