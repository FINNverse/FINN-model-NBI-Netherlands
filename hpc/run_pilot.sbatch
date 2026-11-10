#!/bin/bash
# ===========================================================================
# FINN-model-NBI-Netherlands | run ONE v0.2 pilot (or any calibration script)
# on a SLURM CPU node, single-threaded
#   sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p3_v2_harvest.R
#
# Resources only. Every scientific choice (data set, model, epochs, losses)
# is written in the R script itself. Expects the repository at
# $HOME/FINN-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands and the data at
# $HOME/FINN-zoo/FINN-zoo-data/nbi (the ../../FINN-zoo-data/nbi layout the
# scripts use; see calibration/v02/CLUSTER.md), FINN 0.2.0 (0848b60) in
# $HOME/R/lib-finn-0848b60 (its torch lacks optim_ignite_adam; the scripts
# fall back to optim_adam) and the remaining packages in $R_LIBS_USER.
# ===========================================================================
#SBATCH --job-name=finn-nbi-pilot
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/pilot_%j.out
#SBATCH --error=logs/pilot_%j.err
set -euo pipefail
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/site-library}"
export LD_LIBRARY_PATH="$R_LIBS_USER/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$HOME/FINN-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch hpc/run_pilot.sbatch calibration/v02/<script>.R}"
echo "[job] $SCRIPT on $(hostname), 1 thread, $(date)"
Rscript "$SCRIPT"
echo "[job] done $(date)"
