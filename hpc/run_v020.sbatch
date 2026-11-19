#!/bin/bash
# ===========================================================================
# FINN-model-NBI-Netherlands | run ONE v020-defaults script on a SLURM CPU core
#   sbatch --time=24:00:00 --mem=8G hpc/run_v020.sbatch calibration/v020/fit.R
# Resources only (one CPU, one torch thread; time and memory on the sbatch
# line). Every scientific choice is inside the R script. Expects the repository
# at $HOME/finn-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands and the data at
# $HOME/finn-zoo/FINN-zoo-data/nbi (the ../../FINN-zoo-data/nbi layout the
# scripts use), FINN release-0.2.0 at b6244e5 in $HOME/R/lib-finn-release-b6244e5
# (FINN only; the script puts it first on .libPaths()), torch 0.17 in
# $HOME/R/site-library-torch017 and the remaining packages in $HOME/R/site-library.
# ===========================================================================
#SBATCH --job-name=zoo-nl-v020
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/v020_%j.out
#SBATCH --error=logs/v020_%j.err
set -euo pipefail
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
# torch 0.17 (the release library holds FINN only) before the general site library,
# and its libtorch on LD_LIBRARY_PATH: the libtorch of the older torch in the site
# library fails with "Lantern is not loaded" (the France finding of 2026-11-10).
export R_LIBS_USER="$HOME/R/site-library-torch017:$HOME/R/site-library"
export LD_LIBRARY_PATH="$HOME/R/site-library-torch017/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd "${SLURM_SUBMIT_DIR:-$HOME/finn-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch hpc/run_v020.sbatch calibration/v020/<script>.R}"
echo "[v020] $SCRIPT on $(hostname), $OMP_NUM_THREADS thread(s), $(date)"
Rscript "$SCRIPT"
echo "[v020] done $(date)"
