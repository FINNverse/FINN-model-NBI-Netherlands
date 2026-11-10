# Running the v0.2 pilots (and the v0.2 fit) on the cluster

Written 2026-11-10 for the coordinator, who syncs and submits. Nothing here is a scientific
setting; those are in the scripts.

## Layout expected on the cluster
```
$HOME/FINN-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands/   this repository (git checkout or rsync)
$HOME/FINN-zoo/FINN-zoo-data/nbi/prepared/pairs-v2/          pilot/, train/, holdout_blocked/   (all pilots P1-P7, the fit)
$HOME/FINN-zoo/FINN-zoo-data/nbi/prepared/pairs-v1/          pilot/, holdout_blocked/           (P0 only; already run on the laptop)
$HOME/FINN-zoo/FINN-zoo-data/nbi/fits/                       written by the scripts (v02-pilots/<name>/, pairs-v2-process/)
$HOME/R/lib-finn-0848b60                                     FINN 0.2.0 (main 0848b60); its torch is old: the scripts fall back to optim_adam
$HOME/R/lib-finn-zoo-v02                                     branch build zoo-v02-mortfloor (FINN#51 / PR #52): needed by pilot_p6_floor.R ONLY,
                                                             which puts it first in .libPaths(); the coordinator installs it
```
The scripts address the data as `../../FINN-zoo-data/nbi` relative to the repository root and
set `.libPaths()` themselves; the sbatch wrapper `cd`s to the repository. Other packages
needed: data.table, torch, ggplot2 + hexbin (validation only), jsonlite (fit only). No
WorldClim rasters and no raw data are needed for pilots or fits (only for the preparation
scripts, which ran on the laptop).

## Files to sync (data; small)
| what | path under FINN-zoo-data/nbi | size |
|---|---|---|
| pairs-v2 pilot + train + holdout | `prepared/pairs-v2/{pilot,train,holdout_blocked}/` | ~8 MB |
| species table (card copy) | `prepared/pairs-v2/species.csv` | 1 KB |
Each set folder holds obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv.

## Pilots: one job each, one CPU, 1 thread, <= 2 h, 8 GB
```
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p1_baseline.R     # v0.1 settings on pairs-v2
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p2_cosine.R       # + cosine lr
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p3_theta_cap.R    # + theta frozen, K init 10
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p4_pertree.R      # + per-tree growth / mortality losses (tree_obs_dt)
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p5_harvest.R      # + harvest covariate in mortality
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p7_combo.R        # P2 + P3 + P5 together (no-floor twin of P6)
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p6_floor.R        # P7 + mortality_floor 1e-6 / poisson draw; needs ~/R/lib-finn-zoo-v02
```
P0 (the v0.1 model on pairs-v1 at the pilot budget) already ran on the laptop; its row is in
`pilot_table.csv`. Expected per pilot: 300 sites x 400 epochs at ~2 s/epoch on a laptop core
(P0: 13 min); a cluster core is roughly half as fast, so ~25-30 min plus ~3 min scoring;
< 2 GB RAM. Each pilot writes `fits/v02-pilots/<name>/{model.pt, loss_history.csv}` and
appends (or replaces) its rows in `calibration/v02/pilot_table.csv` and `pilot_strata.csv`
(sync those two files and `fits/v02-pilots/*/loss_history.csv` back; model.pt files are not
needed back). The pilots are independent and may run in parallel. A laptop chain
(`calibration/v02/run_pilots.sh`, armed behind a `< 7 R processes` waiter) may run the same
pilots if the laptop frees up; rows are keyed by pilot name, so a second run replaces the
first (same seed; identical up to the torch build / optimizer fallback).

## The v0.2 production fit (after the pilot decision; not before)
```
sbatch hpc/run_pilot.sbatch calibration/v02/fit.R          # ~1,179 sites; see fit.R header for the expected s/epoch
sbatch --dependency=afterok:<jobid> hpc/run_pilot.sbatch calibration/v02/validation.R
```
Sync back: `fits/pairs-v2-process/{model.pt, parameters.pt, loss_history.csv}`,
`calibration/v02/{fit_summary.md, loss_history.csv, report.md, validation/, figures/}`,
`inst/model/weights.json`.
