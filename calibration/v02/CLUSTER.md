# Running the v0.2 pilots (and the v0.2 fit) on the cluster

Written 2026-11-10 for the coordinator, who syncs and submits. Nothing here is a scientific
setting; those are in the scripts.

## Layout expected on the cluster
```
$HOME/FINN-zoo/FINN-zoo-models/FINN-model-NBI-Netherlands/   this repository (git checkout or rsync)
$HOME/FINN-zoo/FINN-zoo-data/nbi/prepared/pairs-v1/          pilot/, holdout_blocked/  (P0 only)
$HOME/FINN-zoo/FINN-zoo-data/nbi/prepared/pairs-v2/          pilot/, train/, holdout_blocked/
$HOME/FINN-zoo/FINN-zoo-data/nbi/fits/                       written by the scripts (v02-pilots/<name>/, pairs-v2-process/)
$HOME/R/lib-finn-0848b60                                     FINN 0.2.0 (main 0848b60); its torch is old: the scripts use optim_adam there
```
The scripts address the data as `../../FINN-zoo-data/nbi` relative to the repository root and
set `.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))`; the sbatch wrapper `cd`s to the
repository. Other packages needed: data.table, torch, ggplot2 (validation only), hexbin
(validation only), jsonlite (fit only). No WorldClim rasters and no raw data are needed for
pilots or fits (only for the preparation scripts, which ran on the laptop).

## Files to sync (data; small)
| what | path under FINN-zoo-data/nbi | size |
|---|---|---|
| pairs-v1 pilot + holdout (P0) | `prepared/pairs-v1/pilot/`, `prepared/pairs-v1/holdout_blocked/` | ~4 MB |
| pairs-v2 pilot + train + holdout | `prepared/pairs-v2/pilot/`, `prepared/pairs-v2/train/`, `prepared/pairs-v2/holdout_blocked/` | ~8 MB |
| species table (card copy) | `prepared/pairs-v2/species.csv` | 1 KB |
Each set folder holds obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv.

## Pilots: one job each, one CPU, 1 thread, <= 2 h, 8 GB
```
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p0_v1_reference.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p1_v2_data.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p2_v2_hygiene.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p3_v2_harvest.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p4_v2_linear.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p5_v2_countspace.R
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_p6_v2_staged.R
```
Expected: 300 sites x 400 epochs at ~1.6 s/epoch on the laptop (M3, one thread); a cluster
core is roughly half as fast, so ~25 min per pilot plus ~3 min scoring; < 2 GB RAM. Each pilot
writes `fits/v02-pilots/<name>/{model.pt, loss_history*.csv}` and appends (or replaces) its
row in `calibration/v02/pilot_table.csv` (sync that file and `fits/v02-pilots/*/loss_history*.csv`
back; model.pt files are not needed back). The pilots are independent and may run in parallel.
A laptop chain (`calibration/v02/run_pilots.sh`, armed behind a `< 7 R processes` waiter) may
run the same pilots if the laptop frees up; the table row is keyed by pilot name, so a second
run simply replaces the first (same seed, same result up to the torch build).

## The v0.2 production fit (after the pilot decision; not before)
```
sbatch hpc/run_pilot.sbatch calibration/v02/fit.R          # ~1,524 sites; see fit.R header for the expected s/epoch
sbatch --dependency=afterok:<jobid> hpc/run_pilot.sbatch calibration/v02/validation.R
```
Sync back: `fits/pairs-v2-process/{model.pt, parameters.pt, loss_history.csv}`,
`calibration/v02/{fit_summary.md, loss_history.csv, report.md, validation/, figures/}`,
`inst/model/weights.json`.
