# calibration/v020: the refit with the FINN 0.2.0 release defaults ("v020-defaults")

One joint `fit()` of the Dutch NFI model (NBI remeasurement pairs 4-12 years apart,
data `pairs-v2`, 14 taxa) with every FINN 0.2.0 default left at its default, as
decided for the whole zoo on 2026-11-13. The released v0.1.0 and v0.2.0 files and the v0.2 recipe fit folder
(`fits/pairs-v2-recipe`, completed on the cluster 2026-11-10) are untouched; this is a
new fit folder.

## Set-up

- **FINN**: release candidate `release-0.2.0` at b6244e5 (PR FINNverse/FINN#68) in
  `~/R/lib-finn-release-b6244e5` (laptop and cluster; the only library these scripts load).
- **Calibration**: ONE joint fit of competition, growth, mortality and regeneration
  against all six responses (dbh, ba, trees, growth, mortality, regeneration),
  `loss_weights = "auto"`, `update_step = 1`, `loss_recruit = "simulated"`.
- **Defaults in use** (none passed, read back from the built model and recorded in
  `fit_summary.md`): `regeneration_capped` (cap "all", prior "max_obs_recruits" =
  max(reg / period_length) per ha and year = 229.2, prior_lambda 0.1), theta trained
  with `theta_floor = 0.5`, `mortality_draw = "binomial"`, `regeneration_sampler =
  "refined"`, mortality limits 0.001 / 0.999, `reg_floor = 1e-3`, `trees_offset = 0.01`,
  `rb_gradient = TRUE`, `light_mask = TRUE`, competition height bound log(2.3).
- **Model structure** (as the v0.2 recipe fit, `calibration/v02/recipe/`): U-shaped
  mortality; growth and mortality env = SELU network 20-20 (dropout 0.1) on the
  WorldClim normals `tavg_c + prec_mm` plus one linear per-species weight on the
  harvest flag `harvest` (0/1, on its own scale) outside the network (module
  `env_dnn_linmgmt`, pure-shift check in the script), in growth AND mortality as the
  recipe (v0.2.0 had it in the mortality network only); static climate; recruit
  observation operator `recruit_obs_weight` = the median smaller circle of the training
  pairs (10 m, 0.0314 ha; asserted from `site_table.csv`); recruits enter at 5.9 cm; one
  0.1257 ha patch per plot; `growth_rate = "annual"`; growth and mortality scored per
  tracked tree (`tree_obs_dt`, `period_mode = "per_site"`).
- **Schedule**: 300 epochs at lr 0.01, then a cosine decay to 1e-4 over 200 epochs
  (500 in total), plain Adam (`torch::optim_adam`), batch 100 (the recipe's),
  `clip_norm` 2, checkpoints (parameter snapshot + model file) every 25 epochs, seed 1;
  a 24 h cluster segment, continued by `fit_continue.R` if needed.
- **Data**: `FINN-zoo-data/nbi/prepared/pairs-v2/train` (1,179 pairs), unchanged;
  `holdout_blocked` (583 pairs in 25 km blocks) for the validation. The inventory
  problems of the NBI data and how the prepared data handle them are listed in the
  header of `fit.R` (rare-species rule without an "other" class and without shrubs,
  harvest covariate, total stem loss incl. felling, the weak stand level against
  persistence, intervals 4-12 yr, variable-radius circles and the recruit operator,
  ingrowth vs recruitment, dispersion collapse -> theta floor, fate codes per cycle,
  MFV positions, 1 km cells).

## Differences to the v0.2 recipe (`calibration/v02/recipe/`)

| | v0.2 recipe | v020-defaults |
|---|---|---|
| stages | A: demography with the observed recruits fed in; B: regeneration alone with frozen demography | one joint fit |
| regeneration process | `regeneration_saturation`, cap K per species started at 50, theta frozen at its per-species ML value in count space (0.01-0.1) | `regeneration_capped` (0.2.0 default): one cap K shared by all species, started at and pulled towards max(reg / period_length) (prior_lambda 0.1), theta trained above the floor 0.5 (expected to bind) |
| death draw / mortality limits | gamma draw, floor 0.01 / ceiling 0.99 | binomial draw, 0.001 / 0.999 |
| recruit draw | gaussian sampler | refined sampler |
| trees loss | offset 1e-5 | offset 0.01, `rb_gradient` |
| schedule | cosine 0.005 -> 5e-5 over 300 (A) + 150 (B) epochs | constant 0.01 for 300, cosine to 1e-4 over 200 |
| FINN | main 0848b60 in lib-finn-0848b60 | release-0.2.0 b6244e5 in lib-finn-release-b6244e5 |
| API | `weights`, `lr_scheduler`, `optimizer` of `fit()`; `finn(reg_floor =)` | `loss_weights`, `control = finnControl()`, `loss_recruit`; `finn(growth_rate =)` |

Everything else (data, U-shape, env module, harvest term, recruit operator, per-tree
losses, light mask, batch 100) is the recipe's.

## Files

- `fit.R`: the production fit (segment 1); refuses to start over existing checkpoints.
- `fit_continue.R`: resumes from the newest checkpoint with the schedule continued;
  exits at once when `model.pt` exists.
- `fit_smoke.R`: `fit.R` with the schedule cut to 3 epochs and its own output folder
  (`fits/pairs-v2-v020-smoke`); measures s/epoch and exercises every check.
- `validation.R`: the recipe validation adapted to one stage (holdout report with the
  persistence null, change metrics, strata, metrics by harvest flag, per-tree AUC,
  bare-ground succession, parameter trajectories).
- `../../hpc/run_v020.sbatch`: the resource-only wrapper (one CPU, one thread; the
  release library's torch 0.17 on the library path); `../../hpc/submit_v020.sh`: the
  chain fit -> continue -> validation.

## Running

Laptop (one thread): `Rscript calibration/v020/fit_smoke.R`, then
`Rscript calibration/v020/fit.R` and `Rscript calibration/v020/validation.R` from the
repository root. Cluster: the smoke job, then `bash hpc/submit_v020.sh`.

Outputs in `FINN-zoo-data/nbi/fits/pairs-v2-v020/` (outside git): `model_init.pt`,
`model.pt`, `checkpoints/`, `loss_history.csv`, `lr_schedule.csv`, `species_parameters.csv`
(incl. theta and the U-shape), `regeneration_parameters.csv` (predicted / observed
recruitment per species, K), `management_weights.csv`, `module_checks.csv`,
`grad_norms*.csv`, `cap_trajectory.csv`, `convergence.csv`, `fit_summary.md`,
`report.md`, `validation/`.

Convergence is judged on the parameter drift between checkpoints (`convergence.csv`,
`validation/fig6_parameters.png`), never on the loss alone.

## Smoke (laptop, 2026-11-13, one thread, torch 0.15.1)

3 epochs on the full training set: 11.5 s/epoch at batch 100 (12 steps per epoch),
peak RSS 0.7 GB, every check passed, every output written, checkpoint files ~10 MB
each. Baseline-scaled losses after 3 epochs: dbh 0.30, ba 0.14, trees 0.32, growth 0.97,
mortality 0.98, regeneration 1.02. Predicted / observed recruitment on train 1.53 at the
start (K = prior) -> 0.90; theta 1 -> 0.82-0.88; K 229.2 -> 229.5. 98.0 % of the
recruiting rows are whole counts on their own circle; reg x w quartiles 0.83 / 1.61 /
3.64. `validation.R` and `fit_continue.R` were exercised on the smoke folder.
