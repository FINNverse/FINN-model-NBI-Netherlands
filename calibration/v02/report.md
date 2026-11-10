# Validation report: pairs-v2-process (v0.2)

Written 2026-11-10 10:30 CET by calibration/v02/validation.R; FINN 0.2.0 (main 0848b60), torch 0.15.1. Model: ../../FINN-zoo-data/nbi/fits/pairs-v2-process/model.pt; release object parameters.pt sha256 29c8a4c21ca7918c07bc15a21ca0e3c315f933eddc4781f9038790719f797443.

## Sets

- TEST = pairs-v2/train: 1179 sites (plot x remeasurement pair) the model was fitted on; 630 MFV -> NBI-6, 549 NBI-6 -> NBI-7.
- HOLDOUT = pairs-v2/holdout_blocked: 583 sites in 25 km blocks the model never saw. No coordinates in the model; predicted with each site's recorded harvest flag (281 of them 1).
- Aggregation as fit() scores it, per site over its own interval (4-12 yr): state at year = period_length; growth and stem loss = mean of the annual rates; recruitment = sum.

## Metrics

| set | variable | level | n | rmse | r2 | r | bias | obs_mean | auc |
|---|---|---|---|---|---|---|---|---|---|
| test | dbh (cm, species mean) | species x site |  3000 |   6.7400 |  0.7350 | 0.8740 | -1.88000 |  23.3000 | NA |
| test | basal area (m2 per 0.1257 ha patch) | species x site | 16500 |   0.2120 |  0.9300 | 0.9650 | -0.01250 |   0.2430 | NA |
| test | trees (per 0.1257 ha patch) | species x site | 16500 |   9.4500 |  0.8190 | 0.9070 |  0.33400 |   5.7500 | NA |
| test | growth (rel. dbh increment / yr) | species x site |  2760 |   0.0130 |  0.4990 | 0.7180 | -0.00226 |   0.0202 | NA |
| test | stem loss (annual rate, incl. felling) | species x site |  2900 |   0.2100 |  0.0169 | 0.2500 | -0.03710 |   0.0723 | 0.697 |
| test | recruitment (stems / ha / interval) | species x site | 16500 |  43.3000 | -0.0339 | 0.0880 |  2.59000 |   4.9700 | NA |
| holdout | dbh (cm, species mean) | species x site |  1530 |   6.5200 |  0.7380 | 0.8720 | -1.66000 |  22.1000 | NA |
| holdout | basal area (m2 per 0.1257 ha patch) | species x site |  8160 |   0.2320 |  0.9090 | 0.9540 | -0.00777 |   0.2290 | NA |
| holdout | trees (per 0.1257 ha patch) | species x site |  8160 |   9.6500 |  0.8420 | 0.9180 |  0.32600 |   6.1500 | NA |
| holdout | growth (rel. dbh increment / yr) | species x site |  1400 |   0.0142 |  0.3830 | 0.6300 | -0.00184 |   0.0203 | NA |
| holdout | stem loss (annual rate, incl. felling) | species x site |  1490 |   0.2370 | -0.0198 | 0.1570 | -0.04820 |   0.0852 | 0.655 |
| holdout | recruitment (stems / ha / interval) | species x site |  8160 |  44.0000 | -0.0359 | 0.0680 |  2.74000 |   4.8900 | NA |
| test | basal area (m2 per patch) | stand total |  1180 |   0.7850 |  0.6870 | 0.8450 | -0.17400 |   3.4000 | NA |
| test | trees (per patch) | stand total |  1180 |  37.5000 |  0.6620 | 0.8300 |  4.68000 |  80.6000 | NA |
| test | recruitment (stems / ha / interval) | stand total |  1180 | 170.0000 | -0.0736 | 0.0969 | 36.30000 |  69.6000 | NA |
| test | basal area (m2 per patch) | stand total, PERSISTENCE reference |  1180 |   0.8490 |  0.6330 | 0.8230 | -0.22000 |   3.4000 | NA |
| test | trees (per patch) | stand total, PERSISTENCE reference |  1180 |  51.5000 |  0.3620 | 0.8120 | 13.70000 |  80.6000 | NA |
| test | basal-area CHANGE (m2 per patch) | stand total |  1180 |   0.7850 |  0.0847 | 0.4250 | -0.17400 |   0.2200 | NA |
| test | stem-number CHANGE (per patch) | stand total |  1180 |  37.5000 |  0.4300 | 0.6630 |  4.68000 | -13.7000 | NA |
| holdout | basal area (m2 per patch) | stand total |   583 |   0.8350 |  0.6350 | 0.8140 | -0.10900 |   3.2100 | NA |
| holdout | trees (per patch) | stand total |   583 |  38.4000 |  0.6990 | 0.8460 |  4.56000 |  86.1000 | NA |
| holdout | recruitment (stems / ha / interval) | stand total |   583 | 170.0000 | -0.0591 | 0.1080 | 38.30000 |  68.5000 | NA |
| holdout | basal area (m2 per patch) | stand total, PERSISTENCE reference |   583 |   0.8520 |  0.6190 | 0.8100 | -0.15700 |   3.2100 | NA |
| holdout | trees (per patch) | stand total, PERSISTENCE reference |   583 |  53.0000 |  0.4260 | 0.8250 | 15.40000 |  86.1000 | NA |
| holdout | basal-area CHANGE (m2 per patch) | stand total |   583 |   0.8350 |  0.0068 | 0.3240 | -0.10900 |   0.1570 | NA |
| holdout | stem-number CHANGE (per patch) | stand total |   583 |  38.4000 |  0.4270 | 0.6600 |  4.56000 | -15.4000 | NA |

A negative R2 means worse than predicting the mean. The persistence rows predict the second-visit state by the first-visit state (model-free reference).
"Stem loss" is the fitted mortality target: natural death AND felling (see fit.R).

## Stem-loss AUC (interval loss probability against the tallied fates; each site x species row weighted by its trees at risk)

| set | n_trees_at_risk | n_lost | auc | obs_loss_share | pred_loss_share |
|---|---|---|---|---|---|
| test | 23200 | 5080 | 0.697 | 0.219 | 0.243 |
| holdout | 11800 | 2760 | 0.655 | 0.234 | 0.255 |

## Stratified by the end-census share of species outside S (the agreed rare-species rule: where the error grows with the share, the model stops applying)

| set | oos_bin | n_sites | ba_stand_r2 | ba_persistence_r2 | ba_rmse_ha | ba_bias_ha | trees_stand_r2 | trees_bias_ha | dba_r2 | reg_ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| test | < 1 % | 865 | 0.678 | 0.593 | 6.15 | -1.350 | 0.674 | 32.50 |  0.1490 | 1.42 |
| test | 1-5 % | 218 | 0.763 | 0.771 | 5.91 | -1.770 | 0.730 | 33.30 | -0.1470 | 1.76 |
| test | > 5 % |  96 | 0.296 | 0.355 | 7.65 | -0.902 | 0.216 | 88.90 | -0.1100 | 1.94 |
| holdout | < 1 % | 436 | 0.679 | 0.643 | 6.41 | -0.580 | 0.739 | 43.30 |  0.0746 | 1.53 |
| holdout | 1-5 % | 108 | 0.344 | 0.442 | 7.37 | -2.010 | 0.405 |  5.16 | -0.2380 | 1.60 |
| holdout | > 5 % |  39 | 0.566 | 0.581 | 7.09 | -0.879 | 0.299 | 44.30 | -0.0831 | 1.69 |

## Per species, holdout

| species_name | sites_present | ba_r2 | trees_r2 | growth_r | growth_bias | mort_auc | mort_obs | mort_pred | reg_r | reg_obs | reg_pred |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Betula | 344 | 0.852 |  0.877 | 0.559 | -0.000570 | 0.606 | 0.219 | 0.3450 |  0.09310 | 19.9000 | 28.7000 |
| Pinus sylvestris | 341 | 0.914 |  0.870 | 0.733 | -0.000698 | 0.710 | 0.241 | 0.2620 |  0.07680 |  9.6700 |  3.2100 |
| Quercus robur/petraea | 305 | 0.888 |  0.794 | 0.624 | -0.001320 | 0.627 | 0.184 | 0.2020 | -0.05470 |  5.4800 |  6.9100 |
| Pseudotsuga menziesii | 114 | 0.831 |  0.575 | 0.570 | -0.005370 | 0.597 | 0.288 | 0.1830 |  0.07870 |  9.0300 |  5.5500 |
| Fagus sylvatica |  98 | 0.951 |  0.903 | 0.789 | -0.007260 | 0.754 | 0.132 | 0.0529 |  0.11000 |  3.0700 |  7.2800 |
| Larix kaempferi |  97 | 0.854 |  0.483 | 0.553 | -0.001450 | 0.628 | 0.282 | 0.2420 |  0.01670 |  7.3500 |  5.2700 |
| Quercus rubra |  88 | 0.830 |  0.782 | 0.635 | -0.001630 | 0.687 | 0.289 | 0.2100 |  0.18500 |  5.3400 |  5.3000 |
| Pinus nigra |  47 | 0.943 |  0.959 | 0.279 | -0.002520 | 0.686 | 0.205 | 0.2800 | -0.01640 |  0.4370 |  0.1180 |
| Picea abies |  42 | 0.909 |  0.818 | 0.597 |  0.001240 | 0.636 | 0.362 | 0.3270 |  0.04670 |  5.6700 |  2.1200 |
| Alnus glutinosa |  18 | 0.912 |  0.932 | 0.160 |  0.004850 | 0.821 | 0.158 | 0.2260 | -0.01870 |  0.3280 |  5.8900 |
| Fraxinus excelsior |  15 | 0.698 |  0.695 | 0.342 | -0.008370 | 0.458 | 0.224 | 0.3890 | -0.04350 |  0.7670 |  1.2700 |
| Acer |   9 | 0.803 | -1.660 | 0.594 | -0.019600 | 0.464 | 0.378 | 0.2720 |  0.06190 |  0.9020 | 29.4000 |
| Salix |   7 | 0.723 |  0.841 | 0.406 | -0.013500 | 0.298 | 0.333 | 0.4920 | -0.00431 |  0.5440 |  5.8500 |
| Populus |   6 | 0.603 |  0.882 | 0.815 | -0.001690 | 0.784 | 0.307 | 0.2140 |  0.02310 |  0.0189 |  0.0398 |

## Bare-ground equilibrium (stand totals, last 100 vs previous 100 years)

| label | ba_last100 | ba_prev100 | trees_last100 | change_pct |
|---|---|---|---|---|
| 9.0 degC, 825 mm, harvest 0 | 45.9 | 46.3 | 599 | -0.836 |
| 9.5 degC, 830 mm, harvest 0 | 45.8 | 47.9 | 694 | -4.400 |
| 10.1 degC, 769 mm, harvest 0 | 55.1 | 52.6 | 505 |  4.800 |
| 9.5 degC, 830 mm, harvest 1 | 27.7 | 28.0 | 473 | -1.160 |

## Figures

- figures/fig1_loss_epochs.png: loss per response over epochs
- figures/fig2_obs_pred.png: observed vs predicted, both sets
- figures/fig3_mortality_roc.png: stem-loss ROC
- figures/fig4_map.png: observed and predicted stand basal area, both sets
- figures/fig5_bare_ground.png: succession from bare ground at three climates (harvest 0) and at the median climate with harvest 1
- figures/fig6_parameter_trajectories.png, fig6b_parameter_change.png: parameter trajectories and settling

## Provenance

| step | script | inputs | outputs |
|---|---|---|---|
| raw data | calibration/01_download_nbi_raw.R | Probos NBI-7 public database archive (sha256 in PROVENANCE.md) | raw/probos/unz, PROVENANCE.md |
| extraction | calibration/02_extract_nbi_tables.R | the .accdb (mdbtools) | raw/csv/*.csv, extract_log.csv |
| preparation | calibration/v02/04_prepare_pairs_v2.R | raw CSV tables, WorldClim 2.1 30 s | prepared/pairs-v2/{national,train,holdout_blocked,pilot}, v02/PREPARATION.md, v02/expected_preparation_v2.csv |
| audit + pilots | calibration/AUDIT_v02.md, calibration/v02/pilot_*.R | pairs-v1 / pairs-v2 pilot sets, v0.1 fit | v02/pilot_table.csv |
| fit | calibration/v02/fit.R | pairs-v2/train | fits/pairs-v2-process/{model.pt, parameters.pt}, v02/loss_history.csv, v02/fit_summary.md, inst/model/weights.json |
| validation | calibration/v02/validation.R | model.pt, train, holdout_blocked | v02/report.md, v02/validation/*.csv (incl. metrics_strata.csv), v02/figures/*.png |

## Fit summary (from fit.R)

# Fit summary: pairs-v2-process (v0.2)

Fitted 2026-11-10 10:30 CET by calibration/v02/fit.R; FINN 0.2.0 (main 0848b60), torch 0.15.1, R 4.5.0; 1 threads; optimizer optim_ignite_adam; runtime 0.90 h for 600 epochs (5.4 s/epoch).

Data: pairs-v2/train, 1179 sites, 14 species (incl. other), intervals 4-12 yr. Weights: `auto` (baseline-scaled; no numeric vector, weight_scale `auto` moot).

Release assets: parameters.pt (named list of 34 tensors incl. env scaling), sha256 29c8a4c21ca7918c07bc15a21ca0e3c315f933eddc4781f9038790719f797443; model.pt (whole fitted module), sha256 d82572dff39ddaab6fb7fd1246c8f88da1e0d7aff66a28792eb6a9ee5228bf42. build_model(parameters.pt) reproduced the fitted model bit for bit on 200 sites.

Final losses (baseline-scaled, last epoch):

- dbh: 0.2832
- ba: 0.07508
- trees: 0.1898
- growth: 0.3469
- mortality: 1.089
- regeneration: 0.4347
- unused: 0
- epoch: 600

## PARAMS

```r
List of 24
 $ prep_version : chr "pairs-v2"
 $ train_set    : chr "train"
 $ holdout_set  : chr "holdout_blocked"
 $ patch_size_ha: num 0.126
 $ n_patches    : int 1
 $ env_vars     : chr [1:3] "tavg_c" "prec_mm" "harvest"
 $ mort_bounds  : num [1:2] -10 10
 $ freeze_theta : logi FALSE
 $ recruits_dbh : num 5.9
 $ hidden       : int [1:2] 20 20
 $ dropout      : num 0.1
 $ reg_K_init   : num 50
 $ reg_floor    : num 0.001
 $ epochs       : int 600
 $ batchsize    : int 150
 $ lr           : num 0.01
 $ lr_scheduler : chr "cosine"
 $ threads      : int 1
 $ weights      : chr "auto"
 $ weight_scale : chr "auto"
 $ loss_family  : Named chr [1:6] "mse" "mse" "poisson" "mse" ...
 $ checkpoints  : int 10
 $ device       : chr "cpu"
 $ seed         : int 1
```
