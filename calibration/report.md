# Validation report: pairs-v1-process

Written 2026-11-10 01:31 CET by calibration/validation.R; FINN 0.2.0 (main 0848b60), torch 0.15.1. Model: ../../FINN-zoo-data/nbi/fits/pairs-v1-process/model.pt; release object parameters.pt sha256 19d24a826394ff918f7d9cdfb11f79f2cb6b032df90466ef388df2abd3043a99.

## Sets

- TEST = pairs-v1/train: 1531 sites (plot x remeasurement pair) the model was fitted on; 713 MFV -> NBI-6, 818 NBI-6 -> NBI-7.
- HOLDOUT = pairs-v1/holdout_blocked: 740 sites in 25 km blocks the model never saw. No coordinates in the model, predicted as is.
- Aggregation as fit() scores it, per site over its own interval (4-12 yr): state at year = period_length; growth and stem loss = mean of the annual rates; recruitment = sum.

## Metrics

| set | variable | level | n | rmse | r2 | r | bias | obs_mean | auc |
|---|---|---|---|---|---|---|---|---|---|
| test | dbh (cm, species mean) | species x site |  5130 |   7.4700 |  0.7040 | 0.8590 | -2.34e+00 |  20.2000 | NA |
| test | basal area (m2 per 0.1257 ha patch) | species x site | 30600 |   0.2180 |  0.8950 | 0.9470 | -5.61e-03 |   0.1740 | NA |
| test | trees (per 0.1257 ha patch) | species x site | 30600 |  12.8000 |  0.5990 | 0.7970 |  1.89e+00 |   4.7700 | NA |
| test | growth (rel. dbh increment / yr) | species x site |  4160 |   0.0156 |  0.3570 | 0.6160 | -2.78e-03 |   0.0210 | NA |
| test | stem loss (annual rate, incl. felling) | species x site |  4480 |   0.2650 | -0.0639 | 0.0172 | -5.81e-02 |   0.0968 | 0.573 |
| test | recruitment (stems / ha / interval) | species x site | 30600 |  87.0000 | -0.9750 | 0.0694 |  1.75e+01 |   8.2600 | NA |
| holdout | dbh (cm, species mean) | species x site |  2540 |   6.8000 |  0.7290 | 0.8720 | -2.17e+00 |  19.3000 | NA |
| holdout | basal area (m2 per 0.1257 ha patch) | species x site | 14800 |   0.2220 |  0.8800 | 0.9380 | -2.99e-04 |   0.1640 | NA |
| holdout | trees (per 0.1257 ha patch) | species x site | 14800 |  12.8000 |  0.6110 | 0.8040 |  1.99e+00 |   4.9000 | NA |
| holdout | growth (rel. dbh increment / yr) | species x site |  2050 |   0.0155 |  0.3230 | 0.5790 | -2.03e-03 |   0.0208 | NA |
| holdout | stem loss (annual rate, incl. felling) | species x site |  2230 |   0.2770 | -0.0712 | 0.0223 | -6.59e-02 |   0.1050 | 0.556 |
| holdout | recruitment (stems / ha / interval) | species x site | 14800 |  84.1000 | -1.1800 | 0.0674 |  1.82e+01 |   7.8100 | NA |
| test | basal area (m2 per patch) | stand total |  1530 |   0.9390 |  0.5690 | 0.7780 | -1.12e-01 |   3.4800 | NA |
| test | trees (per patch) | stand total |  1530 |  64.8000 |  0.2130 | 0.7160 |  3.77e+01 |  95.4000 | NA |
| test | recruitment (stems / ha / interval) | stand total |  1530 | 469.0000 | -1.6300 | 0.1000 |  3.50e+02 | 165.0000 | NA |
| test | basal area (m2 per patch) | stand total, PERSISTENCE reference |  1530 |   0.8860 |  0.6160 | 0.8220 | -2.86e-01 |   3.4800 | NA |
| test | trees (per patch) | stand total, PERSISTENCE reference |  1530 |  53.5000 |  0.4630 | 0.7920 |  4.95e+00 |  95.4000 | NA |
| holdout | basal area (m2 per patch) | stand total |   740 |   0.9870 |  0.5030 | 0.7450 | -5.98e-03 |   3.2900 | NA |
| holdout | trees (per patch) | stand total |   740 |  65.8000 |  0.1380 | 0.7050 |  3.97e+01 |  97.9000 | NA |
| holdout | recruitment (stems / ha / interval) | stand total |   740 | 461.0000 | -2.1200 | 0.1040 |  3.65e+02 | 156.0000 | NA |
| holdout | basal area (m2 per patch) | stand total, PERSISTENCE reference |   740 |   0.9240 |  0.5640 | 0.7860 | -2.12e-01 |   3.2900 | NA |
| holdout | trees (per patch) | stand total, PERSISTENCE reference |   740 |  56.7000 |  0.3600 | 0.7670 |  7.94e+00 |  97.9000 | NA |

A negative R2 means worse than predicting the mean. The persistence rows predict the second-visit state by the first-visit state (model-free reference).
"Stem loss" is the fitted mortality target: natural death AND felling (see fit.R).

## Stem-loss AUC (interval loss probability against the tallied fates; each site x species row weighted by its trees at risk)

| set | n_trees_at_risk | n_lost | auc | obs_loss_share | pred_loss_share |
|---|---|---|---|---|---|
| test | 30900 | 6880 | 0.573 | 0.223 | 0.253 |
| holdout | 15200 | 3550 | 0.556 | 0.233 | 0.262 |

## Per species, holdout

| species_name | sites_present | ba_r2 | trees_r2 | growth_r | growth_bias | mort_auc | mort_obs | mort_pred | reg_r | reg_obs | reg_pred |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Betula | 431 | 0.8360 |  0.8530 | 0.5500 | -0.000387 | 0.598 | 0.223 | 0.37400 |  0.063200 | 24.2000 |  35.600 |
| Pinus sylvestris | 411 | 0.8860 |  0.8240 | 0.6890 | -0.001330 | 0.572 | 0.236 | 0.27700 |  0.054900 | 10.6000 |   3.710 |
| Quercus robur/petraea | 405 | 0.8950 |  0.7760 | 0.6280 | -0.001570 | 0.534 | 0.176 | 0.15400 | -0.039800 |  8.9500 |   8.250 |
| other | 185 | 0.7720 | -1.6600 | 0.3170 | -0.002490 | 0.470 | 0.259 | 0.22200 |  0.005770 | 23.7000 | 275.000 |
| Sorbus aucuparia | 142 | 0.7490 |  0.4800 | 0.0725 | -0.008160 | 0.527 | 0.207 | 0.03170 |  0.063200 | 12.4000 |  13.200 |
| Pseudotsuga menziesii | 137 | 0.8010 |  0.5120 | 0.5290 | -0.003070 | 0.522 | 0.292 | 0.29800 |  0.065400 | 10.6000 |   4.790 |
| Fagus sylvatica | 131 | 0.9360 |  0.8690 | 0.7110 | -0.008070 | 0.437 | 0.111 | 0.05180 |  0.062400 |  5.3000 |  19.400 |
| Quercus rubra | 124 | 0.8260 |  0.6150 | 0.5140 |  0.002510 | 0.567 | 0.286 | 0.35400 |  0.101000 |  6.1600 |   4.260 |
| Prunus serotina | 118 | 0.2330 | -0.3140 | 0.4250 | -0.014000 | 0.556 | 0.474 | 0.09450 |  0.049400 | 19.2000 |  18.100 |
| Rhamnus frangula | 116 | 0.0119 |  0.0915 | 0.3570 | -0.012600 | 0.559 | 0.399 | 0.00806 |  0.131000 | 10.3000 |   6.200 |
| Larix kaempferi | 107 | 0.7520 |  0.4680 | 0.6360 |  0.002880 | 0.572 | 0.283 | 0.25000 |  0.065400 |  7.7600 |   6.900 |
| Pinus nigra |  52 | 0.9320 |  0.9480 | 0.5110 | -0.001950 | 0.579 | 0.204 | 0.20700 |  0.012400 |  0.8600 |   0.485 |
| Picea abies |  51 | 0.8990 |  0.8350 | 0.6020 |  0.001630 | 0.596 | 0.341 | 0.32500 |  0.043100 |  4.6900 |   2.530 |
| Alnus glutinosa |  33 | 0.9350 |  0.8130 | 0.2790 |  0.004200 | 0.747 | 0.166 | 0.25300 |  0.006020 |  0.5620 |  21.900 |
| Fraxinus excelsior |  29 | 0.7230 |  0.7170 | 0.4970 | -0.003590 | 0.338 | 0.222 | 0.42700 | -0.006330 |  1.0800 |   2.670 |
| Acer |  23 | 0.7420 | -1.4900 | 0.5230 | -0.011900 | 0.317 | 0.281 | 0.48700 |  0.013400 |  2.7600 |  82.600 |
| Corylus avellana |  22 | 0.2770 |  0.3960 | 0.3720 | -0.007590 | 0.666 | 0.149 | 0.01450 |  0.023500 |  5.9500 |  12.800 |
| Populus |  12 | 0.4580 |  0.7380 | 0.4900 | -0.003210 | 0.567 | 0.337 | 0.35500 | -0.000368 |  0.0149 |   0.196 |
| Salix |   9 | 0.8670 |  0.6790 | 0.6030 | -0.006780 | 0.424 | 0.322 | 0.27400 |  0.040700 |  0.5480 |   1.480 |
| Robinia pseudoacacia |   4 | 0.9880 |  0.9930 | 0.5410 | -0.004220 | 0.273 | 0.200 | 0.09830 | -0.022800 |  0.5970 |   0.838 |

## Bare-ground equilibrium (stand totals, last 100 vs previous 100 years)

| label | ba_last100 | ba_prev100 | trees_last100 | change_pct |
|---|---|---|---|---|
| 9.0 degC, 827 mm | 127 | 124 | 6750 |  1.96 |
| 9.5 degC, 824 mm | 120 | 118 | 6560 |  1.37 |
| 10.1 degC, 766 mm | 120 | 126 | 5800 | -5.06 |

## Figures

- figures/fig1_loss_epochs.png: loss per response over epochs
- figures/fig2_obs_pred.png: observed vs predicted, both sets
- figures/fig3_mortality_roc.png: stem-loss ROC
- figures/fig4_map.png: observed and predicted stand basal area, both sets
- figures/fig5_bare_ground.png: succession from bare ground at three climates
- figures/fig6_parameter_trajectories.png, fig6b_parameter_change.png: parameter trajectories and settling

## Provenance

| step | script | inputs | outputs |
|---|---|---|---|
| raw data | calibration/01_download_nbi_raw.R | Probos NBI-7 public database archive (sha256 in PROVENANCE.md) | raw/probos/unz, PROVENANCE.md |
| extraction | calibration/02_extract_nbi_tables.R | the .accdb (mdbtools) | raw/csv/*.csv, extract_log.csv |
| preparation | calibration/03_prepare_pairs_v1.R | raw CSV tables, WorldClim 2.1 30 s | prepared/pairs-v1/{national,train,holdout_blocked,pilot}, PREPARATION.md, expected_preparation.csv |
| fit | calibration/fit.R | pairs-v1/train | fits/pairs-v1-process/{model.pt, parameters.pt}, loss_history.csv, fit_summary.md, inst/model/weights.json |
| validation | calibration/validation.R | model.pt, train, holdout_blocked | report.md, validation/*.csv, figures/*.png |

## Fit summary (from fit.R)

# Fit summary: pairs-v1-process

Fitted 2026-11-10 01:31 CET by calibration/fit.R; FINN 0.2.0 (main 0848b60), torch 0.15.1, R 4.5.0; 1 threads; runtime 1.45 h for 500 epochs (10.4 s/epoch).

Data: pairs-v1/train, 1531 sites, 20 species (incl. other), intervals 4-12 yr. Weights: `auto` (baseline-scaled; no numeric vector, weight_scale `auto` moot).

Release assets: parameters.pt (named list of 28 tensors incl. env scaling), sha256 19d24a826394ff918f7d9cdfb11f79f2cb6b032df90466ef388df2abd3043a99; model.pt (whole fitted module), sha256 60290b1b772268c1d9d9d138889b7953e82a7a0c956e34b3b280bd260c2ce25c. build_model(parameters.pt) reproduced the fitted model bit for bit on 200 sites.

Final losses (baseline-scaled, last epoch):

- dbh: 0.2537
- ba: 0.1821
- trees: 0.7566
- growth: 0.4126
- mortality: 0.9244
- regeneration: 0.529
- unused: 0
- epoch: 500

## PARAMS

```r
List of 21
 $ prep_version : chr "pairs-v1"
 $ train_set    : chr "train"
 $ holdout_set  : chr "holdout_blocked"
 $ patch_size_ha: num 0.126
 $ n_patches    : int 1
 $ env_vars     : chr [1:2] "tavg_c" "prec_mm"
 $ recruits_dbh : num 6.4
 $ hidden       : int [1:2] 20 20
 $ dropout      : num 0.1
 $ reg_K_init   : num 50
 $ reg_floor    : num 0.001
 $ epochs       : int 500
 $ batchsize    : int 150
 $ lr           : num 0.01
 $ threads      : int 1
 $ weights      : chr "auto"
 $ weight_scale : chr "auto"
 $ loss_family  : Named chr [1:6] "mse" "mse" "poisson" "mse" ...
 $ checkpoints  : int 10
 $ device       : chr "cpu"
 $ seed         : int 1
```
