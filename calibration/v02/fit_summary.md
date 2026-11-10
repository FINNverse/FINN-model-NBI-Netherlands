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
