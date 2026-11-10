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
