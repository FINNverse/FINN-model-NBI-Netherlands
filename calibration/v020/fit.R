#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v020/fit.R
##
## WHAT   The zoo refit of the Dutch NFI model (NBI remeasurement pairs 4-12
##        years apart, data pairs-v2, 14 taxa) with the FINN 0.2.0 release
##        defaults ("v020-defaults", the set-up shared by every zoo model):
##        ONE joint fit() of competition, growth, mortality and regeneration
##        against all six responses (dbh, ba, trees, growth, mortality,
##        regeneration), loss_weights "auto", update_step 1. It replaces the two-stage v0.2 recipe
##        (calibration/v02/recipe, fit folder fits/pairs-v2-recipe, completed on
##        the cluster 2026-11-10: stage A with the observed recruits fed in,
##        stage B regeneration_saturation with frozen demography and a frozen
##        ML theta of 0.01-0.1 per species; predicted / observed recruitment on
##        train 4.99 after stage B). The released v0.1.0 and v0.2.0 files and
##        the recipe folder are untouched; this is a new fit folder.
## WHY    The zoo uses the 0.2.0 default set-up (decided on 2026-11-13):
##        every model of the zoo is refitted with it, then the zoo report is
##        updated. The dataset-specific inventory handling stays as prepared.
## FINN   Release candidate release-0.2.0 at b6244e5 (PR FINNverse/FINN#68),
##        installed in ~/R/lib-finn-release-b6244e5 (laptop and cluster).
##        Every 0.2.0 default is left at its default and recorded in
##        fit_summary.md from the built model and formals(): regeneration_capped
##        (cap "all", prior "max_obs_recruits" = max(reg / period_length) per ha
##        and year, prior_lambda 0.1), theta trained with theta_floor 0.5,
##        mortality_draw "binomial", regeneration_sampler "refined", mortality
##        limits 0.001 / 0.999, reg_floor 1e-3, trees_offset 0.01, rb_gradient
##        TRUE, light_mask TRUE, competition height bound log(2.3). No custom
##        floors, no gamma / poisson draw, no weight_scale "raw", no frozen
##        stage, no fed recruits.
## MODEL STRUCTURE (as the v0.2 recipe fit, the only change being the
##        regeneration process):
##   - U-shaped mortality `mortality_ushape` (sen in [0, sen_max] on the logit at
##     dbh 1 m; growth enters as growth / 0.01);
##   - growth and mortality env = module env_dnn_linmgmt: SELU network hidden
##     20-20, dropout 0.1, on the WorldClim normals (tavg_c, prec_mm) PLUS one
##     linear per-species weight on the harvest flag (0/1: harvest activity
##     recorded at the second visit) OUTSIDE the network, so forcing harvest = 0
##     or 1 shifts the growth / mortality logit by exactly beta_s (pure-shift
##     check below); the flag on its 0/1 scale (env_autoscale identity);
##   - static climate (WorldClim 2.1 normals, env rows 1..12 identical per pair);
##   - regeneration = FINN::regeneration_capped with a LINEAR env response on
##     tavg_c + prec_mm + harvest (the 0.2.0 default process instead of the
##     recipe's frozen regeneration_saturation);
##   - recruit_obs_weight = the median smaller circle of the training pairs
##     (10 m, 0.0314 ha), recruits_dbh 5.9 cm, one 0.1257 ha patch per plot,
##     growth_rate "annual", growth and mortality scored per tracked tree
##     (tree_obs_dt, period_mode per_site).
## INVENTORY PROBLEMS of the NBI data and how this fit handles them (all in the
##        prepared data pairs-v2, unchanged here; calibration/AUDIT_v02.md,
##        calibration/v02/PREPARATION.md):
##   1. Rare species / the old "other" class (audit 6, 22): v1's "other" was 38 %
##      shrubs and made 20k of the 23k holdout stem surplus (recruitment 2.3x).
##      pairs-v2 has NO lumped class and no shrubs (NBI group ST out of the
##      population): S = the 14 taxa with the largest national first-visit basal
##      area (genus-level NBI codes Quercus robur/petraea, Betula, Acer, Populus,
##      Salix kept as the NBI records them), a pair is kept when out-of-S
##      species hold < 2 % of its START-census basal area (1,762 of 2,271 pairs,
##      77.6 %), out-of-S trees dropped; the END-census share (oos_ba_share_t1)
##      is a validation stratum, never a filter.
##   2. Harvest covariate (audit 4, pilot P5/P8): harvest is recorded at every
##      second visit (NBI-7 oogstactiviteit, NBI-6 soort kap; felled share 22 %
##      vs < 1.5 % of the first-visit trees). `harvest` (0/1, 510 of 1,179
##      training pairs = 1) is the linear management term in growth and
##      mortality and a linear term of the regeneration env; the harvest
##      covariate was the one pilot change that moved the stem-loss AUC
##      (0.555 -> 0.635). The tool simulates with harvest = 0 (none) or 1.
##   3. Mortality target = total stem loss INCLUDING felled trees (BWI / Spain
##      convention; harvest and death are not separable per tree, censoring
##      biases): `died` = 1 in tree_obs_dt for dead and felled; censored = not
##      found, now < 5 cm, outside forest. Clearcuts (basal-area loss >= 50 %,
##      116 pairs) dropped in the preparation (audit 18). There is no felling
##      mark in tree_obs_dt, so validation.R has no natural-death AUC.
##   4. Stand level weak against persistence (v0.2.0: stand BA R2 0.635 vs
##      persistence 0.619; recipe 0.636 vs 0.619; stems +230/ha): the v0.2.0 card
##      calls the model a short-horizon stand projector. Handled here only
##      through the set-up (all six responses jointly, the recruitment cap);
##      validation.R reports the persistence null, the change over the interval
##      and the strata so that the comparison is explicit.
##   5. Variable interval 4-12 yr and pair type (audit 15): site = plot x
##      remeasurement pair (MFV -> NBI-6 7-12 yr, NBI-6 -> NBI-7 4-9 yr); obs at
##      year = period_length, env rows 1..12 per site; FINN's per-site
##      aggregation path; rates annualised in the preparation; recruitment
##      summed over the window, the cap prior is max(reg / period_length).
##   6. Variable-radius circle 5-20 m (audit 13): every tallied tree is expanded
##      to its per-ha weight on the 0.1257 ha patch (largest remainder per site x
##      species; a 5 m plot is inflated x16); sample sizes are tallied. Recruits
##      are counted inside the smaller circle of the pair, so one recruit stands
##      for 0.005-0.126 ha while FINN's recruit_obs_weight is one number: w =
##      the MEDIAN smaller-circle area of the training pairs (0.0314 ha, asserted
##      from site_table.csv); the share of recruiting rows whole on their own
##      circle and the quantiles of reg x w are reported.
##   7. Ingrowth vs recruitment (audit 2): 51 % of v1's "new" trees had a
##      back-modelled start dbh >= 5 cm; pairs-v2 keeps a new tree as a recruit
##      only if dbh - g75_species x interval < 5 cm (4,802 vs 9,079 recruits).
##   8. Recruit dispersion collapse (audit 7): v1's theta fell to 0.07 for 11
##      species; the recipe's per-species ML theta in count space was 0.01-0.1.
##      Here theta is trained above FINN's default floor 0.5 (the 0.2.0 guard
##      rail; it will bind for most species); the fitted values and the
##      predicted / observed recruitment are reported, and `recruitment` (FINN
##      help) warns that a floored theta does not calibrate the recruitment level.
##   9. Fate codes per cycle (audit 5): class 4 = "dead" in MFV / NBI-6 but
##      "suppressed, alive" in NBI-7; handled in the preparation.
##  10. MFV tree positions mostly unrecorded (audit 16): the "inside the circle"
##      filter cannot act on MFV stands; recruits of MFV -> NBI-6 pairs located
##      with the NBI-6 distance (no visible artefact).
##  11. Coordinates are the 1 km INSPIRE cell only: no coordinates as predictors;
##      holdout = 25 km blocks (nearest training pair median 7 km).
## SCHEDULE 300 epochs at lr 0.01, then a cosine decay to 1e-4 over 200 epochs
##        (500 in total), plain Adam, batch 100 (the recipe's), checkpoints
##        (parameter snapshot + model file) every 25 epochs; a 24 h cluster
##        segment, continued by fit_continue.R if needed. Convergence is judged
##        on the parameter drift between checkpoints (convergence.csv,
##        validation.R fig6), never on the loss.
## DATA   ../../FINN-zoo-data/nbi/prepared/pairs-v2/train (+ prepared/pairs-v2/site_table.csv)
##        (holdout_blocked for validation.R)
## RUN    Rscript calibration/v020/fit.R          (from the repo root, one thread)
##        cluster: bash hpc/submit_v020.sh (fit, fit_continue, validation, chained)
## OUTPUT FIT_DIR = ../../FINN-zoo-data/nbi/fits/pairs-v2-v020 (outside git):
##        model_init.pt, model.pt, checkpoints/epoch_<n>model.pt, loss_history.csv,
##        lr_schedule.csv, species_parameters.csv, regeneration_parameters.csv,
##        management_weights.csv, module_checks.csv, grad_norms.csv,
##        grad_norms_by_epoch.csv, cap_trajectory.csv, convergence.csv, fit_summary.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-release-b6244e5"), .libPaths()))   # FINN release-0.2.0 at b6244e5 (the ONE library of the v020 refits)
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v2-v020")

PARAMS <- list(
  ## ---- per model: what changes from data set to data set ----
  prep_version  = "pairs-v2",              # 14 taxa, no "other", no shrubs, start-census site rule x = 2 %, intervals 4-12 yr (MFV -> NBI-6, NBI-6 -> NBI-7)
  train_set     = "train",                 # 1,179 training pairs in the retained 25 km blocks
  holdout_set   = "holdout_blocked",       # not read here; recorded for validation.R
  env_vars      = c("tavg_c", "prec_mm", "harvest"),   # WorldClim 2.1 normals (bio1, bio12) + the harvest flag of the interval (0/1)
  linear_var    = "harvest",               # outside the network: one linear per-species weight; the climate network sees the other columns
  mgmt_levels   = c(0, 1),                 # the values the covariate takes; the pure-shift check forces each of them
  patch_size_ha = pi * 20^2 / 1e4,         # one FINN patch = the 20 m NBI circle (0.1257 ha); every tallied tree of the variable-radius circle expanded to its per-ha weight
  n_patches     = 1L,
  env_years     = 12L,                     # env rows 1..12 per pair (the longest interval); obs at year = period_length (4-12)
  recruits_dbh  = 5.9,                     # recruits enter just above the 5 cm tally threshold (median dbh of the pairs-v2 recruits)
  growth_rate   = "annual",                # `growth` in obs_dt and tree_obs_dt is the annualised relative dbh increment (PREPARATION.md)
  recruit_obs_weight_ha = pi * 10^2 / 1e4, # the area one observed recruit stands for: the MEDIAN smaller circle of the training pairs (10 m, 0.0314 ha; circles 5-20 m, FINN takes one number; header 6)
  w_tol         = 1e-4,                    # |median smaller-circle area of the training pairs - recruit_obs_weight_ha| allowed (asserted from site_table.csv)
  check_rows    = 200L,                    # env rows of the module checks (forward, gradients, pure shift)
  check_tol     = 1e-5,                    # |output change - beta * delta| / (1 + |output|) allowed in the pure-shift check (float32 rounding)
  ## ---- the v020 set-up (same in every model) ----
  hidden        = c(20L, 20L),             # climate network of growth and mortality (SELU), as the recipe
  dropout       = 0.1,                     # active while training, off in simulateForest (eval mode)
  beta_init     = 0,                       # management weight start
  sen_max       = 5,                       # U-shape cap: at most +5 on the mortality logit at dbh 1 m
  sen_init      = 0.05,                    # U-shape start: +0.05 on the logit at 1 m
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01 (FINN 0.2.0's own default, asserted)
  n_quantiles   = 10L,                     # FINN's default height-quantile competition, written out
  light_mask    = TRUE,                    # 0.2.0's masked competition (FINN's default, written out)
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # FINN's default families, all six responses scored
  loss_weights  = "auto",                  # each term / its intercept-only baseline (FINN's default)
  loss_recruit  = "simulated",             # the regeneration process recruits and is estimated with the rest (FINN's default; no fed recruits)
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # per tracked tree on the observed stand (the recipe); needs period_mode per_site
  clip_norm     = 2,                       # FINN's default per-group clip budget, written out
  update_step   = 1L,                      # FINN's default
  epochs_const  = 300L,                    # epochs at constant lr, then the cosine decay
  anneal_epochs = 200L,                    # epochs of the cosine decay after them
  lr            = 0.01,                    # the constant lr = the start of the cosine
  lr_min        = 1e-4,                    # the end of the cosine (1 % of the start)
  checkpoints   = 25L,                     # a parameter snapshot and a model file (resume point for fit_continue.R) every 25 epochs
  batchsize     = 100L,                    # the recipe's batch: 12 steps per epoch on 1,179 pairs
  threads       = 1L,                      # one torch thread (one CPU per fit job)
  device        = "cpu",
  seed          = 1L                       # the recipe's seed
)
EPOCHS <- PARAMS$epochs_const + PARAMS$anneal_epochs
torch_set_num_threads(PARAMS$threads)
if (file.exists(file.path(FIT_DIR, "model.pt"))) {
  message("model.pt exists: the fit has finished"); quit(save = "no")
}
dir.create(file.path(FIT_DIR, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
if (length(list.files(file.path(FIT_DIR, "checkpoints"), pattern = "model\\.pt$", recursive = TRUE)))
  stop("checkpoints exist already: resume with fit_continue.R instead of starting over")
## the installed FINN must be the release-0.2.0 API (defaults are read, not set)
stopifnot(packageVersion("FINN") == "0.2.0",
          identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale),
          all(c("loss_weights", "loss_recruit", "control", "rb_gradient", "trees_offset", "weight_scale") %in% names(formals(fit))),
          !any(c("weights", "recruits", "lr_scheduler", "optimizer", "clip_norm", "deaths") %in% names(formals(fit))),
          all(c("cap", "prior", "prior_lambda", "theta_floor", "optimizeTheta", "mortality_draw", "mortality_floor",
                "mortality_ceiling", "regeneration_sampler", "light_mask", "reg_floor") %in% names(formals(createProcess))),
          all(c("growth_rate", "recruit_obs_weight") %in% names(formals(finn))),
          is.function(FINN::regeneration_capped), is.function(FINN::finnControl))
DEFAULTS <- list(trees_offset = eval(formals(fit)$trees_offset), rb_gradient = eval(formals(fit)$rb_gradient))
PS    <- PARAMS$patch_size_ha
w_obs <- PARAMS$recruit_obs_weight_ha

## ---- the learning-rate schedule over the whole run --------------------------------
## lr_of(t), t = epochs completed; FINN steps the scheduler once per epoch, so epoch
## t + 1 runs at lr_of(t). torch::lr_lambda multiplies the base lr by lr_lambda(e),
## e = epochs stepped in THIS fit() call (fit_continue.R sets k_done).
lr_of <- function(t) ifelse(t < PARAMS$epochs_const, PARAMS$lr,
  PARAMS$lr_min + (PARAMS$lr - PARAMS$lr_min) * (1 + cos(pi * (t - PARAMS$epochs_const) / PARAMS$anneal_epochs)) / 2)
k_done <- 0L
schedule <- function(optimizer) torch::lr_lambda(optimizer, lr_lambda = function(e) lr_of(k_done + e) / PARAMS$lr)
SCHED <- data.table(epoch = seq_len(EPOCHS))[, lr := lr_of(epoch - 1L)]

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv"))[, .(siteID, patchID, species, dbh, trees, treeName)]
trees <- fread(file.path(train_dir, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
stopifnot(max(obs$species) == n_species, all(obs$year == obs$period_length), max(env$year) == PARAMS$env_years,
          all(env[, .N, by = siteID]$N == PARAMS$env_years), uniqueN(init$patchID) == PARAMS$n_patches,
          all(env[[PARAMS$linear_var]] %in% PARAMS$mgmt_levels), !anyNA(obs$reg), all(obs$reg >= 0))
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
## the recruit observation space: the pair's own smaller circle (site_table.csv: r0, r1 in m,
## area0_m2, area1_m2 = the forest part of the circle) against the one weight FINN takes
site_tab <- merge(fread(file.path(train_dir, "siteID_dt.csv"))[, .(siteID, siteName)],
                  fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "site_table.csv"))[, .(siteName, r0, r1, area0_m2, area1_m2)], by = "siteName")
site_tab[, a_min_ha := pmin(area0_m2, area1_m2) / 1e4]
stopifnot(nrow(site_tab) == uniqueN(obs$siteID), abs(median(site_tab$a_min_ha) - w_obs) < PARAMS$w_tol)
oo <- merge(obs[reg > 0, .(siteID, species, reg)], site_tab[, .(siteID, a_min_ha)], by = "siteID")
whole_share <- oo[, mean(abs(reg * a_min_ha - round(reg * a_min_ha)) < 1e-3)]
y_q <- quantile(obs[reg > 0, reg * w_obs], c(0, 0.25, 0.5, 0.75, 1))
K_prior_expected <- max(obs$reg / obs$period_length)   # what prior = "max_obs_recruits" resolves to (per ha per year)
mgmt_share <- 100 * mean(env[year == 1][[PARAMS$linear_var]])
message(sprintf("train: %d pairs, %d species, intervals %d..%d yr, %d tracked trees (%d lost); harvest on %.1f %% of the pairs; reg > 0 in %d of %d rows; max reg / period %.1f per ha per year",
                uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(trees), sum(trees$died), mgmt_share,
                sum(obs$reg > 0), nrow(obs), K_prior_expected))
message(sprintf("recruit observation space: reg x %.4f ha (median smaller circle; circles %d-%d m, area %.4f-%.4f ha); %.1f %% of the recruiting rows are whole counts on their own circle; reg x w on recruiting rows: min %.2f, quartiles %.2f / %.2f / %.2f, max %.1f",
                w_obs, min(site_tab[, pmin(r0, r1)]), max(site_tab[, pmin(r0, r1)]), min(site_tab$a_min_ha), max(site_tab$a_min_ha), 100 * whole_share, y_q[1], y_q[2], y_q[3], y_q[4], y_q[5]))
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
## scaling per column: climate "auto" (z-scaled), the harvest flag "identity" (beta = the shift of a harvested pair)
autoscale_spec <- c(setNames(rep("auto", length(setdiff(PARAMS$env_vars, PARAMS$linear_var))), setdiff(PARAMS$env_vars, PARAMS$linear_var)),
                    setNames("identity", PARAMS$linear_var))

## ---- 2. the env module of growth and mortality: MLP(climate) + linear management ----
## FINN's env design matrix is stats::model.matrix(formula, env): the intercept
## column first, then the formula's terms in order. The management column's
## position is read from it, not assumed. The module: input [..., n_inputs] (FINN
## passes [sites, years, n] or [sites, n]), output [..., n_species]. mlp = FINN's
## own network layout (linear, SELU, dropout after every hidden activation, linear
## output with bias) on every column except management (the constant intercept
## column stays in, as in FINN's networks); beta = one linear weight per species
## on the management column, added outside the network.
X_names <- colnames(model.matrix(env_formula, env[1:5, PARAMS$env_vars, with = FALSE]))
mgmt_col <- match(PARAMS$linear_var, X_names)
stopifnot(identical(X_names, c("(Intercept)", PARAMS$env_vars)), !is.na(mgmt_col))
env_dnn_linmgmt <- torch::nn_module("env_dnn_linmgmt",
  initialize = function(n_inputs, n_species, mgmt_col, hidden, dropout, beta_init) {
    self$mgmt_col  <- as.integer(mgmt_col)
    self$clim_cols <- setdiff(seq_len(n_inputs), self$mgmt_col)
    layers <- list(); n_in <- length(self$clim_cols)
    for (h in hidden) {
      layers <- c(layers, torch::nn_linear(n_in, h), torch::nn_selu())
      if (dropout > 0) layers <- c(layers, torch::nn_dropout(p = dropout))
      n_in <- h
    }
    self$mlp  <- do.call(torch::nn_sequential, c(layers, torch::nn_linear(n_in, n_species)))
    self$beta <- torch::nn_parameter(torch::torch_tensor(beta_init, dtype = torch::torch_float32()))
  },
  forward = function(x) {
    self$mlp(x[.., self$clim_cols]) + x[.., self$mgmt_col, drop = FALSE] * self$beta
  })
new_env_module <- function() env_dnn_linmgmt(n_inputs = length(X_names), n_species = n_species, mgmt_col = mgmt_col,
                                             hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                             beta_init = rep(PARAMS$beta_init, n_species))
## the checks on check_rows real design-matrix rows (climate z-scaled as FINN will
## see it): (1) forward pass in both call forms FINN uses; (2) gradients reach every
## MLP tensor and beta; (3) forcing management to each of its levels changes the
## output by exactly beta_s * delta on every row (eval mode = without dropout)
scaling_chk <- compute_env_scaling(env, autoscale_spec)
X_chk <- model.matrix(env_formula, apply_env_scaling(env, scaling_chk)[round(seq(1, .N, length.out = PARAMS$check_rows)), PARAMS$env_vars, with = FALSE])
stopifnot(nrow(X_chk) %% 20L == 0L)
check_module <- function(nn, label, grads = TRUE) {
  x <- torch::torch_tensor(unname(X_chk), dtype = torch::torch_float32())
  was_training <- nn$training
  out <- list(module = label, rows = nrow(X_chk))
  if (grads) {
    nn$train()
    y3 <- nn(x$view(c(20L, -1L, ncol(X_chk))))
    stopifnot(identical(as.integer(y3$shape), c(20L, nrow(X_chk) %/% 20L, n_species)))
    g <- torch::autograd_grad(nn(x)$pow(2)$sum(), nn$parameters)
    g <- setNames(vapply(g, function(t) as.numeric(t$abs()$sum()), numeric(1)), names(nn$parameters))
    stopifnot(all(is.finite(g)), all(g[names(g) != "beta"] > 0))
    out$grad_mlp_min_abs_sum <- min(g[names(g) != "beta"]); out$grad_beta_abs_sum <- g[["beta"]]
  }
  nn$eval()
  torch::with_no_grad({
    y2 <- as_array(nn(x))
    stopifnot(identical(dim(y2), c(nrow(X_chk), n_species)))
    beta <- as_array(nn$beta)
    dev <- 0
    for (t in PARAMS$mgmt_levels) {
      xt <- X_chk; xt[, mgmt_col] <- t
      yt <- as_array(nn(torch::torch_tensor(unname(xt), dtype = torch::torch_float32())))
      shift <- outer(xt[, mgmt_col] - X_chk[, mgmt_col], beta)
      dev <- max(dev, max(abs(yt - y2 - shift) / (1 + abs(y2))))
    }
  })
  if (was_training) nn$train()
  stopifnot(dev < PARAMS$check_tol)
  out$shift_max_rel_dev <- dev
  out$max_abs_output <- max(abs(y2))
  as.data.table(out)
}

## ---- 3. the model ------------------------------------------------------------------
## U-shaped mortality (the recipe's function; growth_scale and sen_max are buffers
## so that nothing depends on a formal's default)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / self$mort_growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}
sen_raw_init <- rep(qlogis(PARAMS$sen_init / PARAMS$sen_max), n_species)
build_model <- function() finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  growth_rate  = PARAMS$growth_rate,
  recruit_obs_weight = w_obs,              # the inventory observation operator (header)
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE,
                                       n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, NN = new_env_module(),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape, NN = new_env_module(),   # death draw and limits: FINN's defaults
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_capped,   # FINN's default process; cap, prior, sampler, theta, floor: defaults
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
## FINN builds the env networks lazily; a two-site, two-year bare-ground
## simulation builds them so that model_init.pt survives torch_load. The climate
## is z-scaled here as fit() will scale it (the model carries no scaling yet):
## on the raw normals the untrained linear regeneration response overflows.
build_networks <- function(m) {
  e <- apply_env_scaling(env, scaling_chk)[siteID %in% sort(unique(siteID))[1:2] & year <= 2L]
  e[, siteID := match(siteID, sort(unique(siteID)))]
  invisible(simulateForest(m, env = e, init_cohort = NULL, patches = 1L, patch_size = PS, device = "cpu"))
}
FINN.seed(PARAMS$seed)                     # species parameters and the networks are initialised here
model <- build_model()
build_networks(model)
stopifnot(inherits(model$nn_growth, "env_dnn_linmgmt"), inherits(model$nn_mortality, "env_dnn_linmgmt"),
          !any(c("mort_sen_max", "mort_growth_scale") %in% names(model$parameters)),
          all(c("par_theta_recruits_raw", "reg_logK") %in% names(model$parameters)),
          identical(model$reg_cap, "all"), model$reg_logK$numel() == 1L,
          identical(model$theta_floor, 0.5),   # FINN's default floor on theta (header)
          abs(model$recruit_obs_weight - w_obs) < 1e-12)
CHECKS <- rbind(cbind(stage = "init", check_module(model$nn_growth, "growth")),
                cbind(stage = "init", check_module(model$nn_mortality, "mortality")))
print(CHECKS)
## what the installed FINN made of the defaults (recorded in fit_summary.md)
USED <- list(mortality_draw = model$mortality_draw, mortality_floor = model$mortality_floor, mortality_ceiling = model$mortality_ceiling,
             theta_floor = model$theta_floor, theta_optimized = model$theta_optimized, regeneration_sampler = model$regeneration_sampler,
             sample_regeneration = model$sample_regeneration, reg_floor = model$reg_floor, cap = model$reg_cap, prior = model$reg_prior,
             prior_lambda = model$reg_prior_lambda, light_mask = model$light_mask, trees_offset = DEFAULTS$trees_offset, rb_gradient = DEFAULTS$rb_gradient)
message("FINN defaults in use: ", paste(names(USED), vapply(USED, function(x) paste(format(x), collapse = ","), ""), sep = " = ", collapse = "; "))
torch_save(model, file.path(FIT_DIR, "model_init.pt"))   # reg_logK holds FINN's placeholder until fit() resolves the prior
START <- list(theta = as.numeric(model$par_theta_recruits$detach()$cpu()))
model$period_mode <- "per_site"            # FINN scores loss_aggregation = "tree" on its per-site path only (the intervals vary anyway)

## ---- 4. the optimizer: plain Adam that logs the clipped gradient norms ------------
## FINN clips each parameter group (mechanistic, nn, loss) to clip_norm and then
## calls optimizer$step(); the subclass records, per step and parameter, the L2
## norm of the gradient it receives, and the cap K after the step (cap_trajectory.csv),
## then runs Adam's step unchanged.
GRADLOG <- new.env(); GRADLOG$rows <- list(); GRADLOG$K <- numeric(0)
optim_adam_gradlog <- torch::optimizer("optim_adam_gradlog", inherit = torch::optim_adam,
  initialize = function(params, lr, ...) {
    GRADLOG$names <- names(params)
    super$initialize(params, lr = lr, ...)
  },
  step = function(closure = NULL) {
    ps <- self$param_groups[[1]]$params
    GRADLOG$rows[[length(GRADLOG$rows) + 1L]] <- vapply(ps, function(p) tryCatch({
      g <- p$grad
      if (is.null(g) || !g$numel()) 0 else as.numeric(g$norm()$item())
    }, error = function(e) 0), numeric(1))
    out <- super$step(closure)
    GRADLOG$K <- c(GRADLOG$K, exp(as.numeric(ps[[match("reg_logK", GRADLOG$names)]]$detach()$cpu())))
    out
  })

## ---- 5. the joint fit ----------------------------------------------------------------
## nothing is meant to be fixed: every parameter must require gradient (header)
stopifnot(all(vapply(model$parameters, function(p) isTRUE(p$requires_grad), logical(1))))
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env              = env,
    data             = obs,                # all six responses; max(reg / period_length) resolves the cap's prior
    tree_data        = trees,
    loss_aggregation = PARAMS$loss_aggregation,
    init_cohort      = init_cohort,
    epochs           = EPOCHS,
    batchsize        = PARAMS$batchsize,
    update_step      = PARAMS$update_step,
    lr               = PARAMS$lr,          # the base lr the schedule multiplies
    loss_weights     = PARAMS$loss_weights,
    loss_family      = PARAMS$loss_family,
    loss_recruit     = PARAMS$loss_recruit,
    patches          = PARAMS$n_patches,
    patch_size       = PS,
    control          = finnControl(optimizer = optim_adam_gradlog,   # plain Adam + the gradient-norm log
                                   lr_scheduler = schedule,          # constant, then cosine; stepped once per epoch by FINN
                                   clip_norm = PARAMS$clip_norm),
    device           = PARAMS$device,
    env_autoscale    = autoscale_spec,
    folder           = file.path(FIT_DIR, "checkpoints"),
    checkpoints      = PARAMS$checkpoints,
    plot_progress    = FALSE)              # trees_offset, rb_gradient, weight_scale: FINN's defaults (header)
min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lr_end <- model$optimizer$param_groups[[1]]$lr
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history))[, 1:6]
setnames(hist, c("dbh", "ba", "trees", "growth", "mort", "reg"))
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
fwrite(SCHED, file.path(FIT_DIR, "lr_schedule.csv"))
K_prior <- model$reg_prior_value
K_fit   <- exp(as.numeric(model$reg_logK$detach()$cpu()))
stopifnot(nrow(hist) == EPOCHS, abs(lr_end - lr_of(EPOCHS)) < 1e-9,
          abs(K_prior - K_prior_expected) < 1e-6 * K_prior_expected,   # max_obs_recruits = max(reg / period_length), per ha per year
          abs(as.numeric(model$mort_sen_max$cpu()) - PARAMS$sen_max) < 1e-6,
          abs(as.numeric(model$mort_growth_scale$cpu()) - PARAMS$growth_scale) < 1e-9)

## ---- 6. diagnostics: gradient norms and the cap ------------------------------------
G <- as.data.table(do.call(rbind, GRADLOG$rows)); setnames(G, GRADLOG$names)
G[, step := .I]
steps_per_epoch <- nrow(G) / EPOCHS
G[, epoch := as.integer(ceiling(step / steps_per_epoch))]
GL <- melt(G, id.vars = c("step", "epoch"), variable.name = "parameter", value.name = "norm")
GL[, group := fifelse(grepl("^par_loss_", parameter), "loss", fifelse(grepl("^nn_", parameter), "nn", "mechanistic"))]
GL[, regeneration := grepl("regeneration|theta_recruits|reg_logK", parameter)]
GS <- GL[, .(group_norm = sqrt(sum(norm^2)), reg_share = sum(norm[regeneration]^2) / max(sum(norm^2), 1e-30)),
         by = .(step, epoch, group)]
GS[, at_budget := group_norm > 0.999 * PARAMS$clip_norm]
fwrite(GL, file.path(FIT_DIR, "grad_norms.csv"))
GE <- GS[, .(group_norm_mean = mean(group_norm), share_clipped = mean(at_budget), reg_share_mean = mean(reg_share)),
         by = .(epoch, group)]
fwrite(GE, file.path(FIT_DIR, "grad_norms_by_epoch.csv"))
KT <- data.table(step = seq_along(GRADLOG$K), K = GRADLOG$K)[, epoch := as.integer(ceiling(step / steps_per_epoch))][
  , .(K_end = K[.N], K_min = min(K), K_max = max(K)), by = epoch][, K_prior := K_prior]
fwrite(KT, file.path(FIT_DIR, "cap_trajectory.csv"))

## ---- 7. checks and tables -------------------------------------------------------------
## the trained modules reloaded from model.pt: still the module class, still a pure shift
chk <- torch_load(file.path(FIT_DIR, "model.pt"))
stopifnot(inherits(chk$nn_growth, "env_dnn_linmgmt"), inherits(chk$nn_mortality, "env_dnn_linmgmt"))
CHECKS <- rbind(CHECKS, cbind(stage = "trained, reloaded", check_module(chk$nn_growth, "growth", grads = FALSE)),
                cbind(stage = "trained, reloaded", check_module(chk$nn_mortality, "mortality", grads = FALSE)), fill = TRUE)
fwrite(CHECKS, file.path(FIT_DIR, "module_checks.csv"))
rm(chk)
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(model[[block]]$detach()$cpu())), nm)
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")),
             spar("par_regeneration", "reg_light"))
SPP[, `:=`(sen = PARAMS$sen_max * plogis(as.numeric(model$mort_sen_raw$detach()$cpu())),
           theta_start = START$theta, theta_fit = as.numeric(model$par_theta_recruits$detach()$cpu()))]
SPP[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & mort_size < 0, -mort_size / (2 * sen) * 100, Inf)]
fwrite(SPP, file.path(FIT_DIR, "species_parameters.csv"))
MGMT <- merge(data.table(species = seq_len(n_species), growth_beta = as_array(model$nn_growth$beta$detach()),
                         mortality_beta = as_array(model$nn_mortality$beta$detach())), species_dt, by = "species")
fwrite(MGMT, file.path(FIT_DIR, "management_weights.csv"))   # the shift of log growth / the mortality logit, harvest 0 -> 1
## one simulation of train over fit()'s windows (observation at year p <- sim years
## 1..p): predicted / observed recruitment, and the vanished stem predictions
## (simulated 0, or below FINN's trees_offset, against a positive count)
sim_check <- function(m) {
  m$eval()
  FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = init_cohort, patches = PARAMS$n_patches, patch_size = PS,
                                    device = "cpu", batchsize = PARAMS$batchsize)$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  W <- unique(obs[, .(siteID, obs_year = year, period_length)])[
    , .(year = (obs_year - period_length + 1L):obs_year), by = .(siteID, obs_year)]
  R <- merge(W, S[, .(siteID, year, species, r_mean_ha)], by = c("siteID", "year"))[, .(pred = sum(r_mean_ha)), by = species]
  R <- merge(R, obs[, .(obs = sum(reg)), by = species], by = "species")
  Z <- merge(obs[trees > 0, .(siteID, year, species)], S[, .(siteID, year, species, trees_m = trees)],
             by = c("siteID", "year", "species"), all.x = TRUE)[is.na(trees_m), trees_m := 0]
  list(total = sum(R$pred) / sum(R$obs), species = R[, .(species, ratio = pred / obs)],
       zero = Z[, .(n_pos = .N, n_zero = sum(trees_m == 0), n_below_offset = sum(trees_m < DEFAULTS$trees_offset))])
}
m_start <- torch_load(file.path(FIT_DIR, "model_init.pt"))
invisible(torch::with_no_grad(m_start$reg_logK$fill_(log(K_prior))))   # the start fit() used: K = the resolved prior
m_start$env_scaling <- model$env_scaling   # and the env scaling fit() learned from the same data before its first step
ck_start <- sim_check(m_start); rm(m_start)
ck_end   <- sim_check(model)
REGT <- merge(merge(SPP[, .(species, species_name, reg_light, theta_start, theta_fit)],
                    ck_start$species[, .(species, rec_ratio_start = ratio)], by = "species"),
              ck_end$species[, .(species, rec_ratio_fit = ratio)], by = "species")
REGT <- merge(REGT, obs[, .(obs_max_reg_annual = max(reg / period_length)), by = species], by = "species")
REGT[, `:=`(K_fit = K_fit, K_prior = K_prior, K_over_obs_max_annual = K_fit / obs_max_reg_annual)]
fwrite(REGT, file.path(FIT_DIR, "regeneration_parameters.csv"))
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence.csv"))

## ---- 8. what was fitted --------------------------------------------------------
commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                                    error = function(e) character(0)))
if (!length(commit)) commit <- "unknown (no git checkout)"
sum_l <- function(rows) { rows <- unique(pmin(pmax(rows, 1L), nrow(hist)))
  round(colMeans(as.matrix(hist[rows, .(dbh, ba, trees, growth, mort, reg)])), 4) }
fmt <- function(v) paste(sprintf("%s %.4g", names(v), v), collapse = ", ")
writeLines(c(
  sprintf("# FINN-model-NBI-Netherlands, v020-defaults joint fit (%s)", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v020/fit.R` (FINN %s from library `%s`, release-0.2.0; torch %s; repository %s): %.1f min for %d epochs (%.1f s/epoch), update_step %d, %d thread(s) on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), packageVersion("FINN"), basename(dirname(find.package("FINN"))),
          packageVersion("torch"), commit, min_fit, EPOCHS, 60 * min_fit / EPOCHS, PARAMS$update_step, PARAMS$threads, PARAMS$device), "",
  sprintf("FINN defaults in use (none passed): %s.", paste(names(USED), vapply(USED, function(x) paste(format(x), collapse = ","), ""), sep = " = ", collapse = "; ")), "",
  sprintf("Data: `prepared/%s/%s` (%d pairs, %d species, intervals %d..%d yr, %d tracked trees of which %d lost (dead or felled); harvest recorded on %.1f %% of the pairs). Data convention kept: total stem loss incl. felling is the loss target; not found / now < 5 cm / outside forest censored. All six responses scored; recruitment through recruit_obs_weight = %.4f ha (the median smaller circle of the training pairs; circles %d-%d m; %.1f %% of the recruiting rows whole on their own circle; reg x w quartiles %.2f / %.2f / %.2f). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(trees), sum(trees$died), mgmt_share,
          w_obs, min(site_tab[, pmin(r0, r1)]), max(site_tab[, pmin(r0, r1)]), 100 * whole_share, y_q[2], y_q[3], y_q[4], PARAMS$holdout_set), "",
  sprintf("Plain joint fit; FINN default initialisation (seed %d), theta trained above FINN's default floor %g, no warm start, no L2 other than the cap prior, no fed recruits. Growth and mortality env `%s` through env_dnn_linmgmt (hidden %s, dropout %g, + linear `%s` on its 0/1 scale, start %g, in growth AND mortality as the recipe); U-shaped mortality (sen in [0, %g], start %g at 1 m; growth / %g); light_mask %s, %d quantiles; regeneration FINN::regeneration_capped with a linear env response on the same formula; growth and mortality loss_aggregation = tree (scored on the observed tracked trees, period_mode per_site); growth_rate %s; recruits_dbh %g; %d patch of %.4f ha.",
          PARAMS$seed, model$theta_floor, paste(deparse(env_formula), collapse = ""), paste(PARAMS$hidden, collapse = "-"), PARAMS$dropout,
          PARAMS$linear_var, PARAMS$beta_init, PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, PARAMS$light_mask, PARAMS$n_quantiles,
          PARAMS$growth_rate, PARAMS$recruits_dbh, PARAMS$n_patches, PS), "",
  sprintf("Module checks (`module_checks.csv`, %d rows): forward pass in both call forms ok; gradients reach every MLP tensor and beta at the start; forcing %s to each of %s moves the output by beta_s x delta with max relative deviation %.2g (start) and %.2g (trained, reloaded) = a pure shift.",
          nrow(X_chk), PARAMS$linear_var, paste(PARAMS$mgmt_levels, collapse = "/"), max(CHECKS[stage == "init", shift_max_rel_dev]), max(CHECKS[stage != "init", shift_max_rel_dev])), "",
  sprintf("Management weights (harvest 0 -> 1 shift; `management_weights.csv`): log growth %s; mortality logit %s.",
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$growth_beta), collapse = ", "),
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$mortality_beta), collapse = ", ")), "",
  sprintf("Recruitment cap: prior (max_obs_recruits, per ha per year) %.1f = start K; fitted K %.1f per ha per year (K / prior %.3f).", K_prior, K_fit, K_fit / K_prior), "",
  sprintf("Weights \"auto\": baselines %s; weights %s.", fmt(model$loss_baseline), fmt(model$loss_weights)), "",
  sprintf("Schedule (`lr_schedule.csv`): epochs 1-%d at lr %g, epochs %d-%d cosine %g -> %g; the optimizer ended at lr %.3g. Plain Adam (torch::optim_adam), clip_norm %g, batch %d (%d steps per epoch), checkpoints every %d epochs.",
          PARAMS$epochs_const, PARAMS$lr, PARAMS$epochs_const + 1L, EPOCHS, PARAMS$lr, PARAMS$lr_min, lr_end,
          PARAMS$clip_norm, PARAMS$batchsize, as.integer(round(steps_per_epoch)), PARAMS$checkpoints), "",
  sprintf("Baseline-scaled losses (mean of 10 epochs): start %s; end %s.", fmt(sum_l(1:10)), fmt(sum_l((nrow(hist) - 9):nrow(hist)))), "",
  sprintf("Predicted / observed recruitment on train over the fit windows: %.2f at the start (K = prior), %.2f after the fit (`regeneration_parameters.csv`). Theta: start %s, fit %s (floor %g; the recipe's ML values were 0.01-0.1, so the floor is expected to bind).",
          ck_start$total, ck_end$total, paste(signif(range(START$theta), 3), collapse = "-"), paste(signif(range(SPP$theta_fit), 3), collapse = "-"), model$theta_floor), "",
  sprintf("Vanished stem predictions on train (positive observed count, simulated 0 / below trees_offset %g): start %d / %d, end %d / %d of %d.",
          DEFAULTS$trees_offset, ck_start$zero$n_zero, ck_start$zero$n_below_offset, ck_end$zero$n_zero, ck_end$zero$n_below_offset, ck_end$zero$n_pos), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm): %s (`species_parameters.csv`).",
          paste(sprintf("%.3f", range(SPP$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SPP$species_name, SPP$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Clipped gradients (`grad_norms_by_epoch.csv`, last 10 epochs): %s.",
          paste(GE[epoch > EPOCHS - 10, .(s = sprintf("%s clipped %.0f %% of steps, regeneration share %.2f", group[1], 100 * mean(share_clipped), mean(reg_share_mean))), by = group]$s, collapse = "; ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs; `convergence.csv`): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Validation report: `Rscript calibration/v020/validation.R`. Convergence: judge by the parameter drift between checkpoints, not by the loss."
), file.path(FIT_DIR, "fit_summary.md"))
message(sprintf("done: %.1f min (%.1f s/epoch); final losses: %s", min_fit, 60 * min_fit / EPOCHS,
                paste(sprintf("%s=%.3f", names(hist)[1:6], as.numeric(hist[.N, 1:6])), collapse = " ")))
