#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v02/recipe/fit_stageA.R
##
## WHAT   Stage A of the zoo v0.2 recipe for the Dutch NFI (NBI remeasurement
##        pairs 4-12 years apart, data pairs-v2, 14 taxa): competition, growth
##        and mortality are fitted with the OBSERVED recruitment fed into the
##        simulation; the regeneration process is fitted alone in stage B
##        (fit_stageB.R). The management variant of the recipe
##        (FINN-model-NFI-Sweden/calibration/v02/recipe/fit_stageA.R, the
##        template): the environment response of growth and mortality is
##          env(x) = MLP(climate) + beta_s * harvest
##        a SELU network on the WorldClim normals PLUS one linear per-species
##        weight on the harvest flag (0/1: harvest activity recorded at the
##        second visit) OUTSIDE the network (custom torch module
##        `env_dnn_linmgmt` through createProcess(NN = )), so the network cannot
##        learn a harvest x climate interaction and forcing harvest = 0 or 1 in a
##        simulation shifts the growth / mortality logit by exactly beta_s. The
##        flag keeps its 0/1 scale (env_autoscale "identity").
##        DIFFERENCE TO THE v0.2.0 FIT (calibration/v02/fit.R, harvest in the
##        mortality network only): the recipe puts the management term into
##        growth AND mortality (the module is the env response of both; the
##        Sweden template does the same), and linearly into the regeneration env
##        of stage B. Stated here because the question was open.
## WHY    The recipe (FINN-zoo-data/RECIPE_v02.md, decided 2026-11-10): U-shaped
##        mortality, observed recruitment fed as whole trees, all losses, growth
##        and mortality per tracked tree, weights "auto", cosine lr, one thread,
##        light_mask TRUE, default floor and gamma draw, FINN 0.2.0. The NBI
##        pilots (calibration/AUDIT_v02.md section 24) found that theta frozen
##        at 1 hurt (P3); the recipe freezes theta at its per-species ML value
##        in count space instead (stage B), which is a different thing, and the
##        stage-B summary reports the vector. The refit under the recipe is a
##        new fit folder (fits/pairs-v2-recipe); the released v0.1.0 and the
##        v0.2.0 files are untouched.
## DATA CONVENTION (pairs-v2, kept as prepared): the mortality target is total
##        stem loss incl. felled trees (`died` = 1 in tree_obs_dt for dead and
##        felled; censored = not found, now < 5 cm, outside forest); `harvest`
##        is the covariate that separates plots with recorded harvest from the
##        others, so the tool simulates with harvest = 0 (none) or 1.
## DATA   ../../FINN-zoo-data/nbi/prepared/pairs-v2/train
## NEEDS  FINN 0.2.0 (FINNverse/FINN main, commit 0848b60) in ~/R/lib-finn-0848b60, torch
## RUN    Rscript calibration/v02/recipe/fit_stageA.R      (from the repo root)
##        cluster: hpc/submit_recipe.sh (stage A, then stage B and validation chained)
## OUTPUT <FIT_DIR>/model_stageA.pt, loss_history_stageA.csv, species_parameters_stageA.csv,
##        mortality_size_curve_stageA.csv, management_weights_stageA.csv, module_checks_stageA.csv,
##        convergence_stageA.csv, fit_summary_stageA.md
##        FIT_DIR = ../../FINN-zoo-data/nbi/fits/pairs-v2-recipe (outside git)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60 (as every script of this repository)
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v2-recipe")
dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  ## ---- per model: what changes from data set to data set ----
  prep_version  = "pairs-v2",              # 14 taxa, no "other", no shrubs, start-census site rule x = 2 %, intervals 4-12 yr (MFV -> NBI-6, NBI-6 -> NBI-7)
  train_set     = "train",                 # 1,179 training pairs in the retained 25 km blocks
  env_vars      = c("tavg_c", "prec_mm", "harvest"),   # WorldClim 2.1 normals (bio1, bio12) + the harvest flag of the interval (0/1)
  linear_var    = "harvest",               # outside the network: one linear per-species weight; the climate network sees the other columns
  mgmt_levels   = c(0, 1),                 # the values the covariate takes; the pure-shift check forces each of them
  patch_size_ha = pi * 20^2 / 1e4,         # one FINN patch = the 20 m NBI circle (0.1257 ha); every tallied tree of the variable-radius circle expanded to its per-ha weight
  n_patches     = 1L,
  recruits_dbh  = 5.9,                     # recruits enter just above the 5 cm tally threshold (the v0.2.0 value: median dbh of the pairs-v2 recruits)
  check_rows    = 200L,                    # env rows of the module checks (forward, gradients, pure shift)
  check_tol     = 1e-5,                    # |output change - beta * delta| / (1 + |output|) allowed in the pure-shift check (float32 rounding)
  ## ---- the recipe (same in every model) ----
  hidden        = c(20L, 20L),             # climate network of growth and mortality (SELU)
  dropout       = 0.1,                     # active while training, off in simulateForest (eval mode)
  sen_max       = 5,                       # U-shape cap: at most +5 on the mortality logit at dbh 1 m
  sen_init      = 0.05,                    # U-shape start: +0.05 on the logit at 1 m
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01 (FINN 0.2.0's default)
  reg_floor     = 1e-3,                    # FINN's additive recruitment floor; unused by the observed recruitment, kept for stage B
  reg_obs_eps   = 1e-8,                    # recruits/ha/yr standing in for a year without observed recruits (log scale)
  light_mask    = TRUE,
  n_quantiles   = 10L,
  epochs        = 300L,                    # cosine decay to lr/100 over the epochs
  batchsize     = 100L,
  lr            = 0.005,
  lr_scheduler  = "cosine",
  weights       = "auto",                  # each loss scaled by its intercept-only baseline
  weight_scale  = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "none"),   # all state and rate terms; recruitment is observed here
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # per tracked tree; needs period_mode per_site (intervals 4-12 yr)
  checkpoints   = 10L,
  threads       = 1L,
  device        = "cpu",
  seed          = 1L
)
torch_set_num_threads(PARAMS$threads)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))
PS <- PARAMS$patch_size_ha

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv"))[, .(siteID, patchID, species, dbh, trees, treeName)]
trees <- fread(file.path(train_dir, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
stopifnot(max(obs$species) == n_species, all(obs$year == obs$period_length), uniqueN(init$patchID) == PARAMS$n_patches,
          all(env[[PARAMS$linear_var]] %in% PARAMS$mgmt_levels))
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
Tmax <- max(env$year)
message(sprintf("train: %d sites, %d species, intervals %d..%d yr, Tmax %d, %d tracked trees (%d lost), %.1f %% plots with recorded harvest",
                uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), Tmax, nrow(trees), sum(trees$died),
                100 * mean(env[year == 1][[PARAMS$linear_var]])))

## The observed recruitment as env columns reg_obs_<s>: whole trees per patch on the
## env rows 1..p of the interval, as evenly as integers allow (see the Spain reference)
reg_cols <- paste0("reg_obs_", seq_len(n_species))
sched <- obs[, .(siteID, species, p = period_length, n = round(pmax(reg, 0) * PS))]
sched <- sched[, .(year = seq_len(Tmax)), by = .(siteID, species, p, n)]
sched[, count := fifelse(year <= p, floor(n * year / p) - floor(n * (year - 1) / p), 0)]
stopifnot(sched[, sum(count), by = .(siteID, species)][, all(V1 == sched[year == 1, n])])
sched[, col := paste0("reg_obs_", species)]
REG <- dcast(sched, siteID + year ~ col, value.var = "count")
setcolorder(REG, c("siteID", "year", reg_cols))
for (cc in reg_cols) set(REG, j = cc, value = log(REG[[cc]] / PS + PARAMS$reg_obs_eps))
env_obsreg <- merge(env, REG, by = c("siteID", "year"))
setorder(env_obsreg, siteID, year)
stopifnot(nrow(env_obsreg) == nrow(env))
n_fed <- sched[year == 1, sum(n)]
message(sprintf("observed recruits fed: %d trees on %d plots (%.1f %% of the site x species cells recruit)",
                n_fed, sched[year == 1 & n > 0, uniqueN(siteID)], 100 * mean(obs$reg > 0)))

## ---- 2. the custom process functions (the recipe's code blocks) --------------------
## U-shaped mortality (identical to the Spain reference)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}
sen_raw_init <- rep(qlogis(PARAMS$sen_init / PARAMS$sen_max), n_species)
## observed recruitment (identical to the Spain reference)
regeneration_observed <- function(species, parReg, pred, light) {
  torch_exp(pred)[, NULL]$`repeat`(c(1, species$shape[2], 1))
}
obsreg_formula <- as.formula(paste("~ 0 +", paste(reg_cols, collapse = " + ")))
env_formula    <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
## scaling per column: climate "auto" (z-scaled), the harvest flag and the
## recruitment columns "identity" (so beta is the shift of a harvested plot)
autoscale_spec <- c(setNames(rep("auto", length(setdiff(PARAMS$env_vars, PARAMS$linear_var))), setdiff(PARAMS$env_vars, PARAMS$linear_var)),
                    setNames("identity", PARAMS$linear_var), setNames(rep("identity", n_species), reg_cols))

## ---- 2b. the env module: MLP(climate) + linear management --------------------------
## FINN's env design matrix is stats::model.matrix(formula, env): the intercept
## column first, then the formula's terms in order. The management column's
## position is read from it, not assumed. The module: input [..., n_inputs] (FINN
## passes [sites, years, n] or [sites, n]), output [..., n_species]. mlp = FINN's
## own network layout (build_NN: linear, SELU, dropout after every hidden
## activation, linear output with bias) on every column except management (the
## constant intercept column stays in, as in FINN's networks); beta = one linear
## weight per species on the management column, added outside the network. The
## column indices are plain R fields (kept by torch_save with the module).
X_names <- colnames(model.matrix(env_formula, env[1:5, PARAMS$env_vars, with = FALSE]))
mgmt_col <- match(PARAMS$linear_var, X_names)
stopifnot(identical(X_names, c("(Intercept)", PARAMS$env_vars)), !is.na(mgmt_col))
message(sprintf("env design matrix: %s; management = column %d", paste(X_names, collapse = ", "), mgmt_col))
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
                                             hidden = PARAMS$hidden, dropout = PARAMS$dropout, beta_init = rep(0, n_species))
## the checks on check_rows real design-matrix rows (climate z-scaled as FINN will
## see it): (1) forward pass in both call forms FINN uses; (2) gradients reach every
## MLP tensor and beta; (3) forcing management to each of its levels changes the
## output by exactly beta_s * delta on every row (eval mode = without dropout)
scaling_chk <- compute_env_scaling(env, autoscale_spec[names(autoscale_spec) %in% PARAMS$env_vars])
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

## ---- 3. the training model, written out ---------------------------------------
FINN.seed(PARAMS$seed)                     # the modules' MLPs are initialised here
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE,
                                       n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth,
                                       NN = new_env_module(),   # section 2b; FINN uses it in place of its own network
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape,
                                       NN = new_env_module(),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE))),
  regeneration_process = createProcess(obsreg_formula, regeneration_observed,
                                       initEnv = list(diag(n_species)),
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       sample_regeneration = FALSE)
)
model$period_mode <- "per_site"            # the per-tree losses are scored on FINN's per-site path only
stopifnot(inherits(model$process_growth$NN, "env_dnn_linmgmt"), inherits(model$process_mortality$NN, "env_dnn_linmgmt"))
CHECKS <- rbind(cbind(stage = "init", check_module(model$process_growth$NN, "growth")),
                cbind(stage = "init", check_module(model$process_mortality$NN, "mortality")))
print(CHECKS)

## ---- 4. fit ----------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env           = env_obsreg,
    data          = obs,
    tree_data     = trees,
    init_cohort   = init_cohort,
    epochs        = PARAMS$epochs,
    batchsize     = PARAMS$batchsize,
    lr            = PARAMS$lr,
    lr_scheduler  = PARAMS$lr_scheduler,
    optimizer     = torch::optim_adam,
    weights       = PARAMS$weights,
    weight_scale  = PARAMS$weight_scale,
    loss_family   = PARAMS$loss_family,
    loss_aggregation = PARAMS$loss_aggregation,
    patches       = PARAMS$n_patches,
    patch_size    = PS,
    env_autoscale = autoscale_spec,
    device        = PARAMS$device,
    folder        = NULL,
    checkpoints   = PARAMS$checkpoints,
    plot_progress = FALSE)
fit_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## save first, summarise later
torch_save(model, file.path(FIT_DIR, "model_stageA.pt"))
loss_hist <- as.data.table(do.call(rbind, model$history))
setnames(loss_hist, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(loss_hist))])
loss_hist[, epoch := .I]
fwrite(loss_hist, file.path(FIT_DIR, "loss_history_stageA.csv"))

## ---- 5. checks and what was fitted ---------------------------------------------------
## the trained modules reloaded from model_stageA.pt: still the module class, still a pure shift
chk <- torch_load(file.path(FIT_DIR, "model_stageA.pt"))
stopifnot(inherits(chk$nn_growth, "env_dnn_linmgmt"), inherits(chk$nn_mortality, "env_dnn_linmgmt"))
CHECKS <- rbind(CHECKS, cbind(stage = "trained, reloaded", check_module(chk$nn_growth, "growth", grads = FALSE)),
                cbind(stage = "trained, reloaded", check_module(chk$nn_mortality, "mortality", grads = FALSE)), fill = TRUE)
fwrite(CHECKS, file.path(FIT_DIR, "module_checks_stageA.csv"))
MGMT <- merge(species_dt, data.table(species = seq_len(n_species), growth_beta = as_array(chk$nn_growth$beta), mortality_beta = as_array(chk$nn_mortality$beta)), by = "species")
fwrite(MGMT, file.path(FIT_DIR, "management_weights_stageA.csv"))
rm(chk)
## the fed recruitment must equal the observed trees over fit()'s window (output
## year k is driven by env row k; the observation at year p sums rows 1..p)
model$eval()
FINN.seed(PARAMS$seed)
S <- as.data.table(simulateForest(model, env = env_obsreg, init_cohort = init_cohort, patches = PARAMS$n_patches,
                                  patch_size = PS, device = "cpu")$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
W <- unique(obs[, .(siteID, period_length)])[, .(year = seq_len(period_length)), by = siteID]
fed <- merge(W, S[, .(siteID, year, species, r_mean_ha)], by = c("siteID", "year"))[, .(fed = sum(r_mean_ha) * PS), by = .(siteID, species)]
fed <- merge(fed, sched[year == 1, .(siteID, species, n)], by = c("siteID", "species"))
fed_ratio <- sum(fed$fed) / sum(fed$n)
stopifnot(abs(fed_ratio - 1) < 1e-3, max(abs(fed$fed - fed$n)) < 1e-2)
## the fitted species parameters and the U-shape per species
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(model[[block]]$detach()$cpu())), nm)
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")))
fwrite(SPP, file.path(FIT_DIR, "species_parameters_stageA.csv"))
sen_raw <- as.numeric(model$mort_sen_raw$detach()$cpu())
SEN <- data.table(species = seq_len(n_species), p_dbh = SPP$mort_size, sen = PARAMS$sen_max * plogis(sen_raw))
SEN[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & p_dbh < 0, -p_dbh / (2 * sen) * 100, Inf)]
SEN <- merge(species_dt, SEN, by = "species")
fwrite(SEN, file.path(FIT_DIR, "mortality_size_curve_stageA.csv"))
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence_stageA.csv"))

writeLines(c(
  sprintf("# FINN-model-NBI-Netherlands, zoo v0.2 recipe, stage A (%s)", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v02/recipe/fit_stageA.R` in %.1f min (%.1f s per epoch) on %s with %d thread(s): FINN %s from `%s`, torch %s, optimizer optim_adam, lr %g with %s decay.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$device, PARAMS$threads,
          as.character(packageVersion("FINN")), basename(dirname(find.package("FINN"))), as.character(packageVersion("torch")), PARAMS$lr, PARAMS$lr_scheduler), "",
  sprintf("Data: `prepared/%s/%s` (%d sites, %d species, intervals %d..%d yr, %d tracked trees of which %d lost (dead or felled); %.1f %% plots with recorded harvest). Data convention kept: total stem loss incl. felling is the loss target; not found / now < 5 cm / outside forest censored.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(trees), sum(trees$died),
          100 * mean(env[year == 1][[PARAMS$linear_var]])), "",
  sprintf("Model: competition (light_mask %s, %d quantiles); growth FINN::growth and mortality `mortality_ushape` (sen in [0, %g], start %g at 1 m; growth / %g), each through the module env_dnn_linmgmt = SELU network hidden = c(%s), dropout %g on every design-matrix column of `%s` except `%s` (column %d) + a linear per-species weight on `%s` (0/1 scale), started at 0, in growth AND mortality (v0.2.0 had harvest in the mortality network only); regeneration = the observed recruits fed as whole trees per patch and year (`regeneration_observed`, %d trees on %d plots; fed / observed over the fit windows %.4f).",
          PARAMS$light_mask, PARAMS$n_quantiles, PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, paste(PARAMS$hidden, collapse = ", "), PARAMS$dropout,
          paste(deparse(env_formula), collapse = ""), PARAMS$linear_var, mgmt_col, PARAMS$linear_var, n_fed, sched[year == 1 & n > 0, uniqueN(siteID)], fed_ratio), "",
  sprintf("Module checks (`module_checks_stageA.csv`, %d rows): forward pass in both call forms ok; gradients reach every MLP tensor and beta at the start; forcing %s to each of %s moves the output by beta_s x delta with max relative deviation %.2g (start) and %.2g (trained, reloaded) = a pure shift.",
          nrow(X_chk), PARAMS$linear_var, paste(PARAMS$mgmt_levels, collapse = "/"), max(CHECKS[stage == "init", shift_max_rel_dev]), max(CHECKS[stage != "init", shift_max_rel_dev])), "",
  sprintf("Management weights (harvest 0 -> 1 shift; `management_weights_stageA.csv`): log growth %s; mortality logit %s.",
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$growth_beta), collapse = ", "),
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$mortality_beta), collapse = ", ")), "",
  sprintf("Fit: %d epochs, batch %d, weights %s (weight_scale %s), losses %s, growth and mortality per tracked tree (period_mode per_site), env_autoscale climate auto / %s and reg_obs identity, %d patch of %.4f ha, recruits_dbh %g, seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s=%s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = " "), PARAMS$linear_var,
          PARAMS$n_patches, PS, PARAMS$recruits_dbh, PARAMS$seed), "",
  sprintf("Baseline-scaled losses (1 = no better than the intercept), first and last epoch: %s.",
          paste(sprintf("%s %.3f -> %.3f", names(loss_hist)[1:5], as.numeric(loss_hist[1, 1:5]), as.numeric(loss_hist[.N, 1:5])), collapse = ", ")), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm): %s (`mortality_size_curve_stageA.csv`).",
          paste(sprintf("%.3f", range(SEN$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SEN$species_name, SEN$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Next: `Rscript calibration/v02/recipe/fit_stageB.R` (regeneration), then `validation.R`."
), file.path(FIT_DIR, "fit_summary_stageA.md"))
message(sprintf("stage A done in %.1f min (%.1f s/epoch); final losses: %s", fit_minutes, 60 * fit_minutes / PARAMS$epochs,
                paste(sprintf("%s=%.3f", names(loss_hist)[1:5], as.numeric(loss_hist[.N, 1:5])), collapse = " ")))
