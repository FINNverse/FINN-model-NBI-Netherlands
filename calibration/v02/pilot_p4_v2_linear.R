#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/v02/pilot_p4_v2_linear.R
##
## WHAT   v0.2 pilot P4: P3 (pairs-v2, cosine, theta frozen, mortality bounds,
##        harvest in mortality) with LINEAR environmental responses
##        (hidden = NULL: one linear layer per process instead of the
##        SELU network 20-20).
## HYPOTHESIS the Dutch climate space is ~1.5 degC x 230 mm wide and the data
##        hold 1,524 training sites: a 20-20 network per process has far more
##        freedom than the gradient can identify and absorbs state-loss
##        pressure (BWI: loss competition flipped weakly identified slopes).
##        Linear responses should transfer at least as well to the blocked
##        holdout and settle better.
## RUN    Rscript calibration/v02/pilot_p4_v2_linear.R   (repo root, ~10 min, one thread)
## OUTPUT <DATA_DIR>/fits/v02-pilots/p4_v2_linear/{model.pt, loss_history.csv}; a row in calibration/v02/pilot_table.csv
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))             # generic metric helpers only
torch_set_num_threads(1L)

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
PARAMS <- list(
  label = "p4_v2_linear", prep_version = "pairs-v2", train_set = "pilot", holdout_set = "holdout_blocked",
  patch_size_ha = pi * 20^2 / 1e4, n_patches = 1L, env_vars = c("tavg_c", "prec_mm", "harvest"),
  recruits_dbh = 6.4, hidden = NULL, dropout = 0, reg_K_init = 50, reg_floor = 1e-3,   # hidden = NULL: linear env response
  mort_bounds = c(-3, 3), freeze_theta = TRUE,
  epochs = 400L, batchsize = 100L, lr = 0.01, lr_scheduler = "cosine",
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, seed = 1L)
FIT_DIR <- file.path(DATA_DIR, "fits", "v02-pilots", PARAMS$label); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs <- fread(file.path(d, "obs_dt.csv")); env <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv"))
n_species <- nrow(fread(file.path(d, "species_dt.csv")))
init_cohort <- makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)

env_formula <- ~ tavg_c + prec_mm
mort_formula <- ~ tavg_c + prec_mm + harvest
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(env_formula, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(mort_formula, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       lower = rep(PARAMS$mort_bounds[1], 3), upper = rep(PARAMS$mort_bounds[2], 3)),
  regeneration_process = createProcess(env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species))))
if (PARAMS$freeze_theta) model$par_theta_recruits_raw$requires_grad_(FALSE)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's older torch falls back to optim_adam
optimizer <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
FINN.seed(PARAMS$seed); t0 <- Sys.time()
fit(model, env = env[, c("siteID", "year", PARAMS$env_vars), with = FALSE],
    data = obs[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)], init_cohort = init_cohort,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, optimizer = optimizer, device = "cpu",
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history)); hist[, epoch := .I]; fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
row <- pilot_eval(model, file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set), n_species, PARAMS$patch_size_ha, PARAMS$label, PARAMS$env_vars)
row[, `:=`(data = PARAMS$prep_version, epochs = PARAMS$epochs, fit_min = round(mins, 1), theta_min = round(min(as.numeric(as_array(model$par_theta_recruits))), 4))]
append_pilot_row(row)
