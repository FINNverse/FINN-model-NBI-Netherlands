#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/fit_pilot.R
##
## WHAT   Smoke test of calibration/fit.R: the SAME model on the 300-site pilot
##        set for 3 epochs on one thread. It checks that the prepared tables
##        (variable interval lengths 4-12 yr) go through fit(), that
##        parameters.pt rebuilds the model bit for bit through the package's
##        build_model(), and it times an epoch so the production run can be
##        planned (s/epoch scales with the number of sites). Nothing from this
##        folder is released.
## DIFFERS FROM calibration/fit.R: set "pilot" (300 sites), 3 epochs,
##        batchsize 100, one torch thread, output folder fits/pilot-smoke.
## RUN    Rscript calibration/fit_pilot.R      (from the repository root, ~1 min)
## OUTPUT <DATA_DIR>/fits/pilot-smoke/{model.pt, parameters.pt, loss_history.csv, timing.txt}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)                                          # the laptop is shared with other fits

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pilot-smoke")
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "pairs-v1", train_set = "pilot",
  patch_size_ha = pi * 20^2 / 1e4, n_patches = 1L,       # 0.1257 ha: the 20 m NBI circle as one patch, as prepared
  env_vars      = c("tavg_c", "prec_mm"),
  recruits_dbh  = 6.4, hidden = c(20L, 20L), dropout = 0.1, reg_K_init = 50, reg_floor = 1e-3,
  epochs        = 3L,                # enough to time an epoch and see the loss move
  batchsize     = 100L,              # 3 batches per epoch on 300 sites (every interval length 4-12 yr should be in each batch)
  lr            = 0.01, weights = "auto", weight_scale = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints   = 1L, device = "cpu", seed = 1L
)

## ---- 1. data ---------------------------------------------------------------------
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
train <- list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")),
              init = fread(file.path(d, "init_cohort_dt.csv")), species = fread(file.path(d, "species_dt.csv")))
n_species <- nrow(train$species)
obs <- train$obs[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- train$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(train$init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
message(sprintf("pilot: %d sites, %d species, intervals %d-%d yr, %d obs rows, %d cohorts",
                uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(obs), nrow(train$init)))

## ---- 2. the model, written out (identical to fit.R and R/definition.R) ---------
env_formula <- ~ tavg_c + prec_mm
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(env_formula, FINN::growth,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, FINN::mortality,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(env_formula, FINN::regeneration_saturation,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species)))
)

## ---- 3. fit and time ---------------------------------------------------------------
FINN.seed(PARAMS$seed)
t_start <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr,
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha,
    optimizer = torch::optim_ignite_adam, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
parameters <- lapply(model$parameters, function(p) p$detach()$cpu()$clone())
for (i in seq_len(nrow(model$env_scaling))) {
  v <- model$env_scaling$variable[i]
  parameters[[paste0("env_scaling.", v, ".center")]] <- torch_tensor(model$env_scaling$center[i], dtype = torch_float64())   # double: the model holds them as R doubles
  parameters[[paste0("env_scaling.", v, ".scale")]]  <- torch_tensor(model$env_scaling$scale[i],  dtype = torch_float64())
}
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
hist <- as.data.table(do.call(rbind, model$history))
setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE)
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
print(hist)
cat("parameter tensors saved:\n"); print(vapply(parameters, function(p) paste(dim(p), collapse = "x"), character(1)))

## ---- 4. round trip through the package's build_model() ------------------------------
species_table <- function() utils::read.csv(file.path("inst", "model", "species.csv"))
source(file.path("R", "definition.R"))
rebuilt <- build_model("process", torch_load(file.path(FIT_DIR, "parameters.pt")))
check_ids <- sort(unique(obs$siteID))[seq_len(50L)]
roundtrip <- function(m) {
  m$eval(); FINN.seed(PARAMS$seed)
  as.data.table(simulateForest(m, env = env[siteID %in% check_ids],
                               init_cohort = makeInitCohorts(train$init[siteID %in% check_ids, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                               patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = PARAMS$device, batchsize = 50L)$wide$site)
}
stopifnot(identical(roundtrip(model), roundtrip(rebuilt)))
message("round trip: build_model(parameters.pt) reproduces the fitted model bit for bit on ", length(check_ids), " sites")
timing <- sprintf("pilot: %d sites, %d epochs, 1 thread: %.1f s total, %.1f s/epoch (incl. setup); FINN %s torch %s",
                  uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs, packageVersion("FINN"), packageVersion("torch"))
writeLines(timing, file.path(FIT_DIR, "timing.txt")); message(timing)
