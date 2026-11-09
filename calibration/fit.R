#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/fit.R
##
## WHAT   The production fit of the FINNmodelNBINetherlands "process" variant: a
##        FINN process model (growth, mortality, saturating regeneration,
##        competition) for 19 Dutch tree species / species groups + "other",
##        calibrated to the tree-by-tree remeasurement of the permanent plots of
##        the Dutch national forest inventory (MFV -> NBI-6 and NBI-6 -> NBI-7
##        pairs, 4-12 years apart), prepared data pairs-v1
##        (calibration/03_prepare_pairs_v1.R). Annual timestep, one observation
##        per site at year = period_length (FINN's per-site aggregation windows).
## WHY THIS SHAPE
##        - Mechanistic FINN processes with a small SELU network per process on
##          two climate normals (tavg_c, prec_mm): the network learns the
##          species x climate response and its interactions; no coordinates, so
##          the model transfers to any site with WorldClim normals and can be
##          merged with the other zoo regions (same covariate names as FINN-fia).
##        - FINN::regeneration_saturation (Beverton-Holt cap per species) because
##          a bare-ground run with unbounded recruitment explodes (finn-evaluation).
##        - One patch of 0.1257 ha per site = the 20 m NBI circle; every tallied
##          tree expanded to its per-hectare weight (PREPARATION.md).
##        - The mortality target is total stem LOSS incl. felled trees (the BWI
##          and IFN-Spain decision; natural death and felling cannot be separated
##          tree by tree, censoring would bias the natural rate): the fitted
##          "mortality" is the stem turnover of the managed Dutch forest.
##        - Losses: the FINN 0.2.0 defaults (dbh/ba mse, trees poisson, growth
##          mse, mortality binomial on tallied counts, recruitment nbinom), each
##          scaled by its intercept-only baseline (weights = "auto"). No numeric
##          weight vector is passed, so weight_scale is moot (stated as "auto";
##          FINN#37).
## DATA   <DATA_DIR>/prepared/pairs-v1/train            (fitting; 1,531 sites)
##        <DATA_DIR>/prepared/pairs-v1/holdout_blocked  (validation.R only; 740 sites)
## DIFFERS FROM calibration/fit_pilot.R: the full training set, 500 epochs,
##        batchsize 150 (the pilot fits 300 sites for 3 epochs to time an epoch).
## RUN    Rscript calibration/fit.R      (from the repository root; one thread,
##        ~7 s/epoch for 1,531 sites -> ~1 h; the coordinator's laptop queue or
##        a compute node starts it, never this script's author)
## OUTPUT <DATA_DIR>/fits/pairs-v1-process/model.pt       the whole fitted model (release asset + local use)
##        <DATA_DIR>/fits/pairs-v1-process/parameters.pt  NAMED LIST of parameter tensors
##                                                        (+ env scaling) = the release asset load_model() reads
##        calibration/loss_history.csv, calibration/fit_summary.md, inst/model/weights.json
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")   # zoo layout: data and weights outside git
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v1-process")
CAL_DIR  <- "calibration"
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "pairs-v1", train_set = "train", holdout_set = "holdout_blocked",
  patch_size_ha = pi * 20^2 / 1e4,   # 0.1257 ha: the 20 m NBI circle as one patch, as prepared
  n_patches     = 1L,                # one plot = one patch (taken from the init cohort)
  env_vars      = c("tavg_c", "prec_mm"),   # WorldClim 2.1 bio1 (degC), bio12 (mm/yr); natural units, FINN z-scales
  recruits_dbh  = 6.4,               # median dbh of the observed recruits; tally threshold 5 cm
  hidden        = c(20L, 20L),       # env network per process (SELU), as the other zoo fits
  dropout       = 0.1,
  reg_K_init    = 50,                # Beverton-Holt cap init, stems/ha/yr (observed recruitment ~20/ha/yr mean, ~9 median at the 5 cm threshold; the other zoo fits' init)
  reg_floor     = 1e-3,              # minimum recruitment, stems/ha/yr: lets a species arrive where it is absent
  epochs        = 500L,              # pilot: 1.4 s/epoch per 300 sites on one core -> ~7 s/epoch here, ~1 h; 500 x 11 batches = 5,500 Adam steps; settling is checked in validation (parameter change)
  batchsize     = 150L,              # ~11 batches per epoch over 1,531 sites; every interval length (4-12 yr) is in nearly every batch
  lr            = 0.01,
  threads       = 1L,                # the laptop runs several zoo fits side by side: one torch thread each
  weights       = "auto",            # each loss scaled by its intercept-only baseline
  weight_scale  = "auto",            # no numeric weight vector, so the scale setting is moot; stated for the card
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints   = 10L,               # parameter snapshots every 10 epochs into model$param_history (folder = NULL: no disk cost)
  device        = "cpu",
  seed          = 1L
)

torch_set_num_threads(PARAMS$threads)

## ---- 1. data ---------------------------------------------------------------------
read_set <- function(set) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")),
       init = fread(file.path(d, "init_cohort_dt.csv")), species = fread(file.path(d, "species_dt.csv")))
}
train <- read_set(PARAMS$train_set)
n_species <- nrow(train$species)
obs <- train$obs[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- train$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(train$init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
message(sprintf("train: %d sites, %d species, intervals %d-%d yr, %d obs rows, %d cohorts",
                uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(obs), nrow(train$init)))

## ---- 2. the model, written out (identical to R/definition.R) ---------------------
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

## ---- 3. fit ----------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t_start <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr,
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha,
    optimizer = torch::optim_ignite_adam, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
runtime_h <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))

## save first: the whole model and the release object, a NAMED LIST of parameter
## tensors plus the env standardisation constants the model learned (FINN
## z-scales env internally and must reuse the SAME constants at prediction)
torch_save(model, file.path(FIT_DIR, "model.pt"))
parameters <- lapply(model$parameters, function(p) p$detach()$cpu()$clone())
for (i in seq_len(nrow(model$env_scaling))) {
  v <- model$env_scaling$variable[i]
  parameters[[paste0("env_scaling.", v, ".center")]] <- torch_tensor(model$env_scaling$center[i], dtype = torch_float64())   # double: the model holds them as R doubles
  parameters[[paste0("env_scaling.", v, ".scale")]]  <- torch_tensor(model$env_scaling$scale[i],  dtype = torch_float64())
}
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
## model$history: one unnamed row per epoch in FINN's fixed column order (dbh, ba,
## trees, growth, mortality, regeneration, a seventh unused slot), baseline-scaled
hist <- as.data.table(do.call(rbind, model$history))
setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE)
hist[, epoch := .I]
fwrite(hist, file.path(CAL_DIR, "loss_history.csv"))
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))

## ---- 3b. round trip: the package's build_model() must reproduce this model bit for bit ----
## (name-matched copy of every tensor, env networks included; a mismatch stops build_model)
species_table <- function() utils::read.csv(file.path("inst", "model", "species.csv"))
source(file.path("R", "definition.R"))
rebuilt <- build_model("process", torch_load(file.path(FIT_DIR, "parameters.pt")))
check_ids <- sort(unique(obs$siteID))[seq_len(200L)]
roundtrip <- function(m) {
  m$eval(); FINN.seed(PARAMS$seed)
  as.data.table(simulateForest(m, env = env[siteID %in% check_ids],
                               init_cohort = makeInitCohorts(train$init[siteID %in% check_ids, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                               patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = PARAMS$device, batchsize = 200L)$wide$site)
}
sim_saved <- roundtrip(model); sim_rebuilt <- roundtrip(rebuilt)
stopifnot(identical(sim_saved, sim_rebuilt))
message("round trip: build_model(parameters.pt) reproduces the fitted model bit for bit on ", length(check_ids), " sites")

## ---- 4. release metadata and fit_summary.md -----------------------------------------
## inst/model/weights.json: what load_model() checks. The asset URL follows the
## repository's release layout; the coordinator publishes the tag and the assets.
sha <- unname(tools::sha256sum(file.path(FIT_DIR, "parameters.pt")))
sha_model <- unname(tools::sha256sum(file.path(FIT_DIR, "model.pt")))
release_tag <- "v0.1.0"
## both files ship as release assets: the name-matched parameter dump that
## load_model() uses, and the whole model object (torch_save of the fitted module)
weights <- list(variants = list(process = list(
  tag = release_tag, file = "parameters.pt",
  url = sprintf("https://github.com/FINNverse/FINN-model-NBI-Netherlands/releases/download/%s/parameters.pt", release_tag),
  sha256 = sha,
  model_file = "model.pt",
  model_url = sprintf("https://github.com/FINNverse/FINN-model-NBI-Netherlands/releases/download/%s/model.pt", release_tag),
  model_sha256 = sha_model,
  finn_commit = "0848b60", finn_version = as.character(packageVersion("FINN")), torch_version = as.character(packageVersion("torch")),
  patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha,
  data_version = PARAMS$prep_version, train_set = PARAMS$train_set, epochs = PARAMS$epochs,
  weights = PARAMS$weights, weight_scale = PARAMS$weight_scale,
  fitted = format(Sys.time(), "%Y-%m-%d"))))
writeLines(jsonlite::toJSON(weights, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path("inst", "model", "weights.json"))
last <- hist[.N]
md <- c("# Fit summary: pairs-v1-process", "",
        sprintf("Fitted %s by calibration/fit.R; FINN %s (main 0848b60), torch %s, R %s; %d threads; runtime %.2f h for %d epochs (%.1f s/epoch).",
                format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN"), packageVersion("torch"), getRversion(),
                torch_get_num_threads(), runtime_h, PARAMS$epochs, runtime_h * 3600 / PARAMS$epochs), "",
        sprintf("Data: %s/%s, %d sites, %d species (incl. other), intervals %d-%d yr. Weights: `auto` (baseline-scaled; no numeric vector, weight_scale `auto` moot).",
                PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length)), "",
        sprintf("Release assets: parameters.pt (named list of %d tensors incl. env scaling), sha256 %s; model.pt (whole fitted module), sha256 %s. build_model(parameters.pt) reproduced the fitted model bit for bit on %d sites.",
                length(parameters), sha, sha_model, length(check_ids)), "",
        "Final losses (baseline-scaled, last epoch):", "",
        paste0("- ", names(last), ": ", vapply(last, function(x) format(signif(as.numeric(x), 4)), character(1))), "",
        "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(CAL_DIR, "fit_summary.md"))
message("done: ", FIT_DIR, "  parameters.pt sha256 ", sha)
