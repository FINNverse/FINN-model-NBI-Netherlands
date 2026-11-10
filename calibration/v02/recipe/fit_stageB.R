#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v02/recipe/fit_stageB.R
##
## WHAT   Stage B of the zoo v0.2 recipe for the Dutch NFI (pairs-v2): the
##        regeneration process (FINN::regeneration_saturation: light response x
##        env, Beverton-Holt cap K per species) is fitted with competition, growth
##        and mortality FROZEN at stage A (fit_stageA.R, demography through the
##        module env_dnn_linmgmt = climate network + linear harvest weight) and
##        the recruitment dispersion theta frozen at its maximum-likelihood
##        value. The result is the simulation model of the recipe fit (model.pt +
##        parameters.pt in fits/pairs-v2-recipe; the released v0.1.0 and the
##        v0.2.0 files are untouched). Identical to the Sweden template
##        (FINN-model-NFI-Sweden/calibration/v02/recipe/fit_stageB.R) except for
##        the per-model block of PARAMS and the count-space deviation below; the
##        regeneration env response is linear in the z-scaled climate and the
##        harvest flag (0/1 on its own scale, as in stage A).
## WHY    See the Spain reference and FINN-zoo-data/RECIPE_v02.md: recruitment is
##        identified by the recruitment observations alone; theta per species =
##        intercept-only ML on the observed recruit counts, clamped to
##        [theta_min, theta_max], not trained (the NBI pilot P3 froze theta at 1
##        and hurt: recruitment ratio 2.5; the ML value per species is a
##        different thing and is reported in fit_summary_stageB.md); ALL loss
##        terms, weights "auto" (the stage-B loss rule settled by the Spain
##        divergence test); the frozen networks run without dropout; cosine lr,
##        plain Adam, one thread, light_mask TRUE.
## DEVIATION (documented): the recipe scores the negative binomial on whole
##        recruit counts, recruit_obs_weight = the area one observed recruit
##        stands for, and asserts that reg x w is whole. The NBI plots are
##        variable-radius circles (5-20 m) and pairs-v2 counts the recruits
##        inside the smaller circle of the pair, so one recruit stands for
##        0.005-0.126 ha depending on the plot, while FINN's recruit_obs_weight
##        is one number. w is set to the MEDIAN smaller-circle area of the
##        training pairs (the 10 m circle, 0.0314 ha, computed from
##        site_table.csv and asserted): reg x w is then a whole count on the
##        median plot, 1/4 of a count per recruit on a 5 m plot, 4 per recruit
##        on a 20 m one. The whole-count assertion is replaced by a report of
##        the share of recruiting rows whose reg x (own smaller circle) is whole
##        and of the quantiles of reg x w; theta.ml takes the pseudo-counts as
##        they are.
## DATA   ../../FINN-zoo-data/nbi/prepared/pairs-v2/train  (+ prepared/pairs-v2/site_table.csv) and
##        <FIT_DIR>/model_stageA.pt (written by fit_stageA.R)
## NEEDS  FINN 0.2.0 (commit 0848b60) in ~/R/lib-finn-0848b60, torch, MASS (theta.ml)
## RUN    Rscript calibration/v02/recipe/fit_stageB.R      (from the repo root, after fit_stageA.R)
## OUTPUT <FIT_DIR>/model.pt, parameters.pt, model_init_stageB.pt, loss_history_stageB.csv,
##        regeneration_parameters.csv, convergence_stageB.csv, fit_summary_stageB.md
##        FIT_DIR = ../../FINN-zoo-data/nbi/fits/pairs-v2-recipe (outside git)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v2-recipe")

PARAMS <- list(
  ## ---- per model (must equal fit_stageA.R where both have the entry) ----
  prep_version  = "pairs-v2",
  train_set     = "train",
  holdout_set   = "holdout_blocked",       # not read here; recorded for validation.R
  env_vars      = c("tavg_c", "prec_mm", "harvest"),   # as stage A; harvest enters the regeneration env response linearly too
  linear_var    = "harvest",               # kept on its 0/1 scale (env_autoscale identity), as in stage A
  patch_size_ha = pi * 20^2 / 1e4,         # as stage A
  n_patches     = 1L,
  recruits_dbh  = 5.9,
  recruit_obs_weight_ha = pi * 10^2 / 1e4, # the area one observed recruit stands for: the MEDIAN smaller circle of the training pairs (10 m, 0.0314 ha; the circles span 5-20 m, FINN takes one number; see the header)
  w_tol         = 1e-4,                    # |median smaller-circle area of the training pairs - recruit_obs_weight_ha| allowed (asserted from site_table.csv)
  reg_K_init    = 50,                      # Beverton-Holt cap start, recruits/ha/yr (the 95th percentile of the annual per-ha recruitment of a pair is ~38; per species 5-127 on pairs with recruits)
  ## ---- the recipe ----
  reg_floor     = 1e-3,
  theta_min     = 0.01,
  theta_max     = 100,
  sen_max       = 5,
  growth_scale  = 0.01,
  frozen_dropout = 0,
  light_mask    = TRUE,
  n_quantiles   = 10L,
  epochs        = 150L,
  batchsize     = 100L,
  lr            = 0.005,
  lr_scheduler  = "cosine",
  weights       = "auto",
  weight_scale  = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # ALL terms (the stage-B loss rule)
  loss_aggregation = list(growth = "tree", mortality = "tree"),
  checkpoints   = 10L,
  frozen_tol    = 1e-5,
  threads       = 1L,
  device        = "cpu",
  seed          = 1L
)
torch_set_num_threads(PARAMS$threads)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))
PS <- PARAMS$patch_size_ha
w_obs <- PARAMS$recruit_obs_weight_ha

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv"))[, .(siteID, patchID, species, dbh, trees, treeName)]
trees <- fread(file.path(train_dir, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
autoscale_spec <- c(setNames(rep("auto", length(setdiff(PARAMS$env_vars, PARAMS$linear_var))), setdiff(PARAMS$env_vars, PARAMS$linear_var)),
                    setNames("identity", PARAMS$linear_var))
## the recruit observation space: the plot's own smaller circle (site_table.csv: r0, r1 in m,
## area0_m2, area1_m2 = the forest part of the circle) against the one weight FINN takes
site_tab <- merge(fread(file.path(train_dir, "siteID_dt.csv"))[, .(siteID, siteName)],
                  fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "site_table.csv"))[, .(siteName, r0, r1, area0_m2, area1_m2)], by = "siteName")
site_tab[, a_min_ha := pmin(area0_m2, area1_m2) / 1e4]
stopifnot(nrow(site_tab) == uniqueN(obs$siteID), abs(median(site_tab$a_min_ha) - w_obs) < PARAMS$w_tol)
oo <- merge(obs[reg > 0, .(siteID, species, reg)], site_tab[, .(siteID, a_min_ha)], by = "siteID")
whole_share <- oo[, mean(abs(reg * a_min_ha - round(reg * a_min_ha)) < 1e-3)]
y_q <- quantile(obs[reg > 0, reg * w_obs], c(0, 0.25, 0.5, 0.75, 1))
message(sprintf("recruit observation space: reg x %.4f ha (median smaller circle; circles %d-%d m, area %.4f-%.4f ha); %.1f %% of the recruiting rows are whole counts on their own circle; reg x w on recruiting rows: min %.2f, quartiles %.2f / %.2f / %.2f, max %.1f",
                w_obs, min(site_tab[, pmin(r0, r1)]), max(site_tab[, pmin(r0, r1)]), min(site_tab$a_min_ha), max(site_tab$a_min_ha), 100 * whole_share, y_q[1], y_q[2], y_q[3], y_q[4], y_q[5]))

## ---- 2. stage A, frozen ----------------------------------------------------------
stageA <- torch_load(file.path(FIT_DIR, "model_stageA.pt"))
stopifnot(inherits(stageA$nn_growth, "env_dnn_linmgmt"), inherits(stageA$nn_mortality, "env_dnn_linmgmt"))
nn_pars <- function(nn) lapply(nn$parameters, function(p) as_array(p$detach()$cpu()))
nn_maxdiff <- function(a, b) { stopifnot(length(a) == length(b)); max(mapply(function(x, y) max(abs(x - y)), a, b)) }
for (nn in list(stageA$nn_growth, stageA$nn_mortality)) for (m in nn$modules) if (inherits(m, "nn_dropout")) m$p <- PARAMS$frozen_dropout
drop_p <- unlist(lapply(list(stageA$nn_growth, stageA$nn_mortality), function(nn)
  vapply(Filter(function(m) inherits(m, "nn_dropout"), nn$modules), function(m) m$p, numeric(1))))
stopifnot(length(drop_p) > 0, all(drop_p == PARAMS$frozen_dropout))
frozen <- list(competition = as.matrix(stageA$par_competition$detach()$cpu()), growth = as.matrix(stageA$par_growth$detach()$cpu()),
               mortality = as.matrix(stageA$par_mortality$detach()$cpu()), sen_raw = as.numeric(stageA$mort_sen_raw$detach()$cpu()),
               growth_env = nn_pars(stageA$nn_growth), mortality_env = nn_pars(stageA$nn_mortality))
scaling_A <- as.data.table(stageA$env_scaling)[variable %in% PARAMS$env_vars]
## the stage-A mortality kernel (identical to fit_stageA.R)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}

## ---- 3. the recruitment dispersion: intercept-only ML per species, frozen ----------
theta_ml <- vapply(seq_len(n_species), function(s) {
  y <- obs[species == s, reg * w_obs]
  if (sum(y) == 0) return(PARAMS$theta_max)
  th <- suppressWarnings(tryCatch(as.numeric(MASS::theta.ml(y, mu = mean(y), limit = 100)), error = function(e) NA_real_))
  if (!is.finite(th)) PARAMS$theta_max else min(max(th, PARAMS$theta_min), PARAMS$theta_max)
}, numeric(1))
## warm start of the regeneration env layer: intercept = log of the species' mean
## annual recruitment per ha (the climate columns are z-scaled and harvest is 0/1, so 0
## slopes = the mean of an unharvested plot at the mean climate)
mean_annual_reg <- obs[, .(r = sum(reg) / sum(period_length)), by = species][order(species), r]
W0 <- cbind(log(mean_annual_reg + PARAMS$reg_floor), matrix(0, n_species, length(PARAMS$env_vars)))

## ---- 4. the stage-B model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  recruit_obs_weight = w_obs,
  competition_process  = createProcess(~0, FINN::competition, initSpecies = frozen$competition,
                                       optimizeSpecies = FALSE, n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = frozen$growth,
                                       NN = stageA$nn_growth, optimizeSpecies = FALSE, optimizeEnv = FALSE),
  mortality_process    = createProcess(env_formula, mortality_ushape, initSpecies = frozen$mortality,
                                       NN = stageA$nn_mortality, optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       custom_parameters = list(mort_sen_raw = list(init = frozen$sen_raw, trainable = FALSE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_saturation,
                                       initEnv = list(W0), dispersion_parameter = theta_ml,
                                       custom_parameters = list(reg_logK = list(init = rep(log(PARAMS$reg_K_init), n_species), trainable = TRUE)),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
model$period_mode <- "per_site"
model$env_scaling <- as.data.frame(scaling_A)   # stage A's climate scaling (harvest identity); fit() recomputes it from the same env, checked below
e_build <- env[siteID %in% sort(unique(siteID))[1:2] & year <= 2L]
e_build[, siteID := match(siteID, sort(unique(siteID)))]
invisible(simulateForest(model, env = e_build, init_cohort = NULL, patches = 1L, patch_size = PS, device = "cpu"))
for (nn in list(model$nn_growth, model$nn_mortality)) for (p in nn$parameters) invisible(p$requires_grad_(FALSE))
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(nn_maxdiff(nn_pars(model$nn_growth), frozen$growth_env) < 1e-7, nn_maxdiff(nn_pars(model$nn_mortality), frozen$mortality_env) < 1e-7,
          !any(sapply(c(model$nn_growth$parameters, model$nn_mortality$parameters), function(p) p$requires_grad)),
          max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_ml)) < 1e-4 * max(theta_ml),
          max(abs(as.matrix(model$par_mortality$detach()$cpu()) - frozen$mortality)) < 1e-4,
          max(abs(as_array(model$nn_regeneration$parameters[[1]]$detach()) - W0)) < 1e-6)
torch_save(model, file.path(FIT_DIR, "model_init_stageB.pt"))

## ---- 5. fit -------------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env           = env,
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
model$eval()
torch_save(model, file.path(FIT_DIR, "model.pt"))
parameters <- list(
  state_dict  = model$state_dict(),
  env_scaling = model$env_scaling,
  species     = species_dt,
  settings    = list(prep_version = PARAMS$prep_version, patch_size_ha = PS, n_patches = PARAMS$n_patches,
                     recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor, recruit_obs_weight = w_obs,
                     env_formula = deparse(env_formula), linear_var = PARAMS$linear_var, sen_max = PARAMS$sen_max, growth_scale = PARAMS$growth_scale,
                     theta = theta_ml, epochs_stageB = PARAMS$epochs, recipe = "zoo v0.2 (stage A + B, linear management)",
                     finn_version = as.character(packageVersion("FINN")), torch_version = as.character(packageVersion("torch")))
)
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
loss_hist <- as.data.table(do.call(rbind, model$history))
setnames(loss_hist, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(loss_hist))])
loss_hist[, epoch := .I]
fwrite(loss_hist, file.path(FIT_DIR, "loss_history_stageB.csv"))

## ---- 6. checks and what was fitted ------------------------------------------------------
dev <- c(growth_env    = nn_maxdiff(nn_pars(model$nn_growth), frozen$growth_env),
         mortality_env = nn_maxdiff(nn_pars(model$nn_mortality), frozen$mortality_env),
         growth_par    = max(abs(as.matrix(model$par_growth$detach()$cpu()) - frozen$growth)),
         mortality_par = max(abs(as.matrix(model$par_mortality$detach()$cpu()) - frozen$mortality)),
         competition   = max(abs(as.matrix(model$par_competition$detach()$cpu()) - frozen$competition)),
         sen           = max(abs(as.numeric(model$mort_sen_raw$detach()$cpu()) - frozen$sen_raw)),
         theta_rel     = max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_ml) / theta_ml))
scaling_B <- as.data.table(model$env_scaling)[variable %in% PARAMS$env_vars]
dev["env_scaling"] <- max(abs(merge(scaling_A, scaling_B, by = "variable")[, c(center.x - center.y, scale.x - scale.y)]))
message(sprintf("max |change| of the frozen blocks: %s", paste(sprintf("%s %.2g", names(dev), dev), collapse = ", ")))
stopifnot(all(dev < PARAMS$frozen_tol))
chk <- torch_load(file.path(FIT_DIR, "model.pt"))
stopifnot(nn_maxdiff(nn_pars(chk$nn_growth), frozen$growth_env) < 1e-7, nn_maxdiff(nn_pars(chk$nn_mortality), frozen$mortality_env) < 1e-7,
          inherits(chk$nn_growth, "env_dnn_linmgmt"), setequal(names(parameters$state_dict), names(chk$state_dict())))
rm(chk)
W3 <- as_array(model$nn_regeneration$parameters[[1]]$detach()$cpu())
REG <- data.table(species = seq_len(n_species),
                  mean_annual_reg_obs_ha = mean_annual_reg,
                  theta = theta_ml,
                  K_init = PARAMS$reg_K_init, K_fit = exp(as.numeric(model$reg_logK$detach()$cpu())),
                  light = as.numeric(as.matrix(model$par_regeneration$detach()$cpu())[, 1]),
                  `colnames<-`(as.data.table(W3), paste0("env_", c("intercept", PARAMS$env_vars))))
REG <- merge(species_dt, REG, by = "species")
fwrite(REG, file.path(FIT_DIR, "regeneration_parameters.csv"))
rec_ratio <- function(m) {
  m$eval(); FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = init_cohort, patches = PARAMS$n_patches, patch_size = PS, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  W <- unique(obs[, .(siteID, period_length)])[, .(year = seq_len(period_length)), by = siteID]
  sum(merge(W, S[, .(siteID, year, r_mean_ha)], by = c("siteID", "year"))$r_mean_ha) / sum(pmax(obs$reg, 0))
}
ratio_start <- rec_ratio(torch_load(file.path(FIT_DIR, "model_init_stageB.pt")))
ratio_end   <- rec_ratio(model)
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence_stageB.csv"))

writeLines(c(
  sprintf("# FINN-model-NBI-Netherlands, zoo v0.2 recipe, stage B (%s)", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v02/recipe/fit_stageB.R` in %.1f min (%.1f s per epoch) on %s with %d thread(s): FINN %s from `%s`, torch %s, optimizer optim_adam, lr %g with %s decay.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$device, PARAMS$threads,
          as.character(packageVersion("FINN")), basename(dirname(find.package("FINN"))), as.character(packageVersion("torch")), PARAMS$lr, PARAMS$lr_scheduler), "",
  sprintf("Data: `prepared/%s/%s` (%d sites, %d species). Holdout for validation: `%s`. Frozen demography: `model_stageA.pt` (competition, growth, mortality incl. the U-shape and both env_dnn_linmgmt modules, dropout %g while frozen); max |change| of the frozen blocks %s.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, PARAMS$holdout_set, PARAMS$frozen_dropout,
          paste(sprintf("%s %.1g", names(dev), dev), collapse = ", ")), "",
  sprintf("Regeneration: FINN::regeneration_saturation (light x exp(linear env on `%s`, %s on its 0/1 scale), Beverton-Holt cap K per species, start %g/ha/yr, reg_floor %g), intercept warm-started at the observed mean annual recruitment; theta per species = intercept-only ML on the observed recruits in the observation space reg x %.4f ha (the median smaller circle of the variable-radius pairs, circles %d-%d m: %.1f %% of the recruiting rows are whole counts on their own circle, reg x w quartiles %.2f / %.2f / %.2f), clamped to [%g, %g], frozen: %s.",
          paste(deparse(env_formula), collapse = ""), PARAMS$linear_var, PARAMS$reg_K_init, PARAMS$reg_floor, w_obs,
          min(site_tab[, pmin(r0, r1)]), max(site_tab[, pmin(r0, r1)]), 100 * whole_share, y_q[2], y_q[3], y_q[4], PARAMS$theta_min, PARAMS$theta_max,
          paste(sprintf("%s %.3g", REG$species_name, REG$theta), collapse = ", ")), "",
  sprintf("Fit: %d epochs, batch %d, weights %s (weight_scale %s), losses %s, growth and mortality per tracked tree (period_mode per_site), seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s=%s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = " "), PARAMS$seed), "",
  sprintf("Baseline-scaled losses (1 = no better than the intercept), first and last epoch: %s.",
          paste(sprintf("%s %.3f -> %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[1, 1:6]), as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Predicted / observed recruitment on the training plots (fit windows): %.2f at the start, %.2f after stage B.", ratio_start, ratio_end), "",
  sprintf("Fitted cap K (recruits/ha/yr): %s; harvest slope of log recruitment (harvest 0 -> 1): %s (`regeneration_parameters.csv`).",
          paste(sprintf("%s %.0f", REG$species_name, REG$K_fit), collapse = ", "),
          paste(sprintf("%s %.2f", REG$species_name, REG[[paste0("env_", PARAMS$linear_var)]]), collapse = ", ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Validation report: `Rscript calibration/v02/recipe/validation.R`."
), file.path(FIT_DIR, "fit_summary_stageB.md"))
message(sprintf("stage B done in %.1f min (%.1f s/epoch); recruitment ratio %.2f -> %.2f; final losses: %s", fit_minutes, 60 * fit_minutes / PARAMS$epochs,
                ratio_start, ratio_end, paste(sprintf("%s=%.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = " ")))
