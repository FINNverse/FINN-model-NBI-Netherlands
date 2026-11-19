#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v020/fit_continue.R
##
## WHAT   Continues calibration/v020/fit.R (the v020-defaults joint fit of the
##        Dutch NFI model, pairs-v2; see fit.R's header for the set-up, the
##        model structure and the inventory problems) from its newest checkpoint
##        up to the same 500 epochs, with the learning-rate schedule CONTINUED
##        (constant 0.01 for epochs 1-300, cosine 0.01 -> 1e-4 over 301-500),
##        then does what fit.R does after training. Nothing about the model, the
##        data or the objective differs from fit.R. The resume logic is the
##        one of every v020 continuation in the zoo.
## HOW    - Newest checkpoint: checkpoints/epoch_<n>model.pt (epoch n), or
##          checkpoints/resume_<E>/epoch_<n>model.pt (epoch E + n) from an earlier
##          continuation. Epochs after the newest checkpoint are lost (at most 24).
##        - A fresh finn() is built as fit.R builds it, fit(epochs = 0)
##          registers the modules and the data tensors, and load_state_dict()
##          copies every checkpoint tensor; the script stops unless every tensor
##          is identical.
##        - Loss weights: the 1 / baseline vector fit.R's loss_weights = "auto"
##          stored in the checkpoint (model$loss_weights), passed with
##          weight_scale = "raw"; the same objective as fit.R.
##        - THE CAP: FINN resolves prior = "max_obs_recruits" only once per
##          model object, so the set-up fit(epochs = 0) resolves it on the fresh
##          model, load_state_dict() then copies the trained reg_logK, and the
##          continuing fit() keeps it. fit() builds the lr scheduler after the
##          prior step and before the first step; the schedule function below
##          checks there that the cap is still the checkpoint's, and the script
##          stops otherwise.
##        - param_history is taken over from the checkpoint; `history` and the
##          checkpoint numbering restart at every fit() call: the earlier loss
##          history goes to checkpoints/resume_<E>/previous_history.csv and the
##          new checkpoints into that folder (epoch E + n).
##        - Schedule: lr_lambda(e) = lr_of(E + e) / lr.
##        - Adam starts fresh (model$optimizer <- NULL; FINN cannot carry the
##          optimizer state through torch_save), as every continuation;
##          seed = seed + E.
##        - The gradient-norm log and the cap trajectory cover the last segment only.
## RUN    sbatch --time=24:00:00 --mem=8G --dependency=afterany:<previous job> hpc/run_v020.sbatch calibration/v020/fit_continue.R
##        (as often as needed; exits at once when model.pt exists). hpc/submit_v020.sh chains it.
## OUTPUT as fit.R; checkpoints/resume_<E>/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-release-b6244e5"), .libPaths()))   # FINN release-0.2.0 at b6244e5 (as fit.R)
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v2-v020")

PARAMS <- list(                            # identical to fit.R
  prep_version  = "pairs-v2",              # 14 taxa, no "other", no shrubs, start-census site rule x = 2 %, intervals 4-12 yr
  train_set     = "train",                 # 1,179 training pairs in the retained 25 km blocks
  holdout_set   = "holdout_blocked",       # not read here; recorded for validation.R
  env_vars      = c("tavg_c", "prec_mm", "harvest"),   # WorldClim 2.1 normals (bio1, bio12) + the harvest flag of the interval (0/1)
  linear_var    = "harvest",               # outside the network: one linear per-species weight
  mgmt_levels   = c(0, 1),                 # the values the covariate takes
  patch_size_ha = pi * 20^2 / 1e4,         # one FINN patch = the 20 m NBI circle (0.1257 ha)
  n_patches     = 1L,
  env_years     = 12L,                     # env rows 1..12 per pair; obs at year = period_length (4-12)
  recruits_dbh  = 5.9,                     # recruits enter just above the 5 cm tally threshold
  growth_rate   = "annual",                # `growth` is the annualised relative dbh increment
  recruit_obs_weight_ha = pi * 10^2 / 1e4, # the median smaller circle of the training pairs (0.0314 ha)
  w_tol         = 1e-4,                    # |median smaller-circle area of the training pairs - recruit_obs_weight_ha| allowed
  hidden        = c(20L, 20L),             # climate network of growth and mortality (SELU)
  dropout       = 0.1,                     # on while training, off in simulateForest
  beta_init     = 0,                       # management weight start
  sen_max       = 5,                       # U-shape cap
  sen_init      = 0.05,                    # U-shape start
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01
  n_quantiles   = 10L,
  light_mask    = TRUE,
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # FINN's default families
  loss_weights  = "auto",                  # fit.R: FINN's own 1 / intercept-only baseline (here: the stored vector, raw)
  loss_recruit  = "simulated",
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # per tallied tree on the observed stand
  clip_norm     = 2,
  update_step   = 1L,
  epochs_const  = 300L,                    # epochs at constant lr
  anneal_epochs = 200L,                    # epochs of the cosine decay after them
  lr            = 0.01,
  lr_min        = 1e-4,
  checkpoints   = 25L,
  batchsize     = 100L,
  threads       = 1L,
  device        = "cpu",
  seed          = 1L
)
EPOCHS <- PARAMS$epochs_const + PARAMS$anneal_epochs
torch_set_num_threads(PARAMS$threads)
if (file.exists(file.path(FIT_DIR, "model.pt"))) {
  message("model.pt exists: the fit has finished, nothing to continue"); quit(save = "no")
}
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

## ---- the newest checkpoint ---------------------------------------------------
ck_files <- list.files(file.path(FIT_DIR, "checkpoints"), pattern = "^epoch_[0-9]+model\\.pt$", recursive = TRUE)
stopifnot(length(ck_files) > 0)
ck_dt <- data.table(file = ck_files, n = as.integer(sub(".*epoch_([0-9]+)model\\.pt$", "\\1", ck_files)), offset = 0L)
ck_dt[grepl("^resume_[0-9]+/", file), offset := as.integer(sub("^resume_([0-9]+)/.*", "\\1", file))]
ck_dt[, epoch := offset + n]
newest <- ck_dt[which.max(epoch)]
ck <- torch_load(file.path(FIT_DIR, "checkpoints", newest$file))
seg <- as.data.table(do.call(rbind, ck$history))[seq_len(newest$n), 1:6]
setnames(seg, c("dbh", "ba", "trees", "growth", "mort", "reg"))
prev <- if (newest$offset > 0) fread(file.path(FIT_DIR, "checkpoints", sprintf("resume_%d", newest$offset), "previous_history.csv"))[, 1:6] else NULL
hist_before <- rbind(prev, seg)
stopifnot(nrow(hist_before) == newest$epoch, length(ck$param_history) == newest$epoch / PARAMS$checkpoints)
weights_ck <- as.numeric(ck$loss_weights)  # fit.R's "auto": 1 / baseline, six terms
stopifnot(length(weights_ck) == 6L, all(is.finite(weights_ck)), all(weights_ck > 0))
K_prior_ck <- ck$reg_prior_value           # the prior fit.R resolved
k_done    <- newest$epoch
remaining <- EPOCHS - k_done
stopifnot(remaining >= 0L, is.numeric(K_prior_ck), length(K_prior_ck) == 1L)
message(sprintf("continuing from epoch %d (%s): %d of %d epochs to go", newest$epoch, newest$file, remaining, EPOCHS))

## ---- the learning-rate schedule (as fit.R), continued ---------------------------
lr_of <- function(t) ifelse(t < PARAMS$epochs_const, PARAMS$lr,
  PARAMS$lr_min + (PARAMS$lr - PARAMS$lr_min) * (1 + cos(pi * (t - PARAMS$epochs_const) / PARAMS$anneal_epochs)) / 2)
SCHED <- data.table(epoch = seq_len(EPOCHS))[, lr := lr_of(epoch - 1L)]

## ---- 1. data (as fit.R) ----------------------------------------------------------
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
site_tab <- merge(fread(file.path(train_dir, "siteID_dt.csv"))[, .(siteID, siteName)],
                  fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "site_table.csv"))[, .(siteName, r0, r1, area0_m2, area1_m2)], by = "siteName")
site_tab[, a_min_ha := pmin(area0_m2, area1_m2) / 1e4]
stopifnot(nrow(site_tab) == uniqueN(obs$siteID), abs(median(site_tab$a_min_ha) - w_obs) < PARAMS$w_tol)
oo <- merge(obs[reg > 0, .(siteID, species, reg)], site_tab[, .(siteID, a_min_ha)], by = "siteID")
whole_share <- oo[, mean(abs(reg * a_min_ha - round(reg * a_min_ha)) < 1e-3)]
y_q <- quantile(obs[reg > 0, reg * w_obs], c(0, 0.25, 0.5, 0.75, 1))
K_prior_expected <- max(obs$reg / obs$period_length)
stopifnot(abs(K_prior_ck - K_prior_expected) < 1e-6 * K_prior_expected)   # the same data resolve the same prior (per ha per year)
mgmt_share <- 100 * mean(env[year == 1][[PARAMS$linear_var]])
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
autoscale_spec <- c(setNames(rep("auto", length(setdiff(PARAMS$env_vars, PARAMS$linear_var))), setdiff(PARAMS$env_vars, PARAMS$linear_var)),
                    setNames("identity", PARAMS$linear_var))

## ---- 2. the env module (as fit.R) ----------------------------------------------------
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
scaling_chk <- compute_env_scaling(env, autoscale_spec)
X_chk <- model.matrix(env_formula, apply_env_scaling(env, scaling_chk)[round(seq(1, .N, length.out = 200L)), PARAMS$env_vars, with = FALSE])
check_shift <- function(nn, label) {       # the trained module is still a pure shift in management (as fit.R's check, eval mode)
  x <- torch::torch_tensor(unname(X_chk), dtype = torch::torch_float32())
  nn$eval()
  torch::with_no_grad({
    y2 <- as_array(nn(x)); beta <- as_array(nn$beta); dev <- 0
    for (t in PARAMS$mgmt_levels) {
      xt <- X_chk; xt[, mgmt_col] <- t
      yt <- as_array(nn(torch::torch_tensor(unname(xt), dtype = torch::torch_float32())))
      dev <- max(dev, max(abs(yt - y2 - outer(xt[, mgmt_col] - X_chk[, mgmt_col], beta)) / (1 + abs(y2))))
    }
  })
  stopifnot(dev < 1e-5)
  data.table(stage = "trained, reloaded (continued)", module = label, rows = nrow(X_chk), shift_max_rel_dev = dev, max_abs_output = max(abs(y2)))
}

## ---- 3. the model (as fit.R), then the checkpoint's tensors ----------------------------
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
  recruit_obs_weight = w_obs,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE,
                                       n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, NN = new_env_module(),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape, NN = new_env_module(),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_capped,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
FINN.seed(PARAMS$seed)
model <- build_model()
model$period_mode <- "per_site"
fit(model, env = env, data = obs, tree_data = trees, loss_aggregation = PARAMS$loss_aggregation,
    init_cohort = init_cohort, epochs = 0L, batchsize = PARAMS$batchsize, update_step = PARAMS$update_step,
    lr = PARAMS$lr, loss_weights = weights_ck, weight_scale = "raw", loss_family = PARAMS$loss_family, loss_recruit = PARAMS$loss_recruit,
    patches = PARAMS$n_patches, patch_size = PS, control = finnControl(optimizer = torch::optim_adam, clip_norm = PARAMS$clip_norm),
    device = PARAMS$device, env_autoscale = autoscale_spec, folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
stopifnot(setequal(names(model$state_dict()), names(ck$state_dict())),
          identical(names(model$parameters), names(ck$parameters)))
model$load_state_dict(ck$state_dict())
ck_sd <- ck$state_dict(); new_sd <- model$state_dict()
same <- vapply(names(ck_sd), function(k) identical(as_array(ck_sd[[k]]$cpu()), as_array(new_sd[[k]]$cpu())), logical(1))
if (!all(same)) stop("fresh model differs from the checkpoint in: ", paste(names(same)[!same], collapse = ", "))
stopifnot(inherits(model$nn_growth, "env_dnn_linmgmt"), inherits(model$nn_mortality, "env_dnn_linmgmt"),
          identical(model$reg_cap, "all"), identical(model$theta_floor, 0.5),
          all(vapply(model$parameters, function(p) isTRUE(p$requires_grad), logical(1))))   # nothing is fixed (fit.R header)
USED <- list(mortality_draw = model$mortality_draw, mortality_floor = model$mortality_floor, mortality_ceiling = model$mortality_ceiling,
             theta_floor = model$theta_floor, theta_optimized = model$theta_optimized, regeneration_sampler = model$regeneration_sampler,
             sample_regeneration = model$sample_regeneration, reg_floor = model$reg_floor, cap = model$reg_cap, prior = model$reg_prior,
             prior_lambda = model$reg_prior_lambda, light_mask = model$light_mask, trees_offset = DEFAULTS$trees_offset, rb_gradient = DEFAULTS$rb_gradient)
stopifnot(identical(model$mortality_draw, ck$mortality_draw), identical(model$regeneration_sampler, ck$regeneration_sampler),
          identical(model$mortality_floor, ck$mortality_floor), identical(model$mortality_ceiling, ck$mortality_ceiling),
          identical(model$theta_floor, ck$theta_floor))   # the same FINN defaults as the checkpoint
## the env scaling the checkpoint learned must be the one the fresh set-up learned (same data)
stopifnot(isTRUE(all.equal(as.data.frame(model$env_scaling), as.data.frame(ck$env_scaling), check.attributes = FALSE)))
logK_ck <- model$reg_logK$detach()$clone()  # the trained cap; fit() must start from it (header)
model$param_history <- ck$param_history
model$optimizer <- NULL                    # fit() keeps the set-up call's Adam otherwise; NULL makes it build optim_adam_gradlog
message(sprintf("fresh model: %d tensors identical to the checkpoint; cap K %.1f (prior %.1f)",
                length(same), exp(as.numeric(logK_ck$cpu())), K_prior_ck))
rm(ck, ck_sd, new_sd); invisible(gc())

## the scheduler factory: fit() calls it once, after the cap prior step and before
## the first step; it checks that FINN kept the trained cap (header)
CAPRESET <- new.env(); CAPRESET$K_found <- NA_real_
schedule <- function(optimizer) {
  CAPRESET$K_found <- exp(as.numeric(model$reg_logK$detach()$cpu()))
  stopifnot(abs(log(CAPRESET$K_found) - as.numeric(logK_ck$cpu())) < 1e-6)   # FINN kept the trained cap
  torch::lr_lambda(optimizer, lr_lambda = function(e) lr_of(k_done + e) / PARAMS$lr)
}

## ---- 4. the optimizer (as fit.R) ----------------------------------------------------
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
    if (!length(GRADLOG$K))                # the cap the first step starts from = the checkpoint's
      stopifnot(abs(as.numeric(ps[[match("reg_logK", GRADLOG$names)]]$detach()$cpu()) - as.numeric(logK_ck$cpu())) < 1e-7)
    out <- super$step(closure)
    GRADLOG$K <- c(GRADLOG$K, exp(as.numeric(ps[[match("reg_logK", GRADLOG$names)]]$detach()$cpu())))
    out
  })

## ---- 5. continue the fit --------------------------------------------------------------
min_fit <- 0
hist <- copy(hist_before)
seed_used <- NA_integer_
if (remaining > 0) {
  res_dir <- file.path(FIT_DIR, "checkpoints", sprintf("resume_%d", newest$epoch))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(hist_before, file.path(res_dir, "previous_history.csv"))
  seed_used <- PARAMS$seed + newest$epoch
  FINN.seed(seed_used)
  t0 <- Sys.time()
  fit(model,
      env              = env,
      data             = obs,
      tree_data        = trees,
      loss_aggregation = PARAMS$loss_aggregation,
      init_cohort      = init_cohort,
      epochs           = remaining,
      batchsize        = PARAMS$batchsize,
      update_step      = PARAMS$update_step,
      lr               = PARAMS$lr,        # the base lr the schedule multiplies; fresh Adam
      loss_weights     = weights_ck,       # fit.R's 1 / baseline
      weight_scale     = "raw",
      loss_family      = PARAMS$loss_family,
      loss_recruit     = PARAMS$loss_recruit,
      patches          = PARAMS$n_patches,
      patch_size       = PS,
      control          = finnControl(optimizer = optim_adam_gradlog,
                                     lr_scheduler = schedule,        # lr_of(k_done + e); checks the trained cap (header)
                                     clip_norm = PARAMS$clip_norm),
      device           = PARAMS$device,
      env_autoscale    = autoscale_spec,
      folder           = res_dir,
      checkpoints      = PARAMS$checkpoints,
      plot_progress    = FALSE)            # trees_offset, rb_gradient: FINN's defaults
  min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  stopifnot(abs(model$optimizer$param_groups[[1]]$lr - lr_of(EPOCHS)) < 1e-9, is.finite(CAPRESET$K_found))
  message(sprintf("cap at the start of training %.1f, the checkpoint's (prior %.1f)", CAPRESET$K_found, K_prior_ck))
  new <- as.data.table(do.call(rbind, model$history))[, 1:6]
  setnames(new, c("dbh", "ba", "trees", "growth", "mort", "reg"))
  hist <- rbind(hist_before, new)
}
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
fwrite(SCHED, file.path(FIT_DIR, "lr_schedule.csv"))
K_prior <- model$reg_prior_value
K_fit   <- exp(as.numeric(model$reg_logK$detach()$cpu()))
stopifnot(nrow(hist) == EPOCHS, abs(K_prior - K_prior_ck) < 1e-6 * K_prior_ck,
          abs(as.numeric(model$mort_sen_max$cpu()) - PARAMS$sen_max) < 1e-6,
          abs(as.numeric(model$mort_growth_scale$cpu()) - PARAMS$growth_scale) < 1e-9)

## ---- 6. diagnostics of the last segment (as fit.R) ----------------------------------
seg_epochs <- max(remaining, 1L)
steps_per_epoch <- NA_real_
if (length(GRADLOG$rows)) {
  G <- as.data.table(do.call(rbind, GRADLOG$rows)); setnames(G, GRADLOG$names)
  G[, step := .I]
  steps_per_epoch <- nrow(G) / seg_epochs
  G[, epoch := k_done + as.integer(ceiling(step / steps_per_epoch))]
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
  KT <- data.table(step = seq_along(GRADLOG$K), K = GRADLOG$K)[, epoch := k_done + as.integer(ceiling(step / steps_per_epoch))][
    , .(K_end = K[.N], K_min = min(K), K_max = max(K)), by = epoch][, K_prior := K_prior]
  fwrite(KT, file.path(FIT_DIR, "cap_trajectory.csv"))
} else GE <- data.table(epoch = integer(), group = character(), share_clipped = numeric(), reg_share_mean = numeric())

## ---- 7. checks and tables (as fit.R) ---------------------------------------------------
chk <- torch_load(file.path(FIT_DIR, "model.pt"))
stopifnot(inherits(chk$nn_growth, "env_dnn_linmgmt"), inherits(chk$nn_mortality, "env_dnn_linmgmt"))
CHECKS <- rbind(if (file.exists(file.path(FIT_DIR, "module_checks.csv"))) fread(file.path(FIT_DIR, "module_checks.csv")) else NULL,
                check_shift(chk$nn_growth, "growth"), check_shift(chk$nn_mortality, "mortality"), fill = TRUE)
fwrite(CHECKS, file.path(FIT_DIR, "module_checks.csv"))
rm(chk)
m_init <- torch_load(file.path(FIT_DIR, "model_init.pt"))
theta_start <- as.numeric(m_init$par_theta_recruits$detach()$cpu())
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(model[[block]]$detach()$cpu())), nm)
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")),
             spar("par_regeneration", "reg_light"))
SPP[, `:=`(sen = PARAMS$sen_max * plogis(as.numeric(model$mort_sen_raw$detach()$cpu())),
           theta_start = theta_start, theta_fit = as.numeric(model$par_theta_recruits$detach()$cpu()))]
SPP[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & mort_size < 0, -mort_size / (2 * sen) * 100, Inf)]
fwrite(SPP, file.path(FIT_DIR, "species_parameters.csv"))
MGMT <- merge(data.table(species = seq_len(n_species), growth_beta = as_array(model$nn_growth$beta$detach()),
                         mortality_beta = as_array(model$nn_mortality$beta$detach())), species_dt, by = "species")
fwrite(MGMT, file.path(FIT_DIR, "management_weights.csv"))
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
invisible(torch::with_no_grad(m_init$reg_logK$fill_(log(K_prior))))   # the start fit.R's fit() used: K = the resolved prior
m_init$env_scaling <- model$env_scaling    # and the env scaling fit() learned from the same data before its first step
ck_start <- sim_check(m_init); rm(m_init)
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
  sprintf("Fitted by `calibration/v020/fit.R` and continued by `fit_continue.R` to epoch %d; the last segment resumed at epoch %d (FINN %s from library `%s`, release-0.2.0; torch %s; repository %s; %.1f min for %d epochs, update_step %d, seed %s, %d thread(s) on %s). Adam restarted at the resume; the lr schedule continued; FINN kept the trained cap (K at the start of training %.1f).",
          EPOCHS, k_done, packageVersion("FINN"), basename(dirname(find.package("FINN"))), packageVersion("torch"), commit,
          min_fit, remaining, PARAMS$update_step, seed_used, PARAMS$threads, PARAMS$device, CAPRESET$K_found), "",
  sprintf("FINN defaults in use (none passed): %s.", paste(names(USED), vapply(USED, function(x) paste(format(x), collapse = ","), ""), sep = " = ", collapse = "; ")), "",
  sprintf("Data: `prepared/%s/%s` (%d pairs, %d species, intervals %d..%d yr, %d tracked trees of which %d lost (dead or felled); harvest recorded on %.1f %% of the pairs). Data convention kept: total stem loss incl. felling is the loss target; not found / now < 5 cm / outside forest censored. All six responses scored; recruitment through recruit_obs_weight = %.4f ha (the median smaller circle of the training pairs; %.1f %% of the recruiting rows whole on their own circle; reg x w quartiles %.2f / %.2f / %.2f). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, min(obs$period_length), max(obs$period_length), nrow(trees), sum(trees$died), mgmt_share,
          w_obs, 100 * whole_share, y_q[2], y_q[3], y_q[4], PARAMS$holdout_set), "",
  sprintf("Plain joint fit; FINN default initialisation (seed %d), theta trained above FINN's default floor %g, no warm start, no L2 other than the cap prior, no fed recruits. Growth and mortality env `%s` through env_dnn_linmgmt (hidden %s, dropout %g, + linear `%s` on its 0/1 scale, start %g, in growth AND mortality as the recipe); U-shaped mortality (sen in [0, %g], start %g at 1 m; growth / %g); light_mask %s, %d quantiles; regeneration FINN::regeneration_capped with a linear env response on the same formula; growth and mortality loss_aggregation = tree (scored on the observed tracked trees, period_mode per_site); growth_rate %s; recruits_dbh %g; %d patch of %.4f ha.",
          PARAMS$seed, model$theta_floor, paste(deparse(env_formula), collapse = ""), paste(PARAMS$hidden, collapse = "-"), PARAMS$dropout,
          PARAMS$linear_var, PARAMS$beta_init, PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, PARAMS$light_mask, PARAMS$n_quantiles,
          PARAMS$growth_rate, PARAMS$recruits_dbh, PARAMS$n_patches, PS), "",
  sprintf("Module checks (`module_checks.csv`): forcing %s to each of %s moves the trained, reloaded output by beta_s x delta with max relative deviation %.2g = a pure shift.",
          PARAMS$linear_var, paste(PARAMS$mgmt_levels, collapse = "/"), max(CHECKS[grepl("continued", stage), shift_max_rel_dev])), "",
  sprintf("Management weights (harvest 0 -> 1 shift; `management_weights.csv`): log growth %s; mortality logit %s.",
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$growth_beta), collapse = ", "),
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$mortality_beta), collapse = ", ")), "",
  sprintf("Recruitment cap: prior (max_obs_recruits, per ha per year) %.1f = start K; fitted K %.1f per ha per year (K / prior %.3f).", K_prior, K_fit, K_fit / K_prior), "",
  sprintf("Weights (fit.R's \"auto\", stored in the checkpoint, weight_scale raw): %s.", paste(signif(weights_ck, 4), collapse = ", ")), "",
  sprintf("Schedule (`lr_schedule.csv`): epochs 1-%d at lr %g, epochs %d-%d cosine %g -> %g. Plain Adam (torch::optim_adam), clip_norm %g, batch %d, checkpoints every %d epochs.",
          PARAMS$epochs_const, PARAMS$lr, PARAMS$epochs_const + 1L, EPOCHS, PARAMS$lr, PARAMS$lr_min, PARAMS$clip_norm, PARAMS$batchsize, PARAMS$checkpoints), "",
  sprintf("Baseline-scaled losses (mean of 10 epochs): start %s; end %s.", fmt(sum_l(1:10)), fmt(sum_l((nrow(hist) - 9):nrow(hist)))), "",
  sprintf("Predicted / observed recruitment on train over the fit windows: %.2f at the start (K = prior), %.2f after the fit (`regeneration_parameters.csv`). Theta: start %s, fit %s (floor %g).",
          ck_start$total, ck_end$total, paste(signif(range(theta_start), 3), collapse = "-"), paste(signif(range(SPP$theta_fit), 3), collapse = "-"), model$theta_floor), "",
  sprintf("Vanished stem predictions on train (positive observed count, simulated 0 / below trees_offset %g): start %d / %d, end %d / %d of %d.",
          DEFAULTS$trees_offset, ck_start$zero$n_zero, ck_start$zero$n_below_offset, ck_end$zero$n_zero, ck_end$zero$n_below_offset, ck_end$zero$n_pos), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm): %s (`species_parameters.csv`).",
          paste(sprintf("%.3f", range(SPP$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SPP$species_name, SPP$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Clipped gradients, last segment (`grad_norms_by_epoch.csv`, last 10 epochs): %s.",
          paste(GE[epoch > EPOCHS - 10, .(s = sprintf("%s clipped %.0f %% of steps, regeneration share %.2f", group[1], 100 * mean(share_clipped), mean(reg_share_mean))), by = group]$s, collapse = "; ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs; `convergence.csv`): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Validation report: `Rscript calibration/v020/validation.R`. Convergence: judge by the parameter drift between checkpoints, not by the loss."
), file.path(FIT_DIR, "fit_summary.md"))
message(sprintf("done: %.1f min", min_fit))
