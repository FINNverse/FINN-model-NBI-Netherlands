#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v02/recipe/validation.R
##
## WHAT   The validation report of the zoo v0.2 recipe model (fit_stageA.R +
##        fit_stageB.R, fits/pairs-v2-recipe/model.pt). Runs on the saved model,
##        so it can be redone after a plot change without a refit. Two sets:
##          TEST    = pairs the model was fitted on (prepared/pairs-v2/train; all
##                    1,179, fewer than n_test_sites),
##          HOLDOUT = pairs the model never saw, in 25 km blocks that hold no
##                    training pair (prepared/pairs-v2/holdout_blocked, all 583).
##        The model has no coordinates as predictors, so the holdout is predicted
##        exactly like the test set. Every pair's first-visit stand is simulated
##        forward over its own interval (4-12 yr) with its own harvest flag and
##        aggregated the way fit() scores it (states at the second visit, growth
##        and stem loss as the mean annual rate, recruitment as the sum of the
##        annual per-ha rates). Same report as the Sweden template (metrics with
##        the persistence null and the change over the interval, rate calibration
##        ratios, END-census out-of-S strata, AUC per tree and per site,
##        per-species Spearman, loss curves and parameter trajectories of both
##        stages, map, bare-ground succession with harvest = 0, metrics by
##        harvest flag). No natural-death AUC: tree_obs_dt of pairs-v2 carries
##        the fate as `died` (dead or felled) without the felling mark. Nothing
##        is copied into the model-card folder: the recipe fit is a candidate;
##        which fit the card describes is decided later.
## INPUT  <FIT_DIR>/{model.pt, model_stageA.pt, parameters.pt, loss_history_stage{A,B}.csv,
##        fit_summary_stage{A,B}.md}; prepared/pairs-v2/{train,holdout_blocked}
## OUTPUT <FIT_DIR>/report.md, <FIT_DIR>/validation/*.csv, *.png
## RUN    Rscript calibration/v02/recipe/validation.R        (from the repo root, 4 threads)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({
  library(FINN); library(data.table); library(torch); library(ggplot2)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v2-recipe")
VAL      <- file.path(FIT_DIR, "validation")
dir.create(VAL, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "pairs-v2",
  sets          = list(test = "train", holdout = "holdout_blocked"),
  n_test_sites  = 1500L,                   # cap on the training pairs simulated for the report (all 1,179 here)
  subset_seed   = 7L,
  patch_size_ha = pi * 20^2 / 1e4,         # as fitted
  n_patches     = 1L,                      # as fitted: one plot = one patch
  env_vars      = c("tavg_c", "prec_mm", "harvest"),   # as fitted
  mgmt_var      = "harvest",               # the management covariate: metrics by its level, forced to eq_mgmt in the succession runs
  strata_var    = "oos_ba_share_t1",       # END-census out-of-S basal-area share of siteID_dt, binned below
  strata_breaks = c(0.01, 0.05),           # < 1 %, 1-5 %, > 5 % (the strata of PREPARATION.md)
  sim_batch     = 500L,                    # sites per simulateForest batch
  sim_seed      = 1L,
  ## bare-ground succession: pairs nearest to the 6 corners / centre of the
  ## (tavg_c, prec_mm) space of the training set, simulated for eq_years with the
  ## harvest flag at eq_mgmt (0 = no harvest, the near-natural case)
  eq_years      = 500L,
  eq_patches    = 20L,
  eq_mgmt       = 0L,
  eq_quantiles  = list(c(0.1, 0.1), c(0.1, 0.9), c(0.5, 0.5), c(0.9, 0.1), c(0.9, 0.9), c(0.3, 0.7)),
  n_top_species = 5L,
  threads       = 4L
)
torch_set_num_threads(PARAMS$threads)
PS <- PARAMS$patch_size_ha
ba_tree <- function(dbh) pi * (dbh / 200)^2

## ---- 0. model and data -------------------------------------------------------
model <- torch_load(file.path(FIT_DIR, "model.pt"))
model$eval()
stageA <- torch_load(file.path(FIT_DIR, "model_stageA.pt"))
parameters <- torch_load(file.path(FIT_DIR, "parameters.pt"))
sd_model <- model$state_dict()
stopifnot(setequal(names(parameters$state_dict), names(sd_model)),
          all(vapply(names(sd_model), function(n) identical(as.integer(parameters$state_dict[[n]]$shape), as.integer(sd_model[[n]]$shape)), logical(1))),
          max(vapply(names(sd_model), function(n) max(abs(as.numeric(parameters$state_dict[[n]]) - as.numeric(sd_model[[n]]))), numeric(1))) == 0)
message("parameters.pt: every tensor of model.pt present, same shapes and values")
loss_hist <- rbind(fread(file.path(FIT_DIR, "loss_history_stageA.csv"))[, stage := "A: demography"],
                   fread(file.path(FIT_DIR, "loss_history_stageB.csv"))[, stage := "B: regeneration"], fill = TRUE)
species_dt <- fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$sets$test, "species_dt.csv"))
n_species <- nrow(species_dt)

read_set <- function(set, ids = NULL) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  sites <- fread(file.path(d, "siteID_dt.csv"))
  sites[, oos_bin_end := cut(get(PARAMS$strata_var), c(-Inf, PARAMS$strata_breaks, Inf), labels = c("<1%", "1-5%", ">5%"))]
  out <- list(obs = fread(file.path(d, "obs_dt.csv"))[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)],
              env = fread(file.path(d, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE],
              init = fread(file.path(d, "init_cohort_dt.csv"))[, .(siteID, patchID, species, dbh, trees, treeName)],
              trees = fread(file.path(d, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)],
              sites = sites)
  if (!is.null(ids)) {
    key <- data.table(siteID = ids, new = seq_along(ids))
    out <- lapply(out, function(x) {
      x <- merge(x[siteID %in% ids], key, by = "siteID")
      x[, siteID := new][, new := NULL]
      setorder(x, siteID)
    })
  }
  out
}
set.seed(PARAMS$subset_seed)
n_train_all <- nrow(fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$sets$test, "siteID_dt.csv")))
test_ids <- sort(sample(seq_len(n_train_all), min(PARAMS$n_test_sites, n_train_all)))
sets <- list(test = read_set(PARAMS$sets$test, test_ids), holdout = read_set(PARAMS$sets$holdout))

## ---- 1. simulate every plot over its interval and aggregate like fit() ----------
simulate_set <- function(d) {
  ic <- makeInitCohorts(d$init, Nspecies = n_species)
  FINN.seed(PARAMS$sim_seed)
  sim <- simulateForest(model, env = d$env, patches = PARAMS$n_patches, patch_size = PS,
                        init_cohort = ic, device = "cpu", batchsize = PARAMS$sim_batch)$long$site
  per <- unique(d$obs[, .(siteID, period_length)])
  sim <- merge(sim, per, by = "siteID")
  sim <- sim[year <= period_length]
  w <- dcast(sim, siteID + year + species + period_length ~ variable, value.var = "value")
  w[, .(ba = ba[year == period_length], trees = trees[year == period_length], dbh = dbh[year == period_length],
        growth = mean(growth), mort = mean(mort), reg = sum(r_mean_ha)), by = .(siteID, species)]
}
pred <- lapply(sets, simulate_set)

vars <- c("ba", "trees", "dbh", "growth", "mort", "reg")
OP <- rbindlist(lapply(names(sets), function(s) {
  o <- melt(sets[[s]]$obs[, c("siteID", "species", vars), with = FALSE], id.vars = c("siteID", "species"),
            variable.name = "variable", value.name = "obs")
  p <- melt(pred[[s]][, c("siteID", "species", vars), with = FALSE], id.vars = c("siteID", "species"),
            variable.name = "variable", value.name = "pred")
  merge(o, p, by = c("siteID", "species", "variable"))[, set := s]
}))
OP <- merge(OP, species_dt, by = "species")
TOT <- rbindlist(lapply(names(sets), function(s) {
  o <- sets[[s]]$obs[, .(ba = sum(ba), trees = sum(trees), dbh = if (sum(trees) > 0) weighted.mean(dbh, trees, na.rm = TRUE) else NA_real_,
                          reg = sum(reg)), by = siteID]
  p <- pred[[s]][, .(ba = sum(ba), trees = sum(trees), dbh = if (sum(trees) > 0) weighted.mean(dbh, trees, na.rm = TRUE) else NA_real_,
                     reg = sum(reg)), by = siteID]
  q <- sets[[s]]$init[, .(ba = sum(trees * ba_tree(dbh)), trees = sum(trees), dbh = weighted.mean(dbh, trees), reg = NA_real_), by = siteID]
  m <- merge(melt(o, id.vars = "siteID", variable.name = "variable", value.name = "obs"),
             melt(p, id.vars = "siteID", variable.name = "variable", value.name = "pred"), by = c("siteID", "variable"))
  merge(m, melt(q, id.vars = "siteID", variable.name = "variable", value.name = "persist"), by = c("siteID", "variable"), all.x = TRUE)[, set := s]
}))

metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); o <- obs[ok]; p <- pred[ok]
  if (length(o) < 3) return(list(n = length(o), rmse = NA_real_, r2 = NA_real_, r = NA_real_, bias = NA_real_, obs_mean = NA_real_))
  list(n = length(o), rmse = sqrt(mean((o - p)^2)), r2 = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
       r = suppressWarnings(cor(o, p)), bias = mean(p - o), obs_mean = mean(o))
}
OP[variable %in% c("ba", "trees"), `:=`(obs = obs / PS, pred = pred / PS)]
TOT[variable %in% c("ba", "trees"), `:=`(obs = obs / PS, pred = pred / PS, persist = persist / PS)]
M_species <- OP[, metrics(obs, pred), by = .(set, variable)][, level := "species rows (site x species)"]
M_total   <- TOT[, metrics(obs, pred), by = .(set, variable)][, level := "stand totals"]
M_persist <- TOT[variable %in% c("ba", "trees", "dbh"), metrics(obs, persist), by = .(set, variable)][, level := "stand totals: persistence (end = start)"]
M_change  <- TOT[variable %in% c("ba", "trees"), metrics(obs - persist, pred - persist), by = .(set, variable)][, `:=`(level = "stand totals: change over the interval", variable = paste0("d_", variable))]
M <- rbind(M_species, M_total, M_persist, M_change)
spearman <- function(o, p) { ok <- is.finite(o) & is.finite(p); if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(o[ok], p[ok], method = "spearman")) }
M_bysp <- OP[, .(n = sum(is.finite(obs) & is.finite(pred)), spearman = spearman(obs, pred),
                 bias = mean(pred - obs, na.rm = TRUE), obs_mean = mean(obs, na.rm = TRUE)), by = .(set, variable, species_name)]
fwrite(M, file.path(VAL, "metrics.csv")); fwrite(M_bysp, file.path(VAL, "metrics_by_species.csv"))
BIN <- merge(TOT[set == "holdout" & variable %in% c("ba", "trees")], sets$holdout$sites[, .(siteID, oos_bin_end)], by = "siteID")
M_bin <- BIN[, c(metrics(obs, pred), list(rmse_persist = sqrt(mean((obs - persist)^2)))), by = .(variable, oos_bin_end)]
M_bin[, rmse_ratio_model_persist := rmse / rmse_persist]
fwrite(M_bin, file.path(VAL, "metrics_by_oos_bin.csv"))
## holdout metrics by harvest flag (stand totals)
MG <- merge(TOT[set == "holdout" & variable %in% c("ba", "trees")], sets$holdout$sites[, .(siteID, mgmt = get(PARAMS$mgmt_var))], by = "siteID")
M_mgmt <- MG[, c(metrics(obs, pred), list(rmse_persist = sqrt(mean((obs - persist)^2)))), by = .(variable, mgmt)]
fwrite(M_mgmt, file.path(VAL, "metrics_by_mgmt.csv"))

## ---- 2. rates: AUC per tree and per site, calibration ratios --------------------------
auc_weighted <- function(p, pos, neg) {
  o <- order(p); p <- p[o]; pos <- pos[o]; neg <- neg[o]
  cum_neg <- cumsum(neg) - neg
  ties_neg <- neg
  (sum(pos * (cum_neg + 0.5 * ties_neg))) / (sum(pos) * sum(neg))
}
AUC <- rbindlist(lapply(names(sets), function(s) {
  d <- sets[[s]]
  P <- predictTrees(model, tree_data = d$trees, env = d$env, patches = PARAMS$n_patches, patch_size = PS, device = "cpu")
  tr <- copy(d$trees)[, `:=`(growth_pred = P$growth_pred, mort_pred = P$mort_pred)]
  tree_auc <- tr[is.finite(mort_pred), auc_weighted(mort_pred, died, 1 - died)]
  g <- tr[is.finite(growth) & is.finite(growth_pred), metrics(growth, growth_pred)]
  site <- merge(pred[[s]][, .(siteID, species, mort)], d$obs[n_at_risk > 0, .(siteID, species, n_at_risk, n_died, period_length)],
                by = c("siteID", "species"))
  site[, mort_pred := 1 - (1 - mort)^period_length]
  site_auc <- site[is.finite(mort_pred), auc_weighted(mort_pred, n_died, n_at_risk - n_died)]
  fwrite(tr, file.path(VAL, sprintf("tree_predictions_%s.csv", s)))
  data.table(set = s, tree_auc = tree_auc, site_auc_weighted = site_auc, n_trees = nrow(tr), n_died = sum(tr$died),
             tree_growth_r2 = g$r2, tree_growth_rmse = g$rmse,
             deaths_ratio = site[is.finite(mort_pred), sum(mort_pred * n_at_risk) / sum(n_died)],
             growth_ratio = tr[is.finite(growth) & is.finite(growth_pred), sum(growth_pred) / sum(growth)],
             reg_ratio = OP[set == s & variable == "reg", sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE)])
}))
fwrite(AUC, file.path(VAL, "mortality_auc.csv"))
SP <- merge(OP[set == "holdout" & variable == "reg", .(reg_obs = sum(obs, na.rm = TRUE), reg_ratio = round(sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE), 2)), by = species_name],
            merge(pred$holdout[, .(siteID, species, mort)], sets$holdout$obs[n_at_risk > 0, .(siteID, species, n_at_risk, n_died, period_length)], by = c("siteID", "species"))[
              , .(n_died = sum(n_died), deaths_ratio = round(sum((1 - (1 - mort)^period_length) * n_at_risk) / sum(n_died), 2)), by = species][
              , species_name := species_dt$species_name[species]][, .(species_name, n_died, deaths_ratio)], by = "species_name", all.x = TRUE)[order(-reg_obs)]
fwrite(SP, file.path(VAL, "rates_by_species_holdout.csv"))

## ---- 3. figures ------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10))
LH <- melt(loss_hist[, .(stage, epoch, dbh, ba, trees, growth, mort, reg)], id.vars = c("stage", "epoch"))[is.finite(value) & value != 0]
g1 <- ggplot(LH, aes(epoch, value)) + geom_line(colour = "#1b6ca8") + facet_grid(stage ~ variable, scales = "free") +
  labs(y = "loss (baseline-scaled)", title = "Training loss per response: stage A (demography, observed recruitment) and stage B (regeneration)")
ggsave(file.path(VAL, "fig1_loss.png"), g1, width = 12, height = 4.5, dpi = 130, bg = "white")
OPf <- rbind(OP[, .(set, variable, obs, pred, level = "site x species")], TOT[, .(set, variable, obs, pred, level = "stand total")])
OPf[, variable := factor(variable, levels = vars, labels = c("basal area (m2/ha)", "stems (/ha)", "mean dbh (cm)", "growth (rel./yr)", "stem loss (/yr)", "recruits (/ha/interval)"))]
g2 <- ggplot(OPf[is.finite(obs) & is.finite(pred)], aes(obs, pred, colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = 0.15, size = 0.6) +
  facet_grid(level ~ variable, scales = "free") + scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) +
  labs(x = "observed", y = "predicted", title = "Observed vs predicted at the second visit (4-12 years)", colour = NULL) + theme(legend.position = "bottom")
ggsave(file.path(VAL, "fig2_obs_pred.png"), g2, width = 13, height = 5.5, dpi = 130, bg = "white")
roc_dt <- rbindlist(lapply(names(sets), function(s) {
  tr <- fread(file.path(VAL, sprintf("tree_predictions_%s.csv", s)))[is.finite(mort_pred)]
  setorder(tr, -mort_pred)
  tr[, .(tpr = cumsum(died) / sum(died), fpr = cumsum(1 - died) / sum(1 - died), set = s)][seq(1, .N, length.out = min(.N, 2000))]
}))
g3 <- ggplot(roc_dt, aes(fpr, tpr, colour = set)) + geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_path() +
  scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) + coord_equal() +
  labs(x = "false positive rate", y = "true positive rate", title = "Per-tree stem loss (death or felling): ROC",
       subtitle = paste(sprintf("%s AUC %.3f", AUC$set, AUC$tree_auc), collapse = "; "), colour = NULL)
ggsave(file.path(VAL, "fig3_roc.png"), g3, width = 5, height = 5, dpi = 130, bg = "white")
MAP <- rbindlist(lapply(names(sets), function(s) merge(TOT[set == s & variable == "ba"], sets[[s]]$sites[, .(siteID, lon = x, lat = y)], by = "siteID")))
MAPl <- melt(MAP[, .(set, lon, lat, observed = obs, predicted = pred)], id.vars = c("set", "lon", "lat"))
g4 <- ggplot(MAPl, aes(lon, lat, colour = pmin(value, 50))) + geom_point(size = 0.5) + facet_grid(set ~ variable) + coord_quickmap() +
  scale_colour_viridis_c(name = "basal area\n(m2/ha, capped 50)") + labs(title = "Stand basal area at the second visit", x = NULL, y = NULL)
ggsave(file.path(VAL, "fig4_map.png"), g4, width = 7, height = 9, dpi = 130, bg = "white")
train_sites <- sets$test$sites
pick <- unique(rbindlist(lapply(PARAMS$eq_quantiles, function(q) {
  tq <- quantile(train_sites$tavg_c, q[1]); pq <- quantile(train_sites$prec_mm, q[2])
  z <- train_sites[, (scale(tavg_c)[, 1] - (tq - mean(tavg_c)) / sd(tavg_c))^2 + (scale(prec_mm)[, 1] - (pq - mean(prec_mm)) / sd(prec_mm))^2]
  train_sites[which.min(z), .(siteID, tavg_c, prec_mm)]
})), by = "siteID")
eq_env <- pick[, .(year = seq_len(PARAMS$eq_years)), by = .(siteID, tavg_c, prec_mm)][, .(siteID, year, tavg_c, prec_mm)]
eq_env[, (PARAMS$mgmt_var) := PARAMS$eq_mgmt]
FINN.seed(PARAMS$sim_seed)
eq <- simulateForest(model, env = eq_env, patches = PARAMS$eq_patches, patch_size = PS, init_cohort = NULL, device = "cpu")$long$site
eq <- merge(eq[variable %in% c("ba", "trees")], pick, by = "siteID")
eq <- merge(eq, species_dt, by = "species")
eq[, value := value / PS]
eq[, site_lab := sprintf("%.1f C, %.0f mm", tavg_c, prec_mm)]
top <- eq[variable == "ba", .(mx = max(value)), by = .(siteID, species_name)][order(siteID, -mx)][, head(.SD, PARAMS$n_top_species), by = siteID]
eq[, shown := paste(siteID, species_name) %in% paste(top$siteID, top$species_name)]
g5 <- ggplot() +
  geom_line(data = eq[variable == "ba" & !shown], aes(year, value, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = eq[variable == "ba" & shown], aes(year, value, colour = species_name), linewidth = 0.6) +
  geom_line(data = eq[variable == "ba", .(value = sum(value)), by = .(site_lab, year)], aes(year, value), colour = "black", linewidth = 0.5, linetype = 2) +
  facet_wrap(~ site_lab, ncol = 3) + labs(x = "year", y = "basal area (m2/ha)", colour = NULL,
  title = sprintf("Succession from bare ground, %s = %d (dashed: stand total)", PARAMS$mgmt_var, PARAMS$eq_mgmt),
  subtitle = sprintf("%d patches of %.4f ha per site; the %d species with the highest basal area coloured", PARAMS$eq_patches, PS, PARAMS$n_top_species)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(VAL, "fig5_bare_ground.png"), g5, width = 11, height = 7, dpi = 130, bg = "white")
EQ <- eq[variable == "ba", .(ba_total = sum(value)), by = .(siteID, site_lab, year)][
  , .(ba_last100 = mean(ba_total[year > PARAMS$eq_years - 100]), ba_prev100 = mean(ba_total[year > PARAMS$eq_years - 200 & year <= PARAMS$eq_years - 100])), by = .(siteID, site_lab)]
EQ[, rel_change := (ba_last100 - ba_prev100) / ba_prev100]
EQ <- merge(EQ, eq[variable == "trees" & year > PARAMS$eq_years - 100, .(stems_last100 = sum(value) / 100), by = siteID], by = "siteID")
fwrite(EQ, file.path(VAL, "equilibrium.csv"))
traj <- function(ph, n_epochs, stage) {
  if (length(ph) < 2) return(NULL)
  step <- n_epochs / length(ph)
  rbindlist(lapply(seq_along(ph), function(i) {
    s <- ph[[i]]
    rbindlist(lapply(names(s), function(nm) {
      mat <- s[[nm]]
      data.table(stage = stage, epoch = i * step, block = nm, row = as.integer(row(mat)), col = as.integer(col(mat)), value = as.numeric(mat))
    }))
  }))
}
PT <- rbind(traj(stageA$param_history, loss_hist[stage == "A: demography", .N], "A"), traj(model$param_history, loss_hist[stage == "B: regeneration", .N], "B"))
if (!is.null(PT) && nrow(PT)) {
  PT[, kind := fifelse(grepl("^nn_", block), "environment network", "mechanistic")]
  fwrite(PT, file.path(VAL, "parameter_trajectories.csv"))
  setorder(PT, stage, block, row, col, epoch)
  PT[, delta := abs(value - shift(value)), by = .(stage, block, row, col)]
  moved <- PT[!is.na(delta), .(moved = max(delta) > 0), by = .(stage, block)][moved == TRUE]
  mech <- PT[kind == "mechanistic"][moved, on = .(stage, block)]
  mech[, panel := sprintf("%s: %s [, %d]", stage, sub("_unconstrained", "", block), col)]
  mech[block == "reg_logK", value := exp(value)]
  g6a <- ggplot(mech, aes(epoch, value, group = row)) + geom_line(alpha = 0.35, linewidth = 0.3, colour = "#1b6ca8") +
    facet_wrap(~ panel, scales = "free", ncol = 4) + theme_minimal(base_size = 8) +
    labs(x = "epoch", y = "value", title = "Mechanistic per-species parameters over training (one line per species; A = demography, B = regeneration)")
  conv <- PT[!is.na(delta)][moved, on = .(stage, block)][, .(mean_abs_change = mean(delta)), by = .(stage, epoch, block, kind)]
  g6b <- ggplot(conv, aes(epoch, mean_abs_change, colour = block)) + geom_line(linewidth = 0.4) + scale_y_log10() +
    facet_grid(kind ~ stage, scales = "free") + theme_minimal(base_size = 8) +
    theme(legend.text = element_text(size = 6), legend.key.height = unit(7, "pt")) +
    labs(x = "epoch", y = "mean |change| per snapshot (log)", title = "Are the parameters still moving?")
  ggsave(file.path(VAL, "fig6_parameters.png"), gridExtra::grid.arrange(g6a, g6b, ncol = 1, heights = c(1.3, 1)),
         width = 12, height = 11, dpi = 130, bg = "white")
}

## ---- 4. report -------------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "-", formatC(x, digits = d, format = "f"))
Mtab <- M[, .(set, level, variable, n, RMSE = fmt(rmse), R2 = fmt(r2), r = fmt(r), bias = fmt(bias), obs_mean = fmt(obs_mean))]
summ <- c(fread(file.path(FIT_DIR, "fit_summary_stageA.md"), sep = NULL, header = FALSE)$V1, "",
          fread(file.path(FIT_DIR, "fit_summary_stageB.md"), sep = NULL, header = FALSE)$V1)
provenance <- data.table(
  script = c("calibration/01_download_nbi_raw.R", "calibration/02_extract_nbi_tables.R", "calibration/v02/04_prepare_pairs_v2.R",
             "calibration/v02/recipe/fit_stageA.R", "calibration/v02/recipe/fit_stageB.R", "calibration/v02/recipe/validation.R"),
  output = c("raw/ (NBI-7 database archive; calibration/PROVENANCE.md)", "raw/*.csv (extracted tables)",
             sprintf("prepared/%s/{national,train,holdout_blocked,pilot}, PREPARATION.md, site_table.csv, calibration/v02/expected_preparation_v2.csv", PARAMS$prep_version),
             "fits/pairs-v2-recipe/{model_stageA.pt, loss_history_stageA.csv, species_parameters_stageA.csv, management_weights_stageA.csv, fit_summary_stageA.md}",
             "fits/pairs-v2-recipe/{model.pt, parameters.pt, loss_history_stageB.csv, regeneration_parameters.csv, fit_summary_stageB.md}",
             "fits/pairs-v2-recipe/{report.md, validation/}"))
writeLines(c(
  "# Validation report: FINN-model-NBI-Netherlands, zoo v0.2 recipe (stage A + B, linear management)", "",
  sprintf("Written %s by `calibration/v02/recipe/validation.R` (FINN %s, torch %s).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          as.character(packageVersion("FINN")), as.character(packageVersion("torch"))), "",
  "## Sets", "",
  sprintf("- **test**: %d of the %d training pairs (`prepared/%s/train`, random subset, seed %d).", length(test_ids), n_train_all, PARAMS$prep_version, PARAMS$subset_seed),
  sprintf("- **holdout**: %d pairs in 25 km blocks without any training pair (`prepared/%s/holdout_blocked`); the model has no coordinate predictors, so it is predicted exactly like the test set.",
          nrow(sets$holdout$sites), PARAMS$prep_version), "",
  "Each pair's first-visit stand is simulated over its own interval (4-12 years) with its own harvest flag; states are taken at the second visit, growth and stem loss as the mean annual rate, recruitment as the sum of the annual per-ha rates (the aggregation `fit()` uses). Basal area and stems are reported per ha. Persistence (end state = start state) and the change over the interval are reported beside the model. Stem loss is death or felling (the data convention of pairs-v2), so the per-tree AUC scores stem loss; the felled trees are not marked in tree_obs_dt, so there is no natural-death AUC.", "",
  "## Metrics", "", knitr::kable(Mtab), "",
  "## Holdout stand totals by END-census out-of-S basal-area stratum", "",
  knitr::kable(M_bin[, .(variable, oos_bin_end, n, RMSE = fmt(rmse), R2 = fmt(r2), rmse_persist = fmt(rmse_persist), rmse_ratio_model_persist = fmt(rmse_ratio_model_persist))]), "",
  sprintf("## Holdout stand totals by harvest flag (%s)", PARAMS$mgmt_var), "",
  knitr::kable(M_mgmt[, .(variable, mgmt, n, RMSE = fmt(rmse), R2 = fmt(r2), bias = fmt(bias), rmse_persist = fmt(rmse_persist))]), "",
  "## Rates: stem loss and growth per tree, calibration ratios", "",
  knitr::kable(AUC[, .(set, tree_auc = fmt(tree_auc), site_auc_weighted = fmt(site_auc_weighted), n_trees, n_died,
                       tree_growth_r2 = fmt(tree_growth_r2), tree_growth_rmse = fmt(tree_growth_rmse, 4), deaths_ratio = fmt(deaths_ratio, 2), growth_ratio = fmt(growth_ratio, 2), reg_ratio = fmt(reg_ratio, 2))]), "",
  "Per-tree predictions from `FINN::predictTrees` (the kernels on the observed stand); the site AUC uses the interval death probability implied by the simulated annual rate, each site x species row contributing `n_died` positives and `n_at_risk - n_died` negatives; the ratios are predicted / observed sums.", "",
  "## Per species on the holdout (recruitment and death ratios)", "", knitr::kable(SP), "",
  "## Per-species Spearman correlation (holdout)", "",
  knitr::kable(dcast(M_bysp[set == "holdout", .(species_name, variable, spearman = round(spearman, 2))], species_name ~ variable, value.var = "spearman")), "",
  "## Bare-ground succession", "", knitr::kable(EQ[, .(site_lab, ba_last100 = fmt(ba_last100, 1), ba_prev100 = fmt(ba_prev100, 1), rel_change = fmt(rel_change), stems_last100 = fmt(stems_last100, 0))]), "",
  sprintf("%d patches of %.4f ha, %d years from bare ground at six training pairs spanning the climate space, %s = %d; `rel_change` compares the mean stand basal area of the last 100 years with the 100 before.", PARAMS$eq_patches, PS, PARAMS$eq_years, PARAMS$mgmt_var, PARAMS$eq_mgmt), "",
  "## Figures", "",
  "- `validation/fig1_loss.png` loss per response and stage; `fig2_obs_pred.png` observed vs predicted; `fig3_roc.png` per-tree ROC; `fig4_map.png` map of stand basal area; `fig5_bare_ground.png` succession; `fig6_parameters.png` parameter trajectories of both stages.", "",
  "## Fit", "", summ, "",
  "## Provenance", "", knitr::kable(provenance)
), file.path(FIT_DIR, "report.md"))
message("report: ", file.path(FIT_DIR, "report.md"))
