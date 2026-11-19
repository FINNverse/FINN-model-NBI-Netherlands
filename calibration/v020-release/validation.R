#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v020-release/validation.R
##
## VARIANT v020-release (2026-11-19; see calibration/v020-release/fit.R's
##        header): calibration/v020-soil/validation.R for the v020-soil fit
##        refitted unchanged on the released FINN 0.2.0 (6c9156f,
##        ~/R/lib-finn-release-6c9156f, guarded); changed only the FINN
##        library, the fit folder fits/pairs-v3-soil-v020release and the script
##        names. Note: the release's fit() compounds the simulated yearly growth
##        and mortality over the window (FINN ac23615); the report below still
##        aggregates them as the arithmetic mean of the yearly rates (as
##        v020-soil, kept for comparability). One code change: predictTrees()
##        gets env only for the pairs with tracked trees (section 2), as the
##        release's predictTrees() requires equal site sets; under 20b184f the
##        per-tree predictions after such a pair were scored against the next
##        pair's env. The rest of this header (v020-soil's, then v020's) still
##        holds.
##
## VARIANT v020-soil (2026-11-14): a copy of calibration/v020/validation.R for
##        the soil variant (calibration/v020-soil/fit.R, fits/pairs-v3-soil-v020):
##        data version pairs-v3-soil, env_vars + soil_wet and soil_nutrient_poor,
##        FINN library ~/R/lib-finn-release-20b184f (guarded), and the bare-ground
##        succession run at every one of the six climate picks crossed with a 2 x 2
##        soil grid (eq_soil: soil_wet 1 = Gt VII, dry / 4 = Gt III, wet within the
##        rooting zone; soil_nutrient_poor 2 = enkeerd or lime-free clay / 6 =
##        haar- and holtpodzol, dune sand), 24 runs, each from an EMPTY stand (no
##        stems; never FINN's default start of one 1-cm stem per species and
##        patch). The support of every soil corner among the training pairs is
##        printed in the report. Everything else is v020's validation.
##
## WHAT   The validation report of the v020-defaults joint fit (calibration/
##        v020/fit.R + fit_continue.R, fits/pairs-v2-v020/model.pt): the v0.2
##        recipe validation (calibration/v02/recipe/validation.R) adapted to ONE
##        stage (one model file, one loss history, one parameter trajectory; no
##        stage A file, no parameters.pt). Runs on the saved model, so it can be
##        redone after a plot change without a refit. Two sets:
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
##        per-species Spearman, loss curves and parameter trajectories, map,
##        bare-ground succession with harvest = 0, metrics by harvest flag). No
##        natural-death AUC: tree_obs_dt of pairs-v2 carries the fate as `died`
##        (dead or felled) without the felling mark. Nothing is copied into the
##        model-card folder: which fit the card describes is decided later.
## INPUT  <FIT_DIR>/{model.pt, loss_history.csv, fit_summary.md};
##        prepared/pairs-v2/{train,holdout_blocked}
## OUTPUT <FIT_DIR>/report.md, <FIT_DIR>/validation/*.csv, *.png
## RUN    Rscript calibration/v020-release/validation.R   (from the repo root, 4 threads)
##        cluster: sbatch --time=02:00:00 --mem=16G --cpus-per-task=4 hpc/run_v020release.sbatch calibration/v020-release/validation.R
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
FINN_LIB <- path.expand("~/R/lib-finn-release-6c9156f")
.libPaths(c(FINN_LIB, .libPaths()))        # VARIANT v020-release: FINN 0.2.0 as released, 6c9156f (as fit.R)
suppressPackageStartupMessages({
  library(FINN); library(data.table); library(torch); library(ggplot2)
})
stopifnot(normalizePath(dirname(find.package("FINN"))) == normalizePath(FINN_LIB), packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "pairs-v3-soil-v020release")   # VARIANT v020-release: as fit.R
VAL      <- file.path(FIT_DIR, "validation")
dir.create(VAL, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "pairs-v3-soil",         # VARIANT: as fit.R
  sets          = list(test = "train", holdout = "holdout_blocked"),
  n_test_sites  = 1500L,                   # cap on the training pairs simulated for the report (all 1,179 here)
  subset_seed   = 7L,
  patch_size_ha = pi * 20^2 / 1e4,         # as fitted
  n_patches     = 1L,                      # as fitted: one plot = one patch
  env_vars      = c("tavg_c", "prec_mm", "soil_wet", "soil_nutrient_poor", "harvest"),   # as fitted (VARIANT: + soil)
  soil_vars     = c("soil_wet", "soil_nutrient_poor"),   # the soil columns, set per run in the bare-ground grid
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
  ## VARIANT: every climate pick x this soil grid (scores of the harmonised scales; the
  ## dry-poor corner is the modal Dutch forest site, the wet-poor corner the rarest)
  eq_soil       = data.table(soil_wet = c(1, 1, 4, 4), soil_nutrient_poor = c(2, 6, 2, 6),
                             soil_lab = c("dry, rich", "dry, poor", "wet, rich", "wet, poor")),
  eq_soil_tol   = 0.5,                     # |score - corner| counted as support of a soil corner among the training pairs
  n_top_species = 5L,
  threads       = 4L
)
torch_set_num_threads(PARAMS$threads)
PS <- PARAMS$patch_size_ha
ba_tree <- function(dbh) pi * (dbh / 200)^2

## ---- 0. model and data -------------------------------------------------------
model <- torch_load(file.path(FIT_DIR, "model.pt"))
model$eval()
stopifnot(inherits(model$nn_growth, "env_dnn_linmgmt"), inherits(model$nn_mortality, "env_dnn_linmgmt"),
          identical(model$reg_cap, "all"), length(model$param_history) >= 2)
loss_hist <- fread(file.path(FIT_DIR, "loss_history.csv"))[, stage := "joint fit"]
K_fit <- exp(as.numeric(model$reg_logK$detach()$cpu()))
message(sprintf("model.pt: %d epochs, %d parameter snapshots, cap K %.1f (prior %.1f), theta %s",
                nrow(loss_hist), length(model$param_history), K_fit, model$reg_prior_value,
                paste(signif(range(as.numeric(model$par_theta_recruits$detach()$cpu())), 3), collapse = "-")))
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
  ## VARIANT v020-release: env restricted to the pairs that have tracked trees (1 training, 2 holdout pairs have
  ## none). FINN 0.2.0's predictTrees() stops on unequal site sets (review fix B4); before, it indexed the env array
  ## by position, so the trees of every pair after such a gap were scored against the next pair's climate, soil and harvest
  P <- predictTrees(model, tree_data = d$trees, env = d$env[siteID %in% d$trees$siteID], patches = PARAMS$n_patches, patch_size = PS, device = "cpu")
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
  labs(y = "loss (baseline-scaled)", title = "Training loss per response: the joint fit (all six responses, lr 0.01 for 300 epochs, then cosine to 1e-4)")
ggsave(file.path(VAL, "fig1_loss.png"), g1, width = 12, height = 4.5, dpi = 130, bg = "white")
OPf <- rbind(OP[, .(set, variable, obs, pred, level = "site x species")], TOT[, .(set, variable, obs, pred, level = "stand total")])
OPf[, variable := factor(variable, levels = vars, labels = c("basal area (m2/ha)", "stems (/ha)", "mean dbh (cm)", "growth (rel./yr)", "stem loss (/yr)", "recruits (/ha/interval)"))]
fwrite(OPf, file.path(VAL, "obs_pred.csv"))   # the plotted pairs, so the figure can be redrawn without simulating
# every panel its own scale: x and y share the panel's limits (the 1:1 line stays the diagonal), cut at the
# 99.5th percentile of obs and pred so that a few extreme rows do not flatten the cloud; rows beyond sit on the edge
OPp <- OPf[is.finite(obs) & is.finite(pred)]
OPp[, `:=`(lo = min(0, obs, pred), hi = max(quantile(obs, 0.995), quantile(pred, 0.995))), by = .(level, variable)]
OPp[, `:=`(obs = pmin(obs, hi), pred = pmin(pred, hi))]
g2 <- ggplot(OPp, aes(obs, pred, colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = 0.15, size = 0.6) +
  geom_blank(aes(x = lo, y = lo)) + geom_blank(aes(x = hi, y = hi)) + theme(aspect.ratio = 1) +
  facet_wrap(level ~ variable, scales = "free", nrow = 2) + scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) +
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
## VARIANT: every climate pick x every soil corner is one bare-ground run (siteID = run)
pick[, clim := .I]
eq_grid <- CJ(clim = pick$clim, soil = seq_len(nrow(PARAMS$eq_soil)))
eq_grid <- merge(merge(eq_grid, pick[, .(clim, pair_siteID = siteID, tavg_c, prec_mm)], by = "clim"),
                 cbind(soil = seq_len(nrow(PARAMS$eq_soil)), PARAMS$eq_soil), by = "soil")
setorder(eq_grid, clim, soil)
eq_grid[, siteID := .I]
train_soil <- sets$test$env[year == 1L, c("siteID", PARAMS$soil_vars), with = FALSE]
eq_grid[, support_train_pairs := vapply(seq_len(.N), function(i) train_soil[abs(soil_wet - eq_grid$soil_wet[i]) <= PARAMS$eq_soil_tol &
                                                                          abs(soil_nutrient_poor - eq_grid$soil_nutrient_poor[i]) <= PARAMS$eq_soil_tol, .N], integer(1))]
eq_env <- eq_grid[, .(year = seq_len(PARAMS$eq_years)), by = .(siteID, tavg_c, prec_mm, soil_wet, soil_nutrient_poor)]
eq_env[, (PARAMS$mgmt_var) := PARAMS$eq_mgmt]
setcolorder(eq_env, c("siteID", "year", PARAMS$env_vars))
stopifnot(identical(names(eq_env), c("siteID", "year", PARAMS$env_vars)))
## bare ground is bare ground: an EMPTY stand (one cohort per patch with no stems, dbh 0),
## not simulateForest(init_cohort = NULL), which starts from FINN's default stand of one
## 1-cm stem of every species in every patch (on FIA that start let species grow at sites
## where they barely recruit)
n_eq <- uniqueN(eq_env$siteID)
empty_stand <- FINN::CohortMat(dims = c(n_eq, PARAMS$eq_patches, 1L), dbh = array(0, dim = c(n_eq, PARAMS$eq_patches, 1L)),
                               trees = array(0, dim = c(n_eq, PARAMS$eq_patches, 1L)), species = array(1L, dim = c(n_eq, PARAMS$eq_patches, 1L)),
                               sp = n_species)
stopifnot(sum(as.numeric(as_array(torch_tensor(empty_stand$trees)))) == 0)
FINN.seed(PARAMS$sim_seed)
eq <- simulateForest(model, env = eq_env, patches = PARAMS$eq_patches, patch_size = PS, init_cohort = empty_stand, device = "cpu")$long$site
eq <- merge(eq[variable %in% c("ba", "trees")], eq_grid[, .(siteID, clim, soil, tavg_c, prec_mm, soil_wet, soil_nutrient_poor, soil_lab)], by = "siteID")
eq <- merge(eq, species_dt, by = "species")
eq[, value := value / PS]
eq[, clim_lab := sprintf("%.1f C, %.0f mm", tavg_c, prec_mm)]
eq[, soil_lab := sprintf("%s (wet %g, poor %g)", soil_lab, soil_wet, soil_nutrient_poor)]
eq[, site_lab := paste(clim_lab, soil_lab, sep = "; ")]
top <- eq[variable == "ba", .(mx = max(value)), by = .(siteID, species_name)][order(siteID, -mx)][, head(.SD, PARAMS$n_top_species), by = siteID]
eq[, shown := paste(siteID, species_name) %in% paste(top$siteID, top$species_name)]
g5 <- ggplot() +
  geom_line(data = eq[variable == "ba" & !shown], aes(year, value, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = eq[variable == "ba" & shown], aes(year, value, colour = species_name), linewidth = 0.6) +
  geom_line(data = eq[variable == "ba", .(value = sum(value)), by = .(clim_lab, soil_lab, year)], aes(year, value), colour = "black", linewidth = 0.5, linetype = 2) +
  facet_grid(soil_lab ~ clim_lab) + labs(x = "year", y = "basal area (m2/ha)", colour = NULL,
  title = sprintf("Succession from bare ground (an empty stand), %s = %d (dashed: stand total); columns climate, rows soil", PARAMS$mgmt_var, PARAMS$eq_mgmt),
  subtitle = sprintf("%d patches of %.4f ha per run; the %d species with the highest basal area per run coloured", PARAMS$eq_patches, PS, PARAMS$n_top_species)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7), strip.text = element_text(size = 6.5))
ggsave(file.path(VAL, "fig5_bare_ground.png"), g5, width = 17, height = 10, dpi = 130, bg = "white")
ba_last <- eq[variable == "ba" & year > PARAMS$eq_years - 100, .(ba = mean(value)), by = .(siteID, species_name)]
dom <- ba_last[order(siteID, -ba)][, .(dominant = species_name[1], dominant_share = ba[1] / sum(ba), second = species_name[2]), by = siteID]
EQ <- eq[variable == "ba", .(ba_total = sum(value)), by = .(siteID, clim_lab, soil_lab, year)][
  , .(ba_last100 = mean(ba_total[year > PARAMS$eq_years - 100]), ba_prev100 = mean(ba_total[year > PARAMS$eq_years - 200 & year <= PARAMS$eq_years - 100])), by = .(siteID, clim_lab, soil_lab)]
EQ[, rel_change := (ba_last100 - ba_prev100) / ba_prev100]
EQ <- merge(EQ, eq[variable == "trees" & year > PARAMS$eq_years - 100, .(stems_last100 = sum(value) / 100), by = siteID], by = "siteID")
EQ <- merge(merge(EQ, dom, by = "siteID", all.x = TRUE), eq_grid[, .(siteID, support_train_pairs)], by = "siteID")
setorder(EQ, siteID)
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
PT <- traj(model$param_history, nrow(loss_hist), "joint")
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
    labs(x = "epoch", y = "value", title = "Mechanistic per-species parameters over training (one line per species; reg_logK shown as K)")
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
summ <- readLines(file.path(FIT_DIR, "fit_summary.md"))
provenance <- data.table(
  script = c("calibration/01_download_nbi_raw.R", "calibration/02_extract_nbi_tables.R", "calibration/v02/04_prepare_pairs_v2.R",
             "calibration/v020-soil/prepare_pairs_v3_soil.R",
             "calibration/v020-release/fit.R (+ fit_continue.R)", "calibration/v020-release/validation.R"),
  output = c("raw/ (NBI-7 database archive; calibration/PROVENANCE.md)", "raw/*.csv (extracted tables)",
             "prepared/pairs-v2/{national,train,holdout_blocked,pilot}, PREPARATION.md, site_table.csv, calibration/v02/expected_preparation_v2.csv",
             sprintf("prepared/%s/ = pairs-v2 + soil_wet, soil_nutrient_poor (env_dt.csv, soil_dt.csv, soil_plots.csv, PREPARATION.md)", PARAMS$prep_version),
             sprintf("%s/{model_init.pt, model.pt, checkpoints/, loss_history.csv, lr_schedule.csv, species_parameters.csv, regeneration_parameters.csv, management_weights.csv, module_checks.csv, grad_norms*.csv, cap_trajectory.csv, convergence.csv, fit_summary.md}", basename(FIT_DIR)),
             sprintf("%s/{report.md, validation/}", basename(FIT_DIR))))
writeLines(c(
  "# Validation report: FINN-model-NBI-Netherlands, v020-defaults joint fit + soil on the released FINN 0.2.0 (v020-release; linear management, soil_wet and soil_nutrient_poor)", "",
  sprintf("Written %s by `calibration/v020-release/validation.R` (FINN %s from library `%s`, torch %s).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          as.character(packageVersion("FINN")), basename(dirname(find.package("FINN"))), as.character(packageVersion("torch"))), "",
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
  "## Bare-ground succession", "", knitr::kable(EQ[, .(climate = clim_lab, soil = soil_lab, ba_last100 = fmt(ba_last100, 1), ba_prev100 = fmt(ba_prev100, 1), rel_change = fmt(rel_change),
                                                      stems_last100 = fmt(stems_last100, 0), dominant, dominant_share = fmt(dominant_share, 2), second, support_train_pairs)]), "",
  sprintf("%d patches of %.4f ha, %d years from bare ground (an empty stand, not FINN's default start of one 1-cm stem per species and patch) at the climates of six training pairs spanning the climate space, each crossed with %d soil corners (soil_wet x soil_nutrient_poor: %s), %s = %d; `rel_change` compares the mean stand basal area of the last 100 years with the 100 before; `dominant` = the species with the largest mean basal area over the last 100 years; `support_train_pairs` = training pairs whose two soil scores both lie within %g of the corner.",
          PARAMS$eq_patches, PS, PARAMS$eq_years, nrow(PARAMS$eq_soil), paste(sprintf("%g x %g", PARAMS$eq_soil$soil_wet, PARAMS$eq_soil$soil_nutrient_poor), collapse = ", "),
          PARAMS$mgmt_var, PARAMS$eq_mgmt, PARAMS$eq_soil_tol), "",
  "## Figures", "",
  "- `validation/fig1_loss.png` loss per response and stage; `fig2_obs_pred.png` observed vs predicted; `fig3_roc.png` per-tree ROC; `fig4_map.png` map of stand basal area; `fig5_bare_ground.png` succession; `fig6_parameters.png` parameter trajectories.", "",
  "## Fit", "", summ, "",
  "## Provenance", "", knitr::kable(provenance)
), file.path(FIT_DIR, "report.md"))
message("report: ", file.path(FIT_DIR, "report.md"))
