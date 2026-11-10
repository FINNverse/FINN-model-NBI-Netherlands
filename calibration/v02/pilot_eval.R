## FINN-model-NBI-Netherlands | calibration/v02/pilot_eval.R
##
## Generic helpers for the v0.2 pilots (sourced by every calibration/v02/pilot_*.R):
## metric functions and ONE evaluation routine that simulates a fitted model over
## a prepared set, aggregates as fit() scores (state at year = period_length,
## mean annual rates, summed recruitment), and returns one row of holdout
## metrics next to the persistence (no-change) reference, plus the same
## stand-level metrics STRATIFIED by the end-census out-of-S basal-area share
## (bins < 1 %, 1-5 %, > 5 %; the agreed rare-species rule, point 3). Nothing in
## here decides what a pilot fits; it only measures.
##
## pilot_eval(model, set_dir, n_species, patch_size_ha, label, env_vars)
##   -> list(row = one-row data.table, strata = one row per end-census bin)
##   row: pilot, n_sites, ba_sp_r2, dbh_sp_r2, trees_sp_r2, ba_stand_r2, ba_persist_r2,
##        trees_stand_r2, trees_persist_r2, dba_r2, dn_r2, dn_bias_ha, growth_r, loss_auc,
##        loss_ratio (pred/obs lost trees), reg_ratio (pred/obs recruits), eq_ba_ha
##        (bare-ground BA after bare_years at the set's median climate, harvest 0,
##        10 patches), eq_top_species, eq_top_share
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })

metric_r2 <- function(o, p) { ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  sst <- sum((o - mean(o))^2); if (!length(o) || sst <= 0) return(NA_real_); 1 - sum((p - o)^2) / sst }
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg
  sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}

pilot_eval <- function(model, set_dir, n_species, patch_size_ha, label, env_vars, bare_years = 300L, seed = 1L) {
  model$eval()
  obs <- fread(file.path(set_dir, "obs_dt.csv")); env <- fread(file.path(set_dir, "env_dt.csv"))
  init <- fread(file.path(set_dir, "init_cohort_dt.csv")); species <- fread(file.path(set_dir, "species_dt.csv"))
  sites <- fread(file.path(set_dir, "siteID_dt.csv"))
  per <- unique(obs[, .(siteID, period_length)])
  FINN.seed(seed)
  sim <- simulateForest(model, env = env[, c("siteID", "year", env_vars), with = FALSE],
                        init_cohort = makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                        patches = 1L, patch_size = patch_size_ha, device = "cpu")
  S <- as.data.table(sim$wide$site); for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  S <- merge(S, per, by = "siteID")
  A  <- S[year >= 1L & year <= period_length, .(growth_m = mean(growth), mort_m = mean(mort), reg_m = sum(r_mean_ha)), by = .(siteID, species)]
  St <- S[year == period_length, .(siteID, species, dbh_m = dbh, ba_m = ba, trees_m = trees)]
  M <- merge(obs, merge(St, A, by = c("siteID", "species")), by = c("siteID", "species"))
  M[, mort_period_m := 1 - (1 - pmin(pmax(mort_m, 0), 1))^period_length]
  I <- init[, .(ba_0 = sum(dbh2ba(dbh) * trees), trees_0 = sum(trees)), by = .(siteID, species)]
  M <- merge(M, I, by = c("siteID", "species"), all.x = TRUE); M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
  TOT <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), trees = sum(trees), trees_m = sum(trees_m), trees_0 = sum(trees_0), reg = sum(reg), reg_m = sum(reg_m)), by = siteID]
  ## end-census out-of-S share per site (pairs-v2 and later carry it; older versions get one bin)
  if ("oos_ba_share_t1" %in% names(sites)) {
    TOT <- merge(TOT, sites[, .(siteID, oos = oos_ba_share_t1)], by = "siteID")
    TOT[, bin := cut(oos, c(-Inf, 0.01, 0.05, Inf), labels = c("<1%", "1-5%", ">5%"))]
  } else TOT[, bin := factor("all")]
  L <- M[n_at_risk > 0]
  ## bare ground at the set's median climate (harvest 0 if present), 10 patches
  med <- env[, lapply(.SD, median), .SDcols = env_vars]
  if ("harvest" %in% names(med)) med[, harvest := 0]
  env_bare <- data.table(siteID = 1L, year = seq_len(bare_years))[, (env_vars) := med]
  FINN.seed(seed)
  B <- as.data.table(simulateForest(model, env = env_bare, init_cohort = NULL, patches = 10L, patch_size = patch_size_ha, device = "cpu")$wide$site)
  for (cc in c("year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
  Bl <- B[year > bare_years - 50L, .(ba = sum(ba)), by = species]
  top_id <- Bl[which.max(ba), species]
  row <- data.table(pilot = label, n_sites = uniqueN(M$siteID),
             ba_sp_r2 = metric_r2(M$ba, M$ba_m), dbh_sp_r2 = metric_r2(M[trees > 0, dbh], M[trees > 0, dbh_m]), trees_sp_r2 = metric_r2(M$trees, M$trees_m),
             ba_stand_r2 = metric_r2(TOT$ba, TOT$ba_m), ba_persist_r2 = metric_r2(TOT$ba, TOT$ba_0),
             trees_stand_r2 = metric_r2(TOT$trees, TOT$trees_m), trees_persist_r2 = metric_r2(TOT$trees, TOT$trees_0),
             dba_r2 = metric_r2(TOT$ba - TOT$ba_0, TOT$ba_m - TOT$ba_0), dn_r2 = metric_r2(TOT$trees - TOT$trees_0, TOT$trees_m - TOT$trees_0),
             dn_bias_ha = TOT[, mean((trees_m - trees)) / patch_size_ha],
             growth_r = M[growth_n > 0, cor(growth, growth_m)],
             loss_auc = weighted_auc(L$mort_period_m, L$n_died, L$n_at_risk - L$n_died),
             loss_ratio = L[, sum(n_at_risk * mort_period_m) / sum(n_died)],
             reg_ratio = sum(TOT$reg_m) / sum(TOT$reg),
             eq_ba_ha = sum(Bl$ba) / 50 / patch_size_ha, eq_top_species = species[species == top_id, species_name], eq_top_share = Bl[species == top_id, ba] / sum(Bl$ba))
  strata <- TOT[, .(n_sites = .N, ba_stand_r2 = metric_r2(ba, ba_m), ba_persist_r2 = metric_r2(ba, ba_0),
                    ba_rmse_ha = sqrt(mean((ba_m - ba)^2)) / patch_size_ha, ba_bias_ha = mean(ba_m - ba) / patch_size_ha,
                    dn_bias_ha = mean(trees_m - trees) / patch_size_ha, reg_ratio = sum(reg_m) / max(sum(reg), 1e-9)), by = bin][order(bin)]
  strata[, pilot := label]
  list(row = row, strata = strata)
}

append_pilot_row <- function(res, file = file.path("calibration", "v02", "pilot_table.csv"), strata_file = file.path("calibration", "v02", "pilot_strata.csv")) {
  old <- if (file.exists(file)) fread(file) else NULL
  if (!is.null(old)) old <- old[pilot != res$row$pilot]
  fwrite(rbind(old, res$row, fill = TRUE), file)
  olds <- if (file.exists(strata_file)) fread(strata_file) else NULL
  if (!is.null(olds)) olds <- olds[pilot != res$row$pilot]
  fwrite(rbind(olds, res$strata, fill = TRUE), strata_file)
  print(res$row, digits = 3); print(res$strata, digits = 3)
}
