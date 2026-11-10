#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/v02/04_prepare_pairs_v2.R
##
## WHAT   Prepared-data version "pairs-v2": the pairs-v1 preparation
##        (calibration/03_prepare_pairs_v1.R) with the changes of the v0.2 round
##        (calibration/AUDIT_v02.md; the rare-species rule agreed with Max's side,
##        V02_BRIEF.md UPDATE 05:40 UTC):
##          1. NO LUMPED "OTHER" CLASS AND NO SHRUBS. The model population is a
##             species set S (genus-level NBI codes allowed where the NBI records
##             only the genus: Quercus robur/petraea, Betula, Acer, Populus, Salix;
##             the two Pinus nigra varieties are one species). Shrubs (NBI group
##             ST) and the NBI's own group codes ("other broadleaves", ...) are
##             never in S.
##          2. SITE RULE ON THE START CENSUS ONLY: a pair is kept when species
##             outside S hold < x % of the basal area at the FIRST visit. The end
##             census is never used for the selection (that would remove the
##             model's failures from the holdout); the out-of-S share at the
##             second visit is recorded per site (`oos_ba_share_t1`) so the
##             validation can stratify by it (< 1 %, 1-5 %, > 5 %).
##          3. TRADE-OFF TABLE before the choice: for |S| = 6..24 (species ranked
##             by national first-visit basal area) and x = 1/2/3/5 %, the pairs
##             kept, training / holdout sites, holdout blocks still populated,
##             climate coverage (5-95 % span of tavg_c and prec_mm of the kept
##             pairs relative to all pairs) and the identifiability columns of the
##             weakest species in S (plots with the species at the first visit,
##             trees at risk, holdout blocks with the species, climate span of
##             its plots). Written to calibration/v02/species_subset_tradeoff.csv;
##             the per-species table of the chosen S to species_identifiability.csv.
##             n (holdout blocks) and the span are COVERAGE choices, not measured
##             identifiability limits.
##          4. OUT-OF-S TREES on a retained site are dropped from the stand, the
##             fates, the state and the recruits (option (a) of the agreed rule;
##             at x <= 2 % the competition error is small).
##          5. RECRUITS BY BACK-MODELLED START DIAMETER: a new tree is a recruit
##             only if dbh - g75_species x interval < 5 cm (g75 = 75 % quantile of
##             the annual dbh increment of survivors below 12 cm at the first
##             visit); larger new trees stay in the state as ingrowth (audit 2).
##          6. HARVEST COVARIATE `harvest` (0/1: harvest activity recorded at the
##             second visit; NBI-7 oogstactiviteit, NBI-6 soort kap) in env_dt,
##             constant over the site's years (audit 4).
##          7. tree_obs_dt.csv: one row per tracked first-visit tree (dbh, annual
##             relative growth, died 0/1) for FINN's per-tree losses.
##        Everything else (pairs, patch, fates, intervals, blocks, pilot) is as in
##        pairs-v1.
## OUTPUT <DATA_DIR>/prepared/pairs-v2/PREPARATION.md, site_table.csv, holdout_separation.csv,
##        species.csv and, in national/, train/, holdout_blocked/, pilot/ (same layout each):
##        obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv;
##        calibration/v02/{PREPARATION.md, expected_preparation_v2.csv, species_subset_tradeoff.csv,
##        species_identifiability.csv} (raw layers and WorldClim rasters as in calibration/PROVENANCE.md)
## RUN    Rscript calibration/v02/04_prepare_pairs_v2.R   (from the repository root, ~3 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages({ library(data.table); library(terra) })

DATA_DIR      <- file.path("..", "..", "FINN-zoo-data", "nbi")    # zoo layout: data outside git
WORLDCLIM_DIR <- file.path("..", "..", "worldclim")                 # WorldClim 2.1 30 s rasters

PARAMS <- list(
  prep_version     = "pairs-v2",
  pairs            = list(c(t0 = "MFV", t1 = "NBI6"), c(t0 = "NBI6", t1 = "NBI7")),   # the two remeasurement pairs
  period_range     = c(3L, 15L),           # whole years between the visits; outside = odd dates, dropped
  min_dbh_cm       = 5,                    # NBI tally threshold (diameter_mm >= 50)
  patch_size_ha    = pi * 20^2 / 1e4,      # 0.1257 ha: the 20 m circle, the largest NBI radius, as the one FINN patch
  ## living / dead tree classes differ between the cycles: class 4 is "dead" in MFV and NBI-6, "suppressed" (alive) in NBI-7
  alive_class      = list(MFV = 1:3, NBI6 = 1:3, NBI7 = 1:4),
  dead_class       = list(MFV = 4:8, NBI6 = 4:8, NBI7 = 5:8),
  removed_class    = 0L,                   # "verwijderd": the tree is gone, oogst_status says why
  status_cut       = c(1L, 11L, 12L),      # felled (stem removed; felled in a deforestation; left standing in a deforestation)
  status_dead      = c(2L, 3L, 8L, 9L),    # dead removed, windthrown removed, fallen outside the plot, lying stem rotten < 5 cm
  status_censored  = c(4L, 13L),           # now below 5 cm; trees outside forest: out of the cohort
  forest_forms     = c("1a.1.1", "1a.1.2", "1a.2", "1c.1", "1c.2", "1c.3", "1c.5", "3"),   # stand types kept at the second visit
  max_ba_loss      = 0.5,                  # drop plots that lost >= 50 % of first-visit basal area (clearcut, conversion)
  small_tree_dbh   = 12,                   # cm: survivors below this at the first visit define the small-tree growth quantile per species
  recruit_growth_q = 0.75,                 # quantile of the small-tree annual dbh increment used to back-model a new tree's start dbh
  dbh_change_cm    = c(-2, 3),             # survivor dbh change below -2 cm or above 3 cm/yr x interval: growth NA, tree kept for mortality
  harvest_na_as    = 0L,                   # a second visit without a harvest record (NBI-7 NA) counts as "no harvest recorded"
  ## --- the species set S and the site rule (agreed rule; trade-off table written for every combination below) ---
  shrub_group      = "ST",                 # NBI species group "Struiken": never in S, their BA counts as out-of-S
  tradeoff_sizes   = 6:24,                 # |S| candidates: the top-N species by national first-visit basal area
  tradeoff_shares  = c(0.01, 0.02, 0.03, 0.05),   # x candidates
  species_set_size = 14L,                  # CHOSEN |S| (see PREPARATION.md and AUDIT_v02.md section 23 for the reasons)
  max_oos_share    = 0.02,                 # CHOSEN x: pairs kept when out-of-S species hold < 2 % of the first-visit basal area
  ## identifiability / coverage thresholds reported per species (coverage choices, not measured limits)
  ident_min_plots = 30L, ident_min_at_risk = 100L, ident_min_holdout_blocks = 3L, ident_min_span = 0.5,
  ## env: WorldClim 2.1 normals 1970-2000, 30 arc-seconds, at the centre of the plot's 1 km INSPIRE cell (EPSG:3035)
  env_rasters = c(tavg_c = "wc2.1_30s_bio/wc2.1_30s_bio_1.tif", prec_mm = "wc2.1_30s_bio/wc2.1_30s_bio_12.tif",
                  temp_season = "wc2.1_30s_bio/wc2.1_30s_bio_4.tif", prec_season = "wc2.1_30s_bio/wc2.1_30s_bio_15.tif",
                  tmax_warm_c = "wc2.1_30s_bio/wc2.1_30s_bio_5.tif", tmin_cold_c = "wc2.1_30s_bio/wc2.1_30s_bio_6.tif",
                  elev_m = "wc2.1_30s_elev.tif"),
  ## spatially blocked holdout: 25 km squares on the INSPIRE grid (EPSG:3035, metres); roles drawn over ALL pairs, before the site rule
  block_km = 25, holdout_block_share = 0.2, block_seed = 21L,
  pilot_n = 300L, pilot_seed = 1L          # a small training subset for smoke tests and pilots
)
PS  <- PARAMS$patch_size_ha
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(DATA_DIR, "raw", "csv")
rd  <- function(t) fread(file.path(RAW, paste0(t, ".csv")), encoding = "Latin-1")
log_lines <- character(0)
note <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s); invisible(s) }
parse_date <- function(x) as.IDate(substr(x, 1, 8), format = "%m/%d/%y")   # "07/24/13 00:00:00" as exported
ba_tree <- function(d) pi * (d / 200)^2

## ---- 1. the three cycles in one tree table and one plot table ----------------
sp_mfv <- rd("lookup_MFV_Boomsoort")                      # numeric MFV codes -> NBI two-letter codes
tm <- rd("data_MFV_boomgegevens")
tm <- merge(tm, sp_mfv[, .(CODE_BOOMSOORT = CODE, code = NBI)], by = "CODE_BOOMSOORT", all.x = TRUE)
trees <- rbind(
  tm[, .(cycle = "MFV", plot = PLOTNUMMER, tree = BOOMNUMMER, code, dbh = DIAMETER_MM / 10, klasse = CODE_BOOMKLASSE,
         status = OOGSTSTATUS, dist_m = AFSTAND / 10, nstam = NstamHA)],
  rd("data_NBI6_boomgegevens")[, .(cycle = "NBI6", plot = plotnummer, tree = boomnummer, code = boomsoort, dbh = diameter_mm / 10,
                                   klasse = boomklasse, status = oogst_status, dist_m = afstand / 10, nstam = NstamHA)],
  rd("data_NBI7_boomgegevens")[, .(cycle = "NBI7", plot = plotnummer, tree = boomnummer, code = boomsoort, dbh = diameter_mm / 10,
                                   klasse = boomklasse, status = oogst_status, dist_m = afstand / 10, nstam = NstamHA)])
trees[, alive := mapply(function(k, cy) k %in% PARAMS$alive_class[[cy]], klasse, cycle)]
trees[, dead  := mapply(function(k, cy) k %in% PARAMS$dead_class[[cy]],  klasse, cycle)]
pm <- rd("data_MFV_plotmetingen"); p6 <- rd("data_NBI6_plotmetingen"); p7 <- rd("data_NBI7_plotmetingen")
dm <- rd("data_MFV_plotdefinitie"); d6 <- rd("data_NBI6_plotdefinitie"); d7 <- rd("data_NBI7_plotdefinitie")
## harvest activity recorded at the visit (used for the SECOND visit of a pair):
## NBI-7 oogstactiviteit 2 = traces of harvest present, 1 = none; NBI-6 soort kap 1 = clearcut, 2 = single-tree felling, 3 = none
plots <- rbind(
  merge(pm[, .(plot = PLOTNUMMER, date = parse_date(OPNAMEDATUM), radius_m = STRAAL_HOUT, form = fifelse(VOLLEDIGE_OPNAME == "J", "MFV_complete", "MFV_incomplete"), harvest = NA_integer_)],
        dm[, .(plot = PLOTNUMMER, inspire = GRID)], by = "plot")[, cycle := "MFV"],
  merge(p6[, .(plot = plotnummer, date = parse_date(opnamedatum), radius_m = straal, form = as.character(verschijningsvorm),
               harvest = fcase(`soort kap` %in% c(1L, 2L), 1L, `soort kap` == 3L, 0L, default = NA_integer_))],
        d6[, .(plot, inspire = grid)], by = "plot")[, cycle := "NBI6"],
  merge(p7[, .(plot = plotnummer, date = parse_date(opnamedatum), radius_m = straal, form = as.character(verschijningsvorm),
               harvest = fcase(oogstactiviteit == 2L, 1L, oogstactiviteit == 1L, 0L, default = NA_integer_))],
        d7[, .(plot = PLOT, inspire = `INSPIRE grid`)], by = "plot")[, cycle := "NBI7"])
plots <- plots[!duplicated(plots[, .(cycle, plot)])]
## plot area (forest part of the circle) from the database's own per-hectare weights of the living trees inside
## the circle (NstamHA = 10,000 / area; uniform within a plot), falling back to the nominal circle
area <- trees[alive == TRUE & is.finite(nstam) & nstam > 0, .(area_m2 = 1e4 / median(nstam)), by = .(cycle, plot)]
plots <- merge(plots, area, by = c("cycle", "plot"), all.x = TRUE)
plots[!is.finite(area_m2) & is.finite(radius_m) & radius_m > 0, area_m2 := pi * radius_m^2]
note("tree rows: MFV %d, NBI-6 %d, NBI-7 %d; plots with a measured circle: MFV %d, NBI-6 %d, NBI-7 %d",
     trees[cycle == "MFV", .N], trees[cycle == "NBI6", .N], trees[cycle == "NBI7", .N],
     plots[cycle == "MFV" & radius_m > 0, .N], plots[cycle == "NBI6" & radius_m > 0, .N], plots[cycle == "NBI7" & radius_m > 0, .N])

## ---- 2. the remeasurement pairs = sites --------------------------------------
perm <- rbind(rd("data_permanent_MFV_NBI6")[, .(plot = plotnummer, t0 = "MFV", t1 = "NBI6")],
              rd("data_permanent_NBI6_NBI7")[, .(plot = plotnummer, t0 = "NBI6", t1 = "NBI7")])
perm <- merge(perm, plots[, .(t0 = cycle, plot, date0 = date, r0 = radius_m, area0 = area_m2, form0 = form, inspire0 = inspire)], by = c("plot", "t0"))
perm <- merge(perm, plots[, .(t1 = cycle, plot, date1 = date, r1 = radius_m, area1 = area_m2, form1 = form, inspire1 = inspire, harvest = harvest)], by = c("plot", "t1"))
note("harvest recorded at the second visit on %d pairs, none on %d, no record on %d (set to %d)",
     perm[harvest == 1L, .N], perm[harvest == 0L, .N], perm[is.na(harvest), .N], PARAMS$harvest_na_as)
perm[is.na(harvest), harvest := PARAMS$harvest_na_as]
note("remeasurement pairs listed by the database: %d (MFV -> NBI-6 %d, NBI-6 -> NBI-7 %d)", nrow(perm), perm[t0 == "MFV", .N], perm[t0 == "NBI6", .N])
perm <- perm[is.finite(r0) & r0 > 0 & is.finite(r1) & r1 > 0 & is.finite(area0) & is.finite(area1)]
note("  both visits with a measured circle: %d", nrow(perm))
perm <- perm[form0 != "MFV_incomplete"];                     note("  first visit complete (MFV flag): %d", nrow(perm))
perm <- perm[form1 %in% PARAMS$forest_forms];                note("  stand type at the second visit kept (%s): %d", paste(PARAMS$forest_forms, collapse = " "), nrow(perm))
perm[, period_length := as.integer(round(as.numeric(date1 - date0) / 365.25))]
perm <- perm[is.finite(period_length) & period_length >= PARAMS$period_range[1] & period_length <= PARAMS$period_range[2]]
note("  interval %d-%d whole years: %d", PARAMS$period_range[1], PARAMS$period_range[2], nrow(perm))
note("  pairs whose 1 km cell differs between the visits (neighbouring cells, second-visit cell used): %d", perm[inspire0 != inspire1, .N])
perm[, siteName := paste(plot, t0, t1, sep = "_")]

## ---- 3. first-visit stand, fates, second-visit stand, recruits (all codes) -----
t0 <- merge(trees, perm[, .(cycle = t0, plot, siteName, t1, r0, area0, period_length)], by = c("cycle", "plot"))
t0 <- t0[alive == TRUE & is.finite(dbh)]
n_small <- t0[dbh < PARAMS$min_dbh_cm, .N]; t0 <- t0[dbh >= PARAMS$min_dbh_cm]
n_out <- t0[dist_m > r0, .N]; t0 <- t0[dist_m <= r0 | is.na(dist_m)]
note("first-visit living trees on the pairs: %d (%d below 5 cm dropped, %d outside the first-visit circle dropped)", nrow(t0) + n_small + n_out, n_small, n_out)
t1 <- merge(trees, perm[, .(cycle = t1, plot, siteName, r1, area1)], by = c("cycle", "plot"))
t0 <- merge(t0, t1[, .(siteName, tree, alive1 = alive, dead1 = dead, klasse1 = klasse, status1 = status, dbh1 = dbh, in_t1 = TRUE)], by = c("siteName", "tree"), all.x = TRUE)
t0[is.na(in_t1), in_t1 := FALSE]
t0[, fate := fcase(
  in_t1 & alive1 & is.finite(dbh1),                                    "survivor",
  in_t1 & dead1,                                                       "dead",
  in_t1 & klasse1 == PARAMS$removed_class & status1 %in% PARAMS$status_dead, "dead",
  in_t1 & klasse1 == PARAMS$removed_class & status1 %in% PARAMS$status_cut,  "cut",
  default = "censored")]
fates <- t0[, .N, by = fate][order(-N)]
note("  fates: %s", paste(sprintf("%s %d", fates$fate, fates$N), collapse = ", "))
s1 <- t1[alive == TRUE & is.finite(dbh) & dbh >= PARAMS$min_dbh_cm & (dist_m <= r1 | is.na(dist_m))]
t0[, ef := 1e4 / area0]; s1[, ef := 1e4 / area1]
ba0 <- t0[, .(ba0_ha = sum(ef * ba_tree(dbh))), by = siteName]
ba1 <- s1[, .(ba1_ha = sum(ef * ba_tree(dbh))), by = siteName]
perm <- merge(merge(perm, ba0, by = "siteName", all.x = TRUE), ba1, by = "siteName", all.x = TRUE)
perm[is.na(ba0_ha), ba0_ha := 0]; perm[is.na(ba1_ha), ba1_ha := 0]
n_notree <- perm[ba0_ha == 0, .N]; perm <- perm[ba0_ha > 0];   note("  sites with a living first-visit tree >= 5 cm: %d (%d dropped)", nrow(perm), n_notree)
n_clear <- perm[ba1_ha < (1 - PARAMS$max_ba_loss) * ba0_ha, .N, by = t0]
perm <- perm[ba1_ha >= (1 - PARAMS$max_ba_loss) * ba0_ha]
note("  basal area loss below %.0f %%: %d (%s clearcut / converted dropped; median basal-area change %+.0f %%)", 100 * PARAMS$max_ba_loss, nrow(perm),
     paste(sprintf("%s -> %d", n_clear$t0, n_clear$N), collapse = ", "), 100 * median(perm$ba1_ha / perm$ba0_ha - 1))
t0 <- t0[siteName %in% perm$siteName]; s1 <- s1[siteName %in% perm$siteName]
## recruits: second-visit living trees without a first-visit row, inside the smaller circle, start dbh back-modelled below 5 cm
rec <- merge(s1, perm[, .(siteName, r0, period_length)], by = "siteName")
rec <- merge(rec, t0[, .(siteName, tree, had_t0 = TRUE)], by = c("siteName", "tree"), all.x = TRUE)[is.na(had_t0)]
rec[, r_min := pmin(r0, r1)]
n_rec_out <- rec[is.finite(dist_m) & dist_m > r_min, .N]
rec <- rec[is.na(dist_m) | dist_m <= r_min]
surv_small <- t0[fate == "survivor" & dbh < PARAMS$small_tree_dbh & is.finite(dbh1), .(inc = (dbh1 - dbh) / period_length, code)]
g75 <- surv_small[, .(g75 = quantile(inc, PARAMS$recruit_growth_q), n = .N), by = code]
g75_pooled <- surv_small[, quantile(inc, PARAMS$recruit_growth_q)]
rec <- merge(rec, g75[n >= 20, .(code, g75)], by = "code", all.x = TRUE)
rec[is.na(g75), g75 := g75_pooled]
rec[, start_dbh := dbh - g75 * period_length]
n_rec_big <- rec[start_dbh >= PARAMS$min_dbh_cm, .N]
rec <- rec[start_dbh < PARAMS$min_dbh_cm]
note("  new trees at the second visit: %d recruits (median dbh %.1f cm); %d beyond the smaller circle dropped, %d with a back-modelled start dbh >= %g cm (species g75 %.2f-%.2f cm/yr, pooled %.2f) = ingrowth, not recruitment",
     nrow(rec), median(rec$dbh), n_rec_out, n_rec_big, PARAMS$min_dbh_cm, g75[n >= 20, min(g75)], g75[n >= 20, max(g75)], g75_pooled)

## ---- 4. environment and blocks for ALL pairs (before the site rule) -----------------
## "1kmN3154E4036" = the ETRS89-LAEA (EPSG:3035) cell with lower-left corner N 3,154,000 / E 4,036,000 m
perm[, inspire := inspire1]
perm[, `:=`(x_laea = as.numeric(sub("^1kmN[0-9]+E([0-9]+)$", "\\1", inspire)) * 1000 + 500,
            y_laea = as.numeric(sub("^1kmN([0-9]+)E[0-9]+$", "\\1", inspire)) * 1000 + 500)]
stopifnot(all(is.finite(perm$x_laea)), all(is.finite(perm$y_laea)))
cells <- unique(perm[, .(inspire, x_laea, y_laea)])
v <- project(vect(cells, geom = c("x_laea", "y_laea"), crs = "EPSG:3035"), "EPSG:4326")
cells[, c("x", "y") := as.data.table(crds(v))]
for (nm in names(PARAMS$env_rasters)) {
  r <- rast(file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]]))
  cells[, (nm) := extract(r, v, ID = FALSE)[[1]]]
}
env_vars <- c(names(PARAMS$env_rasters), "x", "y", "harvest")
cells_ok <- cells[complete.cases(cells)]
note("env: WorldClim extracted for %d of %d cells (%d pairs in %d cells without a land pixel dropped)",
     nrow(cells_ok), nrow(cells), perm[!inspire %in% cells_ok$inspire, .N], nrow(cells) - nrow(cells_ok))
perm <- merge(perm[, -c("x_laea", "y_laea")], cells_ok, by = "inspire")
bm <- PARAMS$block_km * 1000
perm[, block := sprintf("E%03d_N%03d", floor(x_laea / bm), floor(y_laea / bm))]
blocks <- perm[, .(n_pairs = .N), by = block][order(block)]
set.seed(PARAMS$block_seed)
blocks[, role := "train"]
blocks[sample(.N, round(PARAMS$holdout_block_share * .N)), role := "holdout"]
perm <- merge(perm, blocks[, .(block, set = role)], by = "block")
note("blocks over all %d pairs: %d of %d %d-km blocks held out (seed %d); holdout blocks hold %d pairs", nrow(perm), blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_km, PARAMS$block_seed, perm[set == "holdout", .N])
t0 <- t0[siteName %in% perm$siteName]; s1 <- s1[siteName %in% perm$siteName]; rec <- rec[siteName %in% perm$siteName]

## ---- 5. species names, the ranking, the trade-off table --------------------------
sp7 <- rd("lookup_NBI7_Boomsoort")[, .(code, dutch = boomsoort, grp = Boomsoortgroep, sci = trimws(`Wetenschappelijke naam`))]
## eligible for S: a named tree taxon (genus allowed); not the NBI group codes ("Other ..."), not the shrub group, not open-phase / unvisited codes
sp7[, eligible := !(grepl("^[Oo]ther|^Onbezocht", sci) | sci == "" | is.na(sci) | code %in% c("KV", "XX") | grp == PARAMS$shrub_group)]
sp7[code %in% c("CD", "OD"), sci := "Pinus nigra"]
sp7[grepl(",", sci), sci := vapply(strsplit(sci, ",\\s*"), function(w) {
  g <- sub(" .*$", "", w); e <- sub("^\\S+ ", "", w)
  if (length(unique(g)) == 1L) paste(g[1], paste(e, collapse = "/")) else paste(w, collapse = " / ") }, character(1))]
sp7[, sci := sub("\\s+(var\\.|subsp\\.)\\s.*$", "", sci)]
t0 <- merge(t0, sp7[, .(code, sci, eligible)], by = "code", all.x = TRUE)
s1 <- merge(s1, sp7[, .(code, sci, eligible)], by = "code", all.x = TRUE)
rec <- merge(rec, sp7[, .(code, sci, eligible)], by = "code", all.x = TRUE)
for (dt in list(t0, s1, rec)) dt[is.na(eligible), `:=`(eligible = FALSE, sci = NA_character_)]
## national ranking by first-visit basal area (per-hectare weights), eligible taxa only
rank_tab <- t0[eligible == TRUE, .(ba_ha_sum = sum(ef * ba_tree(dbh)), stems = .N, plots = uniqueN(siteName)), by = .(species_name = sci)][order(-ba_ha_sum)]
rank_tab[, rank := .I]
note("eligible taxa among first-visit living trees: %d (shrub group and NBI group codes excluded); top 10 by basal area: %s",
     nrow(rank_tab), paste(rank_tab$species_name[1:10], collapse = ", "))
## per pair x taxon basal area at both visits (out-of-S share for any S follows from the ranking)
ba_t0 <- t0[, .(ba = sum(ef * ba_tree(dbh))), by = .(siteName, sci = fifelse(eligible, sci, "__notS__"))]
ba_t1 <- s1[, .(ba = sum(ef * ba_tree(dbh))), by = .(siteName, sci = fifelse(eligible, sci, "__notS__"))]
oos_share <- function(ba_tab, S) ba_tab[, .(oos = sum(ba[!sci %in% S]) / sum(ba)), by = siteName]
clim_all <- perm[, .(tavg_lo = quantile(tavg_c, 0.05), tavg_hi = quantile(tavg_c, 0.95), prec_lo = quantile(prec_mm, 0.05), prec_hi = quantile(prec_mm, 0.95))]
span <- function(x, lo, hi) if (length(x) < 2) 0 else min(1, (quantile(x, 0.95) - quantile(x, 0.05)) / (hi - lo))
identifiability <- function(S, kept) {
  ## per species in S over the KEPT pairs: plots with the species at the first visit, trees at risk, holdout blocks with it, climate span
  pres <- t0[siteName %in% kept$siteName & sci %in% S & eligible == TRUE]
  pres <- merge(pres, kept[, .(siteName, set, block, tavg_c, prec_mm)], by = "siteName")
  out <- pres[, .(plots = uniqueN(siteName), trees_at_risk = sum(fate %in% c("survivor", "dead", "cut")),
                  holdout_blocks = uniqueN(block[set == "holdout"]),
                  tavg_span = span(tavg_c[!duplicated(siteName)], clim_all$tavg_lo, clim_all$tavg_hi),
                  prec_span = span(prec_mm[!duplicated(siteName)], clim_all$prec_lo, clim_all$prec_hi)), by = .(species_name = sci)]
  merge(data.table(species_name = S), out, by = "species_name", all.x = TRUE)[is.na(plots), `:=`(plots = 0L, trees_at_risk = 0L, holdout_blocks = 0L, tavg_span = 0, prec_span = 0)][]
}
tradeoff <- rbindlist(lapply(PARAMS$tradeoff_sizes, function(N) rbindlist(lapply(PARAMS$tradeoff_shares, function(x) {
  S <- rank_tab$species_name[seq_len(min(N, nrow(rank_tab)))]
  o0 <- oos_share(ba_t0, S)
  kept <- merge(perm, o0[oos < x, .(siteName)], by = "siteName")
  id <- identifiability(S, kept)
  o1 <- oos_share(ba_t1, S)[siteName %in% kept$siteName]
  data.table(S_size = length(S), x_share = x, pairs_kept = nrow(kept), share_of_pairs = round(nrow(kept) / nrow(perm), 3),
             train = kept[set == "train", .N], holdout = kept[set == "holdout", .N],
             holdout_blocks_populated = kept[set == "holdout", uniqueN(block)], holdout_blocks_ge10 = kept[set == "holdout", .N, by = block][N >= 10, .N],
             tavg_coverage = round(span(kept$tavg_c, clim_all$tavg_lo, clim_all$tavg_hi), 3), prec_coverage = round(span(kept$prec_mm, clim_all$prec_lo, clim_all$prec_hi), 3),
             ba_share_of_S_national = round(rank_tab[seq_len(length(S)), sum(ba_ha_sum)] / t0[, sum(ef * ba_tree(dbh))], 3),
             min_plots = min(id$plots), min_trees_at_risk = min(id$trees_at_risk), min_holdout_blocks = min(id$holdout_blocks),
             min_tavg_span = round(min(id$tavg_span), 2), min_prec_span = round(min(id$prec_span), 2),
             species_below_thresholds = id[plots < PARAMS$ident_min_plots | trees_at_risk < PARAMS$ident_min_at_risk | holdout_blocks < PARAMS$ident_min_holdout_blocks |
                                            pmin(tavg_span, prec_span) < PARAMS$ident_min_span, .N],
             end_census_oos_gt5pct = round(mean(o1$oos > 0.05), 3), end_census_oos_1to5pct = round(mean(o1$oos >= 0.01 & o1$oos <= 0.05), 3))
}))))
fwrite(tradeoff, file.path("calibration", "v02", "species_subset_tradeoff.csv"))
note("trade-off table written: %d rows (|S| %d-%d x %s)", nrow(tradeoff), min(PARAMS$tradeoff_sizes), max(PARAMS$tradeoff_sizes), paste(100 * PARAMS$tradeoff_shares, collapse = "/"))

## ---- 6. apply the chosen S and x ----------------------------------------------------
S <- rank_tab$species_name[seq_len(PARAMS$species_set_size)]
o0 <- oos_share(ba_t0, S); o1 <- oos_share(ba_t1, S)
perm <- merge(perm, o0[, .(siteName, oos_ba_share_t0 = oos)], by = "siteName")
perm <- merge(perm, o1[, .(siteName, oos_ba_share_t1 = oos)], by = "siteName", all.x = TRUE)
perm[is.na(oos_ba_share_t1), oos_ba_share_t1 := 0]
n_before <- nrow(perm)
perm <- perm[oos_ba_share_t0 < PARAMS$max_oos_share]
note("site rule: |S| = %d (%s), x = %g %%: %d of %d pairs kept (%.1f %%); end-census out-of-S share on kept pairs: < 1 %% %d, 1-5 %% %d, > 5 %% %d",
     length(S), paste(S, collapse = ", "), 100 * PARAMS$max_oos_share, nrow(perm), n_before, 100 * nrow(perm) / n_before,
     perm[oos_ba_share_t1 < 0.01, .N], perm[oos_ba_share_t1 >= 0.01 & oos_ba_share_t1 <= 0.05, .N], perm[oos_ba_share_t1 > 0.05, .N])
ident <- identifiability(S, perm)
ident[, below_thresholds := plots < PARAMS$ident_min_plots | trees_at_risk < PARAMS$ident_min_at_risk | holdout_blocks < PARAMS$ident_min_holdout_blocks | pmin(tavg_span, prec_span) < PARAMS$ident_min_span]
fwrite(ident, file.path("calibration", "v02", "species_identifiability.csv"))
note("identifiability of S on the kept pairs: plots %d-%d, trees at risk %d-%d, holdout blocks %d-%d, climate span >= %.2f; %d species below the thresholds (%s)",
     min(ident$plots), max(ident$plots), min(ident$trees_at_risk), max(ident$trees_at_risk), min(ident$holdout_blocks), max(ident$holdout_blocks),
     min(pmin(ident$tavg_span, ident$prec_span)), sum(ident$below_thresholds), paste(ident[below_thresholds == TRUE, species_name], collapse = ", "))
## out-of-S trees leave the stand, the fates, the state and the recruits (option (a))
n_drop <- c(t0[siteName %in% perm$siteName & !(eligible & sci %in% S), .N], s1[siteName %in% perm$siteName & !(eligible & sci %in% S), .N], rec[siteName %in% perm$siteName & !(eligible & sci %in% S), .N])
t0 <- t0[siteName %in% perm$siteName & eligible == TRUE & sci %in% S]
s1 <- s1[siteName %in% perm$siteName & eligible == TRUE & sci %in% S]
rec <- rec[siteName %in% perm$siteName & eligible == TRUE & sci %in% S]
note("  out-of-S trees dropped on the kept pairs: %d first-visit, %d second-visit, %d recruits", n_drop[1], n_drop[2], n_drop[3])
species_dt <- data.table(species_name = sort(S)); species_dt[, species := .I]
for (dt in list(t0, s1, rec)) { dt[, species_name := sci]; dt[, species := species_dt$species[match(species_name, species_dt$species_name)]] }
sites <- perm[, .(siteName, plot, t0, t1, date0, date1, period_length, r0, r1, area0_m2 = area0, area1_m2 = area1, form1, inspire,
                  x_laea, y_laea, ba0_ha, ba1_ha, oos_ba_share_t0, oos_ba_share_t1, block, set, harvest, x, y,
                  tavg_c, prec_mm, temp_season, prec_season, tmax_warm_c, tmin_cold_c, elev_m)]
setorder(sites, plot, t0); sites[, siteID := .I]
id_of <- function(dt) { dt[, siteID := sites$siteID[match(siteName, sites$siteName)]]; dt }
t0 <- id_of(t0); s1 <- id_of(s1); rec <- id_of(rec)

## ---- 7. virtual trees on the 0.1257 ha patch -------------------------------------
largest_remainder <- function(x) {
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  fl <- floor(x); k <- as.integer(round(sum(x)) - sum(fl))
  if (k > 0L) { ord <- order(x - fl, decreasing = TRUE); fl[ord[seq_len(min(k, length(x)))]] <- fl[ord[seq_len(min(k, length(x)))]] + 1 }
  as.integer(fl)
}
t0[, n_virtual := largest_remainder(ef * PS), by = .(siteName, species)]
s1[, n_virtual := largest_remainder(ef * PS), by = .(siteName, species)]
note("expansion: first-visit basal area of S %.1f m2/ha design-weighted vs %.1f m2/ha from the virtual trees; %d of %d first-visit trees and %d of %d second-visit trees expand to zero",
     t0[, sum(ef * ba_tree(dbh))] / uniqueN(t0$siteName), t0[, sum(n_virtual * ba_tree(dbh))] / PS / uniqueN(t0$siteName),
     t0[n_virtual == 0L, .N], nrow(t0), s1[n_virtual == 0L, .N], nrow(s1))

## ---- 8. obs_dt, tree_obs_dt, env_dt, init_cohort_dt ---------------------------------
P_of <- sites[, .(siteID, period_length)]
state <- s1[, .(ba = sum(n_virtual * ba_tree(dbh)), trees = sum(n_virtual), dbh = weighted.mean(dbh, ef)), by = .(siteID, species)]
t0[, ok_growth := fate == "survivor" & is.finite(dbh1) & (dbh1 - dbh) >= PARAMS$dbh_change_cm[1] & (dbh1 - dbh) <= PARAMS$dbh_change_cm[2] * period_length]
t0[, growth_tree := fifelse(ok_growth & dbh1 / dbh > 0, (dbh1 / dbh)^(1 / period_length) - 1, NA_real_)]   # annualised relative dbh increment
rates <- t0[, {
  at_risk <- fate %in% c("survivor", "dead", "cut"); died <- fate %in% c("dead", "cut")
  g <- fifelse(ok_growth, dbh1 / dbh - 1, NA_real_)
  .(n_at_risk = sum(at_risk), n_died_t = sum(died), growth_n = sum(is.finite(g)),
    growth_rel = if (any(is.finite(g))) sum(ef[is.finite(g)] * g[is.finite(g)]) / sum(ef[is.finite(g)]) else NA_real_)
}, by = .(siteID, species)]
recr <- merge(rec, sites[, .(siteID, area1_m2)], by = "siteID")
recr[, forest_share := pmin(1, area1_m2 / (pi * r1^2))]
recr <- recr[, .(reg = .N / (pi * min(r0[1], r1[1])^2 / 1e4 * forest_share[1])), by = .(siteID, species)]
grid <- CJ(siteID = sites$siteID, species = species_dt$species)
obs_dt <- Reduce(function(a, b) merge(a, b, by = c("siteID", "species"), all.x = TRUE), list(grid, state, rates, recr))
stopifnot(!anyDuplicated(obs_dt[, .(siteID, species)]), nrow(obs_dt) == nrow(grid))
obs_dt <- merge(obs_dt, P_of, by = "siteID")
obs_dt[is.na(ba), `:=`(ba = 0, trees = 0L)]
obs_dt[is.na(reg), reg := 0]
obs_dt[is.na(n_at_risk), `:=`(n_at_risk = 0L, n_died_t = 0L, growth_n = 0L)]
obs_dt[, mort := fifelse(n_at_risk > 0, 1 - (1 - n_died_t / n_at_risk)^(1 / period_length), NA_real_)]
obs_dt[, n_died := as.numeric(n_died_t)]
obs_dt[, growth := fifelse(growth_n > 0 & is.finite(growth_rel) & growth_rel > -1, (1 + growth_rel)^(1 / period_length) - 1, NA_real_)]
obs_dt[growth_n == 0L, growth := NA_real_]
obs_dt[, year := period_length]
obs_dt <- merge(obs_dt, species_dt, by = "species")
obs_cols <- c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg", "species_name", "period_length")
obs_dt <- obs_dt[order(siteID, species), ..obs_cols]
## per tracked tree (FINN per-tree losses): first-visit trees that survived, died or were felled; growth NA where implausible
tree_obs_dt <- t0[fate %in% c("survivor", "dead", "cut"), .(siteID, year = period_length, patchID = 1L, species, dbh, growth = growth_tree,
                                                             died = as.integer(fate %in% c("dead", "cut")), period_length, trees = 1L)]
setorder(tree_obs_dt, siteID, species)
note("obs_dt: %d sites x %d species; intervals %d-%d yr; stem loss over the interval %.3f of the trees at risk (annual, pooled %.4f incl. felling); recruitment mean %.0f, median %.0f stems/ha per interval; tree_obs_dt %d tracked trees",
     uniqueN(obs_dt$siteID), nrow(species_dt), min(obs_dt$period_length), max(obs_dt$period_length), obs_dt[, sum(n_died) / sum(n_at_risk)],
     obs_dt[n_at_risk > 0, weighted.mean(mort, n_at_risk)], obs_dt[, sum(reg), by = siteID][, mean(V1)], obs_dt[, sum(reg), by = siteID][, median(V1)], nrow(tree_obs_dt))
Tmax <- max(sites$period_length)
env_dt <- sites[, .(year = seq_len(Tmax)), by = c("siteID", env_vars)][order(siteID, year), c("siteID", "year", env_vars), with = FALSE]
init_cohort_dt <- t0[n_virtual > 0L, .(siteID, patchID = 1L, year = 0L, treeName = paste(plot, tree, sep = "_"), species, species_name, dbh, trees = n_virtual)]
setorder(init_cohort_dt, siteID, species, treeName)
stopifnot(setequal(unique(init_cohort_dt$siteID), sites$siteID), setequal(unique(env_dt$siteID), sites$siteID), setequal(unique(obs_dt$siteID), sites$siteID))

## ---- 9. the split (blocks drawn in section 4), separation, site table, species table ----
set.seed(PARAMS$pilot_seed)
pilot_ids <- sites[set == "train", sample(siteID, PARAMS$pilot_n)]
sites[, pilot := siteID %in% pilot_ids]
setorder(sites, siteID)
note("sets after the site rule: %d training sites, %d holdout sites in %d of %d holdout blocks", sites[set == "train", .N], sites[set == "holdout", .N],
     sites[set == "holdout", uniqueN(block)], blocks[role == "holdout", .N])
tr_xy <- as.matrix(sites[set == "train", .(x_laea, y_laea)]); ho_xy <- as.matrix(sites[set == "holdout", .(x_laea, y_laea)])
nearest <- apply(outer(ho_xy[, 1], tr_xy[, 1], "-")^2 + outer(ho_xy[, 2], tr_xy[, 2], "-")^2, 1, function(d2) sqrt(min(d2)))
sep <- data.table(stat = c("median_km", "min_km", "q10_km", "q90_km"),
                  value = c(median(nearest), min(nearest), quantile(nearest, 0.1), quantile(nearest, 0.9)) / 1000)
note("holdout separation: nearest training site median %.1f km, minimum %.1f km, 10 %% below %.1f km", sep$value[1], sep$value[2], sep$value[3])
fwrite(sep, file.path(OUT, "holdout_separation.csv"))
fwrite(sites, file.path(OUT, "site_table.csv"))
species_out <- merge(species_dt, t0[, .(nbi_code = paste(sort(unique(code)), collapse = "+"), stems_t0 = .N, sites_t0 = uniqueN(siteName)), by = species_name], by = "species_name")
species_out <- merge(species_out, sp7[, .(dutch_name = paste(unique(dutch), collapse = " / ")), by = .(species_name = sci)], by = "species_name", all.x = TRUE)
species_out <- merge(species_out, rank_tab[, .(species_name, rank_ba_national = rank)], by = "species_name")[order(species)]
fwrite(species_out[, .(species_name, species, nbi_code, dutch_name, stems_t0, sites_t0, rank_ba_national)], file.path(OUT, "species.csv"))

## ---- 10. write the sets (siteIDs renumbered 1..n within each set) --------------------
write_set <- function(ids, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ids <- sort(ids); map <- data.table(siteID_national = ids, siteID_new = seq_along(ids))
  renum <- function(dt) { dt <- merge(dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national"); dt[, siteID := siteID_new][, siteID_new := NULL]; dt }
  o <- renum(obs_dt); e <- renum(env_dt); i <- renum(init_cohort_dt); tr <- renum(tree_obs_dt)
  s <- merge(sites[siteID %in% ids, .(siteID, siteName, plot, t0, t1, period_length, inspire, x, y, tavg_c, prec_mm, harvest, oos_ba_share_t0, oos_ba_share_t1, block, set)], map, by.x = "siteID", by.y = "siteID_national")
  s[, siteID_national := siteID][, siteID := siteID_new][, siteID_new := NULL]
  fwrite(o[order(siteID, species), ..obs_cols], file.path(d, "obs_dt.csv"))
  fwrite(e[order(siteID, year), c("siteID", "year", env_vars), with = FALSE], file.path(d, "env_dt.csv"))
  fwrite(i[order(siteID, species, treeName), .(siteID, patchID, year, treeName, species, species_name, dbh, trees)], file.path(d, "init_cohort_dt.csv"))
  fwrite(tr[order(siteID, species), .(siteID, year, patchID, species, dbh, growth, died, period_length, trees)], file.path(d, "tree_obs_dt.csv"))
  fwrite(species_dt[, .(species, species_name)], file.path(d, "species_dt.csv"))
  fwrite(s[order(siteID), .(siteID, siteID_national, siteName, plot, t0, t1, period_length, inspire, x, y, tavg_c, prec_mm, harvest, oos_ba_share_t0, oos_ba_share_t1, block, set)], file.path(d, "siteID_dt.csv"))
  note("set %-16s %5d sites, %6d obs rows, %6d cohorts, %6d tracked trees", name, length(ids), nrow(o), nrow(i), nrow(tr))
}
write_set(sites$siteID, "national")
write_set(sites[set == "train", siteID], "train")
write_set(sites[set == "holdout", siteID], "holdout_blocked")
write_set(sites[pilot == TRUE, siteID], "pilot")

## ---- 11. PREPARATION.md ----------------------------------------------------------------
site_tot <- obs_dt[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
chosen <- tradeoff[S_size == PARAMS$species_set_size & x_share == PARAMS$max_oos_share]
md <- c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by calibration/v02/04_prepare_pairs_v2.R (data.table %s, terra %s).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("data.table"), packageVersion("terra")), "",
  "Source: Database Zevende Nederlandse Bosinventarisatie (NBI-7), Stichting Probos / Wageningen Environmental Research,",
  "archive as served 2026-11-09 (checksums in raw/SHA256SUMS.txt and calibration/PROVENANCE.md); MFV 2001-2005, NBI-6 2012-2013,",
  "NBI-7 2017-2021 tables. Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the centre of the plot's 1 km INSPIRE cell.", "",
  "## The population this version describes (the agreed rare-species rule)", "",
  sprintf("- Species set S (|S| = %d, ranked by national first-visit basal area; genus-level NBI codes kept as the NBI records them): %s.", length(S), paste(S, collapse = ", ")),
  sprintf("- Site rule on the START census only: a plot pair is kept when species outside S (incl. the NBI shrub group) hold < %g %% of the",
          100 * PARAMS$max_oos_share),
  sprintf("  first-visit basal area: %d of %d pairs (%.1f %%), i.e. the model is valid for stands where the modelled species hold > %g %% of the",
          chosen$pairs_kept, n_before, 100 * chosen$share_of_pairs, 100 - 100 * PARAMS$max_oos_share),
  "  basal area; with one NBI plot per 100 ha the share of pairs is also the share of forest area these stands cover.",
  sprintf("- The end-census out-of-S share is NOT used for the selection; it is recorded per site (`oos_ba_share_t1`: < 1 %% on %d sites, 1-5 %% on %d, > 5 %% on %d) for a stratified validation.",
          sites[oos_ba_share_t1 < 0.01, .N], sites[oos_ba_share_t1 >= 0.01 & oos_ba_share_t1 <= 0.05, .N], sites[oos_ba_share_t1 > 0.05, .N]),
  "- Out-of-S trees on kept sites are dropped from the stand, the fates, the state and the recruits (option (a)).",
  sprintf("- No lumped \"other\" class and no shrubs. The trade-off table behind the choice: calibration/v02/species_subset_tradeoff.csv (|S| %d-%d x %s %%);",
          min(PARAMS$tradeoff_sizes), max(PARAMS$tradeoff_sizes), paste(100 * PARAMS$tradeoff_shares, collapse = "/")),
  "  per-species identifiability of S: species_identifiability.csv. The thresholds there (plots, trees at risk, holdout blocks, climate span) are",
  "  coverage choices, not measured identifiability limits.", "",
  "## What the tables mean", "",
  "- One site = one remeasurement PAIR of a permanent NBI plot (MFV -> NBI-6 or NBI-6 -> NBI-7; a plot measured three times is two",
  "  sites), one patch of 0.1257 ha (the 20 m circle). Every tallied tree stands for 10,000 / plot area stems per ha and is expanded",
  "  to round(that x 0.1257) virtual trees (largest remainder within site x species). `ba` (m2), `trees` and `dbh` (cm) in obs_dt and",
  "  `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.1257 for per-hectare values.",
  "- `year` 0 = first visit, `year` = `period_length` = second visit (interval rounded to whole years, 4-12). Env rows years 1..Tmax.",
  "- `growth` = annualised relative dbh increment of survivors, weighted by the plot's expansion; `growth_n` = tallied survivors behind it.",
  "- `n_at_risk` = tallied first-visit trees that survived, died or were felled (censored: not found, now < 5 cm, outside forest);",
  "  `n_died` = dead + FELLED trees (stem loss; harvest is part of the mortality target); `mort` = annual rate.",
  "- `reg` = new trees >= 5 cm whose back-modelled start dbh (dbh - species g75 x interval) is below 5 cm, per hectare over the",
  "  interval, counted inside the smaller of the two circles; larger new trees are ingrowth (in the state, not recruits).",
  "- tree_obs_dt: one row per tracked first-visit tree (siteID, year = period_length, patchID, species, dbh at the first visit,",
  "  growth = annualised relative dbh increment or NA, died 0/1 incl. felled, period_length, trees = 1) for FINN's per-tree losses.",
  "- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), temp_season (bio4), prec_season (bio15, CV %),",
  "  tmax_warm_c (bio5), tmin_cold_c (bio6), elev_m, x = longitude, y = latitude (WGS84) of the 1 km cell centre,",
  "  harvest = 1 if harvest activity was recorded at the second visit (NBI-7 traces of harvest; NBI-6 clearcut or single-tree felling), else 0.", "",
  "## Processing log", "", paste0("- ", log_lines), "",
  "## Summary (national)", "",
  sprintf("- %d sites (%d MFV -> NBI-6, %d NBI-6 -> NBI-7); mean basal area at the second visit %.1f m2/ha, mean stems %.0f/ha,",
          nrow(sites), sites[t0 == "MFV", .N], sites[t0 == "NBI6", .N], mean(site_tot$ba) / PS, mean(site_tot$trees) / PS),
  sprintf("  mean annual relative growth %.4f, pooled annual stem loss %.4f (incl. felling), mean recruitment %.0f stems/ha per interval",
          obs_dt[growth_n > 0, weighted.mean(growth, growth_n)], obs_dt[n_at_risk > 0, weighted.mean(mort, n_at_risk)], mean(site_tot$reg)), "",
  "## Sets", "",
  "- national: every site passing the filters and the site rule", "- train: every such site in the training blocks",
  "- holdout_blocked: every such site in the held-out blocks (never for fitting)",
  sprintf("- pilot: %d training sites drawn with seed %d (smoke tests and pilots)", PARAMS$pilot_n, PARAMS$pilot_seed), "",
  sprintf("Spatial blocking: %d km squares on EPSG:3035, %d of %d blocks held out (seed %d, drawn over all pairs before the site rule). Holdout site to",
          PARAMS$block_km, blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_seed),
  sprintf("nearest training site: median %.1f km, minimum %.1f km, 10 %% below %.1f km, 90 %% below %.1f km.", sep$value[1], sep$value[2], sep$value[3], sep$value[4]), "",
  "Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set,",
  "siteID_national, plot, cycles, period_length, INSPIRE cell, lon/lat, tavg_c, prec_mm, harvest, oos_ba_share_t0/t1, block, set).",
  "Site attributes: site_table.csv; species list with NBI codes, Dutch names and the national BA rank: species.csv.", "",
  "## Species", "", "| species | name | NBI code | Dutch | stems at the first visit | sites | BA rank |", "|---|---|---|---|---|---|---|",
  species_out[, sprintf("| %d | %s | %s | %s | %d | %d | %d |", species, species_name, nbi_code, dutch_name, stems_t0, sites_t0, rank_ba_national)], "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(OUT, "PREPARATION.md"))
writeLines(md, file.path("calibration", "v02", "PREPARATION.md"))

## ---- 12. committed summary: expected_preparation_v2.csv ----------------------------------
per_set <- rbindlist(lapply(c("national", "train", "holdout_blocked", "pilot"), function(s) {
  o <- fread(file.path(OUT, s, "obs_dt.csv")); i <- fread(file.path(OUT, s, "init_cohort_dt.csv")); tr <- fread(file.path(OUT, s, "tree_obs_dt.csv"))
  st <- o[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
  data.table(set = s, sites = as.numeric(uniqueN(o$siteID)), obs_rows = nrow(o), cohorts = nrow(i), virtual_trees_t0 = sum(i$trees), tracked_trees = nrow(tr),
             ba_ha_t1_mean = round(mean(st$ba) / PS, 3), trees_ha_t1_mean = round(mean(st$trees) / PS, 2),
             growth_annual_mean = round(weighted.mean(o$growth, o$growth_n, na.rm = TRUE), 6),
             n_at_risk = sum(o$n_at_risk), n_died = round(sum(o$n_died), 2), mort_annual_pooled = round(o[n_at_risk > 0, weighted.mean(mort, n_at_risk)], 6),
             reg_ha_interval_mean = round(mean(st$reg), 4))
}))
per_species <- obs_dt[, .(set = "national", sites_present_t1 = as.numeric(sum(trees > 0)), ba_sum = round(sum(ba), 3), n_at_risk = as.numeric(sum(n_at_risk)),
                          n_died = round(sum(n_died), 2), reg_sum = round(sum(reg), 2)), by = .(species, species_name)]
expected <- rbind(melt(per_set, id.vars = "set", variable.name = "quantity", value.name = "value")[, `:=`(species = NA_integer_, species_name = NA_character_)],
                  melt(per_species, id.vars = c("set", "species", "species_name"), variable.name = "quantity", value.name = "value"), use.names = TRUE)
fwrite(expected[, .(set, species, species_name, quantity, value)], file.path("calibration", "v02", "expected_preparation_v2.csv"))
message("done: ", OUT, "; calibration/v02/{PREPARATION.md, expected_preparation_v2.csv, species_subset_tradeoff.csv, species_identifiability.csv}")
