#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-NBI-Netherlands | calibration/03_prepare_pairs_v1.R
##
## WHAT   Turns the extracted NBI tables (02_extract_nbi_tables.R) into FINN's
##        calibration tables for the tree-by-tree remeasurement of the Dutch
##        national forest inventory: obs_dt (site x species at the end of the
##        interval: ba, dbh, trees, growth, mortality counts, recruitment),
##        env_dt (WorldClim 2.1 normals per site and year), init_cohort_dt (the
##        first-visit stand), species_dt, siteID_dt, a spatially blocked
##        train / holdout split, PREPARATION.md and expected_preparation.csv.
##        Prepared-data version "pairs-v1".
##
## THE DATA
##   The NBI is a 1 point per 100 ha sample with one plot per 1 km cell. Half
##   of the MFV (2001-2005) and NBI-6 (2012-2013) points were permanent; from
##   NBI-7 (2017-2021) on all are. On a permanent point every tree is numbered
##   and positioned, and every tree of the previous visit has to be reported
##   back (alive, dead, removed), also when it lies outside a shrunken circle.
##   The database lists the plots whose later visit is a true remeasurement:
##   1,237 MFV -> NBI-6 pairs (7-12 years apart) and 1,413 NBI-6 -> NBI-7 pairs
##   (4-9 years). Each PAIR is one FINN site with its own interval length (a plot
##   measured three times gives two sites; FINN scores each on its own).
##
## THE DESIGN AND HOW IT BECOMES A FINN PATCH
##   One circle per plot whose radius (5-20 m) is chosen so that at least 20
##   trees >= 5 cm dbh are inside; the forest part of the circle is the plot
##   area, and the database gives every tree its per-hectare weight NstamHA
##   (= 10,000 / plot area). FINN needs fixed-area patches with integer tree
##   counts, so -- as the IGN-France zoo model does for its concentric circles --
##   the plot is declared ONE patch of 0.1257 ha (the 20 m circle, the largest
##   the NBI uses) and every tallied tree is expanded into round(NstamHA x
##   0.1257) identical virtual trees (largest-remainder rounding within plot x
##   species, so species totals are preserved). Densities and basal areas per
##   hectare are exact; the sample sizes behind the rates (n_at_risk, growth_n)
##   are the TALLIED trees.
##
## TREE FATES (first-visit living trees inside the first-visit circle, matched by
## plot and tree number):
##   survivor  alive at the second visit with a measured dbh -> growth
##   dead      dead class at the second visit (standing dead, broken, lying), or
##             "removed" with harvest status dead / windthrow / rotten -> loss
##   cut       "removed" with harvest status felled -> LOSS as well: felling and
##             natural death cannot be separated tree by tree in a managed forest
##             and censoring the felled trees biases the natural rate (dying trees
##             are felled pre-emptively); the decision of the BWI and IFN-Spain
##             zoo models. Plots that lost half or more of their basal area are
##             dropped as clearcut / converted stands.
##   censored  no second-visit row, "removed" with status "now < 5 cm" / "outside
##             forest", or no status: out of the cohort, not a loss.
##   recruits  second-visit living trees without a first-visit row, >= 5 cm and
##             not larger than 5 cm + 1.5 cm/yr x interval (a larger "new" tree was
##             missed at the first visit), inside the SMALLER of the two radii
##             (a wider second circle would otherwise book its ring as ingrowth).
##
## TIME   One interval per site: year 0 = first visit, year = period_length =
##        second visit (interval rounded to whole years). Growth and mortality
##        are annualised (FINN scores the mean of its annual rates over the
##        window), recruitment is the per-hectare total over the interval (FINN
##        sums). Env rows for years 1..Tmax for every site (static WorldClim
##        normals; the model has no year-to-year climate).
##
## INPUT  <DATA_DIR>/raw/csv/data_{MFV,NBI6,NBI7}_{boomgegevens,plotmetingen,plotdefinitie}.csv,
##        data_permanent_{MFV_NBI6,NBI6_NBI7}.csv, lookup_{MFV,NBI7}_Boomsoort.csv
##        <WORLDCLIM_DIR>/wc2.1_30s_bio/wc2.1_30s_bio_{1,4,5,6,12,15}.tif, wc2.1_30s_elev.tif
## OUTPUT <DATA_DIR>/prepared/pairs-v1/PREPARATION.md, site_table.csv, holdout_separation.csv,
##        species.csv and, in national/, train/, holdout_blocked/, pilot/ (same layout each):
##        obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv;
##        calibration/PREPARATION.md, calibration/expected_preparation.csv,
##        climate section of calibration/PROVENANCE.md
## RUN    Rscript calibration/03_prepare_pairs_v1.R   (from the repository root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages({ library(data.table); library(terra) })

DATA_DIR      <- file.path("..", "..", "FINN-zoo-data", "nbi")    # zoo layout: data outside git
WORLDCLIM_DIR <- file.path("..", "..", "worldclim")                 # WorldClim 2.1 30 s rasters

PARAMS <- list(
  prep_version     = "pairs-v1",
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
  ## dropped: 1a.1.3 open phase (clearcut), 1b.* park / avenue / hedgerow / tree row, 1c.4 landscape planting,
  ## 2a coppice, 2b osier bed, 91-96 no forest / not visited / deforested / not in selection
  max_ba_loss      = 0.5,                  # drop plots that lost >= 50 % of first-visit basal area (clearcut, conversion)
  recruit_growth   = 1.5,                  # cm/yr: a "new" tree above 5 cm + 1.5 x interval was missed at the first visit, not ingrowth
  dbh_change_cm    = c(-2, 3),             # survivor dbh change below -2 cm or above 3 cm/yr x interval: growth NA, tree kept for mortality
  species_quantile = 0.95,                 # species covering 95 % of living tallied first-visit stems; the rest -> "other"
  ## env: WorldClim 2.1 normals 1970-2000, 30 arc-seconds, at the centre of the plot's 1 km INSPIRE cell (EPSG:3035).
  ## Model covariates are tavg_c and prec_mm (the zoo-wide names); the others ride along for the general zoo model.
  env_rasters = c(tavg_c = "wc2.1_30s_bio/wc2.1_30s_bio_1.tif", prec_mm = "wc2.1_30s_bio/wc2.1_30s_bio_12.tif",
                  temp_season = "wc2.1_30s_bio/wc2.1_30s_bio_4.tif", prec_season = "wc2.1_30s_bio/wc2.1_30s_bio_15.tif",
                  tmax_warm_c = "wc2.1_30s_bio/wc2.1_30s_bio_5.tif", tmin_cold_c = "wc2.1_30s_bio/wc2.1_30s_bio_6.tif",
                  elev_m = "wc2.1_30s_elev.tif"),
  ## spatially blocked holdout: 25 km squares on the INSPIRE grid (EPSG:3035, metres); the Netherlands span ~300 x 300 km
  block_km = 25, holdout_block_share = 0.2, block_seed = 21L,
  pilot_n = 300L, pilot_seed = 1L          # a small training subset for smoke tests
)
PS  <- PARAMS$patch_size_ha
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(DATA_DIR, "raw", "csv")
rd  <- function(t) fread(file.path(RAW, paste0(t, ".csv")), encoding = "Latin-1")
log_lines <- character(0)
note <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s); invisible(s) }
parse_date <- function(x) as.IDate(substr(x, 1, 8), format = "%m/%d/%y")   # "07/24/13 00:00:00" as exported

## ---- 1. the three cycles in one tree table and one plot table ----------------
## tree columns: cycle, plot, tree, code (NBI two-letter species code), dbh (cm), klasse, status, dist_m
## plot columns: cycle, plot, date, radius_m, form (stand type), inspire (1 km cell)
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
plots <- rbind(
  merge(pm[, .(plot = PLOTNUMMER, date = parse_date(OPNAMEDATUM), radius_m = STRAAL_HOUT, form = fifelse(VOLLEDIGE_OPNAME == "J", "MFV_complete", "MFV_incomplete"))],
        dm[, .(plot = PLOTNUMMER, inspire = GRID)], by = "plot")[, cycle := "MFV"],
  merge(p6[, .(plot = plotnummer, date = parse_date(opnamedatum), radius_m = straal, form = as.character(verschijningsvorm))],
        d6[, .(plot, inspire = grid)], by = "plot")[, cycle := "NBI6"],
  merge(p7[, .(plot = plotnummer, date = parse_date(opnamedatum), radius_m = straal, form = as.character(verschijningsvorm))],
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
perm <- merge(perm, plots[, .(t1 = cycle, plot, date1 = date, r1 = radius_m, area1 = area_m2, form1 = form, inspire1 = inspire)], by = c("plot", "t1"))
note("remeasurement pairs listed by the database: %d (MFV -> NBI-6 %d, NBI-6 -> NBI-7 %d)", nrow(perm), perm[t0 == "MFV", .N], perm[t0 == "NBI6", .N])
perm <- perm[is.finite(r0) & r0 > 0 & is.finite(r1) & r1 > 0 & is.finite(area0) & is.finite(area1)]
note("  both visits with a measured circle: %d", nrow(perm))
perm <- perm[form0 != "MFV_incomplete"];                     note("  first visit complete (MFV flag): %d", nrow(perm))
perm <- perm[form1 %in% PARAMS$forest_forms];                note("  stand type at the second visit kept (%s): %d", paste(PARAMS$forest_forms, collapse = " "), nrow(perm))
perm[, period_length := as.integer(round(as.numeric(date1 - date0) / 365.25))]
perm <- perm[is.finite(period_length) & period_length >= PARAMS$period_range[1] & period_length <= PARAMS$period_range[2]]
note("  interval %d-%d whole years: %d", PARAMS$period_range[1], PARAMS$period_range[2], nrow(perm))
## the 1 km cell was derived anew in NBI-7: plots near a cell border sit in a neighbouring cell there; the second-visit cell is used
note("  pairs whose 1 km cell differs between the visits (neighbouring cells, second-visit cell used): %d", perm[inspire0 != inspire1, .N])
perm[, siteName := paste(plot, t0, t1, sep = "_")]

## ---- 3. first-visit stand, fates, second-visit stand, recruits ---------------
## first-visit living trees inside the first-visit circle (positions are known for NBI-6/7; MFV distances are mostly
## unrecorded (0) and MFV was the first visit of its design, so its trees are all inside)
t0 <- merge(trees, perm[, .(cycle = t0, plot, siteName, t1, r0, area0, period_length)], by = c("cycle", "plot"))
t0 <- t0[alive == TRUE & is.finite(dbh)]
n_small <- t0[dbh < PARAMS$min_dbh_cm, .N]; t0 <- t0[dbh >= PARAMS$min_dbh_cm]
n_out <- t0[dist_m > r0, .N]; t0 <- t0[dist_m <= r0 | is.na(dist_m)]
note("first-visit living trees on the pairs: %d (%d below 5 cm dropped, %d outside the first-visit circle dropped)", nrow(t0) + n_small + n_out, n_small, n_out)
t1 <- merge(trees, perm[, .(cycle = t1, plot, siteName, r1, area1)], by = c("cycle", "plot"))
## matched within the PAIR (siteName), not by plot alone: a plot measured three times has two second visits
t0 <- merge(t0, t1[, .(siteName, tree, alive1 = alive, dead1 = dead, klasse1 = klasse, status1 = status, dbh1 = dbh, in_t1 = TRUE)], by = c("siteName", "tree"), all.x = TRUE)
t0[is.na(in_t1), in_t1 := FALSE]
t0[, fate := fcase(
  in_t1 & alive1 & is.finite(dbh1),                                    "survivor",
  in_t1 & dead1,                                                       "dead",
  in_t1 & klasse1 == PARAMS$removed_class & status1 %in% PARAMS$status_dead, "dead",
  in_t1 & klasse1 == PARAMS$removed_class & status1 %in% PARAMS$status_cut,  "cut",
  default = "censored")]                                              # no row, status 4 / 13 / none, alive without dbh
fates <- t0[, .N, by = fate][order(-N)]
note("  fates: %s", paste(sprintf("%s %d", fates$fate, fates$N), collapse = ", "))
## second-visit stand: living trees inside the second-visit circle
s1 <- t1[alive == TRUE & is.finite(dbh) & dbh >= PARAMS$min_dbh_cm & (dist_m <= r1 | is.na(dist_m))]
## basal area loss per plot (per-hectare, from the database weights): clearcut / converted stands out
ba0 <- t0[, .(ba0_ha = sum(1e4 / area0 * pi * (dbh / 200)^2)), by = siteName]
ba1 <- s1[, .(ba1_ha = sum(1e4 / area1 * pi * (dbh / 200)^2)), by = siteName]
perm <- merge(merge(perm, ba0, by = "siteName", all.x = TRUE), ba1, by = "siteName", all.x = TRUE)
perm[is.na(ba0_ha), ba0_ha := 0]; perm[is.na(ba1_ha), ba1_ha := 0]
n_notree <- perm[ba0_ha == 0, .N]; perm <- perm[ba0_ha > 0];   note("  sites with a living first-visit tree >= 5 cm: %d (%d dropped)", nrow(perm), n_notree)
n_clear <- perm[ba1_ha < (1 - PARAMS$max_ba_loss) * ba0_ha, .N, by = t0]
perm <- perm[ba1_ha >= (1 - PARAMS$max_ba_loss) * ba0_ha]
note("  basal area loss below %.0f %%: %d (%s clearcut / converted dropped; median basal-area change %+.0f %%)", 100 * PARAMS$max_ba_loss, nrow(perm),
     paste(sprintf("%s -> %d", n_clear$t0, n_clear$N), collapse = ", "), 100 * median(perm$ba1_ha / perm$ba0_ha - 1))
t0 <- t0[siteName %in% perm$siteName]; s1 <- s1[siteName %in% perm$siteName]
## recruits: second-visit living trees without a first-visit row, inside the smaller circle, within the ingrowth range
rec <- merge(s1, perm[, .(siteName, r0, period_length)], by = "siteName")
rec <- merge(rec, t0[, .(siteName, tree, had_t0 = TRUE)], by = c("siteName", "tree"), all.x = TRUE)[is.na(had_t0)]
rec[, r_min := pmin(r0, r1)]
n_rec_out <- rec[is.finite(dist_m) & dist_m > r_min, .N]
rec <- rec[is.na(dist_m) | dist_m <= r_min]
n_rec_big <- rec[dbh > PARAMS$min_dbh_cm + PARAMS$recruit_growth * period_length, .N]
rec <- rec[dbh <= PARAMS$min_dbh_cm + PARAMS$recruit_growth * period_length]
note("  new trees at the second visit: %d kept as recruits (median dbh %.1f cm); %d in the ring beyond the smaller circle and %d too large for ingrowth dropped",
     nrow(rec), median(rec$dbh), n_rec_out, n_rec_big)

## ---- 4. species: NBI codes -> names, 95 % set ------------------------------------
sp7 <- rd("lookup_NBI7_Boomsoort")[, .(code, dutch = boomsoort, sci = `Wetenschappelijke naam`)]
sp7[, sci := trimws(sci)]
## codes that are groups ("Other ...", shrubs, open phase) can only be "other"; two varieties of Pinus nigra are one species
sp7[grepl("^[Oo]ther|^Onbezocht", sci) | sci == "" | code %in% c("KV", "XX", "ST"), sci := NA_character_]
sp7[code %in% c("CD", "OD"), sci := "Pinus nigra"]
## "Quercus robur, Quercus petraea" -> "Quercus robur/petraea"; "Populus alba, Populus canescens" -> "Populus alba/canescens"
sp7[grepl(",", sci), sci := vapply(strsplit(sci, ",\\s*"), function(w) {
  g <- sub(" .*$", "", w); e <- sub("^\\S+ ", "", w)
  if (length(unique(g)) == 1L) paste(g[1], paste(e, collapse = "/")) else paste(w, collapse = " / ") }, character(1))]
sp7[, sci := sub("\\s+(var\\.|subsp\\.)\\s.*$", "", sci)]
tally <- t0[, .(stems_t0 = .N, sites_t0 = uniqueN(siteName)), by = code][order(-stems_t0)]
tally <- merge(tally, sp7, by = "code", all.x = TRUE)[order(-stems_t0)]
tally[, cum_share := cumsum(stems_t0) / sum(stems_t0)]
n_sel <- which(tally$cum_share >= PARAMS$species_quantile)[1]
tally[, selected := seq_len(.N) <= n_sel & !is.na(sci)]
species_map <- tally[, .(code, species_name = fifelse(selected, sci, "other"))]
note("species: %d codes among first-visit living trees; %d named species cover %.1f %% of stems, %d codes lumped into 'other'",
     nrow(tally), tally[selected == TRUE, uniqueN(sci)], 100 * tally[selected == TRUE, sum(stems_t0)] / sum(tally$stems_t0), tally[selected == FALSE, .N])
species_dt <- data.table(species_name = c(sort(unique(species_map[species_name != "other", species_name])), "other"))
species_dt[, species := .I]
to_species <- function(dt) { dt <- merge(dt, species_map, by = "code", all.x = TRUE); dt[is.na(species_name), species_name := "other"]
                             dt[, species := species_dt$species[match(species_name, species_dt$species_name)]]; dt }
t0 <- to_species(t0); s1 <- to_species(s1); rec <- to_species(rec)

## ---- 5. virtual trees on the 0.1257 ha patch -------------------------------------
largest_remainder <- function(x) {      # integers that preserve round(sum(x)) within the group
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  fl <- floor(x); k <- as.integer(round(sum(x)) - sum(fl))
  if (k > 0L) { ord <- order(x - fl, decreasing = TRUE); fl[ord[seq_len(min(k, length(x)))]] <- fl[ord[seq_len(min(k, length(x)))]] + 1 }
  as.integer(fl)
}
t0[, ef := 1e4 / area0]; s1[, ef := 1e4 / area1]
t0[, n_virtual := largest_remainder(ef * PS), by = .(siteName, species)]
s1[, n_virtual := largest_remainder(ef * PS), by = .(siteName, species)]
ba_tree <- function(d) pi * (d / 200)^2
note("expansion: first-visit basal area %.1f m2/ha design-weighted vs %.1f m2/ha from the virtual trees; %d of %d first-visit trees and %d of %d second-visit trees expand to zero",
     t0[, sum(ef * ba_tree(dbh))] / uniqueN(t0$siteName), t0[, sum(n_virtual * ba_tree(dbh))] / PS / uniqueN(t0$siteName),
     t0[n_virtual == 0L, .N], nrow(t0), s1[n_virtual == 0L, .N], nrow(s1))

## ---- 6. environment: WorldClim at the centre of the 1 km INSPIRE cell ---------------
## "1kmN3154E4036" = the ETRS89-LAEA (EPSG:3035) cell with lower-left corner N 3,154,000 / E 4,036,000 m
sites <- perm[, .(siteName, plot, t0, t1, date0, date1, period_length, r0, r1, area0_m2 = area0, area1_m2 = area1, form1, inspire = inspire1, ba0_ha, ba1_ha)]
sites[, `:=`(x_laea = as.numeric(sub("^1kmN[0-9]+E([0-9]+)$", "\\1", inspire)) * 1000 + 500,
             y_laea = as.numeric(sub("^1kmN([0-9]+)E[0-9]+$", "\\1", inspire)) * 1000 + 500)]
stopifnot(all(is.finite(sites$x_laea)), all(is.finite(sites$y_laea)))
cells <- unique(sites[, .(inspire, x_laea, y_laea)])
v <- project(vect(cells, geom = c("x_laea", "y_laea"), crs = "EPSG:3035"), "EPSG:4326")
cells[, c("x", "y") := as.data.table(crds(v))]
for (nm in names(PARAMS$env_rasters)) {
  r <- rast(file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]]))
  cells[, (nm) := extract(r, v, ID = FALSE)[[1]]]
}
env_vars <- c(names(PARAMS$env_rasters), "x", "y")
cells_ok <- cells[complete.cases(cells)]
note("env: WorldClim extracted for %d of %d cells (%d sites in %d cells without a land pixel dropped)",
     nrow(cells_ok), nrow(cells), sites[!inspire %in% cells_ok$inspire, .N], nrow(cells) - nrow(cells_ok))
sites <- merge(sites[, -c("x_laea", "y_laea")], cells_ok, by = "inspire")
setorder(sites, plot, t0)
sites[, siteID := .I]
t0 <- t0[siteName %in% sites$siteName]; s1 <- s1[siteName %in% sites$siteName]; rec <- rec[siteName %in% sites$siteName]
id_of <- function(dt) { dt[, siteID := sites$siteID[match(siteName, sites$siteName)]]; dt }
t0 <- id_of(t0); s1 <- id_of(s1); rec <- id_of(rec)

## ---- 7. obs_dt: the second-visit observation per site x species ------------------
P_of <- sites[, .(siteID, period_length)]
state <- s1[, .(ba = sum(n_virtual * ba_tree(dbh)), trees = sum(n_virtual), dbh = weighted.mean(dbh, ef)), by = .(siteID, species)]
rates <- t0[, {
  at_risk <- fate %in% c("survivor", "dead", "cut"); died <- fate %in% c("dead", "cut"); surv <- fate == "survivor"
  ok <- surv & is.finite(dbh1) & (dbh1 - dbh) >= PARAMS$dbh_change_cm[1] & (dbh1 - dbh) <= PARAMS$dbh_change_cm[2] * period_length
  g <- fifelse(ok, dbh1 / dbh - 1, NA_real_)
  .(n_at_risk = sum(at_risk), n_died_t = sum(died), growth_n = sum(is.finite(g)),
    growth_rel = if (any(is.finite(g))) sum(ef[is.finite(g)] * g[is.finite(g)]) / sum(ef[is.finite(g)]) else NA_real_)
}, by = .(siteID, species)]
recr <- merge(rec, sites[, .(siteID, area1_m2)], by = "siteID")      # rec already carries r0 and r1
recr[, forest_share := pmin(1, area1_m2 / (pi * r1^2))]                   # the circle's forest part, applied to the smaller circle
recr <- recr[, .(reg = .N / (pi * min(r0[1], r1[1])^2 / 1e4 * forest_share[1])), by = .(siteID, species)]   # [1]: one value per group
grid <- CJ(siteID = sites$siteID, species = species_dt$species)
obs_dt <- Reduce(function(a, b) merge(a, b, by = c("siteID", "species"), all.x = TRUE), list(grid, state, rates, recr))
stopifnot(!anyDuplicated(obs_dt[, .(siteID, species)]), nrow(obs_dt) == nrow(grid))
obs_dt <- merge(obs_dt, P_of, by = "siteID")
obs_dt[is.na(ba), `:=`(ba = 0, trees = 0L)]
obs_dt[is.na(reg), reg := 0]
obs_dt[is.na(n_at_risk), `:=`(n_at_risk = 0L, n_died_t = 0L, growth_n = 0L)]
## annualise: FINN compares the mean of its annual rates over the interval
obs_dt[, mort := fifelse(n_at_risk > 0, 1 - (1 - n_died_t / n_at_risk)^(1 / period_length), NA_real_)]
obs_dt[, n_died := as.numeric(n_died_t)]
obs_dt[, growth := fifelse(growth_n > 0 & is.finite(growth_rel) & growth_rel > -1, (1 + growth_rel)^(1 / period_length) - 1, NA_real_)]
obs_dt[growth_n == 0L, growth := NA_real_]
obs_dt[, year := period_length]
obs_dt <- merge(obs_dt, species_dt, by = "species")
obs_cols <- c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg", "species_name", "period_length")
obs_dt <- obs_dt[order(siteID, species), ..obs_cols]
note("obs_dt: %d sites x %d species; intervals %d-%d yr (Tmax %d); stem loss over the interval %.3f of the trees at risk (annual rate, pooled: %.4f incl. felling); recruitment mean %.0f, median %.0f stems/ha per interval",
     uniqueN(obs_dt$siteID), nrow(species_dt), min(obs_dt$period_length), max(obs_dt$period_length), max(obs_dt$period_length),
     obs_dt[, sum(n_died) / sum(n_at_risk)], obs_dt[n_at_risk > 0, weighted.mean(mort, n_at_risk)],
     obs_dt[, sum(reg), by = siteID][, mean(V1)], obs_dt[, sum(reg), by = siteID][, median(V1)])

## ---- 8. env rows years 1..Tmax, initial cohorts ------------------------------------
Tmax <- max(sites$period_length)
env_dt <- sites[, .(year = seq_len(Tmax)), by = c("siteID", env_vars)][order(siteID, year), c("siteID", "year", env_vars), with = FALSE]
## one cohort per TALLIED first-visit tree with trees = its virtual copies (same stand as one cohort per copy, cheaper)
init_cohort_dt <- t0[n_virtual > 0L, .(siteID, patchID = 1L, year = 0L, treeName = paste(plot, tree, sep = "_"), species, species_name, dbh, trees = n_virtual)]
setorder(init_cohort_dt, siteID, species, treeName)
stopifnot(setequal(unique(init_cohort_dt$siteID), sites$siteID), setequal(unique(env_dt$siteID), sites$siteID), setequal(unique(obs_dt$siteID), sites$siteID))

## ---- 9. site table and the spatially blocked split ----------------------------------
bm <- PARAMS$block_km * 1000
sites[, block := sprintf("E%03d_N%03d", floor(x_laea / bm), floor(y_laea / bm))]
blocks <- sites[, .(n_sites = .N), by = block][order(block)]
set.seed(PARAMS$block_seed)
blocks[, role := "train"]
blocks[sample(.N, round(PARAMS$holdout_block_share * .N)), role := "holdout"]
sites <- merge(sites, blocks[, .(block, set = role)], by = "block")
set.seed(PARAMS$pilot_seed)
pilot_ids <- sites[set == "train", sample(siteID, PARAMS$pilot_n)]
sites[, pilot := siteID %in% pilot_ids]
setorder(sites, siteID)
note("blocks: %d of %d %d-km blocks held out; %d training sites, %d holdout sites",
     blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_km, sites[set == "train", .N], sites[set == "holdout", .N])
tr_xy <- as.matrix(sites[set == "train", .(x_laea, y_laea)]); ho_xy <- as.matrix(sites[set == "holdout", .(x_laea, y_laea)])
nearest <- apply(outer(ho_xy[, 1], tr_xy[, 1], "-")^2 + outer(ho_xy[, 2], tr_xy[, 2], "-")^2, 1, function(d2) sqrt(min(d2)))
sep <- data.table(stat = c("median_km", "min_km", "q10_km", "q90_km"),
                  value = c(median(nearest), min(nearest), quantile(nearest, 0.1), quantile(nearest, 0.9)) / 1000)
note("holdout separation: nearest training site median %.1f km, minimum %.1f km, 10 %% below %.1f km", sep$value[1], sep$value[2], sep$value[3])
fwrite(sep, file.path(OUT, "holdout_separation.csv"))
site_cols <- c("siteID", "siteName", "plot", "t0", "t1", "date0", "date1", "period_length", "r0", "r1", "area0_m2", "area1_m2", "form1",
               "inspire", "x_laea", "y_laea", env_vars, "ba0_ha", "ba1_ha", "block", "set", "pilot")
fwrite(sites[, ..site_cols], file.path(OUT, "site_table.csv"))
species_out <- merge(species_dt, tally[selected == TRUE, .(species_name = sci, nbi_code = paste(code, collapse = "+"), dutch_name = paste(unique(dutch), collapse = " / "),
                                                            stems_t0 = sum(stems_t0), sites_t0 = sum(sites_t0)), by = .(sci)][, sci := NULL],
                     by = "species_name", all.x = TRUE)[order(species)]
species_out[species_name == "other", `:=`(nbi_code = paste(tally[selected == FALSE, code], collapse = "+"), dutch_name = "all other species and shrubs",
                                          stems_t0 = tally[selected == FALSE, sum(stems_t0)], sites_t0 = NA_integer_)]
fwrite(species_out, file.path(OUT, "species.csv"))

## ---- 10. write the sets (siteIDs renumbered 1..n within each set) --------------------
write_set <- function(ids, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ids <- sort(ids); map <- data.table(siteID_national = ids, siteID_new = seq_along(ids))
  renum <- function(dt) { dt <- merge(dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national"); dt[, siteID := siteID_new][, siteID_new := NULL]; dt }
  o <- renum(obs_dt); e <- renum(env_dt); i <- renum(init_cohort_dt)
  s <- merge(sites[siteID %in% ids, .(siteID, siteName, plot, t0, t1, period_length, inspire, x, y, tavg_c, prec_mm, block, set)], map, by.x = "siteID", by.y = "siteID_national")
  s[, siteID_national := siteID][, siteID := siteID_new][, siteID_new := NULL]
  fwrite(o[order(siteID, species), ..obs_cols], file.path(d, "obs_dt.csv"))
  fwrite(e[order(siteID, year), c("siteID", "year", env_vars), with = FALSE], file.path(d, "env_dt.csv"))
  fwrite(i[order(siteID, species, treeName), .(siteID, patchID, year, treeName, species, species_name, dbh, trees)], file.path(d, "init_cohort_dt.csv"))
  fwrite(species_dt[, .(species, species_name)], file.path(d, "species_dt.csv"))
  fwrite(s[order(siteID), .(siteID, siteID_national, siteName, plot, t0, t1, period_length, inspire, x, y, tavg_c, prec_mm, block, set)], file.path(d, "siteID_dt.csv"))
  note("set %-16s %5d sites, %6d obs rows, %6d cohorts", name, length(ids), nrow(o), nrow(i))
}
write_set(sites$siteID, "national")
write_set(sites[set == "train", siteID], "train")
write_set(sites[set == "holdout", siteID], "holdout_blocked")
write_set(sites[pilot == TRUE, siteID], "pilot")

## ---- 11. PREPARATION.md ----------------------------------------------------------------
site_tot <- obs_dt[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
md <- c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by calibration/03_prepare_pairs_v1.R (data.table %s, terra %s).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("data.table"), packageVersion("terra")), "",
  "Source: Database Zevende Nederlandse Bosinventarisatie (NBI-7), Stichting Probos / Wageningen Environmental Research,",
  "archive as served 2026-11-09 (checksums in raw/SHA256SUMS.txt and calibration/PROVENANCE.md); MFV 2001-2005, NBI-6 2012-2013,",
  "NBI-7 2017-2021 tables. Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the centre of the plot's 1 km INSPIRE cell.", "",
  "## What the tables mean", "",
  "- One site = one remeasurement PAIR of a permanent NBI plot (MFV -> NBI-6 or NBI-6 -> NBI-7; a plot measured three times is two",
  "  sites), one patch of 0.1257 ha (the 20 m circle). Every tallied tree stands for 10,000 / plot area stems per ha and is expanded",
  "  to round(that x 0.1257) virtual trees (largest remainder within site x species). `ba` (m2), `trees` and `dbh` (cm) in obs_dt and",
  "  `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.1257 for per-hectare values.",
  "- `year` 0 = first visit, `year` = `period_length` = second visit (interval rounded to whole years, 4-12). Env rows years 1..Tmax.",
  "- `growth` = annualised relative dbh increment of survivors, weighted by the plot's expansion; `growth_n` = tallied survivors behind it.",
  "- `n_at_risk` = tallied first-visit trees that survived, died or were felled (censored: not found, now < 5 cm, outside forest);",
  "  `n_died` = dead + FELLED trees (stem loss; harvest is part of the mortality target, see the script header); `mort` = annual rate.",
  "- `reg` = new trees >= 5 cm (and <= 5 + 1.5 cm/yr x interval) per hectare over the interval, counted inside the smaller of the two circles.",
  "- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), temp_season (bio4), prec_season (bio15, CV %),",
  "  tmax_warm_c (bio5), tmin_cold_c (bio6), elev_m, x = longitude, y = latitude (WGS84) of the 1 km cell centre.", "",
  "## Processing log", "", paste0("- ", log_lines), "",
  "## Summary (national)", "",
  sprintf("- %d sites (%d MFV -> NBI-6, %d NBI-6 -> NBI-7); mean basal area at the second visit %.1f m2/ha, mean stems %.0f/ha,",
          nrow(sites), sites[t0 == "MFV", .N], sites[t0 == "NBI6", .N], mean(site_tot$ba) / PS, mean(site_tot$trees) / PS),
  sprintf("  mean annual relative growth %.4f, pooled annual stem loss %.4f (incl. felling), mean recruitment %.0f stems/ha per interval",
          obs_dt[growth_n > 0, weighted.mean(growth, growth_n)], obs_dt[n_at_risk > 0, weighted.mean(mort, n_at_risk)], mean(site_tot$reg)), "",
  "## Sets", "",
  "- national: every site passing the filters", "- train: every site in the training blocks",
  "- holdout_blocked: every site in the held-out blocks (never for fitting)",
  sprintf("- pilot: %d training sites drawn with seed %d (smoke tests only)", PARAMS$pilot_n, PARAMS$pilot_seed), "",
  sprintf("Spatial blocking: %d km squares on EPSG:3035, %d of %d blocks held out (seed %d). Holdout site to nearest training site:",
          PARAMS$block_km, blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_seed),
  sprintf("median %.1f km, minimum %.1f km, 10 %% below %.1f km, 90 %% below %.1f km.", sep$value[1], sep$value[2], sep$value[3], sep$value[4]), "",
  "Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set, siteID_national,",
  "plot, cycles, period_length, INSPIRE cell, lon/lat, tavg_c, prec_mm, block, set). Site attributes: site_table.csv; species list",
  "with NBI codes and Dutch names: species.csv.", "",
  "## Species", "", "| species | name | NBI code | Dutch | stems at the first visit | sites |", "|---|---|---|---|---|---|",
  species_out[, sprintf("| %d | %s | %s | %s | %d | %s |", species, species_name, nbi_code, dutch_name, stems_t0, ifelse(is.na(sites_t0), "", sites_t0))], "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(OUT, "PREPARATION.md"))
writeLines(md, file.path("calibration", "PREPARATION.md"))          # the committed copy

## ---- 12. committed records: expected_preparation.csv, climate provenance ----------------
per_set <- rbindlist(lapply(c("national", "train", "holdout_blocked", "pilot"), function(s) {
  o <- fread(file.path(OUT, s, "obs_dt.csv")); i <- fread(file.path(OUT, s, "init_cohort_dt.csv"))
  st <- o[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
  data.table(set = s, sites = uniqueN(o$siteID), obs_rows = nrow(o), cohorts = nrow(i), virtual_trees_t0 = sum(i$trees),
             ba_ha_t1_mean = round(mean(st$ba) / PS, 3), trees_ha_t1_mean = round(mean(st$trees) / PS, 2),
             growth_annual_mean = round(weighted.mean(o$growth, o$growth_n, na.rm = TRUE), 6),
             n_at_risk = sum(o$n_at_risk), n_died = round(sum(o$n_died), 2), mort_annual_pooled = round(o[n_at_risk > 0, weighted.mean(mort, n_at_risk)], 6),
             reg_ha_interval_mean = round(mean(st$reg), 4))
}))
per_species <- obs_dt[, .(set = "national", sites_present_t1 = sum(trees > 0), ba_sum = round(sum(ba), 3), n_at_risk = sum(n_at_risk),
                          n_died = round(sum(n_died), 2), reg_sum = round(sum(reg), 2)), by = .(species, species_name)]
expected <- rbind(melt(per_set, id.vars = "set", variable.name = "quantity", value.name = "value")[, `:=`(species = NA_integer_, species_name = NA_character_)],
                  melt(per_species, id.vars = c("set", "species", "species_name"), variable.name = "quantity", value.name = "value"), use.names = TRUE)
fwrite(expected[, .(set, species, species_name, quantity, value)], file.path("calibration", "expected_preparation.csv"))
prov_file <- file.path("calibration", "PROVENANCE.md")
prov <- readLines(prov_file)
if (any(grepl("^## Climate", prov))) prov <- prov[seq_len(which(grepl("^## Climate", prov))[1] - 1L)]
wc <- rbindlist(lapply(names(PARAMS$env_rasters), function(nm) {
  f <- file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]])
  data.table(variable = nm, file = basename(f), sha256 = unname(tools::sha256sum(f)), size = file.size(f))
}))
prov <- c(prov, "## Climate: WorldClim 2.1 (03_prepare_pairs_v1.R)", "",
          "Source: https://www.worldclim.org/data/worldclim21.html (Fick & Hijmans 2017, doi:10.1002/joc.5086), 1970-2000 normals,",
          "30 arc-seconds, files wc2.1_30s_bio.zip and wc2.1_30s_elev.zip; licence: WorldClim's own terms, free for academic and other non-commercial use, no redistribution without permission (the rasters are not redistributed, only values extracted at the plot locations are included); values extracted at the centres of",
          "the plots' 1 km INSPIRE cells (EPSG:3035) on the date of PREPARATION.md.", "",
          "| variable | file | SHA-256 | size (bytes) |", "|---|---|---|---|",
          wc[, sprintf("| %s | %s | %s | %d |", variable, file, sha256, size)], "")
writeLines(prov, prov_file)
message("done: ", OUT, "; calibration/PREPARATION.md, expected_preparation.csv, PROVENANCE.md")
