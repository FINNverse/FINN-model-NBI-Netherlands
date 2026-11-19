#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-NBI-Netherlands/calibration/v020-soil/prepare_pairs_v3_soil.R
##
## WHAT   Prepared-data version "pairs-v3-soil" = pairs-v2 plus two static soil
##        covariates per plot, on the soil scales shared by the zoo
##        inventories (so that the general model can pool the inventories later):
##          soil_wet           0 = no hydromorphy within 1 m .. 6 = water at or
##                             near the surface most of the year
##          soil_nutrient_poor 1 = eutrophic .. 6 = dystrophic
##        Nothing else changes: every pairs-v2 row and column is kept as it is.
## WHY    The pairs-v2 covariates are two WorldClim normals and the harvest flag;
##        nothing describes the site's soil, while Dutch forests are sorted by
##        exactly that (Scots pine and birch on the dry, poor podzols and drift
##        sands, oak and beech on the richer, loamy sands, alder, ash and poplar
##        on wet and clay soils).
## SOURCE The NBI-7 plot table (data_NBI7_plotmetingen.csv; WUR/Probos
##        intersected the TRUE plot positions with the national maps before the
##        public release, so these values are better than anything a map read at
##        the published 1 km INSPIRE cell could give):
##          Bodemtype_2018  unit code of the Bodemkaart van Nederland 1:50 000
##                          (BRO SGM, version 2018)
##          GWT_2006        groundwater-table class (grondwatertrap, Gt) of the
##                          Bodemkaart, version 2006
##          GWT_2018        Gt of the 2018 version (new sub-classes; INFORMATION ONLY)
##        and the MFV plot table (data_MFV_plotmetingen.csv) column BODEM (7 soil
##        classes, 2001-2005). Code meanings checked against Technisch Document 19A
##        (Stiboka / DLO-Staring Centrum 1995, "Handleiding bodemgeografisch
##        onderzoek", Rubriek AII: legenda Bodemkaart 1:50 000, Tabel A-11 Gt,
##        Tabel A-22 associations, sections 10.1-10.15 legend units, 15.3-15.4
##        lime codes A = kalkrijk, B = kalkarm, C = kalkloos).
## soil_wet  from ONE Gt vintage, GWT_2006 (wetness crosswalk: VIII 0, VII 1, VI 2,
##        "V, IV -> 3-4", III 4, II 5, I 6). Why 2006 and not 2018: it covers more
##        of our plots (86 % vs 66 % of the NBI-7 rows), it lies inside the
##        observation window (MFV 2001-05, NBI-6 2012-13, NBI-7 2017-21) while the
##        2018 map postdates most second visits, and its main classes are the
##        ones the crosswalk is written for; the vintages differ because the
##        groundwater changed, so they are not mixed (2018 is kept as information
##        and its agreement reported). Decisions where the wetness crosswalk leaves a choice:
##          - IV = 3 (GHG > 40 cm: groundwater in the subsoil), V = 3.5 (GHG < 40 cm
##            like III, but GLG > 120 cm instead of 80-120: wet in winter only),
##            the "3-4" of the crosswalk;
##          - sub-class b (IIb, IIIb, Vb: GHG 25-40 cm instead of < 25 cm, Tabel
##            A-11) = the main class - 0.5; u (IVu: GHG 40-80) = the main class;
##          - prefix s (sV, sVI, sVII, sVb: temporary perched water over boulder
##            or Tertiary clay, TD 19A 3.5.3) = the main class: soil_wet counts
##            perched and groundwater alike (stagnic and groundwater wetness on one scale).
## soil_nutrient_poor  from the Bodemkaart unit (nutrient crosswalk column "NL
##        Bodemkaart / MFV BODEM"), the MFV class where the unit is missing.
##        Decisions where the nutrient crosswalk leaves a choice (PARAMS$nutr_rules, first match wins):
##          - the lime code decides within clay and sand: calcareous units (suffix A:
##            Mn35A, Zn50A, Zd20A, Sn13A, pZg20A) = 1 like "Kalkrijke klei / Kalkrijk
##            zand"; non-calcareous clay (suffix B or C; old river clay KR, boulder
##            clay KX, Tertiary clay KT, which carry no lime code) = 2 like "Kalkloze
##            klei"; M/R units without a lime code (MOb75, young salt-marsh clay) = 1;
##          - laarpodzols cHd and cY = 3 like cHn (the same 30-50 cm humic topsoil);
##            the earth soils pZn (gooreerd) and tZd/cZd = 3 like pZg;
##          - peat V and peaty W units = 4 ("Veen"); loess L = 1; limestone K[SKM] = 1;
##          - associations (Tabel A-22): AP petgaten 4, AS stuifzand 6, ABk/ABz
##            beekeerd-/beekdal 3, ABl loess 1, ABv peaty 4, AHc/AHl/AHk/AHz loess
##            and limestone slopes 1, AHt/AHs terrace and flint slopes 3, AMm old
##            Meuse meanders 2, AVo/AVk peat 4;
##          - any other code is NOT guessed: it counts as missing and is filled.
## GAPS   (never silent, never a dropped plot; the source of every value is in
##        soil_plots.csv and <set>/soil_dt.csv):
##        soil_wet:  Gt "-" (not mapped) or plot not in NBI-7 -> mean soil_wet of
##                   the NBI-7 plots (all 3,604, not only ours) with a Gt of the
##                   same Bodemkaart unit when >= min_group hold it, else of the
##                   same unit group (code letters before the digits, e.g. Hn,
##                   Zd, Mn), else of the same MFV BODEM class, else of the 25 km
##                   block (our plots), else the national mean;
##        soil_nutrient_poor: no usable unit -> the MFV BODEM class, else the mean
##                   of the 25 km block, else the national mean.
## SCALING env_dt.csv of pairs-v2 is in natural units and FINN z-scales the climate
##        with env_autoscale "auto"; the soil columns are added the same way, as
##        the unscaled ordinal scores (comparable across inventories).
## INPUT  ../../FINN-zoo-data/nbi/prepared/pairs-v2/ (every set; never modified)
##        ../../FINN-zoo-data/nbi/raw/csv/data_NBI7_plotmetingen.csv, data_MFV_plotmetingen.csv
## OUTPUT ../../FINN-zoo-data/nbi/prepared/pairs-v3-soil/
##          PREPARATION.md, soil_plots.csv (every plot of the national set: codes,
##          scores, sources), soil_species_means.csv, soil_env_correlation.csv,
##          v2_v3_env_check.csv, MD5SUMS.txt (manifest of every file)
##          site_table.csv, species.csv, holdout_separation.csv   byte copies of v2
##          <set>/env_dt.csv     v2 + soil_wet, soil_nutrient_poor (v2 columns identical)
##          <set>/soil_dt.csv    siteID, plot, scores, sources, codes
##          <set>/*              every other file a byte copy of v2
## RUN    Rscript calibration/v020-soil/prepare_pairs_v3_soil.R   (from the repository root, < 1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "nbi")    # zoo layout: data outside git

PARAMS <- list(
  version      = "pairs-v3-soil",
  base_version = "pairs-v2",                # the base; never modified
  raw_nbi7     = "raw/csv/data_NBI7_plotmetingen.csv",   # Bodemtype_2018, GWT_2006, GWT_2018 at the true plot position
  raw_mfv      = "raw/csv/data_MFV_plotmetingen.csv",    # BODEM (7 classes)
  gt_column    = "GWT_2006",                # ONE Gt vintage (header: coverage, inside the observation window)
  gt_info_column = "GWT_2018",              # kept as information, never used for soil_wet
  gt_main_score = c(I = 6, II = 5, III = 4, IV = 3, V = 3.5, VI = 2, VII = 1, VIII = 0),   # wetness crosswalk; IV 3 / V 3.5 = its "3-4"
  gt_b_offset  = -0.5,                      # sub-class b = GHG 25-40 cm instead of < 25 cm: half a class drier
  gt_c_offset_2018 = -0.5,                  # GWT_2018 only (information): Ic / IIc = GHG deeper than the class's a/b, as b
  na_codes     = c("", "-", "#N/A", "NA"),  # "not mapped" in the NBI-7 table
  ## Bodemkaart unit -> soil_nutrient_poor, first matching pattern wins (header; nutrient crosswalk)
  nutr_rules   = data.table(
    pattern = c("^AP$", "^AS$", "^A[Bb]k$", "^ABz$", "^ABl$", "^ABv$", "^AH[clkz]$", "^AH[ts]$", "^AMm$", "^AV[ok]$",
                "^[a-z]*V[a-z]", "^[a-z]*W[a-z]",
                "^c(Hn|Hd|Y)", "^[a-z]*Hn", "^[a-z]*(Hd|Y)",
                "^[a-z]*E[ZKL]",
                "^[a-z]*(Z[a-z]|Sn)[0-9]+A$",
                "^[a-z]*Zb", "^[ptc]Z",
                "^[a-z]*Z[nd]",
                "^[a-z]*L[ndh]", "^[a-z]*K[SKM]", "^[a-z]*K[RXT]", "^[a-z]*MZ",
                "^[a-z]*[MR][A-Za-z]*[0-9]*A$", "^[a-z]*[MR][A-Za-z]*[0-9]*[BC]$", "^[a-z]*[MR]"),
    score   = c(4, 6, 3, 3, 1, 4, 1, 3, 2, 4,
                4, 4,
                3, 5, 6,
                2,
                1,
                3, 3,
                6,
                1, 1, 2, 3,
                1, 2, 1),
    label   = c("AP petgaten (peat cuttings)", "AS stuifzandgronden (drift sand)", "ABk kleiige beekeerdgronden", "ABz zandige beekdalgronden",
                "ABl loessige beekdalgronden", "ABv venige beekdalgronden", "AH loess / limestone slopes", "AH terrace / flint slopes",
                "AMm old Meuse meanders", "AV peat in reclamation",
                "V peat", "W peaty (moerige) soils",
                "cHn/cHd/cY laarpodzol", "Hn veldpodzol", "Hd/Y haar- and holtpodzol",
                "E enkeerd (thick plaggen)",
                "calcareous sand (suffix A)",
                "Zb vorstvaaggrond", "pZg/pZn/tZd/cZd earth soils on sand",
                "Zn/Zd non-calcareous sand",
                "L loess", "K limestone", "KR/KX/KT old clays (no lime code)", "MZ old marine sand",
                "M/R clay, calcareous (A)", "M/R clay, low-lime / lime-free (B/C)", "M/R clay, no lime code")),
  mfv_bodem_score = c("Kalkrijke klei" = 1, "Kalkrijk zand" = 1, "Loess" = 1, "Kalkloze klei" = 2,
                      "Rijk zand" = 3, "Veen" = 4, "Arm zand" = 5),   # nutrient crosswalk column "MFV BODEM"
  min_group    = 10L,                       # plots a unit / group / class / block needs to fill a gap with its mean
  soil_vars    = c("soil_wet", "soil_nutrient_poor"),   # the two env columns the fit may use
  env_vars_check = c("tavg_c", "prec_mm", "temp_season", "prec_season", "tmax_warm_c", "tmin_cold_c", "elev_m", "harvest"),   # existing covariates the soil scores are correlated with
  changed_files = c("env_dt.csv")           # the only pairs-v2 file per set that changes (soil_dt.csv is new)
)
BASE <- file.path(DATA_DIR, "prepared", PARAMS$base_version)
OUT  <- file.path(DATA_DIR, "prepared", PARAMS$version)
stopifnot(dir.exists(BASE), !dir.exists(OUT))   # a version is never modified: delete it to rebuild
sets <- list.dirs(BASE, full.names = FALSE, recursive = FALSE)
message("sets: ", paste(sets, collapse = ", "))
na_to <- function(x) { x <- trimws(as.character(x)); x[x %in% PARAMS$na_codes] <- NA_character_; x }

## ---- 1. soil scores per plot from the raw tables ---------------------------------------
p7 <- fread(file.path(DATA_DIR, PARAMS$raw_nbi7), na.strings = NULL, colClasses = "character")[
  , .(plot = as.integer(plotnummer), unit = na_to(Bodemtype_2018), gt = na_to(get(PARAMS$gt_column)), gt_info = na_to(get(PARAMS$gt_info_column)))]
mfv <- fread(file.path(DATA_DIR, PARAMS$raw_mfv), na.strings = NULL, colClasses = "character")[
  , .(plot = as.integer(PLOTNUMMER), mfv_bodem = na_to(BODEM))]
stopifnot(!anyDuplicated(p7$plot), !anyDuplicated(mfv$plot), !anyNA(p7$plot), !anyNA(mfv$plot))

## Gt -> soil_wet: [qualitative prefix][Roman main class][sub-class letters]
gt_score <- function(code, c_offset = 0) {
  code <- sub("^Iv", "IV", code)                         # GWT_2018 writes "Ivc", "Ivu"
  pre  <- sub("^([a-z]*)([IV]+)([a-z]*)$", "\\1", code)
  main <- sub("^([a-z]*)([IV]+)([a-z]*)$", "\\2", code)
  sub_ <- sub("^([a-z]*)([IV]+)([a-z]*)$", "\\3", code)
  ok   <- !is.na(code) & grepl("^[a-z]*[IV]+[a-z]*$", code) & main %in% names(PARAMS$gt_main_score)
  s <- rep(NA_real_, length(code))
  s[ok] <- PARAMS$gt_main_score[main[ok]] + PARAMS$gt_b_offset * grepl("b", sub_[ok]) + c_offset * grepl("c", sub_[ok])
  s
}
p7[, wet_gt := gt_score(gt)]
p7[, wet_gt_info := gt_score(gt_info, PARAMS$gt_c_offset_2018)]
stopifnot(p7[!is.na(gt), all(is.finite(wet_gt))], p7[!is.na(gt_info), all(is.finite(wet_gt_info))])   # every mapped Gt code is understood
GT_TAB <- p7[, .N, by = .(gt, wet_gt)][order(-wet_gt)]

## Bodemkaart unit -> soil_nutrient_poor
nutr_of_unit <- function(u) {
  s <- rep(NA_real_, length(u)); lab <- rep(NA_character_, length(u))
  for (i in seq_len(nrow(PARAMS$nutr_rules))) {
    hit <- is.na(s) & !is.na(u) & grepl(PARAMS$nutr_rules$pattern[i], u)
    s[hit] <- PARAMS$nutr_rules$score[i]; lab[hit] <- PARAMS$nutr_rules$label[i]
  }
  list(s, lab)
}
p7[, c("nutr_unit", "nutr_rule") := nutr_of_unit(unit)]
p7[, unit_group := sub("^([A-Za-z]+).*$", "\\1", unit)]
UNMATCHED <- p7[!is.na(unit) & is.na(nutr_unit), .N, by = unit][order(-N)]
message("Bodemkaart units without a nutrient rule (all NBI-7 plots; filled like missing units): ",
        if (nrow(UNMATCHED)) paste(sprintf("%s (%d)", UNMATCHED$unit, UNMATCHED$N), collapse = ", ") else "none")
mfv[, nutr_mfv := unname(PARAMS$mfv_bodem_score[mfv_bodem])]
stopifnot(mfv[!is.na(mfv_bodem), all(is.finite(nutr_mfv))])   # every MFV class is in the crosswalk

## ---- 2. every plot of the national set, gaps filled ------------------------------------------
site_tab <- fread(file.path(BASE, "site_table.csv"))
nat <- fread(file.path(BASE, "national", "siteID_dt.csv"))
stopifnot(setequal(nat$siteName, site_tab$siteName))
## the block is used for the gap fill only; two plots carry a different INSPIRE cell in the MFV and
## NBI-6 plot definitions (1 km apart, across a block edge): they take the block of their first site
PL <- site_tab[, .(block = block[1], n_blocks = uniqueN(block)), by = plot]
message(sprintf("plots with two blocks (first site's block used for the gap fill): %s", paste(PL[n_blocks > 1, plot], collapse = ", ")))
stopifnot(!anyDuplicated(PL$plot), all(PL$n_blocks <= 2L))
PL <- merge(PL, p7, by = "plot", all.x = TRUE)
PL <- merge(PL, mfv, by = "plot", all.x = TRUE)
PL[, in_nbi7 := plot %in% p7$plot]
## reference tables for the fills: all raw plots, not only ours (the map relation is not forest-specific)
REF <- merge(p7, mfv, by = "plot", all = TRUE)
group_mean <- function(dt, value, by) dt[is.finite(get(value)) & !is.na(get(by)), .(m = mean(get(value)), n = .N), by = by][n >= PARAMS$min_group]

## soil_wet
PL[, soil_wet := wet_gt][, soil_wet_source := fifelse(is.finite(soil_wet), paste0("Gt ", PARAMS$gt_column), NA_character_)]
fill_levels_wet <- list(list(by = "unit",       ref = REF, label = "mean of the Bodemkaart unit"),
                        list(by = "unit_group", ref = REF, label = "mean of the Bodemkaart unit group"),
                        list(by = "mfv_bodem",  ref = REF, label = "mean of the MFV BODEM class"),
                        list(by = "block",      ref = PL[soil_wet_source == paste0("Gt ", PARAMS$gt_column)], label = "mean of the 25 km block"))
for (f in fill_levels_wet) {
  g <- group_mean(f$ref, if (f$by == "block") "soil_wet" else "wet_gt", f$by)
  PL[g, on = f$by, fill := i.m]
  PL[!is.finite(soil_wet) & is.finite(fill), `:=`(soil_wet = fill, soil_wet_source = f$label)]
  PL[, fill := NULL]
}
PL[!is.finite(soil_wet), `:=`(soil_wet = PL[soil_wet_source == paste0("Gt ", PARAMS$gt_column), mean(soil_wet)], soil_wet_source = "national mean")]

## soil_nutrient_poor
PL[, soil_nutrient_poor := nutr_unit][, soil_nutrient_poor_source := fifelse(is.finite(soil_nutrient_poor), "Bodemkaart unit", NA_character_)]
PL[!is.finite(soil_nutrient_poor) & is.finite(nutr_mfv), `:=`(soil_nutrient_poor = nutr_mfv, soil_nutrient_poor_source = "MFV BODEM class")]
g <- group_mean(PL[soil_nutrient_poor_source %in% c("Bodemkaart unit", "MFV BODEM class")], "soil_nutrient_poor", "block")
PL[g, on = "block", fill := i.m]
PL[!is.finite(soil_nutrient_poor) & is.finite(fill), `:=`(soil_nutrient_poor = fill, soil_nutrient_poor_source = "mean of the 25 km block")]
PL[, fill := NULL]
PL[!is.finite(soil_nutrient_poor), `:=`(soil_nutrient_poor = PL[soil_nutrient_poor_source %in% c("Bodemkaart unit", "MFV BODEM class"), mean(soil_nutrient_poor)],
                                        soil_nutrient_poor_source = "national mean")]
stopifnot(PL[, all(is.finite(soil_wet) & soil_wet >= 0 & soil_wet <= 6)],
          PL[, all(is.finite(soil_nutrient_poor) & soil_nutrient_poor >= 1 & soil_nutrient_poor <= 6)])

## coverage, per plot and per site of the national set
site_plot <- nat[, .(siteID, siteName, plot)]
stopifnot(!anyDuplicated(site_plot$siteID), all(site_plot$plot %in% PL$plot))
cover <- function(dt, unit_label) rbindlist(lapply(PARAMS$soil_vars, function(v)
  dt[, .N, by = .(source = get(paste0(v, "_source")))][, `:=`(variable = v, unit = unit_label, share = N / sum(N))]))
COVER <- rbind(cover(PL, "plots"), cover(merge(site_plot, PL, by = "plot"), "sites (national)"))
print(COVER)

## agreement checks
AGREE_MFV <- PL[!is.na(nutr_unit) & !is.na(nutr_mfv), .(n = .N, unit_score_mean = mean(nutr_unit), mfv_score = nutr_mfv[1]), by = mfv_bodem][order(mfv_score)]
AGREE_GT <- PL[is.finite(wet_gt) & is.finite(wet_gt_info), .(n = .N, r = cor(wet_gt, wet_gt_info), mean_2006 = mean(wet_gt), mean_2018 = mean(wet_gt_info),
                                                           share_equal = mean(wet_gt == wet_gt_info))]
r_mfv <- PL[!is.na(nutr_unit) & !is.na(nutr_mfv), cor(nutr_unit, nutr_mfv)]
print(AGREE_MFV); print(AGREE_GT)

## ---- 3. write the version: every set = v2 + soil ------------------------------------------------
dir.create(OUT, recursive = TRUE)
top_files <- setdiff(list.files(BASE), c(sets, "PREPARATION.md"))
invisible(file.copy(file.path(BASE, top_files), OUT, copy.date = TRUE))
fwrite(PL[order(plot), .(plot, block, in_nbi7, bodemtype_2018 = unit, unit_group, gt_2006 = gt, gt_2018 = gt_info, mfv_bodem,
                         soil_wet, soil_wet_source, soil_nutrient_poor, soil_nutrient_poor_source,
                         wet_gt_2006 = wet_gt, wet_gt_2018_info = wet_gt_info, nutr_unit, nutr_rule, nutr_mfv)],
       file.path(OUT, "soil_plots.csv"))
CHECKS <- list()
for (s in sets) {
  dir.create(file.path(OUT, s))
  files <- list.files(file.path(BASE, s))
  invisible(file.copy(file.path(BASE, s, setdiff(files, PARAMS$changed_files)), file.path(OUT, s), copy.date = TRUE))
  sites <- fread(file.path(BASE, s, "siteID_dt.csv"))
  site_soil <- merge(sites[, .(siteID, siteName, plot)],
                     PL[, .(plot, soil_wet, soil_nutrient_poor, soil_wet_source, soil_nutrient_poor_source,
                            bodemtype_2018 = unit, gt_2006 = gt, mfv_bodem)], by = "plot", all.x = TRUE)
  stopifnot(nrow(site_soil) == nrow(sites), !anyDuplicated(site_soil$siteID), !anyNA(site_soil[, PARAMS$soil_vars, with = FALSE]))
  setcolorder(site_soil, c("siteID", "siteName", "plot"))
  fwrite(site_soil[order(siteID)], file.path(OUT, s, "soil_dt.csv"))
  env <- fread(file.path(BASE, s, "env_dt.csv"))
  stopifnot(!any(PARAMS$soil_vars %in% names(env)))
  new <- copy(env)                                       # an update join keeps v2's row order
  new[site_soil, on = "siteID", (PARAMS$soil_vars) := mget(paste0("i.", PARAMS$soil_vars))]
  stopifnot(nrow(new) == nrow(env), !anyNA(new[, PARAMS$soil_vars, with = FALSE]),
            new[, .(k = uniqueN(soil_wet) + uniqueN(soil_nutrient_poor)), by = siteID][, all(k == 2L)])   # static over the years
  fwrite(new, file.path(OUT, s, "env_dt.csv"))
  back <- fread(file.path(OUT, s, "env_dt.csv"))
  same_files <- setdiff(files, PARAMS$changed_files)
  CHECKS[[s]] <- data.table(
    set = s, rows = nrow(back), sites = uniqueN(back$siteID),
    v2_rows_and_columns_identical = isTRUE(all.equal(env, back[, names(env), with = FALSE], check.attributes = FALSE, tolerance = 0)),
    v2_column_order_kept = identical(names(back)[seq_along(names(env))], names(env)),
    other_files_byte_identical = all(unname(tools::md5sum(file.path(BASE, s, same_files))) == unname(tools::md5sum(file.path(OUT, s, same_files)))),
    share_wet_from_gt = mean(site_soil$soil_wet_source == paste0("Gt ", PARAMS$gt_column)),
    share_nutr_from_unit = mean(site_soil$soil_nutrient_poor_source == "Bodemkaart unit"),
    share_nutr_from_unit_or_mfv = mean(site_soil$soil_nutrient_poor_source %in% c("Bodemkaart unit", "MFV BODEM class")))
}
CHECKS <- rbindlist(CHECKS)
print(CHECKS)
stopifnot(all(CHECKS$v2_rows_and_columns_identical), all(CHECKS$v2_column_order_kept), all(CHECKS$other_files_byte_identical),
          all(unname(tools::md5sum(file.path(BASE, top_files))) == unname(tools::md5sum(file.path(OUT, top_files)))))
fwrite(CHECKS, file.path(OUT, "v2_v3_env_check.csv"))

## ---- 4. distributions, correlations, species means (national set) ------------------------------
S <- merge(fread(file.path(OUT, "national", "siteID_dt.csv"))[, .(siteID, set)],
           fread(file.path(OUT, "national", "env_dt.csv"))[year == 1L], by = "siteID")
DIST <- rbindlist(lapply(PARAMS$soil_vars, function(v) S[, .(variable = v, mean = mean(get(v)), sd = sd(get(v)), min = min(get(v)), q10 = quantile(get(v), 0.1),
                                                             q50 = median(get(v)), q90 = quantile(get(v), 0.9), max = max(get(v))), by = set]))
HIST <- rbindlist(lapply(PARAMS$soil_vars, function(v) S[, .N, by = .(value = round(get(v) * 2) / 2)][order(value)][, `:=`(variable = v, share = N / sum(N))]))
XTAB <- S[, .N, by = .(soil_wet = round(soil_wet), soil_nutrient_poor = round(soil_nutrient_poor))]
CORR <- rbindlist(lapply(PARAMS$soil_vars, function(v) rbindlist(lapply(c(setdiff(PARAMS$soil_vars, v), PARAMS$env_vars_check), function(e)
  data.table(soil = v, covariate = e, pearson = cor(S[[v]], S[[e]]), spearman = cor(S[[v]], S[[e]], method = "spearman"))))))
CORR <- CORR[!(soil == "soil_nutrient_poor" & covariate == "soil_wet")]
fwrite(CORR, file.path(OUT, "soil_env_correlation.csv"))
print(DIST); print(CORR)
init <- fread(file.path(OUT, "national", "init_cohort_dt.csv"))
spn  <- fread(file.path(OUT, "national", "species_dt.csv"))
SPM  <- merge(init[, .(siteID, species, ba = trees * pi * (dbh / 200)^2)], S[, c("siteID", PARAMS$soil_vars), with = FALSE], by = "siteID")[
  , .(ba_share = sum(ba), plots = uniqueN(siteID), soil_wet = weighted.mean(soil_wet, ba), soil_nutrient_poor = weighted.mean(soil_nutrient_poor, ba)), by = species]
SPM[, ba_share := ba_share / sum(ba_share)]
SPM <- merge(spn, SPM, by = "species")[order(-soil_wet)]
fwrite(SPM, file.path(OUT, "soil_species_means.csv"))
print(SPM)

## ---- 5. PREPARATION.md and the manifest ------------------------------------------------------------
md_table <- function(x) {
  x <- copy(as.data.table(x))
  for (cc in names(x)) if (is.integer(x[[cc]])) x[, (cc) := as.character(get(cc))] else if (is.numeric(x[[cc]])) x[, (cc) := formatC(get(cc), digits = 3, format = "g")]
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
commit <- suppressWarnings(tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
if (!length(commit)) commit <- "unknown"
rules_used <- PL[!is.na(nutr_rule), .N, by = .(rule = nutr_rule, score = nutr_unit)][order(score, -N)]
prep <- c(
  sprintf("# Prepared data: %s", PARAMS$version), "",
  sprintf("Written %s by `calibration/v020-soil/prepare_pairs_v3_soil.R` (FINN-model-NBI-Netherlands at commit `%s`; data.table %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), commit, packageVersion("data.table")), "",
  sprintf("**%s = %s plus two static soil covariates per plot. Nothing else changes**: every %s row and column is kept (`v2_v3_env_check.csv`); `env_dt.csv` gains `soil_wet` and `soil_nutrient_poor`; each set gains `soil_dt.csv`; every other file is a byte copy (`MD5SUMS.txt`). The population, the sets, the blocks and the holdout separation are those of %s (see its PREPARATION.md).",
          PARAMS$version, PARAMS$base_version, PARAMS$base_version, PARAMS$base_version), "",
  "**Why.** No pairs-v2 covariate describes the soil (two WorldClim normals and the harvest flag), while Dutch forest types are sorted by water and nutrients.", "",
  "**Scales** (harmonised across the zoo inventories): `soil_wet` 0 = no hydromorphy within 1 m .. 6 = water at or near the surface most of the year; `soil_nutrient_poor` 1 = eutrophic .. 6 = dystrophic. Stored unscaled (ordinal scores); the fit z-scales them (env_autoscale \"auto\").", "",
  "**Source.** NBI-7 plot table (`data_NBI7_plotmetingen.csv`, Probos / WENR): `Bodemtype_2018` = unit of the Bodemkaart van Nederland 1:50 000 and `GWT_2006` / `GWT_2018` = groundwater-table class (Gt), intersected by WUR/Probos at the TRUE plot positions; MFV plot table `BODEM` (7 classes). Code meanings: Technisch Document 19A (DLO-Staring Centrum 1995), Tabel A-11 (Gt), Tabel A-22 (associations), sections 10 and 15 (legend units, lime codes A kalkrijk / B kalkarm / C kalkloos).", "",
  sprintf("**soil_wet** from ONE Gt vintage, `%s` (covers more plots than 2018; inside the observation window 2001-2021; the vintages differ because the groundwater changed, so they are never mixed). The wetness crosswalk with the choices it leaves open made explicit: I 6, II 5, III 4, V 3.5, IV 3, VI 2, VII 1, VIII 0 (\"V, IV -> 3-4\": IV has GHG > 40 cm = groundwater in the subsoil, V has GHG < 40 cm like III but GLG > 120 cm); sub-class b (GHG 25-40 cm) = main class %+.1f; u = main class; prefix s (temporary perched water) = main class, because the scale counts perched and groundwater alike.", PARAMS$gt_column, PARAMS$gt_b_offset), "",
  md_table(GT_TAB[, .(Gt = fifelse(is.na(gt), "(not mapped)", gt), soil_wet = wet_gt, nbi7_plots = N)]), "",
  sprintf("Agreement with the 2018 vintage on plots holding both (information only): n = %d, r = %.2f, identical score on %.0f %%, mean %.2f (2006) vs %.2f (2018).",
          AGREE_GT$n, AGREE_GT$r, 100 * AGREE_GT$share_equal, AGREE_GT$mean_2006, AGREE_GT$mean_2018), "",
  "**soil_nutrient_poor** from the Bodemkaart unit by the rules below (nutrient crosswalk; first match wins; the lime code decides within clay and sand), the MFV `BODEM` class where the unit is missing or not covered by a rule. Rules used on our plots:", "",
  md_table(rules_used), "",
  sprintf("Units without a rule (not guessed; filled like a missing unit), over all NBI-7 plots: %s.",
          if (nrow(UNMATCHED)) paste(sprintf("%s (%d)", UNMATCHED$unit, UNMATCHED$N), collapse = ", ") else "none"), "",
  sprintf("Agreement of the unit score with the MFV class on plots holding both (r = %.2f):", r_mfv), "",
  md_table(AGREE_MFV), "",
  sprintf("**Gaps** (never a dropped plot; the source of every value in `soil_plots.csv` and `<set>/soil_dt.csv`). soil_wet: Gt not mapped (\"-\") or plot not in NBI-7 -> mean soil_wet of the NBI-7 plots (all, not only ours) of the same Bodemkaart unit with >= %d plots, else of the same unit group (letters before the digits), else of the same MFV BODEM class, else of the 25 km block, else the national mean. soil_nutrient_poor: unit missing -> MFV BODEM class -> mean of the 25 km block (>= %d plots) -> national mean.", PARAMS$min_group, PARAMS$min_group), "",
  md_table(dcast(COVER, variable + unit ~ source, value.var = "share", fill = 0)), "",
  "**Distributions** over the sites of the national set (one row per site):", "",
  md_table(DIST), "",
  md_table(dcast(HIST, value ~ variable, value.var = "share", fill = 0)), "",
  "Cross-tabulation (rounded scores, sites of the national set; rows soil_wet, columns soil_nutrient_poor):", "",
  md_table(dcast(XTAB, soil_wet ~ soil_nutrient_poor, value.var = "N", fill = 0)), "",
  "**Correlation with the existing covariates** (sites of the national set; `soil_env_correlation.csv`):", "",
  md_table(CORR), "",
  "**Species** (basal-area-weighted mean score of the first-visit stands, national set; `soil_species_means.csv`):", "",
  md_table(SPM[, .(species_name, ba_share, plots, soil_wet, soil_nutrient_poor)]), "",
  "## Checks", "", md_table(CHECKS), "",
  "## PARAMS", "",
  sprintf("- `%s` = `%s`", setdiff(names(PARAMS), "nutr_rules"), vapply(PARAMS[setdiff(names(PARAMS), "nutr_rules")], function(x)
    if (!is.null(names(x)) && is.numeric(x)) paste(names(x), x, sep = " = ", collapse = ", ") else paste(x, collapse = ", "), character(1))),
  "- `nutr_rules` (pattern -> score):", "",
  md_table(PARAMS$nutr_rules)
)
writeLines(prep, file.path(OUT, "PREPARATION.md"))
all_files <- setdiff(list.files(OUT, recursive = TRUE), "MD5SUMS.txt")
writeLines(sprintf("%s  %s", unname(tools::md5sum(file.path(OUT, all_files))), all_files), file.path(OUT, "MD5SUMS.txt"))
message("done: ", OUT)
