# Prepared data: pairs-v2

Written 2026-11-10 07:01 CET by calibration/v02/04_prepare_pairs_v2.R (data.table 1.18.6.1, terra 1.8.60).

Source: Database Zevende Nederlandse Bosinventarisatie (NBI-7), Stichting Probos / Wageningen Environmental Research,
archive as served 2026-11-09 (checksums in raw/SHA256SUMS.txt and calibration/PROVENANCE.md); MFV 2001-2005, NBI-6 2012-2013,
NBI-7 2017-2021 tables. Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the centre of the plot's 1 km INSPIRE cell.

## The population this version describes (the agreed rare-species rule)

- Species set S (|S| = 14, ranked by national first-visit basal area; genus-level NBI codes kept as the NBI records them): Pinus sylvestris, Quercus robur/petraea, Betula, Pseudotsuga menziesii, Fagus sylvatica, Pinus nigra, Larix kaempferi, Quercus rubra, Picea abies, Fraxinus excelsior, Alnus glutinosa, Populus, Acer, Salix.
- Site rule on the START census only: a plot pair is kept when species outside S (incl. the NBI shrub group) hold < 2 % of the
  first-visit basal area: 1762 of 2271 pairs (77.6 %), i.e. the model is valid for stands where the modelled species hold > 98 % of the
  basal area; with one NBI plot per 100 ha the share of pairs is also the share of forest area these stands cover.
- The end-census out-of-S share is NOT used for the selection; it is recorded per site (`oos_ba_share_t1`: < 1 % on 1301 sites, 1-5 % on 326, > 5 % on 135) for a stratified validation.
- Out-of-S trees on kept sites are dropped from the stand, the fates, the state and the recruits (option (a)).
- No lumped "other" class and no shrubs. The trade-off table behind the choice: calibration/v02/species_subset_tradeoff.csv (|S| 6-24 x 1/2/3/5 %);
  per-species identifiability of S: species_identifiability.csv. The thresholds there (plots, trees at risk, holdout blocks, climate span) are
  coverage choices, not measured identifiability limits.

## What the tables mean

- One site = one remeasurement PAIR of a permanent NBI plot (MFV -> NBI-6 or NBI-6 -> NBI-7; a plot measured three times is two
  sites), one patch of 0.1257 ha (the 20 m circle). Every tallied tree stands for 10,000 / plot area stems per ha and is expanded
  to round(that x 0.1257) virtual trees (largest remainder within site x species). `ba` (m2), `trees` and `dbh` (cm) in obs_dt and
  `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.1257 for per-hectare values.
- `year` 0 = first visit, `year` = `period_length` = second visit (interval rounded to whole years, 4-12). Env rows years 1..Tmax.
- `growth` = annualised relative dbh increment of survivors, weighted by the plot's expansion; `growth_n` = tallied survivors behind it.
- `n_at_risk` = tallied first-visit trees that survived, died or were felled (censored: not found, now < 5 cm, outside forest);
  `n_died` = dead + FELLED trees (stem loss; harvest is part of the mortality target); `mort` = annual rate.
- `reg` = new trees >= 5 cm whose back-modelled start dbh (dbh - species g75 x interval) is below 5 cm, per hectare over the
  interval, counted inside the smaller of the two circles; larger new trees are ingrowth (in the state, not recruits).
- tree_obs_dt: one row per tracked first-visit tree (siteID, year = period_length, patchID, species, dbh at the first visit,
  growth = annualised relative dbh increment or NA, died 0/1 incl. felled, period_length, trees = 1) for FINN's per-tree losses.
- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), temp_season (bio4), prec_season (bio15, CV %),
  tmax_warm_c (bio5), tmin_cold_c (bio6), elev_m, x = longitude, y = latitude (WGS84) of the 1 km cell centre,
  harvest = 1 if harvest activity was recorded at the second visit (NBI-7 traces of harvest; NBI-6 clearcut or single-tree felling), else 0.

## Processing log

- tree rows: MFV 69848, NBI-6 79880, NBI-7 78371; plots with a measured circle: MFV 2963, NBI-6 3191, NBI-7 3173
- harvest recorded at the second visit on 1177 pairs, none on 1463, no record on 9 (set to 0)
- remeasurement pairs listed by the database: 2649 (MFV -> NBI-6 1236, NBI-6 -> NBI-7 1413)
-   both visits with a measured circle: 2648
-   first visit complete (MFV flag): 2648
-   stand type at the second visit kept (1a.1.1 1a.1.2 1a.2 1c.1 1c.2 1c.3 1c.5 3): 2426
-   interval 3-15 whole years: 2426
-   pairs whose 1 km cell differs between the visits (neighbouring cells, second-visit cell used): 181
- first-visit living trees on the pairs: 50837 (0 below 5 cm dropped, 1289 outside the first-visit circle dropped)
-   fates: survivor 36448, cut 6164, dead 5928, censored 1008
-   sites with a living first-visit tree >= 5 cm: 2389 (37 dropped)
-   basal area loss below 50 %: 2273 (NBI6 -> 56, MFV -> 60 clearcut / converted dropped; median basal-area change +10 %)
-   new trees at the second visit: 6578 recruits (median dbh 5.9 cm); 3246 beyond the smaller circle dropped, 3090 with a back-modelled start dbh >= 5 cm (species g75 0.15-0.67 cm/yr, pooled 0.37) = ingrowth, not recruitment
- env: WorldClim extracted for 1433 of 1434 cells (2 pairs in 1 cells without a land pixel dropped)
- blocks over all 2271 pairs: 15 of 76 25-km blocks held out (seed 21); holdout blocks hold 740 pairs
- eligible taxa among first-visit living trees: 37 (shrub group and NBI group codes excluded); top 10 by basal area: Pinus sylvestris, Quercus robur/petraea, Betula, Pseudotsuga menziesii, Fagus sylvatica, Pinus nigra, Larix kaempferi, Quercus rubra, Picea abies, Fraxinus excelsior
- trade-off table written: 76 rows (|S| 6-24 x 1/2/3/5)
- site rule: |S| = 14 (Pinus sylvestris, Quercus robur/petraea, Betula, Pseudotsuga menziesii, Fagus sylvatica, Pinus nigra, Larix kaempferi, Quercus rubra, Picea abies, Fraxinus excelsior, Alnus glutinosa, Populus, Acer, Salix), x = 2 %: 1762 of 2271 pairs kept (77.6 %); end-census out-of-S share on kept pairs: < 1 % 1301, 1-5 % 326, > 5 % 135
- identifiability of S on the kept pairs: plots 35-911, trees at risk 264-10018, holdout blocks 5-13, climate span >= 0.64; 0 species below the thresholds ()
-   out-of-S trees dropped on the kept pairs: 619 first-visit, 3558 second-visit, 1976 recruits
- expansion: first-visit basal area of S 24.9 m2/ha design-weighted vs 25.0 m2/ha from the virtual trees; 0 of 35500 first-visit trees and 0 of 31949 second-visit trees expand to zero
- obs_dt: 1762 sites x 14 species; intervals 4-12 yr; stem loss over the interval 0.224 of the trees at risk (annual, pooled 0.0508 incl. felling); recruitment mean 69, median 0 stems/ha per interval; tree_obs_dt 34963 tracked trees
- sets after the site rule: 1179 training sites, 583 holdout sites in 15 of 15 holdout blocks
- holdout separation: nearest training site median 7.0 km, minimum 1.0 km, 10 % below 2.2 km
- set national          1762 sites,  24668 obs rows,  35500 cohorts,  34963 tracked trees
- set train             1179 sites,  16506 obs rows,  23529 cohorts,  23167 tracked trees
- set holdout_blocked    583 sites,   8162 obs rows,  11971 cohorts,  11796 tracked trees
- set pilot              300 sites,   4200 obs rows,   6033 cohorts,   5939 tracked trees

## Summary (national)

- 1762 sites (945 MFV -> NBI-6, 817 NBI-6 -> NBI-7); mean basal area at the second visit 26.6 m2/ha, mean stems 656/ha,
  mean annual relative growth 0.0181, pooled annual stem loss 0.0508 (incl. felling), mean recruitment 69 stems/ha per interval

## Sets

- national: every site passing the filters and the site rule
- train: every such site in the training blocks
- holdout_blocked: every such site in the held-out blocks (never for fitting)
- pilot: 300 training sites drawn with seed 1 (smoke tests and pilots)

Spatial blocking: 25 km squares on EPSG:3035, 15 of 76 blocks held out (seed 21, drawn over all pairs before the site rule). Holdout site to
nearest training site: median 7.0 km, minimum 1.0 km, 10 % below 2.2 km, 90 % below 12.0 km.

Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set,
siteID_national, plot, cycles, period_length, INSPIRE cell, lon/lat, tavg_c, prec_mm, harvest, oos_ba_share_t0/t1, block, set).
Site attributes: site_table.csv; species list with NBI codes, Dutch names and the national BA rank: species.csv.

## Species

| species | name | NBI code | Dutch | stems at the first visit | sites | BA rank |
|---|---|---|---|---|---|---|
| 1 | Acer | ED | Esdoorn | 526 | 71 | 13 |
| 2 | Alnus glutinosa | ZE | Zwarte els | 710 | 94 | 11 |
| 3 | Betula | BE | Berk | 5542 | 897 | 3 |
| 4 | Fagus sylvatica | BU | Beuk | 1979 | 271 | 5 |
| 5 | Fraxinus excelsior | ES | Es | 690 | 62 | 10 |
| 6 | Larix kaempferi | JL | Japanse lariks | 1704 | 237 | 7 |
| 7 | Picea abies | FS | Fijnspar | 1277 | 149 | 9 |
| 8 | Pinus nigra | CD+OD | Corsicaanse den / Oostenrijkse den | 1603 | 130 | 6 |
| 9 | Pinus sylvestris | GD | Grove den | 10154 | 911 | 1 |
| 10 | Populus | PO | Populier | 508 | 47 | 12 |
| 11 | Pseudotsuga menziesii | DG | Douglas | 2360 | 299 | 4 |
| 12 | Quercus robur/petraea | EI | Inlandse eik | 6175 | 894 | 2 |
| 13 | Quercus rubra | AE | Amerikaanse eik | 1985 | 312 | 8 |
| 14 | Salix | WI | Wilg | 287 | 35 | 14 |

## PARAMS

```r
List of 32
 $ prep_version            : chr "pairs-v2"
 $ pairs                   :List of 2
  ..$ : Named chr [1:2] "MFV" "NBI6"
  ..$ : Named chr [1:2] "NBI6" "NBI7"
 $ period_range            : int [1:2] 3 15
 $ min_dbh_cm              : num 5
 $ patch_size_ha           : num 0.126
 $ alive_class             :List of 3
  ..$ MFV : int [1:3] 1 2 3
  ..$ NBI6: int [1:3] 1 2 3
  ..$ NBI7: int [1:4] 1 2 3 4
 $ dead_class              :List of 3
  ..$ MFV : int [1:5] 4 5 6 7 8
  ..$ NBI6: int [1:5] 4 5 6 7 8
  ..$ NBI7: int [1:4] 5 6 7 8
 $ removed_class           : int 0
 $ status_cut              : int [1:3] 1 11 12
 $ status_dead             : int [1:4] 2 3 8 9
 $ status_censored         : int [1:2] 4 13
 $ forest_forms            : chr [1:8] "1a.1.1" "1a.1.2" "1a.2" "1c.1" ...
 $ max_ba_loss             : num 0.5
 $ small_tree_dbh          : num 12
 $ recruit_growth_q        : num 0.75
 $ dbh_change_cm           : num [1:2] -2 3
 $ harvest_na_as           : int 0
 $ shrub_group             : chr "ST"
 $ tradeoff_sizes          : int [1:19] 6 7 8 9 10 11 12 13 14 15 ...
 $ tradeoff_shares         : num [1:4] 0.01 0.02 0.03 0.05
 $ species_set_size        : int 14
 $ max_oos_share           : num 0.02
 $ ident_min_plots         : int 30
 $ ident_min_at_risk       : int 100
 $ ident_min_holdout_blocks: int 3
 $ ident_min_span          : num 0.5
 $ env_rasters             : Named chr [1:7] "wc2.1_30s_bio/wc2.1_30s_bio_1.tif" "wc2.1_30s_bio/wc2.1_30s_bio_12.tif" "wc2.1_30s_bio/wc2.1_30s_bio_4.tif" "wc2.1_30s_bio/wc2.1_30s_bio_15.tif" ...
 $ block_km                : num 25
 $ holdout_block_share     : num 0.2
 $ block_seed              : int 21
 $ pilot_n                 : int 300
 $ pilot_seed              : int 1
```
