# Prepared data: pairs-v2

Written 2026-11-10 06:08 CET by calibration/v02/04_prepare_pairs_v2.R (data.table 1.18.6.1, terra 1.8.60).

pairs-v2 = pairs-v1 with (1) the shrub group ST removed from the population, (2) recruits defined by a back-modelled
start diameter < 5 cm (species 75 % small-tree growth quantile), (3) a `harvest` 0/1 covariate (harvest activity recorded
at the second visit) in env_dt and siteID_dt. See the script header and calibration/AUDIT_v02.md for the reasons.

Source: Database Zevende Nederlandse Bosinventarisatie (NBI-7), Stichting Probos / Wageningen Environmental Research,
archive as served 2026-11-09 (checksums in raw/SHA256SUMS.txt and calibration/PROVENANCE.md); MFV 2001-2005, NBI-6 2012-2013,
NBI-7 2017-2021 tables. Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the centre of the plot's 1 km INSPIRE cell.

## What the tables mean

- One site = one remeasurement PAIR of a permanent NBI plot (MFV -> NBI-6 or NBI-6 -> NBI-7; a plot measured three times is two
  sites), one patch of 0.1257 ha (the 20 m circle). Every tallied tree stands for 10,000 / plot area stems per ha and is expanded
  to round(that x 0.1257) virtual trees (largest remainder within site x species). `ba` (m2), `trees` and `dbh` (cm) in obs_dt and
  `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.1257 for per-hectare values.
- `year` 0 = first visit, `year` = `period_length` = second visit (interval rounded to whole years, 4-12). Env rows years 1..Tmax.
- `growth` = annualised relative dbh increment of survivors, weighted by the plot's expansion; `growth_n` = tallied survivors behind it.
- `n_at_risk` = tallied first-visit trees that survived, died or were felled (censored: not found, now < 5 cm, outside forest);
  `n_died` = dead + FELLED trees (stem loss; harvest is part of the mortality target, see the script header); `mort` = annual rate.
- `reg` = new trees >= 5 cm whose back-modelled start dbh (dbh - species g75 x interval) is below 5 cm, per hectare over the
  interval, counted inside the smaller of the two circles; larger new trees are ingrowth / missed trees (in the state, not recruits).
- Shrubs (NBI species group ST) are not part of the population: not in the stand, the fates, the state or the recruits.
- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), temp_season (bio4), prec_season (bio15, CV %),
  tmax_warm_c (bio5), tmin_cold_c (bio6), elev_m, x = longitude, y = latitude (WGS84) of the 1 km cell centre,
  harvest = 1 if harvest activity was recorded at the second visit (NBI-7 traces of harvest; NBI-6 clearcut or single-tree felling), else 0.

## Processing log

- v2: 15576 rows of the shrub group ST removed from the tree tables (BW DK GW HU HZ KM KN LB LG MD OW SD ST SW TX VB VK VL)
- tree rows: MFV 69848, NBI-6 73595, NBI-7 69080; plots with a measured circle: MFV 2963, NBI-6 3191, NBI-7 3173
- v2: harvest recorded at the second visit on 1177 pairs, none on 1463, no record on 9 (set to 0)
- remeasurement pairs listed by the database: 2649 (MFV -> NBI-6 1236, NBI-6 -> NBI-7 1413)
-   both visits with a measured circle: 2648
-   first visit complete (MFV flag): 2648
-   stand type at the second visit kept (1a.1.1 1a.1.2 1a.2 1c.1 1c.2 1c.3 1c.5 3): 2426
-   interval 3-15 whole years: 2426
-   pairs whose 1 km cell differs between the visits (neighbouring cells, second-visit cell used): 181
- first-visit living trees on the pairs: 48544 (0 below 5 cm dropped, 1283 outside the first-visit circle dropped)
-   fates: survivor 34859, cut 6025, dead 5393, censored 984
-   sites with a living first-visit tree >= 5 cm: 2383 (43 dropped)
-   basal area loss below 50 %: 2264 (NBI6 -> 57, MFV -> 62 clearcut / converted dropped; median basal-area change +10 %)
-   new trees at the second visit: 4802 kept as recruits (median dbh 6.1 cm); 3076 in the ring beyond the smaller circle dropped, 2066 with a back-modelled start dbh >= 5 cm (species g75 0.18-0.67 cm/yr, pooled 0.39) counted as ingrowth, not recruitment
- species: 42 codes among first-visit living trees; 14 named species cover 95.5 % of stems, 27 codes lumped into 'other'
- expansion: first-visit basal area 24.9 m2/ha design-weighted vs 25.0 m2/ha from the virtual trees; 0 of 44934 first-visit trees and 0 of 41357 second-visit trees expand to zero
- env: WorldClim extracted for 1426 of 1427 cells (2 sites in 1 cells without a land pixel dropped)
- obs_dt: 2262 sites x 15 species; intervals 4-12 yr (Tmax 12); stem loss over the interval 0.223 of the trees at risk (annual rate, pooled: 0.0556 incl. felling); recruitment mean 86, median 19 stems/ha per interval
- blocks: 15 of 76 25-km blocks held out; 1524 training sites, 738 holdout sites
- holdout separation: nearest training site median 7.0 km, minimum 1.0 km, 10 % below 2.2 km
- set national          2262 sites,  33930 obs rows,  44885 cohorts
- set train             1524 sites,  22860 obs rows,  30037 cohorts
- set holdout_blocked    738 sites,  11070 obs rows,  14848 cohorts
- set pilot              300 sites,   4500 obs rows,   6066 cohorts

## Summary (national)

- 2262 sites (1067 MFV -> NBI-6, 1195 NBI-6 -> NBI-7); mean basal area at the second visit 26.8 m2/ha, mean stems 685/ha,
  mean annual relative growth 0.0186, pooled annual stem loss 0.0556 (incl. felling), mean recruitment 86 stems/ha per interval

## Sets

- national: every site passing the filters
- train: every site in the training blocks
- holdout_blocked: every site in the held-out blocks (never for fitting)
- pilot: 300 training sites drawn with seed 1 (smoke tests only)

Spatial blocking: 25 km squares on EPSG:3035, 15 of 76 blocks held out (seed 21). Holdout site to nearest training site:
median 7.0 km, minimum 1.0 km, 10 % below 2.2 km, 90 % below 12.0 km.

Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set, siteID_national,
plot, cycles, period_length, INSPIRE cell, lon/lat, tavg_c, prec_mm, block, set). Site attributes: site_table.csv; species list
with NBI codes and Dutch names: species.csv.

## Species

| species | name | NBI code | Dutch | stems at the first visit | sites |
|---|---|---|---|---|---|
| 1 | Acer | ED | Esdoorn | 756 | 124 |
| 2 | Alnus glutinosa | ZE | Zwarte els | 1051 | 148 |
| 3 | Betula | BE | Berk | 6876 | 1160 |
| 4 | Fagus sylvatica | BU | Beuk | 2349 | 355 |
| 5 | Fraxinus excelsior | ES | Es | 1009 | 114 |
| 6 | Larix kaempferi | JL | Japanse lariks | 1969 | 282 |
| 7 | Picea abies | FS | Fijnspar | 1512 | 191 |
| 8 | Pinus nigra | CD+OD | Corsicaanse den / Oostenrijkse den | 1728 | 147 |
| 9 | Pinus sylvestris | GD | Grove den | 11184 | 1094 |
| 10 | Populus | PO | Populier | 583 | 67 |
| 11 | Prunus serotina | AV | Amerikaanse vogelkers | 666 | 167 |
| 12 | Pseudotsuga menziesii | DG | Douglas | 2718 | 362 |
| 13 | Quercus robur/petraea | EI | Inlandse eik | 8106 | 1219 |
| 14 | Quercus rubra | AE | Amerikaanse eik | 2402 | 403 |
| 15 | other | WI+AC+IL+HB+AG+LI+TK+SS+TS+SA+UL+ZD+RP+IE+ZK+CH+WD+PC+GE+OS+JE+NE+AB+SO+EL+PK+TH | all other species and shrubs | 2025 |  |

## PARAMS

```r
List of 25
 $ prep_version       : chr "pairs-v2"
 $ shrub_group        : chr "ST"
 $ small_tree_dbh     : num 12
 $ recruit_growth_q   : num 0.75
 $ harvest_na_as      : int 0
 $ pairs              :List of 2
  ..$ : Named chr [1:2] "MFV" "NBI6"
  ..$ : Named chr [1:2] "NBI6" "NBI7"
 $ period_range       : int [1:2] 3 15
 $ min_dbh_cm         : num 5
 $ patch_size_ha      : num 0.126
 $ alive_class        :List of 3
  ..$ MFV : int [1:3] 1 2 3
  ..$ NBI6: int [1:3] 1 2 3
  ..$ NBI7: int [1:4] 1 2 3 4
 $ dead_class         :List of 3
  ..$ MFV : int [1:5] 4 5 6 7 8
  ..$ NBI6: int [1:5] 4 5 6 7 8
  ..$ NBI7: int [1:4] 5 6 7 8
 $ removed_class      : int 0
 $ status_cut         : int [1:3] 1 11 12
 $ status_dead        : int [1:4] 2 3 8 9
 $ status_censored    : int [1:2] 4 13
 $ forest_forms       : chr [1:8] "1a.1.1" "1a.1.2" "1a.2" "1c.1" ...
 $ max_ba_loss        : num 0.5
 $ dbh_change_cm      : num [1:2] -2 3
 $ species_quantile   : num 0.95
 $ env_rasters        : Named chr [1:7] "wc2.1_30s_bio/wc2.1_30s_bio_1.tif" "wc2.1_30s_bio/wc2.1_30s_bio_12.tif" "wc2.1_30s_bio/wc2.1_30s_bio_4.tif" "wc2.1_30s_bio/wc2.1_30s_bio_15.tif" ...
 $ block_km           : num 25
 $ holdout_block_share: num 0.2
 $ block_seed         : int 21
 $ pilot_n            : int 300
 $ pilot_seed         : int 1
```
