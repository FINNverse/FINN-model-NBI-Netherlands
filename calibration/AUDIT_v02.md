# Audit of the NBI data and the v0.1 fit (zoo v0.2 round, 2026-11-10)

Scope: the pairs-v1 preparation (`calibration/03_prepare_pairs_v1.R`), the v0.1 fit
(`fits/pairs-v1-process`, `calibration/report.md`) and the data traps listed in the v0.2
brief plus the NBI-specific ones. Evidence comes from scratch scripts run on the raw CSV
tables and on `fits/pairs-v1-process/validation/predictions.csv` (holdout = 740 blocked
sites); numbers are quoted as computed. Sections 1-12 follow the brief's trap list,
13-20 the NBI-specific checks, 21 the v0.1 fit diagnostics, 22 the candidates.

## 1. Env rows vs observation windows (trap 1): PASS
env_dt holds years 1..12 for every site, obs sits at `year = period_length` (4-12).
FINN's loop compares `Result[, i, ]` at loop index i = year_sequence with the observation,
and `simulateForest()` output year 1 already differs from the init state (site 1: 96 -> 103
stems, BA 3.18 -> 3.14 m2 per patch), so output year k is the state after k env rows and
the observation after p rows is compared with the state after p steps. No extra row 0.
The holdout `pilot_eval`/validation aggregation uses the same convention
(`S[year == period_length]`, rates over `1..period_length`).

## 2. Ingrowth vs recruitment (trap 2): FAIL in v1, fixed in pairs-v2
v1 called every second-visit tree without a first-visit row a recruit if its dbh was
<= 5 + 1.5 cm/yr x interval (dropped 10 % as "too large"). Back-modelling the start
diameter with the species' survivor growth of small trees (< 12 cm at the first visit;
median 0.18 cm/yr pooled, 0.13-0.36 per species; 75 % quantile 0.35, 90 % 0.59):
- 51 % of the v1 recruits had a start dbh >= 5 cm at the median growth, 23 % at the
  90 % quantile (recruit dbh quantiles: 25 % 5.5, 50 % 6.4, 75 % 8.4, 90 % 15.6, 95 % 29.4 cm).
- Every tree >= 5 cm inside the first-visit circle had to be tallied, so a "new" tree with
  a back-modelled start dbh >= 5 cm is a missed tree or an unusually fast grower, not
  recruitment.
pairs-v2 keeps a new tree as a recruit only if dbh - g75_species x interval < 5 cm
(generous towards recruits): 4,802 recruits instead of 9,079 trees; 2,066 "new" trees
become ingrowth (in the state, not in `reg`). Mean recruitment falls from 162 to 86
stems/ha/interval (median 65 -> 19).

## 3. Sample sizes vs expansion (trap 3): PASS
`n_at_risk`, `n_died`, `growth_n` count TALLIED trees (median n_at_risk per site x
species row 5; quartiles 2-10); only the states are expanded (NstamHA x 0.1257 ha,
largest-remainder integers; design BA 25.0 vs 25.1 m2/ha from the virtual trees).
Side effect to keep in mind: the Poisson loss on `trees` acts on the EXPANDED counts
(a 5 m circle inflates each tree x16), so small circles weigh more in that term than
their information warrants (see 21.4).

## 4. Harvest (trap 4): target kept, covariate added in pairs-v2
Felling is coded per tree (`boomklasse 0` + `oogst_status 1/11/12`), death and
windthrow separately (classes 5-8; status 2/3/8/9). Of the first-visit living trees on
the NBI-6 -> NBI-7 pairs, 12 % were felled and 12 % died in ~6.4 years; pooled annual
stem loss incl. felling 6.0 %/yr (NBI counts every stem >= 5 cm). The target stays
total stem loss (the user's rule; censoring felled trees would bias natural mortality
low). A plot-level harvest record exists at BOTH second visits:
- NBI-7 `oogstactiviteit` (2 = traces of harvest): felled share 21.6 % vs 1.4 % without
  (813 vs 386 plots, 86 % vs 8 % of plots with any felled tree);
- NBI-6 `soort kap` (1 clearcut, 2 single-tree felling, 3 none): 51 % / 22 % / 0.6 %
  (9 / 578 / 482 plots);
- MFV `KAP` is empty (MFV is never a second visit).
pairs-v2 writes `harvest` (0/1) into env_dt; 1,177 of 2,649 pairs are 1, 9 NBI-7 pairs
without a record are set to 0. As in BWI, the flag is realised harvest of the same period
(a leak for prediction, a dial in simulation); the pilots test it in the mortality
formula only.

## 5. Fates and status codes per cycle (trap 5): PASS (handled in v1)
Tree class 4 = "dead" in MFV and NBI-6, "suppressed" (alive) in NBI-7: coded per cycle
(`alive_class`, `dead_class` lists). 894 NBI-6-living trees are class 4 at NBI-7 and
counted as survivors; NBI-6 class-4 trees are dead at the start of a NBI-6 -> NBI-7 pair
and excluded from the stand. Censored (no row / status 4 "now < 5 cm" / 13 "outside
forest"): 1,008 of 49.5k first-visit trees (2 %).

## 6. Lumped "other" (trap 6): FAIL in v1, changed in pairs-v2
v1 "other" = 35 codes, 2,246 first-visit stems, of which 38 % belong to the NBI shrub
group ST (hawthorn 1,288 stems over the three tables, Amelanchier 1,040, elder 955,
holly 736, bird cherry 540, ...); the named set also held three ST codes (Corylus,
Rhamnus frangula, Sorbus aucuparia). 23 % of the v1 recruits were ST codes. In the
v0.1 holdout, "other" received 275 predicted vs 24 observed recruits/ha/interval and
25,693 predicted vs 5,112 observed stems (20k of the 23k stem surplus); in the free run it
took over (120-130 m2/ha). pairs-v2 removes the whole ST group from the population
(15,576 tree rows), which also removes the three named shrub classes; "other" is now 27
tree codes (field maple, hornbeam, elm, lime, sweet chestnut, Tsuga, Abies grandis, Sitka,
...), 2,025 stems. The species set becomes 14 named + other (95.5 % of stems).

## 7. Recruitment dispersion (trap 7): FAIL in v1
`par_theta_recruits` in the v0.1 fit: 0.066 for 11 of 20 species (raw 9.6-15.6 ->
theta = 1/softplus), 0.4-1.3 for the others; `par_theta_recruits_raw` was the one
block still moving at epoch 500 (0.17 mean abs change per snapshot). The collapsed
species are NOT the over-predicted ones (those are Acer, Alnus, Fagus, other with theta
0.5-1.3), so the collapse did not cause the surplus here, but it removes the penalty on
over-prediction. Pilots P2-P6 freeze theta at 1; P5 moves the likelihood to circle-count
space (`recruit_obs_weight` = median circle area 0.0314 ha).

## 8. Steep regeneration climate response (trap 8): not the mechanism here
The climate space is 8.9-10.4 degC x 712-940 mm; the fitted K (Beverton-Holt cap) is
0.2-800 stems/ha/yr per species and binds for the over-predicted species (Acer K 11.5,
pred 83/ha/interval = at the cap; other K 34.5). The surplus comes from the light
response: shade-tolerant classes (Acer, Alnus, Fagus, other) recruit under closed
canopies where the data show almost nothing, while light-demanders are under-predicted
(Pinus sylvestris 3.7 vs 10.6, Douglas 4.8 vs 10.6 /ha/interval).

## 9. Loss competition (trap 9): visible
Mortality loss flat at its baseline for 500 epochs (0.92-1.2); mortality species
parameters drifted to +-9 against the +-10 default bounds (BWI: big trees immortal);
the stem-count loss rose from 0.4 to 0.5 after epoch 250 while dbh/BA/growth stayed flat.
Pilots: bounds [-3, 3] (P2+), two stages (P6).

## 10. Convergence (trap 10): not settled
Mean absolute parameter change per 10-epoch snapshot stays at 0.03-0.07 (noise floor)
from snapshot 10 on; theta and reg_logK still moving. Constant lr 0.01. Pilots use
cosine decay.

## 11. Evaluation (trap 11): extended
v0.1 report had persistence for the states and the weighted AUC. Added for v0.2:
change metrics (dBA R2 -0.21, dN R2 -0.37 on the holdout), per interval length and per
pair type (section 15), recruitment ratio per species, bare-ground BA and "other" share.

## 12. FINN interface (trap 12): PASS
`period_length` present and finite everywhere; per-site windows ("period_length varies
between sites") work with mini-batches that lack a rare interval length (tested with 20-site
batches); patches constant (1); models rebuilt through `build_model()`.

## 13. Variable-radius circle (NBI-specific): understood, one side effect
Radius 5-20 m for >= 20 trees >= 5 cm; `NstamHA` = 10,000 / forest area of the circle,
exactly the nominal circle for 93 % of the plots, smaller where part of the circle is
not forest. Checked: median NstamHA per radius equals 10,000 / (pi r^2) to the digit.
The expansion to the 0.1257 ha patch inflates a 5 m plot x16; the holdout stem surplus
is +411 stems/ha on 14-20 m circles vs -35 on <= 7 m circles, i.e. it is a mature-stand
(closed canopy) effect, not a small-circle artefact. The Poisson `trees` term on inflated
counts remains a concern (21.4).

## 14. Trees beyond the current circle: handled
Permanent-plot trees are reported back even outside a shrunken circle (2,584 NBI-6
trees lie beyond the NBI-7 radius); the end state uses only trees inside the second-visit
circle, the cohort uses first-visit trees inside the first-visit circle. Distances are in
dm in all tables (checked against the radii).

## 15. Interval spread 4-12 yr and pair type
Holdout stand changes by interval (obs vs pred, /ha): dBA tracks the observed increase
(e.g. 8 yr: +1.78 vs +1.94 m2/ha; 12 yr: +4.0 vs +2.5), dN is wrong in sign at every
length (obs -10 to -103, pred +183 to +443) and grows with the interval, as recruitment
does (pred 328 -> 791 vs obs 125 -> 195). MFV -> NBI-6 (7-12 yr) and NBI-6 -> NBI-7
(4-9 yr) behave alike (dBA 2.28 vs 2.30; 1.13 vs 1.03). No interval-specific bias beyond
recruitment.

## 16. MFV positions missing
MFV tree distances are 0 for most trees; the "inside the circle" filter therefore cannot
act on MFV stands (MFV was the first visit of the design, so its trees are inside by
construction) and recruits of MFV -> NBI-6 pairs are located with the NBI-6 distance.
Recruitment per pair type is similar (168 vs 145 /ha/interval observed), so no visible
artefact.

## 17. Genus-level codes
EI (Quercus robur/petraea), BE (Betula), ED (Acer), PO (Populus), WI (Salix) are groups
in the NBI; kept as groups (nothing in the data separates them). For the general zoo
model they cannot map to species-level taxa (open question).

## 18. Clearcut filter
Sites with a basal-area loss >= 50 % (116 of 2,383) were dropped; the felled share of
first-visit trees among the retained NBI-6 -> NBI-7 sites rises steeply with BA loss
(by BA-change class, see scratch audit3), so heavy thinnings are inside the target. Kept
at 50 % (a stand that lost more is a regeneration cut, not stand dynamics).

## 19. Growth: fine
Holdout growth obs vs pred by size class 0.029/0.028, 0.019/0.021, 0.014/0.013,
0.012/0.010, 0.010/0.008 (relative dbh increment per year): a slight under-prediction
for big trees; per species within 10-15 % except Sorbus / Rhamnus (shrubs, gone in v2).

## 20. Stem loss: level right, discrimination absent
Predicted loss share 0.26 vs observed 0.23 over the interval; by size class 0.32/0.26
(< 10 cm), 0.28/0.26, 0.24/0.22, 0.18/0.17, 0.15/0.14 (pred/obs). AUC 0.556 because
felling (half the loss) is not predictable from size and climate; the `harvest` flag
(section 4) is the candidate.

## 20b. FINN's simulated-mortality floor (model side; brief update 05:25 UTC)
FINN clamps the annual death probability used to DRAW deaths at [0.01, 0.99]
(`torch_clamp(m, 0.01, 0.99)` in the forward pass at 0848b60); the rate stored for the
site-level binomial loss is the unclamped `m`. For this model that means:
- the FITTED rates (what the `mortality` loss sees, and the `mort` column of
  `simulateForest()`) are not floored; the AUC and the pred/obs loss ratios in section 20
  and in the pilot table are unaffected;
- the simulated STATES (trees, BA, dbh) and every free run lose at least 1 % of every
  cohort per year, whatever the fitted rate. Where it plausibly acts here: the observed
  annual loss incl. felling is 5.6-6.0 %/yr pooled and above 1.5 %/yr for every named
  species (lowest: Fagus 1.5 %/yr, Quercus robur/petraea 2.5 %/yr on the holdout
  intervals), so within the calibration intervals the floor binds only for large trees of
  low-mortality species on unharvested plots (harvest = 0 sites: natural death ~2 %/yr
  over all sizes, less for dbh > 40 cm). In equilibrium runs it caps the mean residence
  time of a cohort at ~100 years and lowers the equilibrium BA of long-lived species
  (beech, oak) relative to the fitted rates. Not compensated in data or weights (the
  brief's instruction); recorded for Max's side.

## 21. v0.1 fit diagnostics
1. Loss history: dbh, BA, growth plateau by epoch 20; mortality flat at 1.0-1.2
   (baseline) throughout; trees drifts up from epoch 250; regeneration still falling at 500.
2. Parameters: mortality species parameters at +-9 (bounds +-10); theta collapsed
   (section 7); reg_logK 0.2-800; competition and growth parameters well inside bounds.
3. Holdout: species-level BA R2 0.88, dbh 0.73, stems 0.61; stand BA R2 0.50 vs
   persistence 0.56; stand stems 0.14 vs 0.36 (+40 stems/patch); growth r 0.58; loss AUC
   0.556; recruitment R2 < 0, ratio 2.3 (other 11.6, Acer 30, Alnus 39, Fagus 3.7).
4. Candidate-specific concern: `trees` Poisson on x4-x16 inflated counts; a pilot with
   `trees = "none"` (BA + dbh as the only state terms, the IFN-Spain choice) is a fallback
   if P1-P6 leave the stem surplus.

## 22. The species set S and the site rule (agreed rare-species rule, V02_BRIEF 05:40 UTC)
No lumped "other" class and no shrubs in v0.2. S = the top-|S| taxa by national first-visit
basal area (genus-level NBI codes as the NBI records them: Quercus robur/petraea, Betula,
Acer, Populus, Salix; CD + OD = Pinus nigra); a pair is kept when species outside S (NBI
group codes and the shrub group ST included) hold < x % of the FIRST-visit basal area; the
end census is never used for the selection. Trade-off (`calibration/v02/species_subset_tradeoff.csv`;
2,271 pairs before the rule; 15 holdout blocks; coverage = 5-95 % span of the kept pairs'
tavg_c / prec_mm relative to all pairs; identifiability columns = the weakest species in S
over the kept pairs: plots with the species at the first visit, trees at risk, holdout blocks
with it, 5-95 % climate span of its plots relative to the overall span):

| S | x % | pairs kept | share | train / holdout | holdout blocks (>= 10 sites) | tavg / prec coverage | BA share of S | min plots | min at risk | min holdout blocks | min span | species below thresholds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 1,343 | 0.59 | 871 / 472 | 14 (10) | 0.95 / 1.00 | 0.84 | 115 | 1,391 | 8 | 0.79 | 0 |
| 8 | 2 | 1,445 | 0.64 | 932 / 513 | 14 (10) | 0.95 / 1.00 | 0.84 | 121 | 1,461 | 8 | 0.80 | 0 |
| 10 | 2 | 1,595 | 0.70 | 1,038 / 557 | 14 (10) | 0.97 / 1.00 | 0.90 | 22 | 340 | 4 | 0.62 | 1 |
| 12 | 1 | 1,569 | 0.69 | 1,038 / 531 | 15 (10) | 0.98 / 1.00 | 0.93 | 35 | 441 | 5 | 0.54 | 0 |
| 12 | 2 | 1,678 | 0.74 | 1,107 / 571 | 15 (10) | 0.98 / 1.00 | 0.93 | 37 | 442 | 6 | 0.54 | 0 |
| 14 | 1 | 1,651 | 0.73 | 1,108 / 543 | 15 (10) | 0.99 / 1.00 | 0.95 | 33 | 252 | 5 | 0.61 | 0 |
| **14** | **2** | **1,762** | **0.78** | **1,179 / 583** | **15 (10)** | **0.99 / 1.00** | **0.95** | **35** | **264** | **5** | **0.64** | **0** |
| 14 | 3 | 1,836 | 0.81 | 1,228 / 608 | 15 (10) | 0.99 / 1.00 | 0.95 | 35 | 264 | 5 | 0.65 | 0 |
| 14 | 5 | 1,917 | 0.84 | 1,279 / 638 | 15 (10) | 0.99 / 1.00 | 0.95 | 36 | 265 | 6 | 0.65 | 0 |
| 16 | 2 | 1,798 | 0.79 | 1,209 / 589 | 15 (11) | 0.99 / 1.00 | 0.96 | 15 | 66 | 1 | 0.56 | 2 |
| 18 | 2 | 1,881 | 0.83 | 1,264 / 617 | 15 (11) | 1.00 / 1.00 | 0.97 | 16 | 71 | 1 | 0.56 | 3 |
| 20 | 2 | 1,910 | 0.84 | 1,283 / 627 | 15 (11) | 1.00 / 1.00 | 0.97 | 14 | 71 | 1 | 0.51 | 5 |
| 24 | 2 | 1,951 | 0.86 | 1,314 / 637 | 15 (12) | 1.00 / 1.00 | 0.98 | 3 | 27 | 0 | 0.05 | 9 |

(thresholds: >= 30 plots, >= 100 trees at risk, >= 3 holdout blocks, climate span >= 0.5 of
the overall span; the full grid |S| = 6..24 x x = 1/2/3/5 % is in the CSV.)

**Choice: |S| = 14, x = 2 %.** Reasons: (i) every species of S clears the coverage thresholds
(weakest: Salix with 35 plots, 264 trees at risk, 7 holdout blocks, span 0.80-0.91; Acer with
71 plots, 5 holdout blocks, span 0.64); the next two taxa in the ranking (rank 15-16) sit in
ONE holdout block with 66 trees at risk, so |S| = 16 adds 36 pairs and two species the holdout
cannot test; (ii) the 15 holdout blocks and 99-100 % of the climate span stay populated from
x = 1 % on, so x is not a coverage question here; x = 1 % costs 111 pairs for no coverage gain,
x = 3-5 % adds 74-155 pairs but lets up to 5 % of the start basal area be dropped silently
(option (a)); 2 % is the smallest threshold whose cost is small (6 % of the pairs vs x = 3 %);
(iii) S covers 95.1 % of the national first-visit basal area. Result: 1,762 of 2,271 pairs
(77.6 %): the model is valid for stands where the 14 modelled taxa hold > 98 % of the basal
area at the start of the interval; with one NBI plot per 100 ha that share of pairs is also
the share of forest area. End-census out-of-S share on the kept pairs: < 1 % on 1,301 sites,
1-5 % on 326, > 5 % on 135 (the strata of the validation). Dropped on the kept pairs: 619
first-visit trees, 3,558 second-visit trees, 1,976 recruits (option (a)). The n (holdout
blocks) and span thresholds are COVERAGE choices, not measured identifiability limits.

Per-species identifiability of S (`species_identifiability.csv`): plots 35 (Salix) - 911
(Pinus sylvestris); trees at risk 264 - 10,018; holdout blocks 5 (Acer) - 13; climate span
0.64 - 1.00.

Option (b) of the agreed rule (out-of-S trees as one fixed shade-only class) is NOT piloted:
FINN's `createProcess(optimizeSpecies = FALSE)` freezes the species parameters of a whole
process, not of one class, and nothing in the interface freezes one species row or masks its
gradient; its parameters would still move through the shading term of the other species'
losses. A per-species freeze is one of the general agent's FINN proposals; until it exists,
option (a) at x = 2 % is used.

## 23. Candidates tested as pilots (calibration/v02/pilot_*.R, 300 training sites, 400 epochs,
## one thread, scored on the full blocked holdout of the pilot's data version, stratified by the
## end-census out-of-S share; tables calibration/v02/pilot_table.csv and pilot_strata.csv)
| pilot | data | the one change | finding |
|---|---|---|---|
| P0 | pairs-v1 | v0.1 settings (reference with "other") | - |
| P1 | pairs-v2 | v0.1 settings = BASELINE of the new data | 22 |
| P2 | pairs-v2 | + cosine lr decay | 10 |
| P3 | pairs-v2 | + theta frozen at 1, K init 10 (observed level) | 7 |
| P4 | pairs-v2 | + per-tree growth / mortality losses (tree_obs_dt) | 3, 13, 20 |
| P5 | pairs-v2 | + `harvest` in the mortality formula | 4, 20 |
| P7 | pairs-v2 | P2 + P3 + P5 together (presumed production shape) | - |
| P6 | pairs-v2 | P7 + `mortality_floor = 1e-6`, `mortality_draw = "poisson"` (branch build) | 20b |

P0 result (pairs-v1 holdout, 740 sites): species-level BA R2 0.86, stand BA R2 0.44 vs
persistence 0.56, stand stems R2 0.51 vs 0.36, dBA R2 -0.36, dN R2 0.22, stem surplus
+64/ha, growth r 0.53, loss AUC 0.565, loss ratio 1.39, recruitment ratio 2.21, bare ground
38.5 m2/ha after 300 yr (Fagus 57 %). Rows P1-P7 and the v0.2 decision follow in section 24
when the pilots have run (laptop budget permitting, or the cluster via calibration/v02/CLUSTER.md).
