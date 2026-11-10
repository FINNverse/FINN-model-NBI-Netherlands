# FINNmodelNBINetherlands 0.2.0

Second fitted release (tag `v0.2.0`), variant `process`. Every change against
0.1.0, with its reason (audit `calibration/AUDIT_v02.md`, pilots
`calibration/v02/pilot_table.csv`):

* **Species set instead of a lumped "other" class** (the rare-species rule agreed
  for the zoo). The model now has 14 taxa (Pinus sylvestris, Quercus
  robur/petraea, Betula, Pseudotsuga menziesii, Fagus sylvatica, Pinus nigra,
  Larix kaempferi, Quercus rubra, Picea abies, Fraxinus excelsior, Alnus
  glutinosa, Populus, Acer, Salix; genus-level NBI codes where the NBI records
  only the genus), no "other" and no shrubs. Reason: in 0.1.0 "other" (38 %
  shrubs) received 11x its observed recruitment, produced 20k of the 23k
  surplus stems on the holdout and took over the free run (120-130 m2/ha).
  **Users must change their inputs: the species index changed (see
  `species_table()`), 20 classes became 14.**
* **Calibration population = stands dominated by the 14 taxa.** A plot pair is
  used only when species outside the set hold < 2 % of the basal area at the
  START of the interval: 1,762 of 2,271 pairs (77.6 % of the pairs = of the
  forest area at one plot per 100 ha); trees of other species on those sites are
  dropped. The end census is never used for the selection; its out-of-set
  share is recorded per site and the validation is stratified by it. Reason
  and trade-off (|S| x threshold against pairs, holdout blocks, climate
  coverage, per-species coverage): `calibration/AUDIT_v02.md` section 22.
* **Recruits by back-modelled start diameter.** A new tree counts as a recruit
  only if its dbh minus the species' 75 % small-tree growth quantile times the
  interval is below the 5 cm tally threshold; the rest is ingrowth (in the
  state, not in recruitment). Reason: half of the 0.1.0 recruits had been
  tally trees at the first visit (audit section 2). Recruits enter the model at
  5.9 cm (0.1.0: 6.4).
* **`harvest` covariate (0/1) in the mortality process**, the third model input
  next to `tavg_c` and `prec_mm` (`inst/model/covariates.csv`): 1 = harvest
  activity recorded at the second visit (NBI-7 traces of harvest, NBI-6
  clearcut / single-tree felling), 0 = none. The mortality target stays total
  stem loss incl. felling. Reason: the felled share of trees is 22 % on
  recorded-harvest plots and < 1.5 % elsewhere; the covariate lifted the
  holdout stem-loss AUC from 0.56 to 0.64 in the pilots. In simulation it is a
  management dial (a regime, not a forecast).
* **Cosine learning-rate decay** (0.01 -> 1e-4 over 600 epochs; 0.1.0: constant
  0.01 for 500). Reason: 0.1.0 never settled (parameter change at the noise
  floor, stem-count loss drifting up); every parameter block now settles to
  < 1e-3 mean change per snapshot (`calibration/v02/figures/fig6b_parameter_change.png`).
* Data version `pairs-v2` (`calibration/v02/04_prepare_pairs_v2.R`,
  `PREPARATION.md`, `expected_preparation_v2.csv`); 1,179 training sites in
  61 of 76 25-km blocks, 583 holdout sites in the other 15; `tree_obs_dt.csv`
  (one row per tracked tree) added to the prepared sets.
* `build_model()` now reads the model structure from the parameter file
  (species count, mortality covariates, bounds, floor), so one definition
  rebuilds the 0.1.0 and the 0.2.0 release files; `parameters.pt` carries
  `settings.*` scalars (recruits_dbh, reg_floor, mortality bounds).
* Tested and NOT adopted (pilots, 300 sites x 400 epochs, blocked holdout):
  freezing the recruitment dispersion theta at 1 with a low cap (recruitment
  ratio 2.5 instead of 1.4: the per-hectare recruitment of 5-20 m circles is
  far more dispersed than theta = 1); per-tree growth / mortality losses (AUC
  0.59 but growth skill lost and stem loss under-predicted); removing FINN's
  1 %/yr floor on the simulated death draw with a Poisson draw (worse on
  every column; recruitment ratio 4.4). The mortality floor therefore stays
  at FINN's default with the gamma draw; the Poisson draw is out, and that
  pilot (P6) counts as a sensitivity run only (Yannek, 2026-11-10). The
  fitted theta stays free (0.11-1.5 per species after 600 epochs, no collapse).
* Fit: `calibration/v02/fit.R`, 600 epochs, 0.90 h on one laptop core (5.4
  s/epoch), FINN 0.2.0 (main 0848b60), torch 0.15.1, `weights = "auto"`
  (baseline-scaled losses; no numeric vector, so `weight_scale` is moot).
* Release assets (GitHub release `v0.2.0`): `parameters.pt` (named list of 34
  tensors incl. env networks, env standardisation and settings; sha256
  `29c8a4c2…7443`) and `model.pt` (whole fitted module; sha256
  `d82572df…bf42`), both in `inst/model/weights.json`. `build_model(parameters.pt)`
  reproduces the fitted model bit for bit on 200 sites. Pinned for simulation:
  one torch thread and `batchsize = 100`.
* Holdout (583 blocked sites, `calibration/v02/report.md`): species-level basal
  area R2 0.91, species-mean dbh R2 0.74, stems R2 0.84; stand basal area R2
  0.635 (persistence 0.619), stand stems R2 0.70 (persistence 0.43); growth r
  0.63; stem-loss AUC 0.655 (predicted loss share 0.255 vs 0.234 observed);
  recruitment ratio 1.5-1.7 (R2 ~0). Stratified by the end-census out-of-set
  share: stand BA R2 0.68 (persistence 0.64) where it is < 1 % (436 sites),
  0.34 (0.44) at 1-5 % (108), 0.57 (0.58) above 5 % (39): the model applies
  where the modelled taxa keep dominating.
* Bare ground (500 yr, 15 patches): 46-55 m2/ha without harvest, 28 m2/ha
  with the harvest dial at 1, converged (< 5 % change over the last 200 yr);
  composition dominated by Pseudotsuga menziesii (33-47 %) with Fagus, Larix,
  Quercus rubra, Acer and Salix, while Pinus sylvestris (3-5 %), Quercus
  robur/petraea (<= 2 %) and Betula (< 1 %) hardly establish. A limitation,
  stated on the card: the free run is not a succession model for the Dutch
  pine-oak-birch forest.

# FINNmodelNBINetherlands 0.1.0

First fitted release (tag `v0.1.0`), variant `process`.

* Data: the public database of the Dutch national forest inventory (Probos /
  Wageningen Environmental Research; archive as served 2026-11-09), prepared
  data `pairs-v1` (2,271 remeasurement pairs of 1,337 permanent plots: 1,069
  MFV 2001-2005 -> NBI-6 2012-2013 and 1,202 NBI-6 -> NBI-7 2017-2021, tree by
  tree, 4-12 years apart; 19 species / species groups + other; WorldClim 2.1
  `tavg_c` and `prec_mm`), 1,531 training sites in 61 of 76 25-km blocks, 740
  holdout sites in the other blocks.
* Model: FINN process model (growth, mortality, saturating regeneration,
  competition), annual timestep, each process with a two-layer SELU network
  (20, 20) on the two climate normals and per-species mechanistic parameters;
  one patch of 0.1257 ha (the 20 m NBI circle). The mortality target is total
  stem loss, natural death and felling together.
* Fit: `calibration/fit.R`, 500 epochs, 1.45 h on one core, FINN 0.2.0 (main
  0848b60), torch 0.15.1, `weights = "auto"` (baseline-scaled losses; no
  numeric weight vector, so `weight_scale` is moot).
* Release assets (GitHub release `v0.1.0`): `parameters.pt` (named list of 28
  tensors incl. env networks and env standardisation, what `load_model()`
  rebuilds from) and `model.pt` (the whole fitted module for
  `torch::torch_load()` under the same FINN version); sha256 of both in
  `inst/model/weights.json`. `build_model(parameters.pt)` reproduces the fitted
  model bit for bit. Pinned for simulation: one torch thread and
  `batchsize = 100`; a different batchsize changes the random draws but not the
  model.
* Holdout (740 blocked sites, `calibration/report.md`): species-level basal
  area R2 0.88, species-mean dbh R2 0.73, stems R2 0.61; stand basal area R2
  0.50 (persistence of the first visit 0.56), stand stems R2 0.14 (persistence
  0.36; the model adds about 40 % stems); growth R2 0.32 (r 0.58); stem-loss
  AUC 0.56 (predicted loss share 0.26 vs 0.23 observed); recruitment R2 < 0
  (over-predicted 2.3x, mostly in the "other" class).

Known limitations:

* The states are predicted no better than persistence; the rates are weak.
  Fit for short projections from an observed stand only, not for succession:
  the bare-ground run converges at 120-130 m2/ha in the "other" class with the
  dominant species absent (see the model card's verdict).
* "Mortality" is stem turnover of a managed forest (natural death and felling
  cannot be separated tree by tree in the NBI); plots that lost half or more of
  their basal area between the visits (clearcuts) were not used.
* FINN's simulated mortality has a floor of 0.01 per model step
  (FINNverse/FINN#51), here 1 % per year, which the fit did not see; with a
  stem-loss target of about 10 % per year it rarely binds, but it caps how low
  the turnover of any species can go in simulation.
* The calibration range is the Dutch climate (roughly 9-11 degC, 700-950 mm);
  outside it the climate networks extrapolate without constraint.

Licences: code and fitted parameters GPL (>= 3). The NBI database is publicly
available without a data-use agreement; use is permitted with correct
attribution (acknowledgement of the Dutch Ministry of Agriculture, Fisheries,
Food Security and Nature, LVVN, and citation of Schelhaas et al. 2022,
*Zevende Nederlandse Bosinventarisatie: methoden en resultaten*, WOt-rapport
142). This repository does not redistribute the data; `calibration/` downloads
them from the source and the fitted parameters are a derived model. Climate
normals: WorldClim 2.1 (Fick & Hijmans 2017, doi:10.1002/joc.5086), free for
academic and other non-commercial use; the rasters are not redistributed, only
values extracted at the plot locations are included. Details in
`calibration/PROVENANCE.md`.

# FINNmodelNBINetherlands (development)

* Repository created from the FINN model template.
