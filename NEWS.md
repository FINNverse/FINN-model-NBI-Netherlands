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
