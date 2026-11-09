# FINNmodelNBINetherlands (development)

* Calibration pipeline for the first release (variant "process"): the public
  database of the Dutch national forest inventory (Probos / Wageningen
  Environmental Research; archive as served 2026-11-09), prepared data
  `pairs-v1` (2,271 remeasurement pairs of 1,337 permanent plots: 1,069 MFV ->
  NBI-6 and 1,202 NBI-6 -> NBI-7, 4-12 years apart; 19 species / species groups +
  other; WorldClim 2.1 `tavg_c` and `prec_mm`), `fit.R` on 1,531 training
  sites in 61 of 76 25-km blocks, `validation.R` on the 740 holdout sites.
  FINN 0.2.0 (main 0848b60).
* `build_model()` rebuilds the model from the released named list of parameter
  tensors (incl. the env standardisation constants); `inst/model/species.csv`
  and `covariates.csv` describe the inputs.
* No fitted release yet: the production fit is queued; `inst/model/weights.json`
  is filled by `calibration/fit.R` (tag `v0.1.0` planned) and the holdout
  metrics will be listed here with the release.
