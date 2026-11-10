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
* Production fit done 2026-11-10 (500 epochs, 1.45 h on one core); `inst/model/weights.json`
  carries the SHA-256 of both release assets (`parameters.pt`
  `19d24a82…3a99`, `model.pt` `60290b1b…ce25c`); tag `v0.1.0` planned.
* Holdout (740 blocked sites, `calibration/report.md`): species-level basal
  area R² 0.88, dbh R² 0.73, stems R² 0.61; stand basal area R² 0.50 (persistence
  0.56), stand stems R² 0.14 (persistence 0.36, model +40 % stems); growth r 0.58;
  stem-loss AUC 0.56; recruitment R² < 0 (over-predicted 2.3x, "other" class).
  Bare-ground run converges at 120-130 m²/ha in the "other" class with the dominant
  species absent. Fit for short projections from an observed stand only; not for
  succession (see the model card's verdict).
