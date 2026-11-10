First release of FINNmodelNBINetherlands: a fitted FINN forest model for the Netherlands, calibrated to the tree-by-tree remeasurement of the permanent plots of the Dutch National Forest Inventory (NBI: MFV 2001-2005, NBI-6 2012-2013, NBI-7 2017-2021) and packaged as an R simulation tool.

Variant `process`: 19 species / species groups + "other", annual step on one patch of 0.1257 ha (the 20 m NBI circle), FINN process model (growth, mortality, saturating regeneration, competition) with two-layer SELU climate networks on WorldClim 2.1 `tavg_c` and `prec_mm`. The mortality process is fitted to total stem loss, natural death and felling together. Fitted with FINN 0.2.0 (commit 0848b60), torch 0.15.1, one thread, 500 epochs, baseline-scaled losses (`weights = "auto"`).

Install and use:

```r
remotes::install_github("FINNverse/FINN@0848b60")
remotes::install_github("FINNverse/FINN-model-NBI-Netherlands@v0.1.0")
library(FINNmodelNBINetherlands)
torch::torch_set_num_threads(1L)
m <- load_model()   # downloads parameters.pt once, checks its sha256, rebuilds the model with the installed FINN
```

Simulations are pinned to one torch thread and `batchsize = 100`; a different batchsize changes the random draws but not the model.

Assets:

| file | what | sha256 |
|---|---|---|
| parameters.pt | named list of 28 tensors (env networks and env standardisation included); `load_model()` rebuilds from it | 19d24a826394ff918f7d9cdfb11f79f2cb6b032df90466ef388df2abd3043a99 |
| model.pt | the whole fitted module for `torch::torch_load()` under FINN 0848b60 | 60290b1b772268c1d9d9d138889b7953e82a7a0c956e34b3b280bd260c2ce25c |

Holdout skill (740 remeasurement pairs in 25-km blocks the fit never saw; the persistence of the first visit in brackets): species-level basal area R2 0.88, stems R2 0.61; stand basal area R2 0.50 (0.56), stand stems R2 0.14 (0.36; the model adds about 40 % stems); growth r 0.58; stem-loss AUC 0.56; recruitment R2 below zero (over-predicted 2.3x, mostly in "other").

Limitations: the states are predicted no better than persistence and the rates are weak, so the model is for short projections from an observed stand, not for succession (the bare-ground run converges at 120-130 m2/ha in the "other" class); "mortality" is the stem turnover of a managed forest; FINN's simulated mortality has a floor of 0.01 per annual step (FINNverse/FINN#51), which the fit did not see. Details in NEWS.md and `calibration/report.md`.

Data: `pairs-v1`, rebuilt from the public NBI database by `calibration/` (provenance and checksums in `calibration/PROVENANCE.md`). Licences: code and parameters GPL (>= 3). The NBI database is publicly available without a data-use agreement; use is permitted with correct attribution (acknowledgement of the Dutch Ministry of Agriculture, Fisheries, Food Security and Nature, LVVN, and citation of Schelhaas et al. 2022, Zevende Nederlandse Bosinventarisatie: methoden en resultaten, WOt-rapport 142); this repository does not redistribute the data and the fitted parameters are a derived model. WorldClim 2.1: CC BY-SA 4.0.
