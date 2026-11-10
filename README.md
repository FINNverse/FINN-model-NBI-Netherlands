# FINN-model-NBI-Netherlands

A fitted [FINN](https://github.com/FINNverse/FINN) forest model calibrated to the Dutch National Forest Inventory
(Nederlandse Bosinventarisatie, NBI; public database by Stichting Probos and Wageningen Environmental Research),
packaged as a versioned simulation tool (R package `FINNmodelNBINetherlands`). The model card is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Current release: v0.1.0** (see NEWS.md). The fitted parameters are release assets of this repository;
`load_model()` downloads them once and checks their sha256.

```r
remotes::install_github("FINNverse/FINN@0848b60")                   # FINN 0.2.0, the version the model was fitted with
remotes::install_github("FINNverse/FINN-model-NBI-Netherlands@v0.1.0")
library(FINNmodelNBINetherlands)
torch::torch_set_num_threads(1L)   # pinned: the fit and the validation ran single-threaded
m = load_model()          # downloads the released parameters once, checks sha256, rebuilds with the installed FINN
species_table()           # 19 species / species groups + "other": index, name, NBI code, Dutch name
model_info()              # release tag, files, sha256, FINN commit, patch settings, pinned threads and batchsize
vignette("quickstart", package = "FINNmodelNBINetherlands")
```

Simulations are run with one torch thread and an explicit `batchsize` (100 sites per `simulateForest()` batch,
the value pinned in `model_info()`; `simulate_model()` uses it by default). A different batchsize changes the
random draws but not the model.

What the model is: a FINN process model (growth, mortality, saturating regeneration, competition) for 19 Dutch
tree species / species groups + "other", annual timestep, two WorldClim 2.1 covariates (`tavg_c`, `prec_mm`; see
`inst/model/covariates.csv` for how to derive them), calibrated on 1,531 remeasurement pairs of permanent NBI
plots (MFV 2001-2005 -> NBI-6 2012-2013 and NBI-6 -> NBI-7 2017-2021, tree by tree, 4-12 years apart) and
validated on 740 pairs in spatially blocked holdout squares. One patch = 0.1257 ha (the 20 m NBI circle).
The mortality process is fitted to total stem loss, natural death and felling together, so it describes the
stem turnover of the managed Dutch forest rather than natural mortality alone.

This repository is a tool, not the code of a manuscript: it is versioned for
users (see NEWS.md) and may differ from the models in our papers. The
calibration pipeline that produced each release is in `calibration/`
(`calibration/README.md`; `run_all.R` runs it end to end; `PROVENANCE.md`
pins every raw layer).

License: GPL (>= 3) for code and fitted parameters. The calibration data keep
the terms of their source: the NBI database is publicly available without a
data-use agreement; use is permitted with correct attribution (acknowledgement
of the Dutch Ministry of Agriculture, Fisheries, Food Security and Nature, LVVN,
and citation of Schelhaas et al. 2022, *Zevende Nederlandse Bosinventarisatie:
methoden en resultaten*, WOt-rapport 142). This repository does not redistribute
the data; `calibration/` downloads them from the source and the fitted
parameters are a derived model. Climate normals: WorldClim 2.1 (Fick & Hijmans
2017, doi:10.1002/joc.5086), free for academic and other non-commercial use; the
rasters are not redistributed, only values extracted at the plot locations are
included.
