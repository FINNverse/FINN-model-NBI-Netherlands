# Calibration pipeline

Self-contained R scripts that produce the released model, run from the
repository root with `Rscript <script>` and no arguments (house rules:
r-pipeline-conventions, finn-standalone-fit). `run_all.R` runs them in order.

| order | script | what | writes |
|---|---|---|---|
| 1 | `01_download_nbi_raw.R` | Probos public database of the Dutch NFI (NBI-7 archive incl. MFV and NBI-6 tables) and its documentation, SHA-256 pinned | `../../FINN-zoo-data/nbi/raw/`, `PROVENANCE.md` |
| 2 | `02_extract_nbi_tables.R` | the needed Access tables to CSV (mdbtools) | `.../raw/csv/` |
| 3 | `03_prepare_pairs_v1.R` | FINN tables for the two remeasurement pairs (MFV -> NBI-6, NBI-6 -> NBI-7), WorldClim covariates, blocked split | `.../prepared/pairs-v1/{national,train,holdout_blocked,pilot}`, `PREPARATION.md`, `expected_preparation.csv`, climate rows of `PROVENANCE.md` |
| 4 | `fit.R` | the production fit (variant "process") | `.../fits/pairs-v1-process/{model.pt, parameters.pt}`, `loss_history.csv`, `fit_summary.md`, `../inst/model/weights.json` |
| 5 | `validation.R` | holdout report | `report.md`, `validation/*.csv`, `figures/*.png` (and copies into the zoo card folder) |
| – | `fit_pilot.R` | smoke test (300 sites, 3 epochs, one thread), round-trip check and epoch timing | `.../fits/pilot-smoke/` |

Data layout: raw and prepared data, `*.pt` and logs live outside git in
`../../FINN-zoo-data/nbi/` relative to the repository root (the scripts'
`DATA_DIR`); WorldClim 2.1 30 s rasters in `../../worldclim/`. Committed
records: `PROVENANCE.md` (every raw layer with URL, terms, date and checksum),
`PREPARATION.md` (what the tables mean, counts, PARAMS), `expected_preparation.csv`
(summary of the prepared tables for diffing a rerun), `fit_summary.md`,
`loss_history.csv`, `report.md`.

Data terms: the NBI database names no licence; Probos permits use with correct
attribution (acknowledgement of the Dutch ministry LVVN and citation of the
NBI-7 report). This repository therefore redistributes no plot-level rows; the
scripts download them from the source, and the fitted parameters are a derived
model.
