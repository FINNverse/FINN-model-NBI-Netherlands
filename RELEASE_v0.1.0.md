# Release checklist: FINNmodelNBINetherlands v0.1.0

Prepared 2026-11-10. Nothing here has been executed; tag, push and release only
after the go. All commands from the repository root.

## Commit to tag

`58522b8` on `main` ("FINNmodelNBINetherlands 0.1.0: version, FINN pinned at
0848b60, ..."). It is the v0.1 package (`R/`, `inst/model/`, `DESCRIPTION`
0.1.0, vignette); `calibration/` (incl. the v0.2 work in `calibration/v02/`),
`hpc/` and this file are excluded from the build by `.Rbuildignore`.
`R/definition.R` was rewritten on `main` for v0.2 as a self-describing
`build_model()` (species count and covariate width read from the parameter
file, optional `settings.*` entries); the v0.1 parameter file carries none of
the v0.2 entries and loads under it with the v0.1 defaults (checked: 24
parameter tensors, simulation runs). `inst/model/` is unchanged since the v0.1
fit.

## Assets (outside git), sha256 verified 2026-11-10 with `shasum -a 256`

| asset | local path (relative to this repo) | sha256 | size |
|---|---|---|---|
| parameters.pt | ../../FINN-zoo-data/nbi/fits/pairs-v1-process/parameters.pt | 19d24a826394ff918f7d9cdfb11f79f2cb6b032df90466ef388df2abd3043a99 | 17,050 B |
| model.pt | ../../FINN-zoo-data/nbi/fits/pairs-v1-process/model.pt | 60290b1b772268c1d9d9d138889b7953e82a7a0c956e34b3b280bd260c2ce25c | 15.7 MB |

The same values are in `inst/model/weights.json`; `load_model()` refuses a
file whose sha256 differs.

## Checks done (2026-11-10, one thread, FINN 0.2.0 = 0848b60, torch 0.15.1, R 4.5.0)

- `R CMD build` with the vignette (chunks `eval = FALSE`: they need the release
  download) and `R CMD INSTALL` into a temporary library: OK.
- `load_model(path = <parameters.pt>)`: sha256 match, model rebuilt (24
  parameter tensors), 3-year bare-ground simulation with the pinned batchsize
  runs. `model.pt` sha256 match.
- Vignette renders.

## Steps

```sh
# 0. sanity: the sha256 of the files to upload equal weights.json
shasum -a 256 ../../FINN-zoo-data/nbi/fits/pairs-v1-process/parameters.pt ../../FINN-zoo-data/nbi/fits/pairs-v1-process/model.pt

# 1. tag the release commit
git tag -a v0.1.0 58522b8 -m "FINNmodelNBINetherlands 0.1.0: first fitted release (variant process), FINN 0.2.0 commit 0848b60"

# 2. push main and the tag
git push origin main
git push origin v0.1.0

# 3. create the GitHub release with both assets
gh release create v0.1.0 --repo FINNverse/FINN-model-NBI-Netherlands \
  --title "FINNmodelNBINetherlands 0.1.0" \
  --notes-file RELEASE_NOTES_v0.1.0.md \
  ../../FINN-zoo-data/nbi/fits/pairs-v1-process/parameters.pt \
  ../../FINN-zoo-data/nbi/fits/pairs-v1-process/model.pt

# 4. verify the uploaded assets byte for byte
tmp=$(mktemp -d) && gh release download v0.1.0 --repo FINNverse/FINN-model-NBI-Netherlands -D "$tmp" && shasum -a 256 "$tmp"/*.pt

# 5. (user's go) make the repository public
gh repo edit FINNverse/FINN-model-NBI-Netherlands --visibility public --accept-visibility-change-consequences
```

## After the release

```r
remotes::install_github("FINNverse/FINN-model-NBI-Netherlands@v0.1.0")
library(FINNmodelNBINetherlands); torch::torch_set_num_threads(1L)
m <- load_model()             # downloads from the release URL in weights.json and checks the sha256
```

Then the FINN-zoo card cites tag `v0.1.0` and these sha256 values.

## Release notes

The text of `RELEASE_NOTES_v0.1.0.md` (used by `--notes-file` above).
