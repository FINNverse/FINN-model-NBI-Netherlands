#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as `calibration/fit.R`, followed by copying the fitted
#' parameters into the new object. Rebuilding (rather than `torch_load()` of a
#' whole model object) means the model always runs with the methods of the
#' installed FINN.
#'
#' Variant `"process"`: FINN process model for the Dutch tree species of
#' [species_table()] (19 species / species groups + "other") on an annual
#' timestep with two WorldClim 2.1 climate normals as covariates, `tavg_c`
#' (annual mean temperature, degC, bio1) and `prec_mm` (annual precipitation,
#' mm, bio12), each process (growth, mortality, saturating regeneration) with a
#' two-layer SELU network (20, 20) on the covariates and per-species
#' mechanistic parameters; competition without covariates. Calibrated on one
#' patch of 0.1257 ha per site (the 20 m NBI circle) to the tree-by-tree
#' remeasurement of the Dutch national forest inventory (NBI), where the
#' mortality target is total stem loss including felled trees.
#'
#' The parameter list must carry, besides the tensors of `model$parameters`,
#' the standardisation constants `env_scaling.<variable>.center` and
#' `env_scaling.<variable>.scale` learned at fitting: FINN z-scales the
#' covariates internally and reuses these constants at prediction.
#'
#' @param variant Name of the model variant (see [model_info()]).
#' @param parameters Named list of parameter tensors, as saved by the
#'   calibration script (`parameters.pt`).
#' @return A fitted `FINN` model.
#' @keywords internal
build_model = function(variant = "process", parameters) {
  if (!identical(variant, "process"))
    stop("Unknown variant '", variant, "'; this package has the variant 'process'.", call. = FALSE)
  n_species = nrow(species_table())                       # 20 = 19 species / groups + "other"
  env_vars  = c("tavg_c", "prec_mm")
  hidden    = c(20L, 20L)
  patch_size = pi * 20^2 / 1e4                            # 0.1257 ha, the 20 m NBI circle

  ## createProcess() evaluates its formula in the calling frame's parent, so the
  ## formula is written out literally in every call (a local variable would not be found)
  model = FINN::finn(
    N_species    = n_species,
    recruits_dbh = 6.4,                                   # median dbh of the observed recruits (tally threshold 5 cm)
    reg_floor    = 1e-3,                                  # stems/ha/yr: minimum recruitment (immigration)
    competition_process  = FINN::createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = FINN::createProcess(~ tavg_c + prec_mm, FINN::growth,
                                               hidden = hidden, dropout = 0.1,
                                               optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = FINN::createProcess(~ tavg_c + prec_mm, FINN::mortality,
                                               hidden = hidden, dropout = 0.1,
                                               optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = FINN::createProcess(~ tavg_c + prec_mm, FINN::regeneration_saturation,
                                               hidden = hidden, dropout = 0.1,
                                               optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                               custom_parameters = list(reg_logK = rep(log(50), n_species)))
  )

  ## the env standardisation learned at fitting (natural units in, z-scores inside)
  for (nm in paste0("env_scaling.", rep(env_vars, each = 2), c(".center", ".scale")))
    if (is.null(parameters[[nm]])) stop("parameters.pt lacks ", nm, call. = FALSE)
  model$env_scaling = data.frame(
    variable = env_vars,
    center   = vapply(env_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".center")]]), numeric(1)),
    scale    = vapply(env_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".scale")]]), numeric(1)),
    stringsAsFactors = FALSE)
  ## FINN builds the env networks lazily on the first simulation (their input
  ## width comes from the env table), so one minimal bare-ground step creates
  ## them before the fitted weights are copied in
  dummy_env = data.frame(siteID = 1L, year = 1L,
                         tavg_c = model$env_scaling$center[1], prec_mm = model$env_scaling$center[2])
  invisible(utils::capture.output(suppressMessages(
    FINN::simulateForest(model, env = dummy_env, init_cohort = NULL, patches = 1L,
                         patch_size = patch_size, device = "cpu"))))

  ## the fitted tensors, matched by NAME; anything missing or unexpected stops here
  own = model$parameters
  missing = setdiff(names(own), names(parameters))
  if (length(missing))
    stop("parameters.pt lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  unexpected = setdiff(names(parameters), c(names(own), paste0("env_scaling.", rep(env_vars, each = 2), c(".center", ".scale"))))
  if (length(unexpected))
    stop("parameters.pt holds tensors this model structure does not have: ", paste(unexpected, collapse = ", "),
         " (another variant or FINN version?)", call. = FALSE)
  torch::with_no_grad({
    for (nm in names(own)) {
      src = parameters[[nm]]
      if (!identical(as.integer(dim(src)), as.integer(dim(own[[nm]]))))
        stop("parameter '", nm, "' has dims ", paste(dim(src), collapse = "x"),
             ", the model expects ", paste(dim(own[[nm]]), collapse = "x"), call. = FALSE)
      own[[nm]]$copy_(src$to(device = own[[nm]]$device, dtype = own[[nm]]$dtype))
    }
  })
  model$eval()
  model
}
