#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as the calibration fit script, followed by copying the fitted
#' parameters into the new object. Rebuilding (rather than `torch_load()` of a
#' whole model object) means the model always runs with the methods of the
#' installed FINN.
#'
#' Variant `"process"`: FINN process model for the species of [species_table()]
#' on an annual timestep with WorldClim 2.1 climate normals as covariates,
#' `tavg_c` (annual mean temperature, degC, bio1) and `prec_mm` (annual
#' precipitation, mm, bio12), each process (growth, mortality, saturating
#' regeneration) with a two-layer SELU network (20, 20) on its covariates and
#' per-species mechanistic parameters; competition without covariates.
#' Calibrated on one patch of 0.1257 ha per site (the 20 m NBI circle) to the
#' tree-by-tree remeasurement of the Dutch national forest inventory (NBI),
#' where the mortality target is total stem loss including felled trees.
#'
#' The released parameter file describes its own structure, so one definition
#' serves every release: the species count is the first dimension of
#' `par_growth_unconstrained` (checked against [species_table()]), the
#' mortality covariates follow from the input width of `nn_mortality.0.weight`
#' (3 = `tavg_c + prec_mm`; 4 = `tavg_c + prec_mm + harvest`, the v0.2 management
#' dial), and optional scalar entries `settings.recruits_dbh`,
#' `settings.reg_floor`, `settings.mortality_lower`, `settings.mortality_upper`
#' and `settings.mortality_floor` carry the constants of that fit (defaults:
#' the v0.1 values 6.4 cm, 1e-3, -10, 10, and FINN's own floor). The env
#' standardisation constants `env_scaling.<variable>.center` / `.scale` learned
#' at fitting are required: FINN z-scales the covariates internally and reuses
#' them at prediction.
#'
#' @param variant Name of the model variant (see [model_info()]).
#' @param parameters Named list of parameter tensors, as saved by the
#'   calibration script (`parameters.pt`).
#' @param version Structure of the parameter file: `"v0.1"` (the structure
#'   written out below, the released v0.1.0 and v0.2.0 files), `"v0.2-recipe"`
#'   (the zoo v0.2 recipe, [build_model_v02()]) or `"auto"` (the default: the
#'   recipe when the file's `settings$recipe` tag says so, v0.1 otherwise).
#' @return A fitted `FINN` model.
#' @keywords internal
build_model = function(variant = "process", parameters, version = c("auto", "v0.1", "v0.2-recipe")) {
  version = match.arg(version)
  ## the structure switch: a zoo v0.2 recipe file describes itself through its
  ## settings$recipe tag (R/recipe_v02.R); everything else is the v0.1 structure below
  if (version == "auto") version = if (is_recipe_v02(parameters)) "v0.2-recipe" else "v0.1"
  if (version == "v0.2-recipe") return(build_model_v02(parameters, variant = if (missing(variant)) NULL else variant))
  if (!identical(variant, "process"))
    stop("Unknown variant '", variant, "'; this package has the variant 'process'.", call. = FALSE)
  scalar = function(nm, default) if (is.null(parameters[[nm]])) default else as.numeric(parameters[[nm]])
  n_species = as.integer(dim(parameters[["par_growth_unconstrained"]])[1])
  if (!identical(n_species, nrow(species_table())))
    stop("parameters.pt holds ", n_species, " species, species.csv lists ", nrow(species_table()),
         " (another release of the package?)", call. = FALSE)
  mort_width = as.integer(dim(parameters[["nn_mortality.0.weight"]])[2])   # intercept + covariates of the mortality network
  if (!mort_width %in% c(3L, 4L)) stop("unexpected mortality network input width ", mort_width, call. = FALSE)
  env_vars  = c("tavg_c", "prec_mm")
  mort_vars = if (mort_width == 4L) c("tavg_c", "prec_mm", "harvest") else env_vars
  hidden    = c(20L, 20L)
  patch_size = pi * 20^2 / 1e4                            # 0.1257 ha, the 20 m NBI circle
  recruits_dbh = scalar("settings.recruits_dbh", 6.4)     # median dbh of the observed recruits (tally threshold 5 cm)
  reg_floor    = scalar("settings.reg_floor", 1e-3)       # stems/ha/yr: minimum recruitment (immigration)
  mort_lower   = scalar("settings.mortality_lower", -10); mort_upper = scalar("settings.mortality_upper", 10)
  mort_floor   = if (is.null(parameters[["settings.mortality_floor"]])) NULL else as.numeric(parameters[["settings.mortality_floor"]])
  mort_args = list(func = FINN::mortality, hidden = hidden, dropout = 0.1, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                   lower = rep(mort_lower, 3), upper = rep(mort_upper, 3))
  if (!is.null(mort_floor)) {
    if (!"mortality_floor" %in% names(formals(FINN::createProcess)))
      stop("this release was fitted with mortality_floor = ", mort_floor, ", which the installed FINN does not support ",
           "(it needs a build with FINNverse/FINN#51).", call. = FALSE)
    mort_args$mortality_floor = mort_floor; mort_args$mortality_draw = "poisson"
  }
  ## createProcess() evaluates its formula in the calling frame's parent, so the
  ## formulas are written out literally (a local variable would not be found)
  mortality_process = if (mort_width == 4L) do.call(FINN::createProcess, c(list(formula = ~ tavg_c + prec_mm + harvest), mort_args))
                      else                   do.call(FINN::createProcess, c(list(formula = ~ tavg_c + prec_mm), mort_args))
  model = FINN::finn(
    N_species    = n_species,
    recruits_dbh = recruits_dbh,
    reg_floor    = reg_floor,
    competition_process  = FINN::createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = FINN::createProcess(~ tavg_c + prec_mm, FINN::growth,
                                               hidden = hidden, dropout = 0.1,
                                               optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = mortality_process,
    regeneration_process = FINN::createProcess(~ tavg_c + prec_mm, FINN::regeneration_saturation,
                                               hidden = hidden, dropout = 0.1,
                                               optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                               custom_parameters = list(reg_logK = rep(log(50), n_species)))
  )

  ## the env standardisation learned at fitting (natural units in, z-scores inside)
  all_vars = unique(c(env_vars, mort_vars))
  for (nm in paste0("env_scaling.", rep(all_vars, each = 2), c(".center", ".scale")))
    if (is.null(parameters[[nm]])) stop("parameters.pt lacks ", nm, call. = FALSE)
  model$env_scaling = data.frame(
    variable = all_vars,
    center   = vapply(all_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".center")]]), numeric(1)),
    scale    = vapply(all_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".scale")]]), numeric(1)),
    stringsAsFactors = FALSE)
  ## FINN builds the env networks lazily on the first simulation (their input
  ## width comes from the env table), so one minimal bare-ground step creates
  ## them before the fitted weights are copied in
  dummy_env = as.data.frame(c(list(siteID = 1L, year = 1L), stats::setNames(as.list(model$env_scaling$center), all_vars)))
  invisible(utils::capture.output(suppressMessages(
    FINN::simulateForest(model, env = dummy_env, init_cohort = NULL, patches = 1L,
                         patch_size = patch_size, device = "cpu"))))

  ## the fitted tensors, matched by NAME; anything missing or unexpected stops here
  own = model$parameters
  missing = setdiff(names(own), names(parameters))
  if (length(missing))
    stop("parameters.pt lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  known_extra = c(paste0("env_scaling.", rep(all_vars, each = 2), c(".center", ".scale")),
                  paste0("settings.", c("recruits_dbh", "reg_floor", "mortality_lower", "mortality_upper", "mortality_floor")))
  unexpected = setdiff(names(parameters), c(names(own), known_extra))
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
