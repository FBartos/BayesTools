#' @title Compute marginal likelihood of a 'JAGS' model
#'
#' @description A wrapper around
#' \link[bridgesampling]{bridge_sampler} that automatically
#' computes likelihood part dependent on the prior distribution
#' and prepares parameter samples. \code{log_posterior} must
#' specify a function that takes a named parameter list, the data, and any
#' additional arguments, and returns the log likelihood of the model part.
#'
#' @param fit model fitted with either \link[runjags]{runjags} posterior
#' samples obtained with \link[rjags]{rjags-package}
#' @param log_posterior function that takes a named list of samples, the data,
#' and additional list of parameters passed as \code{...} as input and
#' returns the log of the unnormalized posterior density of the model part
#' @param data list containing data to fit the model (not including data for the formulas)
#' @param prior_list named list of prior distributions
#' (names correspond to the parameter names) of parameters not specified within the
#' \code{formula_list}. For \code{BayesTools_fit} objects, stored non-formula
#' priors are used when \code{prior_list = NULL}; if fitted formula priors are
#' supplied here, they are ignored with a warning in favor of the stored formula
#' metadata.
#' @param formula_list named list of formulas to be added to the model
#' (names correspond to the parameter name created by each of the formula). For
#' \code{BayesTools_fit} objects with stored formula-design metadata, formula
#' inputs can usually be omitted; if supplied, they are rebuilt only to check
#' consistency with the fitted design.
#' @param formula_data_list named list of data frames containing data for each formula
#' (names of the lists correspond to the parameter name created by each of the formula)
#' @param formula_prior_list named list of named lists of prior distributions
#' (names of the lists correspond to the parameter name created by each of the formula and
#' the names of the prior distribution correspond to the parameter names) of parameters specified
#' within the \code{formula}
#' @param formula_scale_list named list of named lists for standardizing continuous predictors
#' (names of the lists correspond to the parameter name created by each of the formula).
#' Each entry should be a named list where continuous predictors with \code{TRUE} values will
#' be standardized. Defaults to stored fit metadata when available and to
#' \code{NULL} (no standardization) otherwise.
#' @param add_parameters character vector of additional monitored posterior
#' parameter names to include in the bridge-sampling state and in the
#' `parameters` object passed to `log_posterior`. These are for JAGS nodes that
#' are not owned by `prior_list`, including special likelihood parameters and
#' row-shaped external sources. The `parameters` object is a named superset of
#' prior-owned and additional values; user code should index it by name, for
#' example `parameters[["tau"]]`. Parameters already covered by `prior_list`
#' must not be listed here. For Dirichlet priors generated through BayesTools,
#' add the monitored auxiliary `prior_par_eta_*` coordinates when needed, not
#' the normalized simplex coordinates owned by the prior.
#' @param add_bounds list with two named numeric vectors, \code{"lb"} and
#' \code{"ub"}, containing lower and upper bounds for every
#' \code{add_parameters} entry. Must be supplied whenever
#' \code{add_parameters} is non-empty.
#' @param formula_random_prior_list optional named list of `prior_random()`
#' objects for random effects in `formula_list`. Bridge sampling for formula
#' random effects requires the `prior_random()` interface because the
#' stochastic bridge coordinates are the standardized latent effects and
#' correlation primitives. For \code{BayesTools_fit} objects with stored
#' formula-design metadata, this can be omitted unless formula inputs are being
#' supplied for a consistency check.
#' @param maxiter maximum number of iterations for the
#' \link[bridgesampling]{bridge_sampler}
#' @param silent whether the progress should be printed, defaults to \code{TRUE}
#' @param ... additional argument to the \link[bridgesampling]{bridge_sampler}
#' and \code{log_posterior} function
#'
#' @details Row-shaped external random-effect SD sources, such as
#' `random_sd_source("tau", shape = "row")`, must be reconstructable during
#' bridge sampling. Supply them either as posterior columns named
#' `tau[1]`, ..., `tau[n]` with non-negative lower bounds in `add_bounds`, or
#' as `parameter_source("tau", shape = "row", values = function(parameters,
#' data, n_rows) ...)`. The `values` function is evaluated from the named
#' `parameters` object and row-aligned data; it must return finite,
#' non-negative row values on the support of the model.
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#'
#' # define log posterior for bridge sampling
#' log_posterior <- function(parameters, data){
#'   sum(dnorm(data$x, parameters$mu, 1, log = TRUE))
#' }
#'
#' # get marginal likelihoods
#' marglik <- JAGS_bridgesampling(fit, log_posterior, data, priors_list)
#' }
#' @return \code{JAGS_bridgesampling} returns an object of class 'bridge'.
#'
#' @export
JAGS_bridgesampling <- function(fit, log_posterior, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL, formula_scale_list = NULL,
                                add_parameters = NULL, add_bounds = NULL,
                                formula_random_prior_list = NULL,
                                maxiter = 10000, silent = TRUE, ...){

  ### check input
  check_bool(silent, "silent")
  check_int(maxiter, "maxiter", lower = 1)

  formula_context <- .bt_JAGS_bridge_formula_context(
    fit = fit,
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list
  )
  formula_design_list <- formula_context$formula_design_list
  formula_list <- formula_context$formula_list
  formula_data_list <- formula_context$formula_data_list
  formula_prior_list <- formula_context$formula_prior_list

  if(is.null(prior_list)){
    prior_list <- .bt_JAGS_bridge_prior_list_from_fit(
      fit = fit,
      formula_design_list = formula_design_list
    )
  }else{
    prior_list <- .bt_JAGS_bridge_non_formula_prior_list(
      prior_list = prior_list,
      formula_design_list = formula_design_list,
      warn = TRUE
    )
  }
  if(is.null(prior_list)){
    prior_list <- list()
  }

  # extract the posterior distribution
  posterior <- .fit_to_posterior(fit)

  if(length(formula_prior_list) > 0L){
    all_prior_list <- c(prior_list, do.call(c, unname(formula_prior_list)))
  }else{
    all_prior_list <- prior_list
  }

  if(length(all_prior_list) > 0L && any(sapply(all_prior_list, is.prior.discrete)))
    stop("Discrete or spike and slab priors are not supported with bridgesampling.")

  ### extract relevant variables and upper and lower bound
  random_bridge_parameters <- .bt_JAGS_formula_random_bridge_parameters(formula_design_list)
  if(length(random_bridge_parameters$parameters) > 0L){
    bridge_add <- .bt_JAGS_bridge_merge_add_parameters(
      add_parameters = add_parameters,
      add_bounds = add_bounds,
      bridge_parameters = random_bridge_parameters$parameters,
      bridge_bounds = random_bridge_parameters$bounds
    )
    add_parameters <- bridge_add$add_parameters
    add_bounds <- bridge_add$add_bounds
  }
  .bt_JAGS_bridge_validate_add_parameters_not_formula(
    add_parameters = add_parameters,
    formula_design_list = formula_design_list,
    formula_prior_list = formula_prior_list
  )
  .bt_JAGS_bridge_check_random_posterior(posterior, random_bridge_parameters$parameters)
  bridgesampling_posterior <- JAGS_bridgesampling_posterior(posterior = posterior, prior_list = all_prior_list, add_parameters = add_parameters, add_bounds = add_bounds)
  bridgesampling_posterior <- .bt_JAGS_bridge_apply_random_scalar_rho_bounds(
    bridgesampling_posterior = bridgesampling_posterior,
    formula_design_list = formula_design_list
  )
  .bt_JAGS_bridge_check_row_indexed_external_sd_sources(
    formula_design_list = formula_design_list,
    bridgesampling_posterior = bridgesampling_posterior
  )
  if(ncol(bridgesampling_posterior) == 0)
    stop("Bridge sampling cannot proceed without any estimated parameter")

  bridge_prior_evaluator <- .bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)
  bridge_formula_prior_evaluator <- .bt_JAGS_bridge_compile_formula_prior_evaluator(formula_prior_list)
  bridge_formula_random_prior_evaluator <- .bt_JAGS_bridge_compile_formula_random_prior_evaluator(formula_design_list)
  bridge_formula_parameter_evaluator <- .bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = data
  )


  ### define the marglik function
  full_log_posterior <- function(samples.row, data,
                                 bridge_prior_evaluator,
                                 bridge_formula_prior_evaluator,
                                 bridge_formula_random_prior_evaluator,
                                 bridge_formula_parameter_evaluator,
                                 add_parameters,
                                 ...){

    samples.row <- .bt_JAGS_bridge_cache_posterior_row(
      samples.row,
      bridge_formula_random_prior_evaluator$uses_posterior_row
    )

    # Check prior support before reconstructing formula parameters. Bridge
    # proposals can hit bounded-prior edges where random-effect correlations
    # are intentionally invalid and should contribute zero density.
    marglik <- bridge_prior_evaluator$log_prior(samples.row)
    marglik <- marglik + bridge_formula_prior_evaluator$log_prior(samples.row)
    marglik <- marglik + bridge_formula_random_prior_evaluator$log_prior(samples.row)
    if(is.na(marglik)){
      return(-Inf)
    }
    if(!is.finite(marglik)){
      return(marglik)
    }

    # prepare object for holding the parameters, later accessible to the user specified 'log_posterior'
    parameters <- tryCatch({
      parameters <- list()
      parameters <- c(parameters, bridge_prior_evaluator$parameters(samples.row))
      parameters <- c(parameters, bridge_formula_parameter_evaluator$parameters(samples.row, parameters))
      if(length(add_parameters) > 0){
        parameters <- c(parameters, samples.row[add_parameters])
      }
      parameters
    }, BayesTools_marglik_out_of_support = function(e)e)
    if(inherits(parameters, "BayesTools_marglik_out_of_support")){
      return(-Inf)
    }

    marglik   <- marglik + log_posterior(parameters = parameters, data = data, ...)

    return(marglik)
  }


  ### perform bridgesampling
  marglik <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
      samples            = bridgesampling_posterior,
      data               = data,
      log_posterior      = full_log_posterior,
      bridge_prior_evaluator = bridge_prior_evaluator,
      bridge_formula_prior_evaluator = bridge_formula_prior_evaluator,
      bridge_formula_random_prior_evaluator = bridge_formula_random_prior_evaluator,
      bridge_formula_parameter_evaluator = bridge_formula_parameter_evaluator,
      lb                 = attr(bridgesampling_posterior, "lb"),
      ub                 = attr(bridgesampling_posterior, "ub"),
      silent             = silent,
      maxiter            = maxiter,
      add_parameters     = add_parameters,
      ...
    )), error = function(e)e)

  # add a warning attribute and call the warning if not silent
  if(!inherits(marglik, "error") && marglik[["niter"]] > maxiter){
    attr(marglik, "warning") <- "Marginal likelihood could not be estimated within the maximum number of itetations and might be more variable than usual."
    if(!silent)
      warning(attr(marglik, "warning"), immediate. = TRUE)
  }

  return(marglik)
}

.bt_JAGS_bridge_formula_input_supplied <- function(formula_list,
                                                   formula_data_list,
                                                   formula_prior_list,
                                                   formula_scale_list,
                                                   formula_random_prior_list){

  !is.null(formula_list) ||
    !is.null(formula_data_list) ||
    !is.null(formula_prior_list) ||
    !is.null(formula_scale_list) ||
    !is.null(formula_random_prior_list)
}

.bt_JAGS_bridge_formula_context <- function(fit,
                                            formula_list,
                                            formula_data_list,
                                            formula_prior_list,
                                            formula_scale_list,
                                            formula_random_prior_list){

  formula_input_supplied <- .bt_JAGS_bridge_formula_input_supplied(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list
  )
  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(
    attr(fit, "formula_design")
  )
  has_fitted_formula_design <- length(fitted_formula_design) > 0L

  if(!formula_input_supplied){
    if(has_fitted_formula_design){
      return(.bt_JAGS_bridge_formula_context_from_design(
        fitted_formula_design
      ))
    }
    return(.bt_JAGS_bridge_empty_formula_context())
  }

  rebuildability <- tryCatch(.bt_JAGS_bridge_formula_inputs_rebuildability(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list,
    has_fitted_formula_design = has_fitted_formula_design
  ), error = function(e){
    if(has_fitted_formula_design){
      return(list(
        rebuildable = FALSE,
        detail = conditionMessage(e)
      ))
    }
    stop(conditionMessage(e), call. = FALSE)
  })
  if(!isTRUE(rebuildability$rebuildable)){
    if(has_fitted_formula_design){
      .bt_JAGS_bridge_stop_formula_mismatch(rebuildability$detail)
    }
    stop(rebuildability$detail, call. = FALSE)
  }

  rebuilt_context <- tryCatch(.bt_JAGS_bridge_rebuild_formula_context(
    fit = fit,
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list
  ), error = function(e){
    e
  })
  if(inherits(rebuilt_context, "error")){
    if(has_fitted_formula_design){
      .bt_JAGS_bridge_stop_formula_mismatch(
        paste0("formula inputs could not be rebuilt: ", conditionMessage(rebuilt_context))
      )
    }
    stop(conditionMessage(rebuilt_context), call. = FALSE)
  }
  if(has_fitted_formula_design){
    mismatches <- .bt_JAGS_bridge_formula_design_mismatches(
      fitted_formula_design = fitted_formula_design,
      rebuilt_formula_design = rebuilt_context$formula_design_list
    )
    if(length(mismatches) > 0L){
      .bt_JAGS_bridge_stop_formula_mismatch(mismatches)
    }
    fitted_formula_design <- .bt_JAGS_bridge_formula_design_update_source_values(
      fitted_formula_design = fitted_formula_design,
      rebuilt_formula_design = rebuilt_context$formula_design_list
    )
    return(.bt_JAGS_bridge_formula_context_from_design(
      fitted_formula_design,
      formula_data_list = rebuilt_context$formula_data_list
    ))
  }

  rebuilt_context
}

.bt_JAGS_bridge_formula_inputs_rebuildability <- function(formula_list,
                                                          formula_data_list,
                                                          formula_prior_list,
                                                          formula_scale_list,
                                                          formula_random_prior_list,
                                                          has_fitted_formula_design){

  if(is.null(formula_list) || is.null(formula_data_list) || is.null(formula_prior_list)){
    detail <- if(has_fitted_formula_design){
      "formula-related inputs are incomplete"
    }else{
      "When supplying formula-related inputs to JAGS_bridgesampling(), provide 'formula_list', 'formula_data_list', and 'formula_prior_list'."
    }
    return(list(rebuildable = FALSE, detail = detail))
  }

  check_list(formula_list, "formula_list", allow_NULL = FALSE)
  if(has_fitted_formula_design){
    check_list(formula_data_list, "formula_data_list", allow_NULL = FALSE)
    check_list(formula_prior_list, "formula_prior_list", allow_NULL = FALSE)
    if(!.bt_JAGS_bridge_formula_input_names_match(
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      formula_scale_list = formula_scale_list,
      formula_random_prior_list = formula_random_prior_list
    )){
      return(list(
        rebuildable = FALSE,
        detail = "formula input names are incomplete or inconsistent"
      ))
    }
    scale_check_names <- NULL
    random_check_names <- NULL
  }else{
    check_list(formula_data_list, "formula_data_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = FALSE)
    check_list(formula_prior_list, "formula_prior_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = FALSE)
    scale_check_names <- names(formula_list)
    random_check_names <- names(formula_list)
  }

  check_list(formula_scale_list, "formula_scale_list", check_names = scale_check_names, allow_other = FALSE, all_objects = FALSE, allow_NULL = TRUE)
  check_list(formula_random_prior_list, "formula_random_prior_list", check_names = random_check_names, allow_other = FALSE, all_objects = FALSE, allow_NULL = TRUE)
  if(!is.null(formula_random_prior_list)){
    for(parameter in names(formula_random_prior_list)){
      .bt_check_prior_random(formula_random_prior_list[[parameter]])
    }
  }

  for(parameter in names(formula_list)){
    if(.has_random_effects(formula_list[[parameter]]) &&
       (is.null(formula_random_prior_list) || is.null(formula_random_prior_list[[parameter]]))){
      detail <- paste0(
        "JAGS_bridgesampling() requires 'formula_random_prior_list' with a prior_random() object for formula random effects in parameter '",
        parameter,
        "'. Legacy random-effect priors in 'formula_prior_list' are not bridge-sampling ready."
      )
      if(has_fitted_formula_design){
        detail <- paste0(
          "formula random effects for parameter '",
          parameter,
          "' cannot be rebuilt without a matching 'formula_random_prior_list' entry"
        )
      }
      return(list(rebuildable = FALSE, detail = detail))
    }
  }

  list(rebuildable = TRUE, detail = NULL)
}

.bt_JAGS_bridge_rebuild_formula_context <- function(fit,
                                                    formula_list,
                                                    formula_data_list,
                                                    formula_prior_list,
                                                    formula_scale_list,
                                                    formula_random_prior_list){

  if(is.null(formula_scale_list)){
    formula_scale_list <- .JAGS_formula_scale_list_from_fit(
      fit,
      names(formula_list)
    )
  }

  formula_output <- list()
  for(parameter in names(formula_list)){
    formula_output[[parameter]] <- JAGS_formula(
      formula       = formula_list[[parameter]],
      parameter     = parameter,
      data          = formula_data_list[[parameter]],
      prior_list    = formula_prior_list[[parameter]],
      formula_scale = if(!is.null(formula_scale_list)) formula_scale_list[[parameter]] else NULL,
      prior_random  = if(!is.null(formula_random_prior_list)) formula_random_prior_list[[parameter]] else NULL
    )
  }
  formula_data_output <- lapply(names(formula_output), function(parameter){
    .bt_JAGS_bridge_formula_data_with_raw(
      generated_data = formula_output[[parameter]][["data"]],
      raw_data = formula_data_list[[parameter]]
    )
  })
  names(formula_data_output) <- names(formula_output)

  list(
    formula_design_list = lapply(formula_output, function(output) output[["formula_design"]]),
    formula_list        = lapply(formula_output, function(output) output[["formula"]]),
    formula_data_list   = formula_data_output,
    formula_prior_list  = lapply(formula_output, function(output) output[["prior_list"]])
  )
}

.bt_JAGS_bridge_formula_data_with_raw <- function(generated_data, raw_data){

  out <- list()
  out <- .bt_JAGS_marglik_merge_source_data(out, raw_data)
  out <- .bt_JAGS_marglik_merge_source_data(out, generated_data)

  out
}

.bt_JAGS_bridge_empty_formula_context <- function(){

  list(
    formula_design_list = NULL,
    formula_list        = NULL,
    formula_data_list   = NULL,
    formula_prior_list  = list()
  )
}

.bt_JAGS_bridge_formula_context_from_design <- function(formula_design_list,
                                                       formula_data_list = NULL){

  list(
    formula_design_list = formula_design_list,
    formula_list        = .bt_JAGS_bridge_formula_list_from_design(formula_design_list),
    formula_data_list   = formula_data_list,
    formula_prior_list  = .bt_JAGS_bridge_formula_prior_list_from_design(formula_design_list)
  )
}

.bt_JAGS_bridge_prior_list_from_fit <- function(fit, formula_design_list){

  fit_prior_list <- attr(fit, "prior_list")
  .bt_JAGS_bridge_non_formula_prior_list(
    prior_list = fit_prior_list,
    formula_design_list = formula_design_list,
    warn = FALSE
  )
}

.bt_JAGS_bridge_non_formula_prior_list <- function(prior_list,
                                                   formula_design_list,
                                                   warn = FALSE){

  if(is.null(prior_list)){
    return(list())
  }
  if(!is.list(prior_list)){
    return(prior_list)
  }

  formula_prior_names <- .bt_JAGS_bridge_formula_prior_names(formula_design_list)
  if(length(formula_prior_names) == 0L){
    return(prior_list)
  }

  overlap <- intersect(names(prior_list), formula_prior_names)
  if(length(overlap) > 0L){
    if(isTRUE(warn)){
      warning(
        "JAGS_bridgesampling() received formula priors in 'prior_list'; using fitted formula metadata for these priors and ignoring duplicate 'prior_list' entries: ",
        paste(utils::head(overlap, 8L), collapse = ", "),
        if(length(overlap) > 8L) ", ..." else "",
        ".",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    prior_list <- prior_list[setdiff(names(prior_list), overlap)]
  }

  prior_list
}

.bt_JAGS_bridge_formula_prior_names <- function(formula_design_list){

  formula_prior_list <- .bt_JAGS_bridge_formula_prior_list_from_design(
    formula_design_list
  )
  unique(unlist(lapply(formula_prior_list, names), use.names = FALSE))
}

.bt_JAGS_bridge_formula_list_from_design <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(NULL)
  }

  formula_list <- lapply(formula_design_list, function(design){
    if(inherits(design, "BayesTools_formula_design")){
      return(design$formula)
    }
    NULL
  })
  formula_list[!vapply(formula_list, is.null, logical(1))]
}

.bt_JAGS_bridge_formula_prior_list_from_design <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list())
  }

  prior_list <- lapply(formula_design_list, function(design){
    if(inherits(design, "BayesTools_formula_design") &&
       is.list(design$prior_list)){
      return(design$prior_list)
    }
    list()
  })
  prior_list[!vapply(prior_list, function(x) length(x) == 0L, logical(1))]
}

.bt_JAGS_bridge_formula_input_names_match <- function(formula_list,
                                                      formula_data_list,
                                                      formula_prior_list,
                                                      formula_scale_list = NULL,
                                                      formula_random_prior_list = NULL){

  formula_names <- names(formula_list)
  if(is.null(formula_names) || anyNA(formula_names) || any(!nzchar(formula_names))){
    return(FALSE)
  }

  if(!setequal(formula_names, names(formula_data_list)) ||
     !setequal(formula_names, names(formula_prior_list))){
    return(FALSE)
  }
  if(!is.null(formula_scale_list) &&
     any(!names(formula_scale_list) %in% formula_names)){
    return(FALSE)
  }
  if(!is.null(formula_random_prior_list) &&
     any(!names(formula_random_prior_list) %in% formula_names)){
    return(FALSE)
  }

  TRUE
}

.bt_JAGS_bridge_stop_formula_mismatch <- function(mismatches){

  if(length(mismatches) == 0L){
    mismatches <- "unknown formula-design mismatch"
  }

  stop(
    "JAGS_bridgesampling() supplied formula-related inputs do not fully match ",
    "the fitted formula design. First mismatch: ",
    mismatches[[1L]],
    if(length(mismatches) > 1L) paste0(" Additional mismatches: ", length(mismatches) - 1L, ".") else "",
    call. = FALSE
  )
}

.bt_JAGS_bridge_formula_design_mismatches <- function(fitted_formula_design,
                                                      rebuilt_formula_design){

  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(fitted_formula_design)
  rebuilt_formula_design <- .bt_JAGS_bridge_formula_design_list(rebuilt_formula_design)
  mismatches <- character()

  if(!setequal(names(fitted_formula_design), names(rebuilt_formula_design))){
    mismatches <- c(
      mismatches,
      paste0(
        "formula parameter names differ; fitted: ",
        paste(names(fitted_formula_design), collapse = ", "),
        "; supplied: ",
        paste(names(rebuilt_formula_design), collapse = ", ")
      )
    )
  }

  for(parameter in intersect(names(fitted_formula_design), names(rebuilt_formula_design))){
    fitted <- fitted_formula_design[[parameter]]
    rebuilt <- rebuilt_formula_design[[parameter]]
    mismatches <- c(
      mismatches,
      .bt_JAGS_bridge_formula_design_parameter_mismatches(
        parameter = parameter,
        fitted = fitted,
        rebuilt = rebuilt
      )
    )
  }

  mismatches
}

.bt_JAGS_bridge_formula_design_parameter_mismatches <- function(parameter,
                                                               fitted,
                                                               rebuilt){

  mismatches <- character()
  if(!inherits(fitted, "BayesTools_formula_design") ||
     !inherits(rebuilt, "BayesTools_formula_design")){
    return(paste0("formula design metadata for parameter '", parameter, "' are incomplete"))
  }

  if(!.bt_JAGS_bridge_formulas_equal(fitted$formula, rebuilt$formula)){
    mismatches <- c(mismatches, paste0("formula differs for parameter '", parameter, "'"))
  }
  if(!identical(isTRUE(fitted$log_intercept), isTRUE(rebuilt$log_intercept))){
    mismatches <- c(mismatches, paste0("log-intercept metadata differ for parameter '", parameter, "'"))
  }
  if(!identical(dim(fitted$model_matrix), dim(rebuilt$model_matrix)) ||
     !identical(colnames(fitted$model_matrix), colnames(rebuilt$model_matrix))){
    mismatches <- c(mismatches, paste0("fixed-effect model matrix shape or columns differ for parameter '", parameter, "'"))
  }else if(!isTRUE(all.equal(
    unname(fitted$model_matrix),
    unname(rebuilt$model_matrix),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))){
    mismatches <- c(mismatches, paste0("fixed-effect model matrix values differ for parameter '", parameter, "'"))
  }
  if(!identical(fitted$assign, rebuilt$assign)){
    mismatches <- c(mismatches, paste0("fixed-effect model term assignments differ for parameter '", parameter, "'"))
  }
  if(!.bt_JAGS_bridge_metadata_equal(fitted$contrasts, rebuilt$contrasts) ||
     !.bt_JAGS_bridge_metadata_equal(fitted$xlevels, rebuilt$xlevels)){
    mismatches <- c(mismatches, paste0("factor contrasts or levels differ for parameter '", parameter, "'"))
  }
  if(!.bt_JAGS_bridge_metadata_equal(fitted$formula_scale, rebuilt$formula_scale)){
    mismatches <- c(mismatches, paste0("formula scaling metadata differ for parameter '", parameter, "'"))
  }
  if(!identical(names(fitted$prior_list), names(rebuilt$prior_list))){
    mismatches <- c(mismatches, paste0("formula prior names differ for parameter '", parameter, "'"))
  }else{
    for(prior_name in names(fitted$prior_list)){
      if(!.bt_JAGS_bridge_metadata_equal(fitted$prior_list[[prior_name]], rebuilt$prior_list[[prior_name]])){
        mismatches <- c(mismatches, paste0("formula prior metadata differ for parameter '", parameter, "', prior '", prior_name, "'"))
        break
      }
    }
  }

  random_mismatch <- tryCatch({
    .bt_JAGS_bridge_validate_formula_random_design(
      parameter = parameter,
      fitted = fitted,
      rebuilt = rebuilt
    )
    NULL
  }, error = function(e) conditionMessage(e))
  if(!is.null(random_mismatch)){
    mismatches <- c(mismatches, random_mismatch)
  }

  mismatches
}

.bt_JAGS_bridge_formula_design_update_source_values <- function(fitted_formula_design,
                                                                rebuilt_formula_design){

  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(fitted_formula_design)
  rebuilt_formula_design <- .bt_JAGS_bridge_formula_design_list(rebuilt_formula_design)
  for(parameter in intersect(names(fitted_formula_design), names(rebuilt_formula_design))){
    fitted <- fitted_formula_design[[parameter]]
    rebuilt <- rebuilt_formula_design[[parameter]]
    if(!.bt_formula_design_has_random_effects(fitted) ||
       !.bt_formula_design_has_random_effects(rebuilt)){
      next
    }
    fitted_blocks <- .bt_JAGS_bridge_random_block_names(fitted$random_effects)
    rebuilt_blocks <- .bt_JAGS_bridge_random_block_names(rebuilt$random_effects)
    for(block in intersect(fitted_blocks, rebuilt_blocks)){
      fitted_i <- match(block, fitted_blocks)
      rebuilt_i <- match(block, rebuilt_blocks)
      fitted$random_effects[[fitted_i]] <- .bt_JAGS_bridge_random_term_update_source_values(
        fitted = fitted$random_effects[[fitted_i]],
        rebuilt = rebuilt$random_effects[[rebuilt_i]]
      )
    }
    fitted_formula_design[[parameter]] <- fitted
  }

  fitted_formula_design
}

.bt_JAGS_bridge_random_term_update_source_values <- function(fitted, rebuilt){

  if(!.bt_random_effect_has_row_indexed_external_sd(fitted) ||
     !.bt_random_effect_has_row_indexed_external_sd(rebuilt)){
    return(fitted)
  }

  fitted_source <- .bt_random_effect_row_indexed_source(fitted)
  rebuilt_source <- .bt_random_effect_row_indexed_source(rebuilt)
  if(!.bt_JAGS_bridge_parameter_sources_same(fitted_source, rebuilt_source)){
    return(fitted)
  }

  values <- .bt_parameter_source_values_function(rebuilt_source)
  if(is.null(values)){
    return(fitted)
  }

  if(!is.null(fitted$sd_binding)){
    fitted$sd_binding$source <- .bt_JAGS_bridge_parameter_source_set_values(
      source = fitted$sd_binding$source,
      values = values
    )
    if(length(fitted$sd_binding$allocations) > 0L){
      for(allocation_i in seq_along(fitted$sd_binding$allocations)){
        if(!is.null(fitted$sd_binding$allocations[[allocation_i]]$source)){
          fitted$sd_binding$allocations[[allocation_i]]$source <- .bt_JAGS_bridge_parameter_source_set_values(
            source = fitted$sd_binding$allocations[[allocation_i]]$source,
            values = values
          )
        }
      }
    }
  }

  fitted
}

.bt_JAGS_bridge_parameter_sources_same <- function(x, y){

  identical(x$name, y$name) &&
    identical(x$shape, y$shape)
}

.bt_JAGS_bridge_parameter_source_set_values <- function(source, values){

  if(is.null(source)){
    return(NULL)
  }
  .bt_check_parameter_source_values_function(values, "values")

  if(inherits(source, "random_sd_source")){
    .bt_check_random_sd_source(source)
    source$source$values <- values
    return(source)
  }
  .bt_check_parameter_source(source)
  source$values <- values

  source
}

.bt_JAGS_bridge_formulas_equal <- function(x, y){

  isTRUE(all.equal(x, y, check.environment = FALSE))
}

.bt_JAGS_bridge_check_row_indexed_external_sd_sources <- function(formula_design_list,
                                                                  bridgesampling_posterior){

  if(length(formula_design_list) == 0L){
    return(invisible(TRUE))
  }
  if(!is.matrix(bridgesampling_posterior)){
    stop("'bridgesampling_posterior' must be a matrix.", call. = FALSE)
  }

  bridge_parameters <- colnames(bridgesampling_posterior)
  bridge_lb <- attr(bridgesampling_posterior, "lb")
  bridge_ub <- attr(bridgesampling_posterior, "ub")

  for(parameter in names(formula_design_list)){
    design <- formula_design_list[[parameter]]
    if(!.bt_formula_design_has_random_effects(design)){
      next
    }
    for(random_term in design$random_effects){
      if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
        source <- .bt_random_effect_row_indexed_source(random_term)
        source_names <- .bt_parameter_source_row_names(
          source$source,
          nrow(random_term$model_matrix)
        )
        if(.bt_parameter_source_has_values(source)){
          present <- intersect(source_names, bridge_parameters)
          if(length(present) > 0L){
            stop(
              "JAGS_bridgesampling() row-indexed external SD source '",
              .bt_random_effect_external_sd_source_label(random_term),
              "' for random-effect block '",
              random_term$block_name,
              "' provides a parameter_source() values function and also appears ",
              "as bridge parameter column(s). Use one row-source reconstruction ",
              "path only. Conflicting source column(s): ",
              paste0("'", present[seq_len(min(3L, length(present)))], "'",
                     collapse = ", "),
              if(length(present) > 3L) ", ..." else "",
              ".",
              call. = FALSE
            )
          }
          next
        }
        if(all(source_names %in% bridge_parameters)){
          .bt_JAGS_bridge_check_row_indexed_external_sd_bounds(
            random_term = random_term,
            source_names = source_names,
            bridge_lb = bridge_lb,
            bridge_ub = bridge_ub
          )
          next
        }
        missing <- setdiff(source_names, bridge_parameters)
        if(length(missing) > 0L){
          stop(
            "JAGS_bridgesampling() cannot reconstruct row-indexed external SD source '",
            .bt_random_effect_external_sd_source_label(random_term),
            "' for random-effect block '",
            random_term$block_name,
            "'. Provide a parameter_source(..., values = function(parameters, data, n_rows) ...) ",
            "or include the row-wise source column(s) as bridge parameters with appropriate priors or add_parameters/add_bounds. ",
            "Missing source column(s): ",
            paste0("'", missing[seq_len(min(3L, length(missing)))], "'", collapse = ", "),
            if(length(missing) > 3L) ", ..." else "",
            ".",
            call. = FALSE
          )
        }
      }
    }
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_check_row_indexed_external_sd_bounds <- function(random_term,
                                                                 source_names,
                                                                 bridge_lb,
                                                                 bridge_ub){

  if(is.null(bridge_lb) || is.null(bridge_ub) ||
     !all(source_names %in% names(bridge_lb)) ||
     !all(source_names %in% names(bridge_ub))){
    stop(
      "JAGS_bridgesampling() row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for random-effect block '",
      random_term$block_name,
      "' is missing lower or upper bridge bounds.",
      call. = FALSE
    )
  }

  lower <- bridge_lb[source_names]
  upper <- bridge_ub[source_names]
  bad_bounds <- is.na(lower) | is.na(upper) | lower < 0 | lower >= upper
  if(any(bad_bounds)){
    bad_names <- source_names[bad_bounds]
    stop(
      "JAGS_bridgesampling() row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for random-effect block '",
      random_term$block_name,
      "' is included as bridge parameter(s), but its bridge bounds must ",
      "be valid standard-deviation bounds with non-negative lower bounds. ",
      "Invalid source column(s): ",
      paste0("'", bad_names[seq_len(min(3L, length(bad_names)))], "'",
             collapse = ", "),
      if(length(bad_names) > 3L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_JAGS_formula_random_bridge_parameters <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  parameters <- character()
  lb <- numeric()
  ub <- numeric()

  for(parameter in names(formula_design_list)){
    design <- formula_design_list[[parameter]]
    if(!.bt_formula_design_has_random_effects(design)){
      next
    }
    if(!identical(design$random_effects_interface, "prior_random")){
      stop(
        "JAGS_bridgesampling() supports formula random effects only through the 'prior_random()' interface for parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }
    for(random_term in design$random_effects){
      n_groups <- random_term$n_groups
      n_columns <- random_term$n_columns
      z_names <- as.vector(.bt_random_effect_latent_names(
        random_term = random_term,
        n_groups = n_groups,
        n_columns = n_columns
      ))
      parameters <- c(parameters, z_names)
      lb <- c(lb, stats::setNames(rep(-Inf, length(z_names)), z_names))
      ub <- c(ub, stats::setNames(rep( Inf, length(z_names)), z_names))

      if(identical(.bt_JAGS_bridge_random_term_structure(random_term), "us") &&
         n_columns > 1L){
        u_names <- .bt_random_effect_lkj_primitive_names(
          random_term,
          n_columns,
          context = "Bridge sampling random-effect metadata"
        )
        parameters <- c(parameters, u_names)
        lb <- c(lb, stats::setNames(rep(0, length(u_names)), u_names))
        ub <- c(ub, stats::setNames(rep(1, length(u_names)), u_names))
      }
    }
  }

  keep <- !duplicated(parameters)
  parameters <- parameters[keep]
  lb <- lb[parameters]
  ub <- ub[parameters]

  list(
    parameters = parameters,
    bounds = list(lb = lb, ub = ub)
  )
}

.bt_JAGS_bridge_merge_add_parameters <- function(add_parameters, add_bounds,
                                                 bridge_parameters,
                                                 bridge_bounds){

  if((is.null(add_parameters) || length(add_parameters) == 0L) &&
     !is.null(add_bounds)){
    stop("'add_bounds' requires at least one 'add_parameters' entry.", call. = FALSE)
  }

  if(length(bridge_parameters) == 0L){
    return(list(add_parameters = add_parameters, add_bounds = add_bounds))
  }

  user_parameters <- if(is.null(add_parameters)) character() else add_parameters
  if(length(user_parameters) > 0L){
    user_bounds <- .bt_JAGS_bridge_validate_add_bounds(user_parameters, add_bounds)
    user_lb <- user_bounds$lb
    user_ub <- user_bounds$ub
  }else{
    user_lb <- numeric()
    user_ub <- numeric()
  }

  bridge_lb <- bridge_bounds$lb
  bridge_ub <- bridge_bounds$ub
  overlapping_bridge <- intersect(user_parameters, bridge_parameters)
  if(length(overlapping_bridge) > 0L){
    conflicting_bridge <- overlapping_bridge[
      vapply(overlapping_bridge, function(parameter){
        !isTRUE(all.equal(unname(user_lb[[parameter]]), unname(bridge_lb[[parameter]]))) ||
          !isTRUE(all.equal(unname(user_ub[[parameter]]), unname(bridge_ub[[parameter]])))
      }, logical(1))
    ]
    if(length(conflicting_bridge) > 0L){
      stop(
        "User-supplied bounds conflict with automatically inferred formula random-effect bridge bounds for parameter(s): ",
        paste(conflicting_bridge, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  missing_bridge <- setdiff(bridge_parameters, names(user_lb))
  combined_parameters <- unique(c(user_parameters, bridge_parameters))
  combined_lb <- c(user_lb, bridge_lb[missing_bridge])
  combined_ub <- c(user_ub, bridge_ub[missing_bridge])
  combined_lb <- combined_lb[combined_parameters]
  combined_ub <- combined_ub[combined_parameters]

  list(
    add_parameters = combined_parameters,
    add_bounds = list(lb = combined_lb, ub = combined_ub)
  )
}

.bt_JAGS_bridge_validate_add_bounds <- function(add_parameters, add_bounds){

  if(!is.character(add_parameters)){
    stop("'add_parameters' must be a character vector.", call. = FALSE)
  }
  if(anyDuplicated(add_parameters)){
    stop("'add_parameters' must be unique.", call. = FALSE)
  }
  if(!is.list(add_bounds)){
    stop("'add_bounds' must be a list.", call. = FALSE)
  }
  if(length(add_bounds) != 2L || !all(c("lb", "ub") %in% names(add_bounds)) ||
     !all(names(add_bounds) %in% c("lb", "ub"))){
    stop("'add_bounds' must contain lower and upper bounds ('lb' and 'ub').", call. = FALSE)
  }

  lb <- add_bounds[["lb"]]
  ub <- add_bounds[["ub"]]
  if(length(lb) != length(add_parameters) || length(ub) != length(add_parameters)){
    stop("'lb' and 'ub' must have the same length as 'add_parameters'.", call. = FALSE)
  }
  if(!is.numeric(lb) || !is.numeric(ub)){
    stop("'lb' and 'ub' must be numeric vectors.", call. = FALSE)
  }

  lb <- .bt_JAGS_bridge_normalize_bound_names(lb, add_parameters, "lb")
  ub <- .bt_JAGS_bridge_normalize_bound_names(ub, add_parameters, "ub")
  if(any(is.na(lb)) || any(is.na(ub))){
    stop("'add_bounds' must not contain NA values.", call. = FALSE)
  }
  if(any(lb >= ub)){
    stop("Lower bounds in 'add_bounds' must be smaller than upper bounds.", call. = FALSE)
  }

  list(lb = lb, ub = ub)
}

.bt_JAGS_bridge_normalize_bound_names <- function(bounds, add_parameters,
                                                  bound_name){

  bound_names <- names(bounds)
  if(is.null(bound_names)){
    stop(
      "'add_bounds$", bound_name,
      "' names must be unique and match 'add_parameters'.",
      call. = FALSE
    )
  }

  if(any(!nzchar(bound_names)) || anyDuplicated(bound_names)){
    stop(
      "'add_bounds$", bound_name,
      "' names must be unique and match 'add_parameters'.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(add_parameters, bound_names)
  unknown_names <- setdiff(bound_names, add_parameters)
  if(length(missing_names) > 0L || length(unknown_names) > 0L){
    stop(
      "'add_bounds$", bound_name,
      "' names must match 'add_parameters'.",
      call. = FALSE
    )
  }

  bounds[add_parameters]
}

.bt_JAGS_bridge_check_random_posterior <- function(posterior, bridge_parameters){

  if(length(bridge_parameters) == 0L){
    return(invisible(TRUE))
  }
  missing <- setdiff(bridge_parameters, colnames(posterior))
  if(length(missing) == 0L){
    return(invisible(TRUE))
  }

  stop(
    "Bridge sampling for formula random effects requires posterior samples of standardized latent random effects and LKJ primitive coordinates. ",
    "Refit with 'prior_random()' and 'random_monitor(latent = TRUE)' for the affected blocks. Missing parameter(s): ",
    paste(utils::head(missing, 8L), collapse = ", "),
    if(length(missing) > 8L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_JAGS_bridge_apply_random_scalar_rho_bounds <- function(bridgesampling_posterior,
                                                           formula_design_list){

  rho_bridge <- .bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    formula_design_list
  )
  if(length(rho_bridge$parameters) == 0L){
    return(bridgesampling_posterior)
  }

  missing <- setdiff(rho_bridge$parameters, colnames(bridgesampling_posterior))
  if(length(missing) > 0L){
    stop(
      "Bridge sampling for formula random effects requires posterior samples of scalar random-effect correlation coordinates. ",
      "Missing parameter(s): ",
      paste(utils::head(missing, 8L), collapse = ", "),
      if(length(missing) > 8L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  lb <- attr(bridgesampling_posterior, "lb")
  ub <- attr(bridgesampling_posterior, "ub")
  for(parameter in rho_bridge$parameters){
    if(!parameter %in% names(lb) || !parameter %in% names(ub)){
      stop(
        "Bridge sampling scalar random-effect correlation coordinate '",
        parameter,
        "' is missing lower or upper bridge bounds.",
        call. = FALSE
      )
    }
    next_lb <- max(lb[[parameter]], rho_bridge$bounds$lb[[parameter]])
    next_ub <- min(ub[[parameter]], rho_bridge$bounds$ub[[parameter]])
    if(is.na(next_lb) || is.na(next_ub) || next_lb >= next_ub){
      stop(
        "Bridge sampling scalar random-effect correlation bounds conflict for parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }
    lb[[parameter]] <- next_lb
    ub[[parameter]] <- next_ub
  }

  attr(bridgesampling_posterior, "lb") <- lb
  attr(bridgesampling_posterior, "ub") <- ub

  bridgesampling_posterior
}

.bt_JAGS_formula_random_scalar_rho_bridge_parameters <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  parameters <- character()
  lb <- numeric()
  ub <- numeric()

  for(design in formula_design_list){
    if(!.bt_formula_design_has_random_effects(design)){
      next
    }
    for(random_term in design$random_effects){
      rho_parameter <- .bt_JAGS_bridge_scalar_rho_parameter(random_term)
      if(is.null(rho_parameter)){
        next
      }
      parameters <- c(parameters, rho_parameter$parameter)
      lb <- c(lb, stats::setNames(rho_parameter$lower, rho_parameter$parameter))
      ub <- c(ub, stats::setNames(rho_parameter$upper, rho_parameter$parameter))
    }
  }

  if(length(parameters) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  unique_parameters <- unique(parameters)
  out_lb <- vapply(unique_parameters, function(parameter){
    max(lb[parameters == parameter])
  }, numeric(1))
  out_ub <- vapply(unique_parameters, function(parameter){
    min(ub[parameters == parameter])
  }, numeric(1))
  names(out_lb) <- unique_parameters
  names(out_ub) <- unique_parameters

  if(any(out_lb >= out_ub)){
    stop(
      "Bridge sampling scalar random-effect correlation bounds are inconsistent.",
      call. = FALSE
    )
  }

  list(
    parameters = unique_parameters,
    bounds = list(lb = out_lb, ub = out_ub)
  )
}

.bt_JAGS_bridge_scalar_rho_parameter <- function(random_term){

  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(NULL)
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation'.",
      call. = FALSE
    )
  }
  if(!is.null(.bt_random_effect_rho_fixed_sample_metadata(
    correlation,
    random_term
  ))){
    return(NULL)
  }

  parameter <- .bt_JAGS_bridge_scalar_rho_sample_name(correlation, random_term)
  bounds <- .bt_JAGS_bridge_scalar_rho_sample_bounds(correlation, random_term)

  list(
    parameter = parameter,
    lower = bounds[["lower"]],
    upper = bounds[["upper"]]
  )
}

.bt_JAGS_bridge_scalar_rho_sample_name <- function(correlation, random_term){

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term)
  parameter <- if(identical(rho_scale, "rho")){
    correlation$rho_name
  }else{
    correlation$sample_name
  }

  if(!is.character(parameter) || length(parameter) != 1L ||
     is.na(parameter) || !nzchar(parameter)){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar correlation sample name.",
      call. = FALSE
    )
  }

  parameter
}

.bt_JAGS_bridge_scalar_rho_sample_bounds <- function(correlation,
                                                     random_term){

  .bt_random_effect_rho_sample_bounds(
    correlation = correlation,
    random_term = random_term,
    context = "Bridge sampling random-effect metadata"
  )
}

.bt_JAGS_bridge_validate_formula_random_designs <- function(fitted_formula_design,
                                                            rebuilt_formula_design){

  if(is.null(fitted_formula_design) || length(fitted_formula_design) == 0L){
    return(invisible(TRUE))
  }

  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(fitted_formula_design)
  rebuilt_formula_design <- .bt_JAGS_bridge_formula_design_list(rebuilt_formula_design)
  if(length(fitted_formula_design) == 0L){
    return(invisible(TRUE))
  }

  fitted_random_parameters <- names(fitted_formula_design)[
    vapply(fitted_formula_design, .bt_formula_design_has_random_effects, logical(1))
  ]
  rebuilt_random_parameters <- names(rebuilt_formula_design)[
    vapply(rebuilt_formula_design, .bt_formula_design_has_random_effects, logical(1))
  ]
  if(length(fitted_random_parameters) == 0L &&
     length(rebuilt_random_parameters) == 0L){
    return(invisible(TRUE))
  }
  if(!setequal(fitted_random_parameters, rebuilt_random_parameters)){
    stop(
      "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design. ",
      "Random-effect formula parameter(s) differ; fitted: ",
      paste(fitted_random_parameters, collapse = ", "),
      "; rebuilt: ",
      paste(rebuilt_random_parameters, collapse = ", "),
      ". Supply the same formula, data, scaling, and prior_random() metadata used to fit the model.",
      call. = FALSE
    )
  }

  for(parameter in fitted_random_parameters){
    .bt_JAGS_bridge_validate_formula_random_design(
      parameter = parameter,
      fitted = fitted_formula_design[[parameter]],
      rebuilt = rebuilt_formula_design[[parameter]]
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_formula_design_list <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    parameter <- formula_design$parameter
    if(is.null(parameter) || length(parameter) != 1L || !nzchar(parameter)){
      parameter <- ""
    }
    out <- list(formula_design)
    names(out) <- parameter
    return(out)
  }
  if(!is.list(formula_design)){
    return(list())
  }

  design_names <- names(formula_design)
  if(is.null(design_names)){
    design_names <- rep("", length(formula_design))
  }
  for(i in seq_along(formula_design)){
    if(!nzchar(design_names[i]) &&
       inherits(formula_design[[i]], "BayesTools_formula_design") &&
       !is.null(formula_design[[i]]$parameter) &&
       length(formula_design[[i]]$parameter) == 1L){
      design_names[i] <- formula_design[[i]]$parameter
    }
  }
  names(formula_design) <- design_names
  formula_design
}

.bt_JAGS_bridge_validate_formula_random_design <- function(parameter, fitted,
                                                          rebuilt){

  if(!identical(fitted$random_effects_interface, rebuilt$random_effects_interface)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect interface differs"
    )
  }

  fitted_blocks <- .bt_JAGS_bridge_random_block_names(fitted$random_effects)
  rebuilt_blocks <- .bt_JAGS_bridge_random_block_names(rebuilt$random_effects)
  if(!identical(fitted_blocks, rebuilt_blocks)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect block names or order differ"
    )
  }

  for(block in fitted_blocks){
    fitted_term <- fitted$random_effects[[match(block, fitted_blocks)]]
    rebuilt_term <- rebuilt$random_effects[[match(block, rebuilt_blocks)]]
    .bt_JAGS_bridge_validate_random_term(
      parameter = parameter,
      block = block,
      fitted = fitted_term,
      rebuilt = rebuilt_term
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_random_block_names <- function(random_effects){

  if(length(random_effects) == 0L){
    return(character())
  }

  vapply(random_effects, function(random_term){
    block_name <- random_term$block_name
    if(is.null(block_name) || length(block_name) != 1L){
      return("")
    }
    as.character(block_name)
  }, character(1))
}

.bt_JAGS_bridge_validate_random_term <- function(parameter, block, fitted,
                                                rebuilt){

  if(!identical(.bt_JAGS_bridge_random_term_structure(fitted),
                .bt_JAGS_bridge_random_term_structure(rebuilt))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect covariance structure differs",
      block = block
    )
  }
  if(!identical(fitted$group_label, rebuilt$group_label)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect grouping label differs",
      block = block
    )
  }
  if(!identical(as.character(fitted$group_levels),
                as.character(rebuilt$group_levels))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group levels or their order differ",
      block = block
    )
  }
  if(!identical(fitted$n_groups, rebuilt$n_groups) ||
     !identical(fitted$n_columns, rebuilt$n_columns)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group or column counts differ",
      block = block
    )
  }
  if(!identical(fitted$column_names, rebuilt$column_names) ||
     !identical(fitted$raw_column_names, rebuilt$raw_column_names)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect design column names differ",
      block = block
    )
  }
  if(!identical(dim(fitted$model_matrix), dim(rebuilt$model_matrix)) ||
     !identical(colnames(fitted$model_matrix), colnames(rebuilt$model_matrix))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect model matrix shape or columns differ",
      block = block
    )
  }
  if(!isTRUE(all.equal(
    unname(fitted$model_matrix),
    unname(rebuilt$model_matrix),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect model matrix values differ",
      block = block
    )
  }
  if(!identical(as.integer(fitted$group_map),
                as.integer(rebuilt$group_map))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group map differs",
      block = block
    )
  }
  if(!.bt_JAGS_bridge_metadata_equal(
    .bt_JAGS_bridge_scale_metadata(fitted),
    .bt_JAGS_bridge_scale_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect scale/allocation metadata differ",
      block = block
    )
  }
  if(!.bt_JAGS_bridge_metadata_equal(
    .bt_JAGS_bridge_structured_index_metadata(fitted),
    .bt_JAGS_bridge_structured_index_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "structured random-effect index metadata differ",
      block = block
    )
  }
  if(!.bt_JAGS_bridge_metadata_equal(
    .bt_JAGS_bridge_car_metadata(fitted),
    .bt_JAGS_bridge_car_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "CAR random-effect index metadata differ",
      block = block
    )
  }
  if(!.bt_JAGS_bridge_metadata_equal(
    .bt_JAGS_bridge_correlation_metadata(fitted),
    .bt_JAGS_bridge_correlation_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect correlation metadata differ",
      block = block
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_random_term_structure <- function(random_term){

  .bt_random_effect_structure(
    random_term,
    context = "Bridge sampling random-effect metadata"
  )
}

.bt_JAGS_bridge_scale_metadata <- function(random_term){

  list(
    sd_parameter_names = random_term$sd_parameter_names,
    sd_binding = .bt_JAGS_bridge_sd_binding_metadata(
      random_term$sd_binding,
      n_columns = random_term$n_columns
    ),
    homogeneous_sd = .bt_random_effect_homogeneous_sd_metadata(
      random_term,
      context = "Bridge sampling random-effect metadata"
    )
  )
}

.bt_JAGS_bridge_sd_binding_metadata <- function(binding, n_columns = NULL){

  if(is.null(binding)){
    return(NULL)
  }

  .bt_check_random_sd_binding(binding)
  if(isTRUE(binding$true_allocation) &&
     length(binding$allocations) > 0L &&
     identical(binding$allocations[[1L]]$target, "sd_component")){
    .bt_check_random_sd_component_binding(
      binding = binding,
      n_columns = n_columns,
      context = "Bridge sampling random-effect metadata"
    )
  }
  out <- binding[intersect(
    names(binding),
    c(
      "application", "factors", "factors_by_column",
      "true_allocation", "allocations", "sd_component_names",
      "sd_component_terms", "sd_component_index_by_column"
    )
  )]
  out$source <- .bt_JAGS_bridge_parameter_source_metadata(binding$source)
  out$sources_by_column <- lapply(
    binding$sources_by_column,
    .bt_JAGS_bridge_parameter_source_metadata
  )
  if(!is.null(out$factors)){
    out$factors <- lapply(out$factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }
  if(!is.null(out$factors_by_column)){
    out$factors_by_column <- lapply(out$factors_by_column, function(factors){
      lapply(factors, .bt_JAGS_bridge_allocation_factor_metadata)
    })
  }
  if(!is.null(out$allocations)){
    out$allocations <- lapply(out$allocations, .bt_JAGS_bridge_allocation_metadata)
  }

  out
}

.bt_JAGS_bridge_allocation_metadata <- function(allocation){

  if(is.null(allocation)){
    return(NULL)
  }

  allocation <- allocation[intersect(
    names(allocation),
    c(
      "label", "terms", "index", "target", "scale",
      "parent", "source", "factors", "parent_factors", "n_targets", "total_name",
      "weight_name", "leaf_names", "leaf_terms",
      "leaf_index_by_column", "source_node", "weight_suffix",
      "total_suffix", "sd_component_names", "sd_component_terms",
      "sd_component_index_by_column"
    )
  )]
  if(!is.null(allocation$source)){
    allocation$source <- .bt_JAGS_bridge_parameter_source_metadata(allocation$source)
  }
  if(!is.null(allocation$factors)){
    allocation$factors <- lapply(allocation$factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }
  if(!is.null(allocation$parent_factors)){
    allocation$parent_factors <- lapply(allocation$parent_factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }

  allocation
}

.bt_JAGS_bridge_allocation_factor_metadata <- function(factor){

  factor[intersect(
    names(factor),
    c("weight_name", "index", "scale", "n_targets")
  )]
}

.bt_JAGS_bridge_parameter_source_metadata <- function(source){

  if(is.null(source)){
    return(NULL)
  }
  if(!is.list(source)){
    return(source)
  }

  source <- source[intersect(
    names(source),
    c(
      "name", "shape", "kind", "owned",
      "total_name", "total_suffix", "source"
    )
  )]
  if(!is.null(source$source)){
    source$source <- .bt_JAGS_bridge_parameter_source_metadata(source$source)
  }

  source
}

.bt_JAGS_bridge_structured_index_metadata <- function(random_term){

  index <- random_term$structured_index
  if(is.null(index)){
    return(NULL)
  }

  index[intersect(names(index), c("variables", "name", "label", "structure"))]
}

.bt_JAGS_bridge_car_metadata <- function(random_term){

  car <- random_term$car
  if(is.null(car)){
    return(NULL)
  }

  car[intersect(names(car), c("time_variable", "time_values", "distance_matrix"))]
}

.bt_JAGS_bridge_correlation_metadata <- function(random_term){

  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation)){
    return(NULL)
  }

  correlation[intersect(
    names(correlation),
    c(
      "type", "eta", "primitive_names", "primitive_bounds",
      "structure", "rho_name", "sample_name", "prior_name", "parameter_name",
      "sample_fixed", "rho_scale", "bounds", "cholesky_name", "correlation_name",
      "time_variable", "time_values", "distance_matrix"
    )
  )]
}

.bt_JAGS_bridge_metadata_equal <- function(x, y){

  isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE))
}

.bt_JAGS_bridge_random_design_mismatch <- function(parameter, detail,
                                                   block = NULL){

  stop(
    "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design for parameter '",
    parameter,
    "'",
    if(!is.null(block)) paste0(", block '", block, "'") else "",
    ": ",
    detail,
    ". Supply the same formula, data, scaling, and prior_random() metadata used to fit the model.",
    call. = FALSE
  )
}

.JAGS_formula_scale_list_from_fit <- function(fit, formula_parameters){

  formula_scale <- attr(fit, "formula_scale")
  if(is.null(formula_scale) || length(formula_scale) == 0L){
    return(NULL)
  }

  scale_list <- vector("list", length(formula_parameters))
  names(scale_list) <- formula_parameters

  for(parameter in intersect(formula_parameters, names(formula_scale))){
    parameter_scale <- formula_scale[[parameter]]
    if(is.null(parameter_scale) || length(parameter_scale) == 0L){
      next
    }

    scaled_terms <- names(parameter_scale)
    if(is.null(scaled_terms) || length(scaled_terms) == 0L){
      next
    }

    parameter_prefix <- paste0(parameter, "_")
    predictor_terms  <- ifelse(
      startsWith(scaled_terms, parameter_prefix),
      substring(scaled_terms, nchar(parameter_prefix) + 1L),
      scaled_terms
    )

    scale_list[[parameter]] <- as.list(stats::setNames(rep(TRUE, length(predictor_terms)), predictor_terms))
  }

  scale_list <- scale_list[!vapply(scale_list, is.null, logical(1))]
  if(length(scale_list) == 0L){
    return(NULL)
  }

  return(scale_list)
}

.fit_to_posterior <- function(fit){

  ### check the input and split it on posterior and data
  if(inherits(fit, "runjags")){

    # get posterior and merge chains
    posterior <- .extract_posterior_samples(fit, as_list = FALSE)

  }else if(is.list(fit) & all(sapply(fit, inherits, what = "mcarray"))){

    # rjags model with rjags::jags.samples
    # merge chains
    posterior <- do.call(cbind, lapply(names(fit), function(par){
      if(dim(fit[[par]])[1] > 1){
        samples <- do.call(rbind, lapply(1:(dim(fit[[par]]))[3],  function(chain)t(fit[[par]][,,chain])))
        colnames(samples) <- paste0(attr(fit[[par]], "varname"), "[",1:ncol(samples),"]")
      }else{
        samples <- matrix(do.call(c, lapply(1:(dim(fit[[par]]))[3],  function(chain)fit[[par]][,,chain])), ncol = 1)
        colnames(samples) <- attr(fit[[par]], "varname")
      }
      return(samples)
    }))

  }else if(inherits(fit, "mcmc.list")){

    # rjags model with rjags::coda.samples or samples extracted via coda::as.mcmc.list
    # merge chains
    posterior <- do.call(rbind, fit)

  }else if (inherits(fit, "mcmc") && length(dim(fit)) == 2) {

    # rjags model with samples extracted via coda::as.mcmc
    return(fit)

  } else {

    stop("the method is not implemented for this output")

  }

  return(posterior)
}


#' @title Create a 'bridgesampling' object
#'
#' @description prepares a 'bridgesampling' object with a given
#' log marginal likelihood.
#'
#' @param logml log marginal likelihood. Defaults to \code{-Inf}.
#'
#'
#' @return \code{JAGS_bridgesampling} returns an object of class 'bridge'.
#'
#' @export
bridgesampling_object <- function(logml = -Inf){

  marglik        <- list()
  marglik$logml  <- logml
  class(marglik) <- "bridge"

  return(marglik)
}


#' @title Prepare 'JAGS' posterior for 'bridgesampling'
#'
#' @description prepares posterior distribution for 'bridgesampling'
#' by removing unnecessary parameters and attaching lower and upper
#' bounds of parameters based on a list of prior distributions.
#'
#' @param posterior matrix of mcmc samples from the posterior distribution.
#' @param prior_list named list of prior distributions. Names correspond to
#' parameter names owned by BayesTools priors and determine the posterior
#' columns retained for bridge sampling together with their bounds.
#' @param add_parameters character vector of additional monitored posterior
#' parameter names to retain in addition to the prior-owned parameters. These
#' parameters are passed through as named entries and must be indexed by name
#' by downstream likelihood code.
#' @param add_bounds list with two named numeric vectors, \code{"lb"} and
#' \code{"ub"}, containing lower and upper bounds for every
#' \code{add_parameters} entry.
#'
#' @return \code{JAGS_bridgesampling_posterior} returns a matrix of
#' posterior samples with 'lb' and 'ub' attributes carrying the
#' lower and upper boundaries.
#'
#' @export
JAGS_bridgesampling_posterior <- function(posterior, prior_list, add_parameters = NULL, add_bounds = NULL){

  # check the input
  if(!is.matrix(posterior))
    stop("'posterior' must be a matrix")
  if(is.null(prior_list)){
    prior_list <- list()
  }
  if(!is.null(prior_list)){
    if(!is.list(prior_list))
      stop("'prior_list' must be a list.")
    if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
      stop("'prior_list' must be a list of priors.")
  }
  if(is.null(add_parameters) && !is.null(add_bounds)){
    stop("'add_bounds' requires at least one 'add_parameters' entry.", call. = FALSE)
  }
  if(!is.null(add_parameters)){
    if(!is.character(add_parameters))
      stop("'add_parameters' must be a character vector.")
    if(length(add_parameters) == 0L){
      if(!is.null(add_bounds)){
        stop("'add_bounds' requires at least one 'add_parameters' entry.", call. = FALSE)
      }
      add_parameters <- NULL
    }else{
      add_bounds <- .bt_JAGS_bridge_validate_add_bounds(add_parameters, add_bounds)
    }
  }

  # these are not generally possible because the component indicators are discrete and bridgesampling
  # package cannot currently deal with them
  if(length(prior_list) > 0L && any(sapply(prior_list, is.prior.spike_and_slab)))
    stop("Marginal likelihood computation for spike and slab priors is not implemented.")
  if(length(prior_list) > 0L && any(sapply(prior_list, is.prior.mixture))){
    .JAGS_marglik_stop_unsupported_mixture(prior_list[[which(sapply(prior_list, is.prior.mixture))[1L]]])
  }

  # get information about the specified parameters
  parameters_names <- .JAGS_bridgesampling_posterior_info(prior_list)
  owned_parameter_names <- .bt_JAGS_bridge_owned_parameter_names(prior_list)

  # add the user defined parameters
  if(!is.null(add_parameters)){
    overlapping_parameters <- intersect(
      add_parameters,
      unique(c(parameters_names, owned_parameter_names))
    )
    if(length(overlapping_parameters) > 0L){
      stop(
        "'add_parameters' contains BayesTools-owned parameter(s) already covered by 'prior_list': ",
        paste(overlapping_parameters, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    .bt_JAGS_bridge_validate_add_parameters_not_prior_dirichlet(
      add_parameters = add_parameters,
      prior_list = prior_list
    )
    parameters_names_lb <- c(attr(parameters_names, "lb"), add_bounds[["lb"]])
    parameters_names_ub <- c(attr(parameters_names, "ub"), add_bounds[["ub"]])
    parameters_names <- c(parameters_names, add_parameters)
    attr(parameters_names, "lb") <- parameters_names_lb
    attr(parameters_names, "ub") <- parameters_names_ub
  }


  # check that all parameter names exist in the posterior
  if(!all(parameters_names %in% colnames(posterior)))
    stop("'posterior' does not contain all of the parameters corresponding to the 'prior_list' and the 'add_parameter' argument.")

  posterior <- posterior[,parameters_names, drop = FALSE]
  attr(posterior, "lb") <- attr(parameters_names, "lb")
  attr(posterior, "ub") <- attr(parameters_names, "ub")

  return(posterior)
}

.bt_JAGS_bridge_validate_add_parameters_not_prior_dirichlet <- function(
    add_parameters, prior_list){

  if(length(add_parameters) == 0L || length(prior_list) == 0L){
    return(invisible(TRUE))
  }

  normalized_dirichlet <- character()
  for(parameter_name in names(prior_list)){
    prior <- prior_list[[parameter_name]]
    if(is.prior.simplex(prior) && identical(prior$distribution, "dirichlet")){
      K <- prior$parameters[["K"]]
      normalized_dirichlet <- c(
        normalized_dirichlet,
        paste0(parameter_name, "[", seq_len(K), "]")
      )
    }
  }

  overlapping <- intersect(add_parameters, normalized_dirichlet)
  if(length(overlapping) > 0L){
    stop(
      "'add_parameters' must not contain normalized Dirichlet coordinate(s) ",
      "already represented by auxiliary bridge parameters from 'prior_list': ",
      paste(overlapping, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_validate_add_parameters_not_formula <- function(
    add_parameters, formula_design_list = NULL, formula_prior_list = NULL){

  if(is.null(add_parameters) || length(add_parameters) == 0L){
    return(invisible(TRUE))
  }

  formula_parameter_names <- unique(c(
    names(formula_design_list),
    names(formula_prior_list)
  ))
  formula_parameter_names <- formula_parameter_names[
    !is.na(formula_parameter_names) & nzchar(formula_parameter_names)
  ]
  overlapping <- intersect(add_parameters, formula_parameter_names)
  if(length(overlapping) > 0L){
    stop(
      "'add_parameters' contains BayesTools-owned formula parameter(s): ",
      paste(overlapping, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_owned_parameter_names <- function(prior_list){

  if(length(prior_list) == 0L){
    return(character())
  }

  if(!is.list(prior_list)){
    stop("'prior_list' must be a list.")
  }
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior))){
    stop("'prior_list' must be a list of priors.")
  }
  .check_prior_list_unique_names(prior_list)

  owned <- character()
  for(i in seq_along(prior_list)){
    prior <- prior_list[[i]]
    parameter_name <- names(prior_list)[i]
    if(is.prior.weightfunction(prior)){
      owned <- c(owned, .bt_JAGS_bridge_owned_weightfunction_names(prior))
    }else if(is_prior_phacking(prior)){
      owned <- c(owned, .bt_JAGS_bridge_owned_phacking_names(prior))
    }else if(is_prior_bias(prior)){
      spec <- selection_backend_spec(prior)
      owned <- c(owned, spec$monitor)
      if(!is.null(prior$selection)){
        owned <- c(owned, .bt_JAGS_bridge_owned_weightfunction_names(prior$selection))
      }
      if(!is.null(prior$phacking)){
        owned <- c(owned, .bt_JAGS_bridge_owned_phacking_names(prior$phacking))
      }
    }else if(is.prior.mixture(prior)){
      .JAGS_marglik_stop_unsupported_mixture(prior)
    }else if(is.prior.PET(prior)){
      owned <- c(owned, "PET")
    }else if(is.prior.PEESE(prior)){
      owned <- c(owned, "PEESE")
    }else if(is.prior.factor(prior) ||
             is.prior.vector(prior) ||
             is.prior.simple(prior)){
      owned <- c(owned, .bt_JAGS_bridge_owned_prior_parameter_names(
        prior = prior,
        parameter_name = parameter_name
      ))
    }
  }

  unique(owned)
}

.bt_JAGS_bridge_owned_weightfunction_names <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior)){
    stop("improper prior provided")
  }

  spec <- selection_backend_spec(prior)
  J <- .weightfunction_n_bins(prior)
  owned <- c(
    spec$monitor,
    spec$step$coefficient_ids,
    "omega",
    paste0("omega[", seq_len(J), "]")
  )
  if(identical(prior$weights$type, "cumulative")){
    owned <- c(
      owned,
      "eta", paste0("eta[", seq_len(J), "]"),
      "std_eta", paste0("std_eta[", seq_len(J), "]")
    )
  }else if(identical(prior$weights$type, "independent") &&
           identical(prior$weights$scale, "log_omega")){
    owned <- c(owned, "log_omega")
    if(J > 1L){
      owned <- c(owned, paste0("log_omega[", 2:J, "]"))
    }
  }

  owned
}

.bt_JAGS_bridge_owned_phacking_names <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior)){
    stop("improper prior provided")
  }

  spec <- selection_backend_spec(prior)
  c(
    spec$monitor,
    spec$step$coefficient_ids,
    "alpha", "pi_null", "beta_null", "phack_kind",
    "phack_z_source", "phack_z_source[1]", "phack_z_source[2]",
    "phack_z_dest", "phack_z_dest[1]", "phack_z_dest[2]"
  )
}

.bt_JAGS_bridge_owned_prior_parameter_names <- function(prior,
                                                        parameter_name){

  .check_prior(prior)
  check_char(parameter_name, "parameter_name")

  owned <- parameter_name
  if(is.prior.factor(prior)){
    K <- .get_prior_factor_levels(prior)
    if(K > 1L){
      owned <- c(owned, paste0(parameter_name, "[", seq_len(K), "]"))
    }
  }else if(is.prior.vector(prior)){
    K <- prior$parameters[["K"]]
    if(is.numeric(K) && length(K) == 1L && !is.na(K) && K > 1L){
      owned <- c(owned, paste0(parameter_name, "[", seq_len(K), "]"))
    }
  }

  owned
}

.JAGS_marglik_stop_unsupported_mixture <- function(prior){

  if(inherits(prior, "prior.bias_mixture")){
    selection_backend_spec(prior)
    stop(
      "Marginal likelihood computation for bias mixture priors is not implemented because bridge sampling does not support discrete bias indicators.",
      call. = FALSE
    )
  }

  stop("Marginal likelihood computation for prior mixture priors is not implemented.", call. = FALSE)
}

.JAGS_bridgesampling_posterior_info                <- function(prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    parameters <- character()
    attr(parameters, "lb") <- numeric()
    attr(parameters, "ub") <- numeric()
    return(parameters)
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)


  # add the resulting parameters
  parameters    <- character()
  parameters_lb <- numeric()
  parameters_ub <- numeric()
  for(i in seq_along(prior_list)){

    add_parameter <- NULL

    if(is.prior.weightfunction(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.weightfunction(prior_list[[i]])

    }else if(is_prior_phacking(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.phacking(prior_list[[i]])

    }else if(is_prior_bias(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.bias(prior_list[[i]])

    }else if(is.prior.mixture(prior_list[[i]])){

      .JAGS_marglik_stop_unsupported_mixture(prior_list[[i]])

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.PP(prior_list[[i]])

    }else if(is.prior.factor(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.factor(prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.vector(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.vector(prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.simple(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.simple(prior_list[[i]], names(prior_list)[i])

    }

    if(!is.null(add_parameter)){
      parameters    <- c(parameters,    add_parameter)
      parameters_lb <- c(parameters_lb, attr(add_parameter, "lb"))
      parameters_ub <- c(parameters_ub, attr(add_parameter, "ub"))
    }
  }

  attr(parameters, "lb") <- parameters_lb
  attr(parameters, "ub") <- parameters_ub

  return(parameters)
}
.JAGS_bridgesampling_posterior_info.simple         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")


  if(prior[["distribution"]] == "invgamma"){
    parameter <- paste0("inv_", parameter_name)
    attr(parameter, "lb") <- prior$truncation[["upper"]]^-1
    attr(parameter, "ub") <- prior$truncation[["lower"]]^-1
  }else if(prior[["distribution"]] == "point"){
    parameter <- NULL
  }else{
    parameter <- parameter_name
    attr(parameter, "lb") <- prior$truncation[["lower"]]
    attr(parameter, "ub") <- prior$truncation[["upper"]]
  }

  names(attr(parameter, "lb")) <- parameter
  names(attr(parameter, "ub")) <- parameter

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.vector         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(prior$parameters[["K"]], "K", lower = 1)
  if(prior[["distribution"]] != "mpoint")
    .check_vector_truncation_unsupported(prior$truncation)

  if(prior[["distribution"]] == "mpoint"){
    parameter <- NULL
  }else if(prior[["distribution"]] == "dirichlet"){
    parameter <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(prior$parameters[["K"]]), "]")
    attr(parameter, "lb") <- rep(0, length(parameter))
    attr(parameter, "ub") <- rep(Inf, length(parameter))
    names(attr(parameter, "lb")) <- parameter
    names(attr(parameter, "ub")) <- parameter
  }else{
    if(prior$parameters[["K"]] == 1){
      parameter <- parameter_name
    }else{
      parameter <- paste0(parameter_name, "[", 1:prior$parameters[["K"]], "]")
    }

    attr(parameter, "lb") <- rep(prior$truncation[["lower"]], prior$parameters[["K"]])
    attr(parameter, "ub") <- rep(prior$truncation[["upper"]], prior$parameters[["K"]])

    names(attr(parameter, "lb")) <- parameter
    names(attr(parameter, "ub")) <- parameter
  }

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.factor         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    if(.get_prior_factor_levels(prior) == 1){

      parameter <- .JAGS_bridgesampling_posterior_info.simple(prior, parameter_name)

    }else{

      parameter    <- NULL
      parameter_lb <- NULL
      parameter_ub <- NULL

      for(i in 1:.get_prior_factor_levels(prior)){

        add_parameter <- .JAGS_bridgesampling_posterior_info.simple(prior, paste0(parameter_name, "[", i, "]"))

        parameter    <- c(parameter,    add_parameter)
        parameter_lb <- c(parameter_lb, attr(add_parameter, "lb"))
        parameter_ub <- c(parameter_ub, attr(add_parameter, "ub"))
      }

      attr(parameter, "lb") <- parameter_lb
      attr(parameter, "ub") <- parameter_ub

    }

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    parameter <- .JAGS_bridgesampling_posterior_info.vector(prior, parameter_name)

  }

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.PP             <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    parameter <- .JAGS_bridgesampling_posterior_info.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    parameter <- .JAGS_bridgesampling_posterior_info.simple(prior, "PEESE")
  }

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.weightfunction <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "cumulative"){

    parameter <- paste0("eta[", seq_len(J), "]")
    attr(parameter, "lb") <- rep(0,   length(parameter))
    attr(parameter, "ub") <- rep(Inf, length(parameter))

  }else if(prior$weights$type == "independent" && prior$weights$scale == "omega"){

    parameter <- if(J > 1L) paste0("omega[", 2:J, "]") else NULL
    attr(parameter, "lb") <- rep(prior$weights$prior$truncation[["lower"]], length(parameter))
    attr(parameter, "ub") <- rep(prior$weights$prior$truncation[["upper"]], length(parameter))

  }else if(prior$weights$type == "independent" && prior$weights$scale == "log_omega"){

    parameter <- if(J > 1L) paste0("log_omega[", 2:J, "]") else NULL
    attr(parameter, "lb") <- rep(prior$weights$prior$truncation[["lower"]], length(parameter))
    attr(parameter, "ub") <- rep(prior$weights$prior$truncation[["upper"]], length(parameter))

  }else if(prior$weights$type == "fixed"){

    parameter <- NULL

  }

  if(!is.null(parameter)){
    names(attr(parameter, "lb")) <- parameter
    names(attr(parameter, "ub")) <- parameter
  }

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.phacking <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  parameter <- .JAGS_bridgesampling_posterior_info.simple(prior$alpha, "alpha")

  return(parameter)
}
.JAGS_bridgesampling_posterior_info.bias <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)

  parameter <- NULL
  parameter_lb <- NULL
  parameter_ub <- NULL

  if(!is.null(prior$selection)){
    selection_parameter <- .JAGS_bridgesampling_posterior_info.weightfunction(prior$selection)
    parameter <- c(parameter, selection_parameter)
    parameter_lb <- c(parameter_lb, attr(selection_parameter, "lb"))
    parameter_ub <- c(parameter_ub, attr(selection_parameter, "ub"))
  }
  if(!is.null(prior$phacking)){
    phacking_parameter <- .JAGS_bridgesampling_posterior_info.phacking(prior$phacking)
    parameter <- c(parameter, phacking_parameter)
    parameter_lb <- c(parameter_lb, attr(phacking_parameter, "lb"))
    parameter_ub <- c(parameter_ub, attr(phacking_parameter, "ub"))
  }

  attr(parameter, "lb") <- parameter_lb
  attr(parameter, "ub") <- parameter_ub
  return(parameter)
}
# .JAGS_bridgesampling_posterior_info.spike_and_slab <- function(prior, parameter_name){
#
#   .check_prior(prior)
#   if(!is.prior.spike_and_slab(prior))
#     stop("improper prior provided")
#   check_char(parameter_name, "parameter_name")
#
#   if(!is.prior.point(prior[["inclusion"]])){
#
#     parameter_variable  <- .JAGS_bridgesampling_posterior_info.simple(prior[["variable"]],  paste0(parameter_name, "_variable"))
#     parameter_inclusion <- .JAGS_bridgesampling_posterior_info.simple(prior[["inclusion"]], paste0(parameter_name, "_inclusion"))
#
#     parameter <- c(parameter_variable, parameter_inclusion)
#
#     attr(parameter, "lb") <- c(attr(parameter_variable, "lb"), attr(parameter_inclusion, "lb"))
#     attr(parameter, "ub") <- c(attr(parameter_variable, "ub"), attr(parameter_inclusion, "ub"))
#
#     names(attr(parameter, "lb")) <- c(names(attr(parameter_variable, "lb")), names(attr(parameter_inclusion, "lb")))
#     names(attr(parameter, "ub")) <- c(names(attr(parameter_variable, "ub")), names(attr(parameter_inclusion, "ub")))
#
#   }else{
#     parameter  <- .JAGS_bridgesampling_posterior_info.simple(prior[["variable"]],  paste0(parameter_name, "_variable"))
#   }
#
#
#   return(parameter)
# }

#' @title Compute marginal likelihood for 'JAGS' priors
#'
#' @description Computes marginal likelihood for the
#' prior part of a 'JAGS' model within 'bridgesampling'
#' function
#'
#' @param samples samples provided by bridgesampling
#' function
#'
#' @inheritParams JAGS_bridgesampling
#'
#' @return \code{JAGS_marglik_priors} returns a numeric value
#' of likelihood evaluated at the current posterior sample.
#'
#' @export JAGS_marglik_priors
#' @export JAGS_marglik_priors_formula
#' @name JAGS_marglik_priors
NULL

#' @rdname JAGS_marglik_priors
JAGS_marglik_priors                <- function(samples, prior_list){

  # return zero log prior contribution in case that no prior was specified
  if(length(prior_list) == 0){
    return(0)
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)


  # add the resulting parameters
  marglik <- 0
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.weightfunction(samples, prior_list[[i]])

    }else if(is_prior_phacking(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.phacking(samples, prior_list[[i]])

    }else if(is_prior_bias(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.bias(samples, prior_list[[i]])

    }else if(is.prior.mixture(prior_list[[i]])){

      .JAGS_marglik_stop_unsupported_mixture(prior_list[[i]])

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.PP(samples, prior_list[[i]])

    }else if(is.prior.factor(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.factor(samples, prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.vector(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.vector(samples, prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.simple(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.simple(samples, prior_list[[i]], names(prior_list)[i])

    }
  }

  return(marglik)
}


.JAGS_marglik_priors.simple         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(prior[["distribution"]] == "invgamma"){

    inv_name <- paste0("inv_", parameter_name)
    inv_value <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = inv_name,
      missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters."
    )
    if(is.null(inv_value)){
      return(-Inf)
    }
    sampling_prior <- prior(
      "distribution" = "gamma",
      "parameters"   = list("shape" = prior$parameters[["shape"]], "rate" = prior$parameters[["scale"]]),
      "truncation"   = list("lower" = prior$truncation[["upper"]]^-1, "upper" = prior$truncation[["lower"]]^-1))
    marglik <- lpdf(sampling_prior, inv_value)

  }else if(prior[["distribution"]] == "point"){

    marglik <- 0

  }else{

    marglik <- lpdf(prior, samples[[ parameter_name ]])

  }

  return(marglik)
}
.JAGS_marglik_priors.vector         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(prior$parameters[["K"]], "K", lower = 1)
  if(prior[["distribution"]] != "mpoint")
    .check_vector_truncation_unsupported(prior$truncation)

  if(prior[["distribution"]] == "mpoint"){
    marglik <- 0
  }else if(prior[["distribution"]] == "dirichlet"){
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(prior$parameters[["K"]]), "]")
    eta <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = eta_names,
      missing_message = "'samples' does not contain all monitored Dirichlet prior parameters."
    )
    if(is.null(eta)){
      return(-Inf)
    }
    marglik <- sum(stats::dgamma(
      eta,
      shape = prior$parameters[["alpha"]],
      rate = 1,
      log = TRUE
    ))
  }else if(prior$parameters[["K"]] == 1){
    marglik <- lpdf(prior, samples[[ parameter_name ]])
  }else{
    marglik <- lpdf(prior, samples[ paste0(parameter_name, "[", 1:prior$parameters[["K"]], "]") ])
  }

  return(marglik)
}
.JAGS_marglik_priors.factor         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    if(.get_prior_factor_levels(prior) == 1){

      marglik <- .JAGS_marglik_priors.simple(samples, prior, parameter_name)

    }else{

      marglik <- sum(sapply(1:.get_prior_factor_levels(prior), function(i) .JAGS_marglik_priors.simple(samples, prior, paste0(parameter_name, "[", i, "]"))))

    }

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    marglik <- .JAGS_marglik_priors.vector(samples, prior, parameter_name)

  }

  return(marglik)
}
.JAGS_marglik_priors.PP             <- function(samples, prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    marglik <- .JAGS_marglik_priors.simple(samples, prior, "PET")
  }else if(is.prior.PEESE(prior)){
    marglik <- .JAGS_marglik_priors.simple(samples, prior, "PEESE")
  }

  return(marglik)
}
.JAGS_marglik_priors.weightfunction <- function(samples, prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){

    marglik <- 0

  }else if(prior$weights$type == "cumulative"){

    eta_names <- paste0("eta[", seq_len(J), "]")
    eta <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = eta_names,
      missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters."
    )
    if(is.null(eta)){
      return(-Inf)
    }
    marglik <- sum(stats::dgamma(eta, shape = prior$weights$alpha, rate = 1, log = TRUE))

  }else if(prior$weights$type == "independent"){

    if(J == 1L){
      marglik <- 0
    }else if(prior$weights$scale == "omega"){
      marglik <- sum(mlpdf(prior$weights$prior, samples[paste0("omega[", 2:J, "]")]))
    }else if(prior$weights$scale == "log_omega"){
      marglik <- sum(mlpdf(prior$weights$prior, samples[paste0("log_omega[", 2:J, "]")]))
    }

  }

  return(marglik)
}
.JAGS_marglik_priors.phacking <- function(samples, prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  .JAGS_marglik_priors.simple(samples, prior$alpha, "alpha")
}
.JAGS_marglik_priors.bias <- function(samples, prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)

  marglik <- 0
  if(!is.null(prior$selection)){
    marglik <- marglik + .JAGS_marglik_priors.weightfunction(samples, prior$selection)
  }
  if(!is.null(prior$phacking)){
    marglik <- marglik + .JAGS_marglik_priors.phacking(samples, prior$phacking)
  }

  return(marglik)
}
# .JAGS_marglik_priors.spike_and_slab <- function(samples, prior, parameter_name){
#
#   .check_prior(prior)
#   if(!is.prior.spike_and_slab(prior))
#     stop("improper prior provided")
#   check_char(parameter_name, "parameter_name")
#
#   marglik <- 0
#   if(!is.prior.point(prior[["inclusion"]])){
#     marglik <- marglik + .JAGS_marglik_priors.simple(samples, prior[["inclusion"]], paste0(parameter_name, "_inclusion"))
#   }
#
#   if(samples[[ paste0(if(prior[["variable"]][["distribution"]] == "invgamma") "inv_" else "", parameter_name, "_variable") ]] != 0){
#     marglik <- marglik + .JAGS_marglik_priors.simple(samples, prior[["variable"]], paste0(parameter_name, "_variable"))
#   }
#
#   return(marglik)
# }

#' @rdname JAGS_marglik_priors
JAGS_marglik_priors_formula <- function(samples, formula_prior_list){

  marglik <- 0

  for(parameter in names(formula_prior_list)){
    marglik <- marglik + JAGS_marglik_priors(samples, formula_prior_list[[parameter]])
  }

  return(marglik)
}

.bt_JAGS_marglik_priors_formula_random <- function(samples, formula_design_list){

  if(length(formula_design_list) == 0L){
    return(0)
  }

  marglik <- 0
  for(design in formula_design_list){
    if(!.bt_formula_design_has_random_effects(design)){
      next
    }
    for(random_term in design$random_effects){
      contribution <- .bt_JAGS_marglik_random_effect_prior(samples, random_term)
      if(is.na(contribution)){
        return(-Inf)
      }
      marglik <- marglik + contribution
      if(is.na(marglik)){
        return(-Inf)
      }
    }
  }

  marglik
}

.bt_JAGS_marglik_random_effect_prior <- function(samples, random_term){

  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  ))
  if(!all(z_names %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  z_values <- samples[z_names]
  if(any(is.na(z_values))){
    return(-Inf)
  }
  marglik <- sum(stats::dnorm(z_values, mean = 0, sd = 1, log = TRUE))
  if(is.na(marglik)){
    return(-Inf)
  }

  scalar_rho_support <- .bt_JAGS_marglik_random_effect_scalar_rho_support(
    samples = samples,
    random_term = random_term
  )
  if(!is.finite(scalar_rho_support)){
    return(scalar_rho_support)
  }
  marglik <- marglik + scalar_rho_support

  if(identical(.bt_JAGS_bridge_random_term_structure(random_term), "us") &&
     n_columns > 1L){
    u_names <- .bt_random_effect_lkj_primitive_names(
      random_term,
      n_columns,
      context = "Bridge sampling random-effect metadata"
    )
    if(!all(u_names %in% names(samples))){
      stop(
        "Bridge samples are missing LKJ primitive coordinates for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    u_values <- unname(samples[u_names])
    if(any(is.na(u_values) | u_values <= 0 | u_values >= 1)){
      return(-Inf)
    }
    correlation <- .bt_random_effect_correlation_metadata(
      random_term,
      structure = "us",
      context = "Bridge sampling random-effect metadata"
    )
    eta <- correlation$eta
    if(!is.numeric(eta) || length(eta) != 1L || is.na(eta)){
      stop(
        "Bridge sampling random-effect metadata",
        .bt_random_effect_metadata_block_detail(random_term),
        " is missing canonical 'random_term$correlation$eta'.",
        call. = FALSE
      )
    }
    marglik <- marglik + .bt_lkj_cholesky_cpc_u_log_prior(
      u_values,
      K = n_columns,
      eta = eta
    )
    if(is.na(marglik)){
      return(-Inf)
    }
  }

  marglik
}

.bt_JAGS_marglik_random_effect_scalar_rho_support <- function(samples,
                                                              random_term){

  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(0)
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation'.",
      call. = FALSE
    )
  }

  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
    missing = "error",
    out_of_support = "null",
    sample_space = TRUE,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(rho) || any(is.na(rho) | !is.finite(rho))){
    return(-Inf)
  }

  R <- .bt_random_effect_structured_correlation_matrix(
    structure = structure,
    K = random_term$n_columns,
    rho = rho[1L],
    distance_matrix = if(identical(structure, "car")) correlation$distance_matrix else NULL
  )
  chol_ok <- tryCatch({
    chol(R)
    TRUE
  }, error = function(e) FALSE)
  if(!chol_ok){
    return(-Inf)
  }

  0
}

#' @title Extract parameters for 'JAGS' priors
#'
#' @description Extracts transformed parameters from the
#' prior part of a 'JAGS' model inside of a 'bridgesampling'
#' function (returns them as a named list)
#'
#' @param samples samples provided by bridgesampling
#' function
#' @param prior_list_parameters named list of prior distributions on model parameters
#' (not specified within the formula but that might scale the formula parameters)
#' @param formula_design_list optional formula-design metadata produced by
#' \code{JAGS_formula()}, used internally to reconstruct formula random effects
#' for bridge sampling.
#' @param model_data optional data passed to row-wise external parameter source
#' reconstruction functions during formula random-effect bridge sampling.
#'
#' @return \code{JAGS_marglik_parameters} returns a named list
#' of (transformed) posterior samples.
#'
#' @inheritParams JAGS_bridgesampling
#' @export JAGS_marglik_parameters
#' @export JAGS_marglik_parameters_formula
#' @name JAGS_marglik_parameters
NULL

#' @rdname JAGS_marglik_parameters
JAGS_marglik_parameters                <- function(samples, prior_list){

  # return empty list in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # add the resulting parameters
  parameters <- list()
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.weightfunction(samples, prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.phacking(samples, prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.bias(samples, prior_list[[i]]))

    }else if(is.prior.mixture(prior_list[[i]])){

      .JAGS_marglik_stop_unsupported_mixture(prior_list[[i]])

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.PP(samples, prior_list[[i]]))

    }else if(is.prior.factor(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.factor(samples, prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.vector(samples, prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.simple(samples, prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(parameters)
}


.JAGS_marglik_parameters.simple         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")


  parameter <- list()
  if(prior[["distribution"]] == "invgamma"){
    inv_value <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = paste0("inv_", parameter_name),
      missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters.",
      signal = TRUE
    )
    parameter[[parameter_name]] <- inv_value^-1
  }else if(prior[["distribution"]] == "point"){
    parameter[[parameter_name]] <- prior$parameters[["location"]]
  }else{
    parameter[[parameter_name]] <- samples[[ parameter_name ]]
  }

  return(parameter)
}
.JAGS_marglik_parameter_values          <- function(samples, prior, parameter_names){

  if(is.prior.point(prior)){
    return(rep(prior$parameters[["location"]], length(parameter_names)))
  }

  sample_names <- parameter_names
  if(prior[["distribution"]] == "invgamma"){
    sample_names <- paste0("inv_", parameter_names)
  }

  if(!all(sample_names %in% names(samples))){
    stop("'samples' does not contain all monitored formula prior parameters.", call. = FALSE)
  }

  if(prior[["distribution"]] == "invgamma"){
    values <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = sample_names,
      missing_message = "'samples' does not contain all monitored formula prior parameters.",
      signal = TRUE
    )
  }else{
    values <- unname(unlist(samples[sample_names], use.names = FALSE))
  }
  if(prior[["distribution"]] == "invgamma"){
    values <- values^-1
  }

  return(values)
}
.JAGS_marglik_parameters.vector         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")


  parameter <- list()
  if(prior[["distribution"]] == "dirichlet"){
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(prior$parameters[["K"]]), "]")
    eta <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = eta_names,
      missing_message = "'samples' does not contain all monitored Dirichlet prior parameters.",
      signal = TRUE
    )
    parameter[[parameter_name]] <- eta / sum(eta)
    return(parameter)
  }

  if(prior$parameters[["K"]] == 1){
    parameter_monitor_name <- parameter_name
  }else{
    parameter_monitor_name <- paste0(parameter_name, "[", 1:prior$parameters[["K"]], "]")
  }

  if(prior[["distribution"]] == "mpoint"){
    parameter[[parameter_name]] <- rep(prior$parameters[["location"]], length(parameter_monitor_name))
  }else{
    parameter[[parameter_name]] <- samples[ parameter_monitor_name ]
  }

  return(parameter)
}
.JAGS_marglik_parameters.factor         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")


  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    parameter <- list()
    if(.get_prior_factor_levels(prior) == 1){
      parameter_names <- parameter_name
    }else{
      parameter_names <- paste0(parameter_name, "[", 1:.get_prior_factor_levels(prior), "]")
    }
    parameter[[parameter_name]] <- .JAGS_marglik_parameter_values(samples, prior, parameter_names)

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    parameter <- .JAGS_marglik_parameters.vector(samples, prior, parameter_name)

  }


  return(parameter)
}
.JAGS_marglik_parameters.PP             <- function(samples, prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    parameter <- .JAGS_marglik_parameters.simple(samples, prior, "PET")
  }else if(is.prior.PEESE(prior)){
    parameter <- .JAGS_marglik_parameters.simple(samples, prior, "PEESE")
  }

  return(parameter)
}
.JAGS_marglik_parameters.weightfunction <- function(samples, prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  parameter <- list()
  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "cumulative"){

    eta     <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = paste0("eta[", seq_len(J), "]"),
      missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters.",
      signal = TRUE
    )
    std_eta <- eta / sum(eta)
    omega <- unname(rev(cumsum(rev(std_eta))))

  }else if(prior$weights$type == "independent"){

    omega <- rep(1, J)
    if(J > 1L){
      if(prior$weights$scale == "omega"){
        omega[2:J] <- samples[paste0("omega[", 2:J, "]")]
      }else if(prior$weights$scale == "log_omega"){
        omega[2:J] <- exp(samples[paste0("log_omega[", 2:J, "]")])
      }
    }
  }else if(prior$weights$type == "fixed"){

    omega <- unname(prior$weights$omega)

  }

  expansion <- .weightfunction_mapping_expansion(prior, force_one_sided = TRUE)
  parameter[["omega"]] <- unname(omega[expansion$index])

  return(parameter)
}

.bt_JAGS_marglik_positive_auxiliary_values <- function(samples,
                                                       parameter_names,
                                                       missing_message,
                                                       signal = FALSE){

  if(!all(parameter_names %in% names(samples))){
    stop(missing_message, call. = FALSE)
  }

  values <- unname(unlist(samples[parameter_names], use.names = FALSE))
  invalid <- !is.finite(values) | values <= 0
  if(any(invalid)){
    if(isTRUE(signal)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support positive auxiliary coordinate '",
        parameter_names[which(invalid)[1L]],
        "'."
      )
    }
    return(NULL)
  }

  values
}
.JAGS_marglik_parameters.phacking <- function(samples, prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  alpha <- .JAGS_marglik_parameters.simple(samples, prior$alpha, "alpha")[["alpha"]]
  constants <- phack_backend_constants(prior$form, prior$source, prior$destination, target = prior$target)
  list(
    alpha     = alpha,
    pi_null   = alpha * constants$pi_null_per_alpha,
    beta_null = alpha * constants$beta_null_per_alpha
  )
}
.JAGS_marglik_parameters.bias <- function(samples, prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)

  parameter <- list()
  if(!is.null(prior$selection)){
    parameter <- c(parameter, .JAGS_marglik_parameters.weightfunction(samples, prior$selection))
  }
  if(!is.null(prior$phacking)){
    parameter <- c(parameter, .JAGS_marglik_parameters.phacking(samples, prior$phacking))
  }

  return(parameter)
}

.bt_JAGS_bridge_compile_prior_list_evaluator <- function(prior_list){

  if(length(prior_list) == 0L){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list()
    ))
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)

  evaluators <- vector("list", length(prior_list))
  for(i in seq_along(prior_list)){
    evaluators[[i]] <- .bt_JAGS_bridge_compile_prior_evaluator(
      prior_list[[i]],
      names(prior_list)[i]
    )
  }

  list(
    log_prior = function(samples){
      marglik <- 0
      for(evaluator in evaluators){
        marglik <- marglik + evaluator$log_prior(samples)
      }
      marglik
    },
    parameters = function(samples){
      parameters <- list()
      for(evaluator in evaluators){
        parameters <- c(parameters, evaluator$parameters(samples))
      }
      parameters
    }
  )
}

.bt_JAGS_bridge_compile_prior_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(is.prior.weightfunction(prior_object)){
    return(.bt_JAGS_bridge_compile_weightfunction_evaluator(prior_object))
  }else if(is.prior.none(prior_object)){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list()
    ))
  }else if(is_prior_phacking(prior_object)){
    return(.bt_JAGS_bridge_compile_phacking_evaluator(prior_object))
  }else if(is_prior_bias(prior_object)){
    return(.bt_JAGS_bridge_compile_bias_evaluator(prior_object))
  }else if(is.prior.mixture(prior_object)){
    .JAGS_marglik_stop_unsupported_mixture(prior_object)
  }else if(is.prior.PET(prior_object) | is.prior.PEESE(prior_object)){
    return(.bt_JAGS_bridge_compile_PP_evaluator(prior_object))
  }else if(is.prior.factor(prior_object)){
    return(.bt_JAGS_bridge_compile_factor_evaluator(prior_object, parameter_name))
  }else if(is.prior.vector(prior_object)){
    return(.bt_JAGS_bridge_compile_vector_evaluator(prior_object, parameter_name))
  }else if(is.prior.simple(prior_object)){
    return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, parameter_name))
  }

  stop("Unsupported prior object.", call. = FALSE)
}

.bt_JAGS_bridge_compile_simple_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  distribution <- prior_object[["distribution"]]

  if(identical(distribution, "invgamma")){
    inv_name <- paste0("inv_", parameter_name)
    sampling_prior <- prior(
      "distribution" = "gamma",
      "parameters"   = list(
        "shape" = prior_object$parameters[["shape"]],
        "rate"  = prior_object$parameters[["scale"]]
      ),
      "truncation"   = list(
        "lower" = prior_object$truncation[["upper"]]^-1,
        "upper" = prior_object$truncation[["lower"]]^-1
      )
    )
    return(list(
      log_prior = function(samples){
        inv_value <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = inv_name,
          missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters."
        )
        if(is.null(inv_value)){
          return(-Inf)
        }
        lpdf(sampling_prior, inv_value)
      },
      parameters = function(samples){
        inv_value <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = inv_name,
          missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters.",
          signal = TRUE
        )
        parameter <- list()
        parameter[[parameter_name]] <- inv_value^-1
        parameter
      }
    ))
  }

  if(identical(distribution, "point")){
    location <- prior_object$parameters[["location"]]
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- location
        parameter
      }
    ))
  }

  list(
    log_prior = function(samples){
      lpdf(prior_object, samples[[parameter_name]])
    },
    parameters = function(samples){
      parameter <- list()
      parameter[[parameter_name]] <- samples[[parameter_name]]
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_parameter_values <- function(prior_object, parameter_names){

  force(prior_object)
  force(parameter_names)

  if(is.prior.point(prior_object)){
    location <- prior_object$parameters[["location"]]
    return(function(samples) rep(location, length(parameter_names)))
  }

  sample_names <- parameter_names
  if(identical(prior_object[["distribution"]], "invgamma")){
    sample_names <- paste0("inv_", parameter_names)
    return(function(samples){
      values <- .bt_JAGS_marglik_positive_auxiliary_values(
        samples = samples,
        parameter_names = sample_names,
        missing_message = "'samples' does not contain all monitored formula prior parameters.",
        signal = TRUE
      )
      values^-1
    })
  }

  function(samples){
    if(!all(sample_names %in% names(samples))){
      stop("'samples' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
    unname(unlist(samples[sample_names], use.names = FALSE))
  }
}

.bt_JAGS_bridge_compile_vector_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(identical(prior_object[["distribution"]], "dirichlet")){
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(prior_object$parameters[["K"]]), "]")
    alpha <- prior_object$parameters[["alpha"]]
    return(list(
      log_prior = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored Dirichlet prior parameters."
        )
        if(is.null(eta)){
          return(-Inf)
        }
        sum(stats::dgamma(eta, shape = alpha, rate = 1, log = TRUE))
      },
      parameters = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored Dirichlet prior parameters.",
          signal = TRUE
        )
        parameter <- list()
        parameter[[parameter_name]] <- eta / sum(eta)
        parameter
      }
    ))
  }

  if(prior_object$parameters[["K"]] == 1){
    parameter_monitor_names <- parameter_name
  }else{
    parameter_monitor_names <- paste0(parameter_name, "[", seq_len(prior_object$parameters[["K"]]), "]")
  }

  if(identical(prior_object[["distribution"]], "mpoint")){
    location <- prior_object$parameters[["location"]]
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- rep(location, length(parameter_monitor_names))
        parameter
      }
    ))
  }

  list(
    log_prior = function(samples){
      if(length(parameter_monitor_names) == 1L){
        lpdf(prior_object, samples[[parameter_monitor_names]])
      }else{
        lpdf(prior_object, samples[parameter_monitor_names])
      }
    },
    parameters = function(samples){
      parameter <- list()
      parameter[[parameter_name]] <- samples[parameter_monitor_names]
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_factor_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(is.prior.treatment(prior_object) | is.prior.independent(prior_object)){
    levels <- .get_prior_factor_levels(prior_object)
    if(levels == 1L){
      return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, parameter_name))
    }

    parameter_names <- paste0(parameter_name, "[", seq_len(levels), "]")
    density_evaluators <- lapply(parameter_names, function(name){
      .bt_JAGS_bridge_compile_simple_evaluator(prior_object, name)$log_prior
    })
    parameter_values <- .bt_JAGS_bridge_compile_parameter_values(prior_object, parameter_names)
    return(list(
      log_prior = function(samples){
        marglik <- 0
        for(evaluator in density_evaluators){
          marglik <- marglik + evaluator(samples)
        }
        marglik
      },
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- parameter_values(samples)
        parameter
      }
    ))
  }

  factor_prior <- prior_object
  factor_prior$parameters[["K"]] <- .get_prior_factor_levels(prior_object)
  .bt_JAGS_bridge_compile_vector_evaluator(factor_prior, parameter_name)
}

.bt_JAGS_bridge_compile_PP_evaluator <- function(prior_object){

  force(prior_object)

  if(is.prior.PET(prior_object)){
    return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, "PET"))
  }
  .bt_JAGS_bridge_compile_simple_evaluator(prior_object, "PEESE")
}

.bt_JAGS_bridge_compile_weightfunction_evaluator <- function(prior_object){

  force(prior_object)

  J <- .weightfunction_n_bins(prior_object)
  expansion <- .weightfunction_mapping_expansion(prior_object, force_one_sided = TRUE)

  if(prior_object$weights$type == "fixed"){
    omega_fixed <- unname(prior_object$weights$omega[expansion$index])
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(omega = omega_fixed)
    ))
  }

  if(prior_object$weights$type == "cumulative"){
    eta_names <- paste0("eta[", seq_len(J), "]")
    alpha <- prior_object$weights$alpha
    return(list(
      log_prior = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters."
        )
        if(is.null(eta)){
          return(-Inf)
        }
        sum(stats::dgamma(eta, shape = alpha, rate = 1, log = TRUE))
      },
      parameters = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters.",
          signal = TRUE
        )
        std_eta <- eta / sum(eta)
        omega <- unname(rev(cumsum(rev(std_eta))))
        list(omega = unname(omega[expansion$index]))
      }
    ))
  }

  omega <- rep(1, J)
  if(J == 1L){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(omega = unname(omega[expansion$index]))
    ))
  }

  if(prior_object$weights$scale == "omega"){
    omega_names <- paste0("omega[", 2:J, "]")
    weight_prior <- prior_object$weights$prior
    return(list(
      log_prior = function(samples){
        sum(mlpdf(weight_prior, samples[omega_names]))
      },
      parameters = function(samples){
        omega[2:J] <- samples[omega_names]
        list(omega = unname(omega[expansion$index]))
      }
    ))
  }

  log_omega_names <- paste0("log_omega[", 2:J, "]")
  weight_prior <- prior_object$weights$prior
  list(
    log_prior = function(samples){
      sum(mlpdf(weight_prior, samples[log_omega_names]))
    },
    parameters = function(samples){
      omega[2:J] <- exp(samples[log_omega_names])
      list(omega = unname(omega[expansion$index]))
    }
  )
}

.bt_JAGS_bridge_compile_phacking_evaluator <- function(prior_object){

  force(prior_object)

  alpha_evaluator <- .bt_JAGS_bridge_compile_simple_evaluator(prior_object$alpha, "alpha")
  constants <- phack_backend_constants(
    prior_object$form,
    prior_object$source,
    prior_object$destination,
    target = prior_object$target
  )

  list(
    log_prior = alpha_evaluator$log_prior,
    parameters = function(samples){
      alpha <- alpha_evaluator$parameters(samples)[["alpha"]]
      list(
        alpha     = alpha,
        pi_null   = alpha * constants$pi_null_per_alpha,
        beta_null = alpha * constants$beta_null_per_alpha
      )
    }
  )
}

.bt_JAGS_bridge_compile_bias_evaluator <- function(prior_object){

  force(prior_object)

  selection_backend_spec(prior_object)

  selection_evaluator <- if(!is.null(prior_object$selection)){
    .bt_JAGS_bridge_compile_weightfunction_evaluator(prior_object$selection)
  }else{
    NULL
  }
  phacking_evaluator <- if(!is.null(prior_object$phacking)){
    .bt_JAGS_bridge_compile_phacking_evaluator(prior_object$phacking)
  }else{
    NULL
  }

  list(
    log_prior = function(samples){
      marglik <- 0
      if(!is.null(selection_evaluator)){
        marglik <- marglik + selection_evaluator$log_prior(samples)
      }
      if(!is.null(phacking_evaluator)){
        marglik <- marglik + phacking_evaluator$log_prior(samples)
      }
      marglik
    },
    parameters = function(samples){
      parameter <- list()
      if(!is.null(selection_evaluator)){
        parameter <- c(parameter, selection_evaluator$parameters(samples))
      }
      if(!is.null(phacking_evaluator)){
        parameter <- c(parameter, phacking_evaluator$parameters(samples))
      }
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_formula_prior_evaluator <- function(formula_prior_list){

  if(length(formula_prior_list) == 0L){
    return(list(log_prior = function(samples) 0))
  }

  evaluators <- lapply(formula_prior_list, .bt_JAGS_bridge_compile_prior_list_evaluator)
  list(
    log_prior = function(samples){
      marglik <- 0
      for(evaluator in evaluators){
        marglik <- marglik + evaluator$log_prior(samples)
      }
      marglik
    }
  )
}

.bt_JAGS_bridge_compile_formula_random_prior_evaluator <- function(formula_design_list){

  evaluators <- list()
  if(length(formula_design_list) > 0L){
    for(design in formula_design_list){
      if(.bt_formula_design_has_random_effects(design)){
        for(random_term in design$random_effects){
          evaluators[[length(evaluators) + 1L]] <-
            .bt_JAGS_bridge_compile_random_effect_prior_evaluator(random_term)
        }
      }
    }
  }

  list(
    uses_posterior_row = length(evaluators) > 0L,
    log_prior = function(samples){
      marglik <- 0
      for(evaluator in evaluators){
        contribution <- evaluator$log_prior(samples)
        if(is.null(contribution)){
          return(-Inf)
        }
        marglik <- marglik + contribution
      }
      marglik
    }
  )
}

.bt_JAGS_bridge_compile_random_effect_prior_evaluator <- function(random_term){

  force(random_term)

  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  block_name <- random_term$block_name
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  ))
  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  scalar_rho_support_evaluator <- .bt_JAGS_bridge_compile_random_effect_scalar_rho_support(
    random_term = random_term,
    structure = structure
  )
  lkj_evaluator <- .bt_JAGS_bridge_compile_random_effect_lkj_prior(
    random_term = random_term,
    structure = structure,
    n_columns = n_columns
  )

  list(
    log_prior = function(samples){
      if(!all(z_names %in% names(samples))){
        stop(
          "Bridge samples are missing standardized latent random effects for block '",
          block_name,
          "'.",
          call. = FALSE
        )
      }

      z_values <- samples[z_names]
      if(any(is.na(z_values))){
        return(-Inf)
      }
      marglik <- sum(stats::dnorm(z_values, mean = 0, sd = 1, log = TRUE))
      if(is.na(marglik)){
        return(-Inf)
      }

      scalar_rho_support <- scalar_rho_support_evaluator(samples)
      if(!is.finite(scalar_rho_support)){
        return(scalar_rho_support)
      }
      marglik <- marglik + scalar_rho_support

      if(!is.null(lkj_evaluator)){
        marglik <- marglik + lkj_evaluator(samples)
        if(is.na(marglik)){
          return(-Inf)
        }
      }

      marglik
    }
  )
}

.bt_JAGS_bridge_compile_random_effect_scalar_rho_support <- function(random_term,
                                                                    structure){

  force(random_term)
  force(structure)

  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(function(samples) 0)
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation'.",
      call. = FALSE
    )
  }

  K <- random_term$n_columns
  distance_matrix <- if(identical(structure, "car")) correlation$distance_matrix else NULL

  function(samples){
    rho <- .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
      missing = "error",
      out_of_support = "null",
      sample_space = TRUE,
      context = "Bridge sampling random-effect metadata"
    )
    if(is.null(rho) || any(is.na(rho) | !is.finite(rho))){
      return(-Inf)
    }

    R <- .bt_random_effect_structured_correlation_matrix(
      structure = structure,
      K = K,
      rho = rho[1L],
      distance_matrix = distance_matrix
    )
    chol_ok <- tryCatch({
      chol(R)
      TRUE
    }, error = function(e) FALSE)
    if(!chol_ok){
      return(-Inf)
    }

    0
  }
}

.bt_JAGS_bridge_compile_random_effect_lkj_prior <- function(random_term,
                                                            structure,
                                                            n_columns){

  force(random_term)
  force(structure)
  force(n_columns)

  if(!identical(structure, "us") || n_columns <= 1L){
    return(NULL)
  }

  block_name <- random_term$block_name
  u_names <- .bt_random_effect_lkj_primitive_names(
    random_term,
    n_columns,
    context = "Bridge sampling random-effect metadata"
  )
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = "us",
    context = "Bridge sampling random-effect metadata"
  )
  eta <- correlation$eta
  if(!is.numeric(eta) || length(eta) != 1L || is.na(eta)){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical 'random_term$correlation$eta'.",
      call. = FALSE
    )
  }

  function(samples){
    if(!all(u_names %in% names(samples))){
      stop(
        "Bridge samples are missing LKJ primitive coordinates for block '",
        block_name,
        "'.",
        call. = FALSE
      )
    }
    u_values <- unname(samples[u_names])
    if(any(is.na(u_values) | u_values <= 0 | u_values >= 1)){
      return(-Inf)
    }

    .bt_lkj_cholesky_cpc_u_log_prior(
      u_values,
      K = n_columns,
      eta = eta
    )
  }
}

.bt_JAGS_bridge_cache_posterior_row <- function(samples, cache){

  if(isTRUE(cache)){
    posterior_row <- matrix(
      unname(samples),
      nrow = 1L,
      dimnames = list(NULL, names(samples))
    )
    attr(
      posterior_row,
      "BayesTools_random_effect_dirichlet_draw_cache"
    ) <- new.env(parent = emptyenv())
    attr(samples, "BayesTools_marglik_posterior_row") <- posterior_row
  }

  samples
}

.bt_JAGS_bridge_compile_formula_parameter_evaluator <- function(formula_list,
                                                                formula_data_list,
                                                                formula_prior_list,
                                                                formula_design_list,
                                                                model_data){

  if(length(formula_prior_list) == 0L){
    return(list(parameters = function(samples, prior_list_parameters) list()))
  }

  fixed_plans <- list()
  random_plans <- list()

  for(parameter in names(formula_prior_list)){
    formula_parameter <- if(!is.null(formula_list)) formula_list[[parameter]] else NULL
    formula_data <- if(!is.null(formula_data_list)) formula_data_list[[parameter]] else NULL
    log_intercept <- if(!is.null(formula_parameter)){
      isTRUE(attr(formula_parameter, "log(intercept)"))
    }else{
      FALSE
    }
    parameter_prior_list <- formula_prior_list[[parameter]]
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL

    if(.bt_formula_design_has_random_effects(design)){
      design <- .bt_JAGS_bridge_prepare_random_effect_allocation_design(design)
      fixed_prior_list <- .bt_JAGS_marglik_formula_fixed_priors(
        parameter_prior_list,
        parameter
      )
      source_data <- .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = formula_data,
        design = design
      )
      random_plans[[parameter]] <- list(
        parameter = parameter,
        design = design,
        formula_prior_list = parameter_prior_list,
        source_data = source_data
      )
    }else{
      fixed_prior_list <- parameter_prior_list
    }

    fixed_plans[[parameter]] <- .bt_JAGS_bridge_compile_formula_fixed_plan(
      parameter = parameter,
      formula_data = formula_data,
      formula_prior_list = fixed_prior_list,
      design = design,
      log_intercept = log_intercept
    )
  }

  list(
    parameters = function(samples, prior_list_parameters){
      parameters <- list()
      for(parameter in names(fixed_plans)){
        parameters[[parameter]] <- fixed_plans[[parameter]]$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }

      if(length(random_plans) > 0L){
        source_base <- .bt_JAGS_bridge_formula_source_base_parameters(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
        for(random_plan in random_plans){
          source_parameters <- .bt_JAGS_bridge_formula_source_parameters(
            source_base = source_base,
            formula_parameters = parameters
          )
          parameters[[random_plan$parameter]] <- parameters[[random_plan$parameter]] +
            .bt_JAGS_marglik_random_effects_value(
              samples = samples,
              design = random_plan$design,
              formula_prior_list = random_plan$formula_prior_list,
              data = random_plan$source_data,
              parameters = source_parameters
            )
        }
      }

      parameters
    }
  )
}

.bt_JAGS_bridge_prepare_random_effect_allocation_design <- function(design){

  if(!.bt_formula_design_has_random_effects(design)){
    return(design)
  }

  for(random_i in seq_along(design$random_effects)){
    design$random_effects[[random_i]] <-
      .bt_JAGS_bridge_prepare_random_effect_allocation_term(
        design$random_effects[[random_i]]
      )
  }

  design
}

.bt_JAGS_bridge_prepare_random_effect_allocation_term <- function(random_term){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(random_term)
  }

  .bt_check_random_sd_binding(binding)
  sd_component_allocation <- isTRUE(binding$true_allocation) &&
    length(binding$allocations) > 0L &&
    is.list(binding$allocations[[1L]]) &&
    identical(binding$allocations[[1L]]$target, "sd_component")

  if(!isTRUE(sd_component_allocation)){
    binding$factors <- .bt_random_effect_allocation_factors_with_plan(
      binding$factors
    )
    if(is.list(binding$factors_by_column) && length(binding$factors_by_column) > 0L){
      for(column in seq_along(binding$factors_by_column)){
        binding$factors_by_column[[column]] <-
          .bt_random_effect_allocation_factors_with_plan(
            binding$factors_by_column[[column]]
          )
      }
    }
  }
  if(is.list(binding$allocations) && length(binding$allocations) > 0L){
    for(allocation_i in seq_along(binding$allocations)){
      allocation <- binding$allocations[[allocation_i]]
      if(is.list(allocation$factors) &&
         !identical(allocation$target, "sd_component")){
        allocation$factors <- .bt_random_effect_allocation_factors_with_plan(
          allocation$factors
        )
      }
      if(is.list(allocation$parent_factors) &&
         !identical(allocation$target, "sd_component")){
        allocation$parent_factors <- .bt_random_effect_allocation_factors_with_plan(
          allocation$parent_factors
        )
      }
      binding$allocations[[allocation_i]] <- allocation
    }
  }

  random_term$sd_binding <- binding
  random_term
}

.bt_JAGS_bridge_compile_formula_fixed_plan <- function(parameter,
                                                       formula_data,
                                                       formula_prior_list,
                                                       design,
                                                       log_intercept){

  force(parameter)
  force(formula_data)
  force(formula_prior_list)
  force(log_intercept)

  if(.bt_JAGS_formula_design_can_reconstruct(design)){
    return(.bt_JAGS_bridge_compile_formula_design_plan(
      design = design,
      formula_prior_list = formula_prior_list,
      log_intercept = log_intercept
    ))
  }

  list(
    value = function(samples, prior_list_parameters){
      .JAGS_marglik_parameters_formula_get(
        samples = samples,
        parameter = parameter,
        formula_data_list = formula_data,
        formula_prior_list = formula_prior_list,
        prior_list_parameters = prior_list_parameters,
        log_intercept = log_intercept
      )
    }
  )
}

.bt_JAGS_bridge_compile_formula_design_plan <- function(design,
                                                        formula_prior_list,
                                                        log_intercept){

  force(design)
  force(formula_prior_list)

  parameter <- design$parameter
  n_rows <- nrow(design$model_matrix)
  intercept_name <- paste0(parameter, "_intercept")
  intercept_plan <- NULL
  if(intercept_name %in% names(formula_prior_list)){
    intercept_plan <- .bt_JAGS_bridge_compile_formula_intercept_plan(
      prior_object = formula_prior_list[[intercept_name]],
      parameter_name = intercept_name,
      log_intercept = isTRUE(log_intercept) || isTRUE(design$log_intercept)
    )
  }

  term_names <- setdiff(names(formula_prior_list), intercept_name)
  term_plans <- lapply(term_names, function(term_name){
    term_prior <- formula_prior_list[[term_name]]
    model_term <- sub(paste0("^", parameter, "_"), "", term_name)
    columns <- .bt_JAGS_formula_design_term_columns(design, model_term)
    .bt_JAGS_bridge_compile_formula_term_plan(
      term_name = term_name,
      term_prior = term_prior,
      term_data = design$model_matrix[, columns, drop = FALSE]
    )
  })

  list(
    value = function(samples, prior_list_parameters){
      output <- rep(0, n_rows)
      if(!is.null(intercept_plan)){
        output <- output + intercept_plan$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }
      for(term_plan in term_plans){
        output <- output + term_plan$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }
      as.vector(output)
    }
  )
}

.bt_JAGS_bridge_compile_formula_intercept_plan <- function(prior_object,
                                                           parameter_name,
                                                           log_intercept){

  force(log_intercept)
  value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
    prior_object = prior_object,
    parameter_names = parameter_name
  )
  multiply_by_evaluator <- .bt_JAGS_bridge_compile_prior_multiply_by(prior_object)

  list(
    value = function(samples, prior_list_parameters){
      value <- value_evaluator(samples)
      if(isTRUE(log_intercept)){
        value <- log(value)
      }
      multiply_by_evaluator(prior_list_parameters) * value
    }
  )
}

.bt_JAGS_bridge_compile_formula_term_plan <- function(term_name,
                                                      term_prior,
                                                      term_data){

  force(term_name)
  force(term_prior)
  force(term_data)

  n_rows <- nrow(term_data)
  multiply_by_evaluator <- .bt_JAGS_bridge_compile_prior_multiply_by(term_prior)

  if(is.prior.point(term_prior) && !is.prior.factor(term_prior)){
    location <- term_prior[["parameters"]][["location"]]
    term_vector <- as.vector(term_data)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) * location * term_vector
      }
    ))
  }

  if(is.prior.point(term_prior) && is.prior.factor(term_prior)){
    location <- term_prior[["parameters"]][["location"]]
    levels <- .get_prior_factor_levels(term_prior)
    if(levels == 1L){
      term_vector <- as.vector(term_data)
      return(list(
        value = function(samples, prior_list_parameters){
          multiply_by_evaluator(prior_list_parameters) * location * term_vector
        }
      ))
    }

    term_values <- rep(location, levels)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          as.vector(term_data %*% term_values)
      }
    ))
  }

  if(is.prior.factor(term_prior)){
    levels <- .get_prior_factor_levels(term_prior)
    if(levels == 1L){
      value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
        prior_object = term_prior,
        parameter_names = term_name
      )
      term_vector <- as.vector(term_data)
      return(list(
        value = function(samples, prior_list_parameters){
          multiply_by_evaluator(prior_list_parameters) *
            value_evaluator(samples) * term_vector
        }
      ))
    }

    parameter_names <- paste0(term_name, "[", seq_len(levels), "]")
    value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
      prior_object = term_prior,
      parameter_names = parameter_names
    )
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          as.vector(term_data %*% value_evaluator(samples))
      }
    ))
  }

  if(is.prior.simple(term_prior)){
    value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
      prior_object = term_prior,
      parameter_names = term_name
    )
    term_vector <- as.vector(term_data)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          value_evaluator(samples) * term_vector
      }
    ))
  }

  list(
    value = function(samples, prior_list_parameters) rep(0, n_rows)
  )
}

.bt_JAGS_bridge_compile_prior_multiply_by <- function(prior_object){

  multiply_by <- attr(prior_object, "multiply_by")
  if(is.null(multiply_by)){
    return(function(prior_list_parameters) 1)
  }
  if(is.numeric(multiply_by)){
    force(multiply_by)
    return(function(prior_list_parameters) multiply_by)
  }

  force(multiply_by)
  function(prior_list_parameters){
    prior_list_parameters[[multiply_by]]
  }
}

.bt_JAGS_bridge_formula_source_base_parameters <- function(samples,
                                                           prior_list_parameters){

  out <- as.list(samples)
  if(length(prior_list_parameters) > 0L){
    out[names(prior_list_parameters)] <- prior_list_parameters
  }

  out
}

.bt_JAGS_bridge_formula_source_parameters <- function(source_base,
                                                      formula_parameters){

  out <- source_base
  if(length(formula_parameters) > 0L){
    out[names(formula_parameters)] <- formula_parameters
  }

  out
}

# .JAGS_marglik_parameters.spike_and_slab <- function(samples, prior, parameter_name){
#
#   .check_prior(prior)
#   if(!is.prior.spike_and_slab(prior))
#     stop("improper prior provided")
#   check_char(parameter_name, "parameter_name")
#
#   parameter <- list()
#   parameter[paste0(parameter_name, "_variable")]  <- .JAGS_marglik_parameters.simple(samples, prior[["variable"]],  paste0(parameter_name, "_variable"))
#   if(!is.prior.point(prior[[parameter_name]][["inclusion"]])){
#     parameter[paste0(parameter_name, "_inclusion")] <- .JAGS_marglik_parameters.simple(samples, prior[["inclusion"]], paste0(parameter_name, "_inclusion"))
#   }
#
#   return(parameter)
# }

#' @rdname JAGS_marglik_parameters
JAGS_marglik_parameters_formula      <- function(samples, formula_list, formula_data_list, formula_prior_list, prior_list_parameters,
                                                 formula_design_list = NULL,
                                                 model_data = NULL){

  # return empty list in case that no prior was specified
  if(length(formula_prior_list) == 0){
    return(list())
  }

  parameters <- list()
  random_parameters <- character()

  for(parameter in names(formula_prior_list)){
    # check for log(intercept) attribute on the formula
    formula_parameter <- if(!is.null(formula_list)) formula_list[[parameter]] else NULL
    log_intercept <- if(!is.null(formula_parameter)) isTRUE(attr(formula_parameter, "log(intercept)")) else FALSE
    parameter_prior_list <- formula_prior_list[[parameter]]
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL
    if(.bt_formula_design_has_random_effects(design)){
      parameter_prior_list <- .bt_JAGS_marglik_formula_fixed_priors(parameter_prior_list, parameter)
      random_parameters <- c(random_parameters, parameter)
    }
    if(.bt_JAGS_formula_design_can_reconstruct(design)){
      parameters[[parameter]] <- .bt_JAGS_marglik_parameters_formula_design(
        samples = samples,
        design = design,
        formula_prior_list = parameter_prior_list,
        prior_list_parameters = prior_list_parameters,
        log_intercept = log_intercept
      )
    }else{
      parameters[[parameter]] <- .JAGS_marglik_parameters_formula_get(samples, parameter, formula_data_list[[parameter]], parameter_prior_list, prior_list_parameters, log_intercept)
    }
  }

  for(parameter in random_parameters){
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL
    if(.bt_formula_design_has_random_effects(design)){
      source_parameters <- .bt_JAGS_marglik_parameter_source_parameters(
        samples = samples,
        prior_list_parameters = prior_list_parameters,
        formula_parameters = parameters
      )
      source_formula_data <- if(!is.null(formula_data_list)){
        formula_data_list[[parameter]]
      }else{
        NULL
      }
      source_data <- .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = source_formula_data,
        design = design
      )
      parameters[[parameter]] <- parameters[[parameter]] +
        .bt_JAGS_marglik_random_effects_value(
          samples = samples,
          design = design,
          formula_prior_list = formula_prior_list[[parameter]],
          data = source_data,
          parameters = source_parameters
        )
    }
  }

  return(parameters)
}

.bt_JAGS_marglik_parameter_source_parameters <- function(samples,
                                                         prior_list_parameters,
                                                         formula_parameters){

  out <- as.list(samples)
  if(length(prior_list_parameters) > 0L){
    out[names(prior_list_parameters)] <- prior_list_parameters
  }
  if(length(formula_parameters) > 0L){
    out[names(formula_parameters)] <- formula_parameters
  }

  out
}

.bt_JAGS_marglik_parameter_source_data <- function(model_data,
                                                   formula_data,
                                                   design){

  out <- list()
  out <- .bt_JAGS_marglik_merge_source_data(out, formula_data)
  out <- .bt_JAGS_marglik_merge_source_data(out, model_data)
  if(inherits(design, "BayesTools_formula_design") &&
     !is.null(design$source_data)){
    out <- .bt_JAGS_marglik_fill_source_data(out, design$source_data)
  }

  out
}

.bt_JAGS_marglik_fill_source_data <- function(out, data){

  if(is.null(data)){
    return(out)
  }
  data_list <- if(is.data.frame(data)){
    as.list(data)
  }else if(is.list(data)){
    data
  }else{
    return(out)
  }
  data_names <- names(data_list)
  if(is.null(data_names)){
    return(out)
  }
  keep <- !is.na(data_names) & nzchar(data_names)
  if(!any(keep)){
    return(out)
  }
  data_list <- data_list[keep]
  missing <- setdiff(names(data_list), names(out))
  out[missing] <- data_list[missing]

  out
}

.bt_JAGS_marglik_merge_source_data <- function(out, data){

  if(is.null(data)){
    return(out)
  }
  data_list <- if(is.data.frame(data)){
    as.list(data)
  }else if(is.list(data)){
    data
  }else{
    return(out)
  }
  data_names <- names(data_list)
  if(is.null(data_names)){
    return(out)
  }
  keep <- !is.na(data_names) & nzchar(data_names)
  if(!any(keep)){
    return(out)
  }
  data_list <- data_list[keep]
  overlap <- intersect(names(out), names(data_list))
  for(name in overlap){
    if(!.bt_JAGS_marglik_source_data_equal(out[[name]], data_list[[name]])){
      stop(
        "JAGS_bridgesampling() row-indexed source reconstruction received ",
        "conflicting data for variable '", name, "'. Supply formula data, fitted ",
        "formula metadata, and model data with matching row-aligned values, or ",
        "remove the duplicate variable from the non-formula data.",
        call. = FALSE
      )
    }
  }
  out[names(data_list)] <- data_list

  out
}

.bt_JAGS_marglik_source_data_equal <- function(x, y){

  if(identical(x, y)){
    return(TRUE)
  }

  isTRUE(all.equal(
    x,
    y,
    check.attributes = FALSE,
    tolerance = 1e-12
  ))
}

.bt_JAGS_formula_design_can_reconstruct <- function(design){

  inherits(design, "BayesTools_formula_design") &&
    !is.null(design$parameter) &&
    !is.null(design$model_matrix) &&
    !is.null(design$assign) &&
    !is.null(design$model_terms)
}

.bt_JAGS_marglik_parameters_formula_design <- function(samples, design,
                                                       formula_prior_list,
                                                       prior_list_parameters,
                                                       log_intercept = FALSE){

  parameter <- design$parameter
  output <- rep(0, nrow(design$model_matrix))
  if(length(formula_prior_list) == 0L){
    return(output)
  }

  intercept_name <- paste0(parameter, "_intercept")
  if(intercept_name %in% names(formula_prior_list)){
    intercept_prior <- formula_prior_list[[intercept_name]]
    intercept_value <- .JAGS_marglik_parameter_values(
      samples,
      intercept_prior,
      intercept_name
    )
    if(isTRUE(log_intercept) || isTRUE(design$log_intercept)){
      intercept_value <- log(intercept_value)
    }
    output <- output + .bt_JAGS_marglik_prior_multiply_by(
      intercept_prior,
      prior_list_parameters
    ) * intercept_value
  }

  remaining_terms <- setdiff(names(formula_prior_list), intercept_name)
  for(term in remaining_terms){
    term_prior <- formula_prior_list[[term]]
    model_term <- sub(paste0("^", parameter, "_"), "", term)
    columns <- .bt_JAGS_formula_design_term_columns(design, model_term)
    term_data <- design$model_matrix[, columns, drop = FALSE]
    multiply_by <- .bt_JAGS_marglik_prior_multiply_by(
      term_prior,
      prior_list_parameters
    )

    if(is.prior.point(term_prior) && !is.prior.factor(term_prior)){
      output <- output + multiply_by * term_prior[["parameters"]][["location"]] * as.vector(term_data)

    }else if(is.prior.point(term_prior) && is.prior.factor(term_prior)){
      if(.get_prior_factor_levels(term_prior) == 1){
        output <- output + multiply_by * term_prior[["parameters"]][["location"]] * as.vector(term_data)
      }else{
        output <- output + multiply_by * as.vector(term_data %*% rep(term_prior[["parameters"]][["location"]], .get_prior_factor_levels(term_prior)))
      }

    }else if(is.prior.factor(term_prior)){
      if(.get_prior_factor_levels(term_prior) == 1){
        term_value <- .JAGS_marglik_parameter_values(samples, term_prior, term)
        output <- output + multiply_by * term_value * as.vector(term_data)
      }else{
        term_names <- paste0(term, "[", 1:.get_prior_factor_levels(term_prior), "]")
        term_values <- .JAGS_marglik_parameter_values(samples, term_prior, term_names)
        output <- output + multiply_by * as.vector(term_data %*% term_values)
      }

    }else if(is.prior.simple(term_prior)){
      term_value <- .JAGS_marglik_parameter_values(samples, term_prior, term)
      output <- output + multiply_by * term_value * as.vector(term_data)
    }
  }

  as.vector(output)
}

.bt_JAGS_formula_design_term_columns <- function(design, model_term){

  term_index <- match(model_term, design$model_terms)
  if(is.na(term_index)){
    stop(
      "Stored formula design for parameter '",
      design$parameter,
      "' is missing model term '",
      model_term,
      "'.",
      call. = FALSE
    )
  }

  columns <- which(design$assign == (term_index - 1L))
  if(length(columns) == 0L){
    stop(
      "Stored formula design for parameter '",
      design$parameter,
      "' is missing model-matrix columns for term '",
      model_term,
      "'.",
      call. = FALSE
    )
  }

  columns
}

.bt_JAGS_marglik_prior_multiply_by <- function(prior, prior_list_parameters){

  multiply_by <- attr(prior, "multiply_by")
  if(is.null(multiply_by)){
    return(1)
  }
  if(is.numeric(multiply_by)){
    return(multiply_by)
  }

  prior_list_parameters[[multiply_by]]
}

.bt_JAGS_marglik_formula_fixed_priors <- function(formula_prior_list, parameter){

  if(length(formula_prior_list) == 0L){
    return(formula_prior_list)
  }

  random_prefix <- paste0(parameter, "__xREx__")
  is_random <- startsWith(names(formula_prior_list), random_prefix) |
    vapply(formula_prior_list, function(prior){
      .bt_is_random_effect_prior(prior, include_summary = FALSE)
    }, logical(1))

  formula_prior_list[!is_random]
}

.bt_JAGS_marglik_random_effects_value <- function(samples, design,
                                                  formula_prior_list,
                                                  data = NULL,
                                                  parameters = NULL){

  output <- rep(0, nrow(design$model_matrix))
  for(random_term in design$random_effects){
    output <- output + .bt_JAGS_marglik_random_effect_value(
      samples = samples,
      random_term = random_term,
      prior_list = formula_prior_list,
      data = data,
      parameters = parameters
    )
  }

  output
}

.bt_JAGS_marglik_random_effect_value <- function(samples, random_term,
                                                 prior_list,
                                                 data = NULL,
                                                 parameters = NULL){

  .bt_JAGS_marglik_check_random_effect_dirichlet_samples(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(.bt_JAGS_marglik_random_effect_row_indexed_value(
      samples = samples,
      random_term = random_term,
      prior_list = prior_list,
      data = data,
      parameters = parameters
    ))
  }

  model_matrix <- random_term$model_matrix
  group_map <- random_term$group_map
  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns

  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }
  z <- matrix(
    unname(samples[as.vector(z_names)]),
    nrow = n_groups,
    ncol = n_columns
  )
  z_draws <- lapply(seq_len(n_columns), function(column){
    matrix(z[, column], nrow = 1L)
  })

  sd_values <- .bt_JAGS_marglik_random_effect_sd_values(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )
  L <- .bt_JAGS_marglik_random_effect_cholesky(
    samples = samples,
    random_term = random_term
  )

  contribution <- .bt_random_effect_contribution_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = matrix(sd_values, nrow = 1L),
    cholesky = array(L, dim = c(1L, n_columns, n_columns))
  )
  as.vector(contribution[, 1L])
}

.bt_JAGS_marglik_check_random_effect_dirichlet_samples <- function(
    samples, random_term, prior_list){

  parameter_names <- .bt_JAGS_marglik_random_effect_allocation_parameters(
    random_term
  )
  for(parameter_name in parameter_names){
    if(!parameter_name %in% names(prior_list)){
      next
    }
    prior <- prior_list[[parameter_name]]
    if(!is.prior.simplex(prior) || !identical(prior$distribution, "dirichlet")){
      next
    }
    K <- prior$parameters[["K"]]
    weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
    if(all(weight_names %in% names(samples)) && all(eta_names %in% names(samples))){
      stop(
        "Bridge samples contain both normalized Dirichlet allocation coordinates ",
        "and auxiliary eta coordinates for parameter '",
        parameter_name,
        "'. Use the auxiliary eta bridge coordinates only.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_JAGS_marglik_random_effect_allocation_parameters <- function(random_term){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(character())
  }
  .bt_check_random_sd_binding(binding)

  out <- character()
  collect_factors <- function(factors){
    if(!is.list(factors) || length(factors) == 0L){
      return(character())
    }
    vapply(factors, function(factor){
      if(is.list(factor) &&
         is.character(factor$weight_name) &&
         length(factor$weight_name) == 1L &&
         !is.na(factor$weight_name)){
        return(factor$weight_name)
      }
      NA_character_
    }, character(1))
  }
  out <- c(out, collect_factors(binding$factors))
  if(is.list(binding$factors_by_column)){
    for(factors in binding$factors_by_column){
      out <- c(out, collect_factors(factors))
    }
  }
  if(is.list(binding$allocations)){
    for(allocation in binding$allocations){
      if(!is.list(allocation)){
        next
      }
      if(is.character(allocation$weight_name) &&
         length(allocation$weight_name) == 1L &&
         !is.na(allocation$weight_name)){
        out <- c(out, allocation$weight_name)
      }
      out <- c(out, collect_factors(allocation$factors))
      out <- c(out, collect_factors(allocation$parent_factors))
    }
  }

  unique(out[!is.na(out) & nzchar(out)])
}

.bt_JAGS_marglik_random_effect_row_indexed_value <- function(samples,
                                                             random_term,
                                                             prior_list,
                                                             data = NULL,
                                                             parameters = NULL){

  model_matrix <- random_term$model_matrix
  group_map <- random_term$group_map
  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns

  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }
  z <- matrix(
    unname(samples[as.vector(z_names)]),
    nrow = n_groups,
    ncol = n_columns
  )
  z_draws <- lapply(seq_len(n_columns), function(column){
    matrix(z[, column], nrow = 1L)
  })

  L <- .bt_JAGS_marglik_random_effect_cholesky(
    samples = samples,
    random_term = random_term
  )
  unit_columns <- .bt_random_effect_column_contributions_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = matrix(1, nrow = 1L, ncol = n_columns),
    cholesky = array(L, dim = c(1L, n_columns, n_columns))
  )

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  source_draws <- .bt_JAGS_marglik_row_indexed_external_sd_source_draws(
    random_term = random_term,
    n_rows = nrow(model_matrix),
    posterior = posterior,
    data = data,
    parameters = parameters,
    context = "Bridge-sampling reconstruction"
  )
  if(any(!is.finite(source_draws) | source_draws < 0)){
    .bt_JAGS_marglik_out_of_support(
      "Bridge samples contain out-of-support row-indexed external SD source values for block '",
      random_term$block_name,
      "'."
    )
  }
  column_allocation_draws <- .bt_JAGS_marglik_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = n_columns
  )
  if(!is.null(column_allocation_draws)){
    if(any(!is.finite(column_allocation_draws) | column_allocation_draws < 0)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support row-indexed external SD allocation values for block '",
        random_term$block_name,
        "'."
      )
    }
    return(as.vector(.bt_random_effect_apply_row_indexed_source_to_unit_columns(
      unit_columns = unit_columns,
      source_draws = source_draws,
      allocation_draws = column_allocation_draws
    )[, 1L]))
  }
  allocation_draws <- .bt_JAGS_marglik_row_indexed_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list
  )
  if(any(!is.finite(allocation_draws) | allocation_draws < 0)){
    .bt_JAGS_marglik_out_of_support(
      "Bridge samples contain out-of-support row-indexed external SD allocation values for block '",
      random_term$block_name,
      "'."
    )
  }

  as.vector(.bt_random_effect_apply_row_indexed_source_to_unit_columns(
    unit_columns = unit_columns,
    source_draws = source_draws,
    allocation_draws = allocation_draws
  )[, 1L])
}

.bt_JAGS_marglik_row_indexed_column_allocation_draws <- function(
    random_term, posterior, prior_list, n_columns){

  tryCatch(
    .bt_random_effect_row_indexed_column_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list,
      n_columns = n_columns
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
}

.bt_JAGS_marglik_row_indexed_allocation_draws <- function(random_term,
                                                          posterior,
                                                          prior_list){

  tryCatch(
    .bt_random_effect_row_indexed_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
}

.bt_JAGS_marglik_random_effect_support_error <- function(error){

  inherits(error, "BayesTools_random_effect_allocation_out_of_support")
}

.bt_JAGS_marglik_row_indexed_external_sd_source_draws <- function(random_term,
                                                                  n_rows,
                                                                  posterior,
                                                                  data = NULL,
                                                                  parameters = NULL,
  context = "Bridge-sampling reconstruction"){

  source <- .bt_random_effect_row_indexed_source(random_term)
  source_parameter <- source$source
  source_label <- .bt_random_sd_source_label(source)
  source_names <- .bt_parameter_source_row_names(source_parameter, n_rows)
  values_function <- .bt_parameter_source_values_function(source_parameter)
  if(is.null(values_function)){
    return(.bt_random_effect_row_indexed_source_draws(
      random_term = random_term,
      n_rows = n_rows,
      posterior = posterior,
      data = data,
      parameters = parameters,
      context = context
    ))
  }
  present <- intersect(source_names, colnames(posterior))
  if(length(present) > 0L){
    stop(
      context, " for source '", source_label,
      "' is ambiguous: the source provides a values function and the posterior ",
      "also contains sampled row-source column(s): ",
      paste0("'", present[seq_len(min(3L, length(present)))], "'",
             collapse = ", "),
      if(length(present) > 3L) ", ..." else "",
      ". Use one row-source reconstruction path only.",
      call. = FALSE
    )
  }
  if(!is.matrix(posterior)){
    stop("'posterior' must be a matrix.", call. = FALSE)
  }
  check_int(n_rows, "n_rows", lower = 1, allow_NA = FALSE)

  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_rows)
  for(draw in seq_len(nrow(posterior))){
    draw_parameters <- .bt_JAGS_marglik_row_indexed_source_parameters(
      posterior = posterior,
      draw = draw,
      parameters = parameters
    )
    values <- tryCatch(
      values_function(
        parameters = draw_parameters,
        data = data,
        n_rows = n_rows
      ),
      error = function(e)e
    )
    if(inherits(values, "error")){
      stop(
        context, " for source '", source_label,
        "' failed: ", conditionMessage(values),
        call. = FALSE
      )
    }
    if(!is.numeric(values)){
      stop(
        context, " for source '", source_label,
        "' must return a numeric vector.",
        call. = FALSE
      )
    }
    values <- as.numeric(values)
    if(length(values) != n_rows){
      stop(
        context, " for source '", source_label,
        "' must return a numeric vector of length ", n_rows, ".",
        call. = FALSE
      )
    }
    if(anyNA(values)){
      .bt_JAGS_marglik_out_of_support(
        context, " for source '", source_label,
        "' returned missing values."
      )
    }
    out[draw, ] <- values
  }

  colnames(out) <- source_names
  out
}

.bt_JAGS_marglik_row_indexed_source_parameters <- function(posterior,
                                                           draw,
                                                           parameters = NULL){

  if(nrow(posterior) == 1L &&
     !is.null(parameters) &&
     is.list(parameters) &&
     all(colnames(posterior) %in% names(parameters))){
    return(parameters)
  }

  .bt_parameter_source_draw_parameters(
    posterior = posterior,
    draw = draw,
    parameters = parameters
  )
}

.bt_JAGS_marglik_random_effect_sd_values <- function(samples, random_term,
                                                     prior_list){

  .bt_JAGS_marglik_check_random_effect_dirichlet_samples(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    .bt_JAGS_marglik_row_indexed_external_sd_stop(random_term)
  }

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  sd_draws <- tryCatch(
    .bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = random_term$n_columns,
      posterior = posterior,
      prior_list = prior_list
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
  if(is.null(sd_draws)){
    .bt_JAGS_marglik_explain_random_effect_sd_missing(samples, random_term, prior_list)
    stop(
      "Random-effect SD metadata are incomplete for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  unname(sd_draws[1L, ])
}

.bt_JAGS_marglik_row_indexed_external_sd_stop <- function(random_term){

  stop(
    "Random-effect block '",
    random_term$block_name,
    "' uses row-indexed external SD source '",
    .bt_random_effect_external_sd_source_label(random_term),
    "'; scalar random-effect SD values are not defined for row-indexed sources.",
    call. = FALSE
  )
}

.bt_JAGS_marglik_explain_random_effect_sd_missing <- function(samples,
                                                              random_term,
                                                              prior_list){

  binding <- random_term$sd_binding
  if(!is.null(binding) && isTRUE(binding$true_allocation) &&
     length(binding$allocations) > 0L){
    allocation <- binding$allocations[[1L]]
    if(!.bt_JAGS_marglik_random_effect_parameter_available(
      samples = samples,
      parameter_name = .bt_random_sd_binding_source_name(allocation$source),
      prior_list = prior_list
    )){
      stop("'posterior' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
    factors <- if(identical(allocation$target, "sd_component")){
      .bt_random_effect_allocation_parent_factors_metadata(allocation)
    }else{
      .bt_random_effect_allocation_factors_metadata(allocation)
    }
    missing_factor <- vapply(factors, function(factor){
      !.bt_JAGS_marglik_random_effect_dirichlet_available(
        samples = samples,
        parameter_name = factor$weight_name,
        prior_list = prior_list
      )
    }, logical(1))
    if(any(missing_factor)){
      stop(
        "Bridge samples are missing Dirichlet allocation coordinates for parameter '",
        factors[[which(missing_factor)[1L]]]$weight_name,
        "'.",
        call. = FALSE
      )
    }
    if(identical(allocation$target, "sd_component") &&
       !.bt_JAGS_marglik_random_effect_dirichlet_available(
         samples = samples,
         parameter_name = allocation$weight_name,
         prior_list = prior_list
       )){
      stop(
        "Bridge samples are missing Dirichlet allocation coordinates for parameter '",
        allocation$weight_name,
        "'.",
        call. = FALSE
      )
    }
  }else{
    sd_names <- random_term$sd_parameter_names
    if(is.null(sd_names) || length(sd_names) != random_term$n_columns || any(is.na(sd_names))){
      return(invisible(FALSE))
    }
    missing_sd <- !vapply(sd_names, function(parameter_name){
      .bt_JAGS_marglik_random_effect_parameter_available(
        samples = samples,
        parameter_name = parameter_name,
        prior_list = prior_list
      )
    }, logical(1))
    if(any(missing_sd)){
      stop("'posterior' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
  }

  invisible(FALSE)
}

.bt_JAGS_marglik_random_effect_parameter_available <- function(samples,
                                                               parameter_name,
                                                               prior_list){

  if(parameter_name %in% names(samples)){
    return(TRUE)
  }

  prior_name <- sub("\\[[0-9]+\\]$", "", parameter_name)
  if(!prior_name %in% names(prior_list)){
    return(FALSE)
  }

  is.prior.point(prior_list[[prior_name]])
}

.bt_JAGS_marglik_random_effect_dirichlet_available <- function(samples,
                                                               parameter_name,
                                                               prior_list){

  if(!parameter_name %in% names(prior_list)){
    return(FALSE)
  }
  prior <- prior_list[[parameter_name]]
  if(!is.prior.simplex(prior) || !identical(prior$distribution, "dirichlet")){
    return(FALSE)
  }

  K <- prior$parameters[["K"]]
  weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
  eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
  all(weight_names %in% names(samples)) || all(eta_names %in% names(samples))
}

.bt_JAGS_marglik_random_effect_cholesky <- function(samples, random_term){

  n_columns <- random_term$n_columns
  L <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
    sample_space = TRUE
  )
  if(is.null(L)){
    if(.bt_JAGS_marglik_random_effect_correlation_sample_available(samples, random_term)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support random-effect correlation coordinates for block '",
        random_term$block_name,
        "'."
      )
    }
    stop(
      "Bridge samples are missing random-effect correlation coordinates for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  L[1L, , ]
}

.bt_JAGS_marglik_random_effect_correlation_sample_available <- function(samples,
                                                                        random_term){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Bridge sampling random-effect metadata"
  )
  if(structure %in% c("diag", "id")){
    return(FALSE)
  }
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation)){
    return(FALSE)
  }
  sample_names <- character()
  if(identical(correlation$type, "rho")){
    sample_names <- c(correlation$rho_name, correlation$sample_name)
  }else if(identical(correlation$type, "lkj")){
    sample_names <- .bt_random_effect_lkj_primitive_names(
      random_term,
      random_term$n_columns,
      context = "Bridge sampling random-effect metadata"
    )
  }
  sample_names <- sample_names[!is.na(sample_names) & nzchar(sample_names)]

  any(sample_names %in% names(samples))
}

.bt_JAGS_marglik_out_of_support <- function(...){

  stop(structure(
    list(message = paste0(...), call = NULL),
    class = c("BayesTools_marglik_out_of_support", "error", "condition")
  ))
}

.bt_JAGS_marglik_random_effect_rho <- function(samples, random_term){

  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
    sample_space = TRUE,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(rho)){
    stop(
      "Bridge samples are missing or invalid scalar correlation coordinates for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  unname(rho[1L])
}

.bt_JAGS_marglik_random_effect_posterior_row <- function(samples){

  cached <- attr(samples, "BayesTools_marglik_posterior_row", exact = TRUE)
  if(!is.null(cached)){
    return(cached)
  }

  matrix(
    unname(samples),
    nrow = 1L,
    dimnames = list(NULL, names(samples))
  )
}

.JAGS_marglik_parameters_formula_get <- function(samples, parameter, formula_data_list, formula_prior_list, prior_list_parameters, log_intercept = FALSE){

  formula_terms            <- names(formula_prior_list)
  names(formula_data_list) <- sub(paste0("^", parameter, "_data_"), paste0(parameter, "_"), names(formula_data_list))

  # start with intercept
  if(sum(formula_terms == paste0(parameter, "_intercept")) == 1){

    # check for scaling factors
    if(!is.null(attr(formula_prior_list[[paste0(parameter, "_intercept")]], "multiply_by"))){
      if(is.numeric(attr(formula_prior_list[[paste0(parameter, "_intercept")]], "multiply_by"))){
        multiply_by <- attr(formula_prior_list[[paste0(parameter, "_intercept")]], "multiply_by")
      }else{
        multiply_by <- prior_list_parameters[[attr(formula_prior_list[[paste0(parameter, "_intercept")]], "multiply_by")]]
      }
    }else{
      multiply_by <- 1
    }

    intercept_prior <- formula_prior_list[[paste0(parameter, "_intercept")]]
    intercept_value <- .JAGS_marglik_parameter_values(samples, intercept_prior, paste0(parameter, "_intercept"))
    # apply log transformation if log(intercept) attribute is set
    if(log_intercept){
      intercept_value <- log(intercept_value)
    }
    output <- multiply_by * rep(intercept_value, formula_data_list[[paste0("N_", parameter)]])

  }else{
    output <- rep(0, formula_data_list[[paste0("N_", parameter)]])
  }

  # add the remaining terms
  remaining_terms <- formula_terms[formula_terms != paste0(parameter, "_intercept")]
  if(length(remaining_terms) > 0){
    for(term in remaining_terms){

      # check for scaling factors
      if(!is.null(attr(formula_prior_list[[term]], "multiply_by"))){
        if(is.numeric(attr(formula_prior_list[[term]], "multiply_by"))){
          multiply_by <- attr(formula_prior_list[[term]], "multiply_by")
        }else{
          multiply_by <- prior_list_parameters[[attr(formula_prior_list[[term]], "multiply_by")]]
        }
      }else{
        multiply_by <- 1
      }


      if(is.prior.point(formula_prior_list[[term]]) && !is.prior.factor(formula_prior_list[[term]])){

        output <- output + multiply_by * formula_prior_list[[term]][["parameters"]][["location"]] * formula_data_list[[term]]

      }else if(is.prior.point(formula_prior_list[[term]]) && is.prior.factor(formula_prior_list[[term]])){

        if(.get_prior_factor_levels(formula_prior_list[[term]]) == 1){
          output <- output + multiply_by * formula_prior_list[[term]][["parameters"]][["location"]] * formula_data_list[[term]]
        }else{
          output <- output + multiply_by * formula_data_list[[term]] %*% rep(formula_prior_list[[term]][["parameters"]][["location"]], .get_prior_factor_levels(formula_prior_list[[term]]))
        }

      }else if(is.prior.factor(formula_prior_list[[term]])){

        if(.get_prior_factor_levels(formula_prior_list[[term]]) == 1){
          term_value <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term)
          output     <- output + multiply_by * term_value * formula_data_list[[term]]
        }else{
          term_names  <- paste0(term,"[", 1:.get_prior_factor_levels(formula_prior_list[[term]]), "]")
          term_values <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term_names)
          output      <- output + multiply_by * formula_data_list[[term]] %*% term_values
        }


      }else if(is.prior.simple(formula_prior_list[[term]])){

        term_value <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term)
        output     <- output + multiply_by * term_value * formula_data_list[[term]]

      }

    }
  }


  return(as.vector(output))
}
