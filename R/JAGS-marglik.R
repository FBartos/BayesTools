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
#' @param bridge_context whether the \code{log_posterior} callback receives a
#' read-only \code{bridge_context} argument. Defaults to \code{FALSE}, preserving
#' the historical \code{log_posterior(parameters, data, ...)} call. The context
#' contains the current independent bridge state and BayesTools-resolved
#' deterministic formula/random-effect nodes; it is not necessarily an
#' original MCMC posterior row.
#' @param formula_random_prior_list optional named list of `prior_random()`
#' objects for random effects in `formula_list`. Bridge sampling for formula
#' random effects requires the `prior_random()` interface because the
#' stochastic bridge coordinates are the standardized latent effects and
#' correlation primitives. For \code{BayesTools_fit} objects with stored
#' formula-design metadata, this can be omitted unless formula inputs are being
#' supplied for a consistency check.
#' @param formula_random_effects_compile_list optional named list of
#' `random_effects_compile()` objects. When formula inputs are supplied for
#' bridge-sampling rebuild/validation, this must match the fitted
#' random-effect compilation policy; otherwise a fitted marginalized model would
#' not be rebuilt as the same model.
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
#' When `bridge_context = TRUE`, the callback receives an object of class
#' `BayesTools_bridge_context` with fields `state`, `state_matrix`, `nodes`,
#' `prior_parameters`, `formula_prior_parameters`, `formula_parameters`,
#' `add_parameters`, `random`, `node_info`, and `metadata`. The `state` field
#' contains the independent bridge coordinates. The `nodes` field is a flat
#' named numeric vector that also includes deterministic BayesTools-resolved
#' nodes such as normalized Dirichlet allocation weights reconstructed from
#' `prior_par_eta_*` coordinates. The `node_info` field records the owner and
#' role of each exposed node. The `random` field contains structured
#' random-effect block state including resolved SDs, allocation weights,
#' correlations, Cholesky factors, covariance matrices where row-invariant, and
#' row-indexed SD source values where applicable. Each random block is keyed by
#' formula parameter and block name, and contains `block_name`, `compile_mode`,
#' `dimensions`, `levels`, `scale`, `allocation`, `correlation`, `covariance`,
#' `latent`, and `nodes` fields.
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
                                formula_random_effects_compile_list = NULL,
                                bridge_context = FALSE,
                                maxiter = 10000, silent = TRUE, ...){

  ### check input
  check_bool(bridge_context, "bridge_context", allow_NA = FALSE)
  check_bool(silent, "silent")
  check_int(maxiter, "maxiter", lower = 1)

  formula_context <- .bt_JAGS_bridge_formula_context(
    fit = fit,
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list
  )
  formula_design_list <- formula_context$formula_design_list
  formula_list <- formula_context$formula_list
  formula_data_list <- formula_context$formula_data_list
  formula_prior_list <- formula_context$formula_prior_list
  .bt_JAGS_bridge_check_no_allocation_inclusion(formula_design_list)

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
                                 bridge_context,
                                 formula_design_list,
                                 formula_data_list,
                                 formula_prior_list,
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
    evaluated <- tryCatch({
      prior_parameters <- bridge_prior_evaluator$parameters(samples.row)
      formula_prior_parameters <- bridge_formula_prior_evaluator$parameters(samples.row)
      formula_parameters <- bridge_formula_parameter_evaluator$parameters(samples.row, prior_parameters)
      parameters <- c(prior_parameters, formula_parameters)
      if(length(add_parameters) > 0){
        parameters <- c(parameters, samples.row[add_parameters])
      }
      context <- NULL
      if(isTRUE(bridge_context)){
        context <- .bt_JAGS_bridge_context(
          samples = samples.row,
          prior_parameters = prior_parameters,
          formula_prior_parameters = formula_prior_parameters,
          formula_parameters = formula_parameters,
          add_parameters = add_parameters,
          formula_design_list = formula_design_list,
          formula_data_list = formula_data_list,
          formula_prior_list = formula_prior_list,
          model_data = data
        )
      }
      list(parameters = parameters, context = context)
    }, BayesTools_marglik_out_of_support = function(e)e)
    if(inherits(evaluated, "BayesTools_marglik_out_of_support")){
      return(-Inf)
    }

    marglik <- marglik + .bt_JAGS_bridge_call_log_posterior(
      log_posterior = log_posterior,
      parameters = evaluated$parameters,
      data = data,
      context = evaluated$context,
      bridge_context = bridge_context,
      ...
    )

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
      bridge_context     = bridge_context,
      formula_design_list = formula_design_list,
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
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
