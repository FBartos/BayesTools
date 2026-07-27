.fit_to_posterior <- function(fit){

  ### check the input and split it on posterior and data
  if(inherits(fit, "runjags")){

    # get posterior and merge chains
    posterior <- .extract_posterior_samples(fit, as_list = FALSE)

  }else if(is.list(fit) &&
           length(fit) > 0L &&
           all(vapply(fit, inherits, logical(1), what = "mcarray"))){

    # rjags model with rjags::jags.samples
    # flatten parameter dimensions and merge chains in coda order
    fit_names <- names(fit)
    converted <- lapply(seq_along(fit), function(i){
      fallback_name <- NULL
      if(!is.null(fit_names) &&
         !is.na(fit_names[i]) &&
         nzchar(fit_names[i])){
        fallback_name <- fit_names[i]
      }
      .bt_fit_mcarray_to_matrix(fit[[i]], fallback_name = fallback_name)
    })

    reference_layout <- c(
      converted[[1L]]$n_iteration,
      converted[[1L]]$n_chain
    )
    matching_layout <- vapply(converted, function(parameter){
      identical(
        c(parameter$n_iteration, parameter$n_chain),
        reference_layout
      )
    }, logical(1))
    if(!all(matching_layout)){
      stop(
        "All 'mcarray' parameters must have matching iteration and chain dimensions.",
        call. = FALSE
      )
    }

    reference_iterations <- converted[[1L]]$iterations
    matching_iterations <- vapply(converted, function(parameter){
      identical(parameter$iterations, reference_iterations)
    }, logical(1))
    if(!all(matching_iterations)){
      stop(
        "All 'mcarray' parameters must have matching iteration metadata.",
        call. = FALSE
      )
    }

    posterior <- do.call(cbind, lapply(converted, `[[`, "samples"))

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

.bt_fit_mcarray_to_matrix <- function(parameter, fallback_name = NULL){

  parameter_dim <- dim(parameter)
  if(is.null(parameter_dim) || length(parameter_dim) < 3L){
    stop(
      "Each 'mcarray' parameter must have at least one parameter dimension followed by iteration and chain dimensions.",
      call. = FALSE
    )
  }
  if(any(parameter_dim <= 0L)){
    stop(
      "All 'mcarray' dimensions must be positive.",
      call. = FALSE
    )
  }

  variable_dim <- parameter_dim[seq_len(length(parameter_dim) - 2L)]
  n_iteration <- parameter_dim[length(parameter_dim) - 1L]
  n_chain <- parameter_dim[length(parameter_dim)]

  varname <- attr(parameter, "varname", exact = TRUE)
  if(!is.character(varname) ||
     length(varname) != 1L ||
     is.na(varname) ||
     !nzchar(varname)){
    varname <- fallback_name
  }
  if(!is.character(varname) ||
     length(varname) != 1L ||
     is.na(varname) ||
     !nzchar(varname)){
    stop(
      "Each 'mcarray' parameter must have a non-empty 'varname' attribute or list name.",
      call. = FALSE
    )
  }

  samples <- t(matrix(
    as.vector(parameter),
    nrow = prod(variable_dim),
    ncol = n_iteration * n_chain
  ))
  colnames(samples) <- .bt_fit_mcarray_column_names(
    varname = varname,
    variable_dim = variable_dim
  )

  list(
    samples = samples,
    n_iteration = n_iteration,
    n_chain = n_chain,
    iterations = attr(parameter, "iterations", exact = TRUE)
  )
}

.bt_fit_mcarray_column_names <- function(varname, variable_dim){

  # Match the ordering used by rjags::coda.samples: the first parameter
  # dimension varies fastest, followed by the remaining dimensions.
  if(prod(variable_dim) == 1L){
    return(varname)
  }

  basename <- varname
  n_dimension <- length(variable_dim)
  lower <- rep.int(1L, n_dimension)
  upper <- variable_dim

  parsed_name <- .bt_fit_mcarray_parse_varname(varname)
  if(!is.null(parsed_name$lower) &&
     !is.null(parsed_name$upper) &&
     length(parsed_name$lower) == length(parsed_name$upper)){
    monitored_dim <- parsed_name$upper - parsed_name$lower + 1
    matching_dim <- isTRUE(all.equal(
      variable_dim[variable_dim != 1L],
      monitored_dim[monitored_dim != 1L],
      check.attributes = FALSE
    ))
    if(matching_dim){
      basename <- parsed_name$name
      lower <- parsed_name$lower
      upper <- parsed_name$upper
      n_dimension <- length(monitored_dim)
    }
  }

  indices <- as.character(seq.int(lower[1L], upper[1L]))
  if(n_dimension > 1L){
    for(i in 2:n_dimension){
      indices <- outer(
        indices,
        seq.int(lower[i], upper[i]),
        FUN = paste,
        sep = ","
      )
    }
  }

  paste0(basename, "[", as.vector(indices), "]")
}

.bt_fit_mcarray_parse_varname <- function(varname){

  parsed_name <- try(parse(text = varname, n = 1L), silent = TRUE)
  if(!is.expression(parsed_name) || length(parsed_name) != 1L){
    return(NULL)
  }

  parsed_name <- parsed_name[[1L]]
  if(is.name(parsed_name)){
    return(list(name = deparse(parsed_name)))
  }
  if(!is.call(parsed_name) ||
     !identical(deparse(parsed_name[[1L]]), "[") ||
     length(parsed_name) <= 2L){
    return(NULL)
  }

  parsed_elements <- vapply(
    parsed_name,
    function(element) paste(deparse(element), collapse = ""),
    character(1)
  )
  if(any(!nzchar(parsed_elements))){
    return(NULL)
  }

  n_dimension <- length(parsed_name) - 2L
  lower <- upper <- numeric(n_dimension)
  for(i in seq_len(n_dimension)){
    index <- parsed_name[[i + 2L]]
    if(is.numeric(index)){
      lower[i] <- upper[i] <- index
    }else if(is.call(index) &&
             length(index) == 3L &&
             identical(deparse(index[[1L]]), ":") &&
             is.numeric(index[[2L]]) &&
             is.numeric(index[[3L]])){
      lower[i] <- index[[2L]]
      upper[i] <- index[[3L]]
    }else{
      return(NULL)
    }
  }
  if(any(upper < lower)){
    return(NULL)
  }

  list(
    name = deparse(parsed_name[[2L]]),
    lower = lower,
    upper = upper
  )
}


#' @title Create a BayesTools marginal-likelihood object
#'
#' @description Creates a marginal-likelihood result from a supplied scalar
#' natural-log marginal likelihood. This is useful for analytically known or
#' deliberately zero-evidence models that were not evaluated by bridge
#' sampling.
#'
#' @param logml scalar natural-log marginal likelihood. Defaults to \code{-Inf},
#' representing zero marginal likelihood.
#'
#'
#' @return A `BayesTools_marglik` object with the `"supplied_scalar"`
#' aggregation rule and no repetition rows.
#'
#' @export
bridgesampling_object <- function(logml = -Inf){

  .bt_marglik_manual_result(logml)
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
  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma support.
  posterior <- .bt_JAGS_bridge_materialize_legacy_invgamma_posterior(
    posterior = posterior,
    prior_list = prior_list
  )

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

# TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma support.
.bt_JAGS_bridge_materialize_legacy_invgamma_posterior <- function(posterior,
                                                                  prior_list){

  if(length(prior_list) == 0L){
    return(posterior)
  }

  for(i in seq_along(prior_list)){
    posterior <- .bt_JAGS_bridge_materialize_legacy_invgamma_prior(
      posterior = posterior,
      prior = prior_list[[i]],
      parameter_name = names(prior_list)[i]
    )
  }

  posterior
}

.bt_JAGS_bridge_materialize_legacy_invgamma_prior <- function(posterior,
                                                              prior,
                                                              parameter_name){

  if(is_prior_bias(prior)){
    if(!is.null(prior$phacking)){
      posterior <- .bt_JAGS_bridge_materialize_legacy_invgamma_prior(
        posterior = posterior,
        prior = prior$phacking,
        parameter_name = "alpha"
      )
    }
    return(posterior)
  }

  if(is_prior_phacking(prior)){
    return(.bt_JAGS_bridge_materialize_legacy_invgamma_prior(
      posterior = posterior,
      prior = prior$alpha,
      parameter_name = "alpha"
    ))
  }

  if(is.prior.PET(prior)){
    parameter_name <- "PET"
  }

  if(is.prior.PEESE(prior)){
    parameter_name <- "PEESE"
  }

  if(!is.prior.simple(prior) || !identical(prior[["distribution"]], "invgamma")){
    return(posterior)
  }

  parameter_names <- if(is.prior.factor(prior)){
    .JAGS_prior_factor_names(parameter_name, prior)
  }else{
    parameter_name
  }

  .bt_JAGS_bridge_materialize_legacy_invgamma_columns(
    posterior = posterior,
    parameter_names = parameter_names
  )
}

.bt_JAGS_bridge_materialize_legacy_invgamma_columns <- function(posterior,
                                                                parameter_names){

  missing <- !parameter_names %in% colnames(posterior)
  if(!any(missing)){
    return(posterior)
  }

  missing_parameter_names <- parameter_names[missing]
  legacy_names <- paste0("inv_", missing_parameter_names)
  if(!all(legacy_names %in% colnames(posterior))){
    return(posterior)
  }

  values <- 1 / posterior[, legacy_names, drop = FALSE]
  colnames(values) <- missing_parameter_names
  cbind(posterior, values)
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
  if(is.prior.ordered(prior)){
    K <- .get_prior_factor_levels(prior)
    if(K > 1L){
      owned <- c(owned, paste0(parameter_name, "[", seq_len(K), "]"))
    }
    for(record in .prior_ordered_dirichlet_records(prior)){
      owned <- c(owned, record$node, paste0(record$node, "[", seq_len(record$dim), "]"))
    }
  }else if(is.prior.factor(prior)){
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
