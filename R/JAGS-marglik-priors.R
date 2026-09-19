#' @title Compute marginal likelihood for 'JAGS' priors
#'
#' @description Computes marginal likelihood for the
#' prior part of a 'JAGS' model within 'bridgesampling'
#' function
#'
#' @param samples samples provided by the bridgesampling function. Supply one
#' named posterior row to `JAGS_marglik_priors()` and a named matrix or data
#' frame to `JAGS_marglik_priors_rows()`. Row-wise samples must have real
#' numeric columns.
#' @param prior_list named list of model-level prior distributions. For
#' `JAGS_marglik_priors_formula()`, optional model priors whose shared ordered
#' allocations have already been included by `JAGS_marglik_priors()`.
#'
#' @details `JAGS_marglik_priors_rows()` preserves joint prior boundaries
#' within each posterior row. In particular, auxiliary-gamma contributions for
#' vector and Dirichlet priors are summed separately within each row.
#'
#' When model and formula priors share an ordered Dirichlet allocation, supply
#' the model priors as `prior_list` to `JAGS_marglik_priors_formula()`. Its
#' contribution then excludes allocations already included by
#' `JAGS_marglik_priors(samples, prior_list)`, so the two log prior contributions
#' can be added without counting shared allocations twice.
#' Evaluate the combined contribution as
#' `JAGS_marglik_priors(samples, prior_list) + JAGS_marglik_priors_formula(samples, formula_prior_list, prior_list)`.
#'
#' @inheritParams JAGS_bridgesampling
#'
#' @return \code{JAGS_marglik_priors} returns a numeric value
#' of likelihood evaluated at the current posterior sample.
#' \code{JAGS_marglik_priors_rows} returns one numeric value for every row of
#' posterior samples. \code{JAGS_marglik_priors_rows_evaluator} returns a
#' function that applies the same row-wise evaluation without recompiling the
#' prior structure on each call.
#'
#' @export JAGS_marglik_priors
#' @export JAGS_marglik_priors_rows
#' @export JAGS_marglik_priors_rows_evaluator
#' @export JAGS_marglik_priors_formula
#' @name JAGS_marglik_priors
NULL

#' @rdname JAGS_marglik_priors
JAGS_marglik_priors                <- function(samples, prior_list){

  .bt_JAGS_marglik_priors(samples, prior_list)
}


.bt_JAGS_marglik_priors <- function(samples, prior_list, emitted_allocations = character()){

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
  ordered_allocation_keys <- unique(emitted_allocations)
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

    }else if(is.prior.ordered(prior_list[[i]])){

      ordered_marglik <- .JAGS_marglik_priors.ordered(
        samples,
        prior_list[[i]],
        names(prior_list)[i],
        emitted_allocations = ordered_allocation_keys
      )
      marglik <- marglik + ordered_marglik[["marglik"]]
      ordered_allocation_keys <- unique(c(
        ordered_allocation_keys,
        ordered_marglik[["allocation_keys"]]
      ))

    }else if(is.prior.factor(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.factor(samples, prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.vector(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.vector(samples, prior_list[[i]], names(prior_list)[i])

    }else if(is.prior.simple(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.simple(samples, prior_list[[i]], names(prior_list)[i])

    }else if(!is.prior.none(prior_list[[i]])){

      stop("Unsupported prior object.", call. = FALSE)

    }
  }

  return(marglik)
}


#' @rdname JAGS_marglik_priors
JAGS_marglik_priors_rows <- function(samples, prior_list){

  evaluator <- JAGS_marglik_priors_rows_evaluator(prior_list)
  return(evaluator(samples))
}


#' @rdname JAGS_marglik_priors
JAGS_marglik_priors_rows_evaluator <- function(prior_list){

  evaluator <- .bt_JAGS_marglik_compile_prior_rows_evaluator(prior_list)
  return(function(samples){
    samples <- .bt_JAGS_marglik_prior_rows(samples)
    return(unname(evaluator(samples)))
  })
}


.bt_JAGS_marglik_prior_rows <- function(samples){

  if(is.null(dim(samples))){
    if(is.null(names(samples)))
      stop("'samples' must contain named posterior samples.", call. = FALSE)
    samples <- matrix(
      samples,
      nrow     = 1L,
      dimnames = list(NULL, names(samples))
    )
  }

  if(!is.matrix(samples) && !is.data.frame(samples))
    stop("'samples' must be a matrix or data frame.", call. = FALSE)
  numeric_columns <- if(is.data.frame(samples)){
    all(vapply(samples, function(column){
      is.numeric(column) && !is.complex(column)
    }, logical(1)))
  }else{
    is.numeric(samples) && !is.complex(samples)
  }
  if(!numeric_columns)
    stop("'samples' must contain only real numeric posterior sample columns.", call. = FALSE)
  if(ncol(samples) > 0L && is.null(colnames(samples)))
    stop("'samples' must contain named posterior samples.", call. = FALSE)

  return(samples)
}


.bt_JAGS_marglik_compile_prior_rows_evaluator <- function(prior_list){

  scalar_evaluator <- .bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)
  if(length(prior_list) == 0L){
    return(function(samples) numeric(nrow(samples)))
  }

  evaluators <- Map(
    .bt_JAGS_marglik_compile_prior_rows_component,
    prior_list,
    names(prior_list)
  )
  if(any(vapply(evaluators, is.null, logical(1)))){
    return(function(samples){
      vapply(seq_len(nrow(samples)), function(i){
        sample_row <- stats::setNames(
          as.numeric(unlist(samples[i, , drop = FALSE], use.names = FALSE)),
          colnames(samples)
        )
        scalar_evaluator$log_prior(sample_row)
      }, numeric(1))
    })
  }

  return(function(samples){
    marglik <- numeric(nrow(samples))
    for(evaluator in evaluators){
      marglik <- marglik + evaluator(samples)
    }
    return(marglik)
  })
}


.bt_JAGS_marglik_compile_prior_rows_component <- function(prior_object,
                                                           parameter_name){

  force(prior_object)
  force(parameter_name)

  if(is.prior.none(prior_object) || is.prior.point(prior_object)){
    return(function(samples) numeric(nrow(samples)))
  }

  if(is.prior.factor(prior_object) &&
     (is.prior.treatment(prior_object) || is.prior.independent(prior_object))){
    parameter_names <- .JAGS_prior_factor_names(parameter_name, prior_object)
    if(identical(prior_object[["distribution"]], "invgamma")){
      return(.bt_JAGS_marglik_compile_invgamma_prior_rows(prior_object, parameter_names))
    }
    log_density <- .prior_simple_lpdf_evaluator(prior_object)
    return(function(samples){
      if(!all(parameter_names %in% colnames(samples)))
        stop("'samples' does not contain all monitored factor prior parameters.", call. = FALSE)

      marglik <- numeric(nrow(samples))
      for(name in parameter_names){
        marglik <- marglik + log_density(samples[, name])
      }
      return(marglik)
    })
  }

  if(is.prior.PET(prior_object)){
    parameter_name <- "PET"
  }else if(is.prior.PEESE(prior_object)){
    parameter_name <- "PEESE"
  }
  is_plain_simple <- is.prior.simple(prior_object) &&
    !is.prior.factor(prior_object)
  if(is_plain_simple &&
     identical(prior_object[["distribution"]], "invgamma")){
    return(.bt_JAGS_marglik_compile_invgamma_prior_rows(prior_object, parameter_name))
  }

  if(is_plain_simple){
    log_density <- .prior_simple_lpdf_evaluator(prior_object)
    return(function(samples){
      if(!parameter_name %in% colnames(samples))
        stop("'samples' does not contain all monitored prior parameters.", call. = FALSE)
      return(log_density(samples[, parameter_name]))
    })
  }

  if(is.prior.vector(prior_object) &&
     identical(prior_object[["distribution"]], "dirichlet")){
    alpha <- prior_object$parameters[["alpha"]]
    eta_names <- paste0(
      .JAGS_prior_dirichlet_eta_name(parameter_name),
      "[", seq_along(alpha), "]"
    )
    return(function(samples){
      if(!all(eta_names %in% colnames(samples)))
        stop("'samples' does not contain all monitored Dirichlet prior parameters.", call. = FALSE)

      marglik <- numeric(nrow(samples))
      invalid <- logical(nrow(samples))
      for(i in seq_along(alpha)){
        eta      <- samples[, eta_names[i]]
        invalid  <- invalid | !is.finite(eta) | eta <= 0
        marglik <- marglik + stats::dgamma(
          eta,
          shape = alpha[i],
          rate  = 1,
          log   = TRUE
        )
      }
      marglik[invalid] <- -Inf
      return(marglik)
    })
  }

  return(NULL)
}


.bt_JAGS_marglik_compile_invgamma_prior_rows <- function(prior_object, parameter_names){

  log_density <- .prior_simple_lpdf_evaluator(prior_object)
  force(parameter_names)
  function(samples){
    marglik <- numeric(nrow(samples))
    for(parameter_name in parameter_names){
      if(parameter_name %in% colnames(samples)){
        values <- samples[, parameter_name]
      }else{
        legacy_name <- paste0("inv_", parameter_name)
        if(!legacy_name %in% colnames(samples))
          stop("'samples' does not contain all monitored inverse-gamma prior parameters.", call. = FALSE)
        values <- samples[, legacy_name]^-1
      }

      supported <- is.finite(values) & values > 0
      contribution <- rep(-Inf, nrow(samples))
      contribution[supported] <- log_density(values[supported])
      marglik <- marglik + contribution
    }
    marglik
  }
}


.JAGS_marglik_priors.ordered        <- function(samples, prior, parameter_name, emitted_allocations = character()){

  .prior_ordered_bridge_check(prior)
  total_names <- .prior_ordered_total_monitor_names(prior, parameter_name)

  marglik <- 0
  if(!is.prior.point(prior$total)){
    if(!all(total_names %in% names(samples))){
      stop("'samples' does not contain all monitored ordered total prior parameters.", call. = FALSE)
    }
    total_values <- unname(unlist(samples[total_names], use.names = FALSE))
    marglik <- marglik + sum(lpdf(prior$total, total_values))
  }

  emitted_now <- character()
  for(record in .prior_ordered_dirichlet_records(prior)){
    if(record$key %in% emitted_allocations || record$key %in% emitted_now){
      next
    }
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(record$node), "[", seq_len(record$dim), "]")
    eta <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = eta_names,
      missing_message = "'samples' does not contain all monitored ordered Dirichlet allocation parameters."
    )
    if(is.null(eta)){
      return(list(marglik = -Inf, allocation_keys = emitted_now))
    }
    marglik <- marglik + sum(stats::dgamma(eta, shape = record$spec$alpha, rate = 1, log = TRUE))
    emitted_now <- c(emitted_now, record$key)
  }

  list(marglik = marglik, allocation_keys = emitted_now)
}

.JAGS_marglik_priors.simple         <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(prior[["distribution"]] == "invgamma"){

    value <- .bt_JAGS_marglik_invgamma_values(
      samples = samples,
      parameter_names = parameter_name,
      missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters."
    )
    if(is.null(value)){
      return(-Inf)
    }
    marglik <- lpdf(prior, value)

  }else if(prior[["distribution"]] == "point"){

    marglik <- 0

  }else{

    if(!parameter_name %in% names(samples)){
      stop("'samples' does not contain all monitored prior parameters.", call. = FALSE)
    }
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
  }else{
    parameter_names <- if(prior$parameters[["K"]] == 1L){
      parameter_name
    }else{
      paste0(parameter_name, "[", seq_len(prior$parameters[["K"]]), "]")
    }
    if(!all(parameter_names %in% names(samples))){
      stop("'samples' does not contain all monitored vector prior parameters.", call. = FALSE)
    }
    marglik <- lpdf(prior, unlist(samples[parameter_names], use.names = FALSE))
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

    if(J == 2L){
      beta_parameters <- .weightfunction_alpha_marginal(
        prior$weights$alpha,
        2L
      )
      omega <- .bt_JAGS_marglik_binary_cumulative_weight(samples)
      if(is.null(omega)){
        return(-Inf)
      }
      marglik <- stats::dbeta(
        omega,
        shape1 = beta_parameters$alpha,
        shape2 = beta_parameters$beta,
        log = TRUE
      )
    }else{
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
    }

  }else if(prior$weights$type == "independent"){

    if(J == 1L){
      marglik <- 0
    }else if(prior$weights$scale == "omega"){
      if(!all(paste0("omega[", 2:J, "]") %in% names(samples))){
        stop("'samples' does not contain all monitored independent weightfunction parameters.", call. = FALSE)
      }
      marglik <- sum(mlpdf(prior$weights$prior, samples[paste0("omega[", 2:J, "]")]))
    }else if(prior$weights$scale == "log_omega"){
      if(!all(paste0("log_omega[", 2:J, "]") %in% names(samples))){
        stop("'samples' does not contain all monitored independent weightfunction parameters.", call. = FALSE)
      }
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

  selection_backend_spec(prior, include_init = FALSE)

  marglik <- 0
  if(!is.null(prior$selection)){
    marglik <- marglik + .JAGS_marglik_priors.weightfunction(samples, prior$selection)
  }
  if(!is.null(prior$phacking)){
    marglik <- marglik + .JAGS_marglik_priors.phacking(samples, prior$phacking)
  }

  return(marglik)
}
#' @rdname JAGS_marglik_priors
JAGS_marglik_priors_formula <- function(samples, formula_prior_list, prior_list = NULL){

  if(length(formula_prior_list) == 0L){
    return(0)
  }

  model_evaluator <- .bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)
  formula_priors <- do.call(c, unname(formula_prior_list))
  .bt_JAGS_marglik_priors(
    samples = samples,
    prior_list = formula_priors,
    emitted_allocations = model_evaluator$allocation_keys
  )
}

.bt_JAGS_random_effect_scalar_rho_support_spec <- function(
    random_term,
    structure,
    context = "Bridge sampling random-effect metadata"){

  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = context
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation'.",
      call. = FALSE
    )
  }

  bounds <- .bt_random_effect_rho_bounds_metadata(
    correlation,
    random_term = random_term,
    context = context
  )
  exact_bounds <- .bt_random_effect_structured_rho_bounds(
    K = random_term$n_columns,
    structure = structure
  )
  bounds_values <- as.numeric(bounds[c("lower", "upper")])
  if(!identical(bounds_values, as.numeric(exact_bounds))){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " contain scalar correlation bounds that do not match the exact ",
      toupper(structure),
      " support.",
      call. = FALSE
    )
  }

  if(identical(structure, "car")){
    time_values <- correlation$time_values
    if(!is.numeric(time_values) || length(time_values) != random_term$n_columns ||
       any(!is.finite(time_values)) || anyDuplicated(time_values) ||
       any(diff(time_values) <= 0)){
      stop(
        context,
        .bt_random_effect_metadata_block_detail(random_term),
        " are missing canonical ordered CAR time coordinates.",
        call. = FALSE
      )
    }
  }

  list(
    correlation = correlation,
    bounds = exact_bounds
  )
}
