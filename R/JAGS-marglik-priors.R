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
  ordered_allocation_keys <- character()
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

    }
  }

  return(marglik)
}


.JAGS_marglik_priors.ordered        <- function(samples, prior, parameter_name, emitted_allocations = character()){

  .prior_ordered_bridge_check(prior)
  total_names <- .prior_ordered_total_monitor_names(prior, parameter_name)

  marglik <- 0
  if(!is.prior.point(prior$total)){
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
#' @rdname JAGS_marglik_priors
JAGS_marglik_priors_formula <- function(samples, formula_prior_list){

  if(length(formula_prior_list) == 0L){
    return(0)
  }

  prior_list <- do.call(c, unname(formula_prior_list))
  JAGS_marglik_priors(samples, prior_list)
}

.bt_JAGS_marglik_priors_formula_random <- function(samples, formula_design_list){

  if(length(formula_design_list) == 0L){
    return(0)
  }
  .bt_JAGS_bridge_check_no_allocation_inclusion(formula_design_list)

  marglik <- 0
  design_names <- names(formula_design_list)
  if(is.null(design_names)){
    design_names <- rep("", length(formula_design_list))
  }
  for(design_i in seq_along(formula_design_list)){
    design <- formula_design_list[[design_i]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      next
    }
    .bt_JAGS_bridge_validate_formula_random_compile(
      parameter = .bt_JAGS_bridge_design_parameter_name(
        design,
        fallback = design_names[[design_i]]
      ),
      design = design,
      label = "stored"
    )
    for(random_term in .bt_formula_design_random_effects(design)){
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

  n_columns <- random_term$n_columns
  sampled_random_effect <- identical(
    .bt_random_effect_term_compile_mode(random_term),
    "sampled"
  )
  marglik <- 0
  if(isTRUE(sampled_random_effect)){
    marglik <- .bt_random_effect_latent_log_density(random_term, samples)
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

  support_spec <- .bt_JAGS_random_effect_scalar_rho_support_spec(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  correlation <- support_spec$correlation
  rho_scale <- .bt_random_effect_rho_scale_metadata(
    correlation,
    random_term = random_term,
    context = "Bridge sampling random-effect metadata"
  )
  source_name <- if(!identical(rho_scale, "rho") &&
                    correlation$sample_name %in% colnames(posterior)){
    correlation$sample_name
  }else if(correlation$rho_name %in% colnames(posterior)){
    correlation$rho_name
  }else{
    correlation$sample_name
  }
  if(source_name %in% colnames(posterior) &&
     any(!is.finite(posterior[, source_name]))){
    return(-Inf)
  }

  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "null",
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(rho) || any(is.na(rho) | !is.finite(rho))){
    return(-Inf)
  }
  if(any(.bt_random_effect_rho_outside_support(
    rho,
    bounds = support_spec$bounds,
    structure = structure
  ))){
    return(-Inf)
  }

  0
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
