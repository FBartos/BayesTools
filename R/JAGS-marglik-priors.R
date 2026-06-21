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

  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  sampled_random_effect <- identical(
    .bt_random_effect_term_compile_mode(random_term),
    "sampled"
  )
  marglik <- 0
  if(isTRUE(sampled_random_effect)){
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

