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

    }else if(is.prior.ordered(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.ordered(samples, prior_list[[i]], names(prior_list)[i]))

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
    value <- .bt_JAGS_marglik_invgamma_values(
      samples = samples,
      parameter_names = parameter_name,
      missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters.",
      signal = TRUE
    )
    parameter[[parameter_name]] <- value
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

  if(prior[["distribution"]] == "invgamma"){
    return(.bt_JAGS_marglik_invgamma_values(
      samples = samples,
      parameter_names = parameter_names,
      missing_message = "'samples' does not contain all monitored formula prior parameters.",
      signal = TRUE
    ))
  }

  sample_names <- parameter_names
  if(!all(sample_names %in% names(samples))){
    stop("'samples' does not contain all monitored formula prior parameters.", call. = FALSE)
  }

  values <- unname(unlist(samples[sample_names], use.names = FALSE))

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
.JAGS_marglik_parameters.ordered       <- function(samples, prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.ordered(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  parameter <- list()
  parameter_names <- .JAGS_prior_factor_names(parameter_name, prior)
  parameter[[parameter_name]] <- .JAGS_marglik_parameter_values(samples, prior, parameter_names)

  parameter
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
.bt_JAGS_marglik_invgamma_values <- function(samples,
                                             parameter_names,
                                             missing_message,
                                             signal = FALSE){

  if(all(parameter_names %in% names(samples))){
    values <- unname(unlist(samples[parameter_names], use.names = FALSE))
    invalid <- !is.finite(values) | values <= 0
    if(any(invalid)){
      if(isTRUE(signal)){
        .bt_JAGS_marglik_out_of_support(
          "Bridge samples contain out-of-support inverse-gamma coordinate '",
          parameter_names[which(invalid)[1L]],
          "'."
        )
      }
      return(NULL)
    }
    return(values)
  }

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma support.
  legacy_names <- paste0("inv_", parameter_names)
  if(all(legacy_names %in% names(samples))){
    legacy_values <- .bt_JAGS_marglik_positive_auxiliary_values(
      samples = samples,
      parameter_names = legacy_names,
      missing_message = missing_message,
      signal = signal
    )
    if(is.null(legacy_values)){
      return(NULL)
    }
    return(legacy_values^-1)
  }

  stop(missing_message, call. = FALSE)
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

