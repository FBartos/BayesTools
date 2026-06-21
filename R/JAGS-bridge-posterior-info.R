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


  if(prior[["distribution"]] == "point"){
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

