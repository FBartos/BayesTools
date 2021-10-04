#' @title Compute marginal likelihood of a 'JAGS' model
#'
#' @description A wrapper around
#' \link[bridgesampling]{bridge_sampler} that automatically
#' computes likelihood part dependent on the prior distribution
#' and prepares parameter samples. \code{log_posterior} must
#' specify a function that takes two arguments - a named list
#' of samples from the prior distributions and the data, and returns
#' log likelihood of the model part.
#'
#' @param fit model fitted with either \link[runjags]{runjags} posterior
#' samples obtained with \link[rjags]{rjags-package}
#' @param data data that were used to fit the model
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param log_posterior function that takes a named list of samples, the data,
#' and additional list of parameters passed as \code{...} as input and
#' returns the log of the unnormalized posterior density of the model part
#' @param add_parameters vector of additional parameter names that should be used
#' in bridgesampling but were not specified in the \code{prior_list}
#' @param add_bounds list with two name vectors (\code{"lb"} and \code{"up"})
#' containing lower and upper bounds of the additional parameters that were not
#' specified in the \code{prior_list}
#' @param maxiter maximum number of iterations for the
#' \link[bridgesampling]{bridge_sampler}
#' @param silent whether the progress should be printed, defaults to \code{TRUE}
#' @param ... additional argument to the \link[bridgesampling]{bridge_sampler}
#' and \code{log_posterior} function
#'
#' @examples
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
#' fit <- JAGS_fit(model_syntax, data, priors_list,
#'                 chains = 1, burnin = 100, sample = 500, seed = 0)
#'
#' # define log posterior for bridge sampling
#' log_posterior <- function(parameters, data){
#'   sum(dnorm(data$x, parameters$mu, 1, log = TRUE))
#' }
#'
#' # get marginal likelihoods
#' marglik <- JAGS_bridgesampling(fit, data, priors_list, log_posterior)
#'
#' @return \code{JAGS_bridgesampling} returns an object of class 'bridge'.
#'
#' @export
JAGS_bridgesampling <- function(fit, data, prior_list, log_posterior,
                                add_parameters = NULL, add_bounds = NULL,
                                maxiter = 10000, silent = TRUE, ...){


  check_bool(silent, "silent")
  check_int(maxiter, "maxiter", lower = 1)


  ### check the input and split it on posterior and data
  if(inherits(fit, "runjags")){

    # get posterior and merge chains
    posterior <- suppressWarnings(coda::as.mcmc(fit))

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

    # rjags model with rjags::coda.samples
    # merge chains
    posterior <- do.call(rbind, fit)

  }else{

    stop("the method is not implemented for this output")

  }


  ### extract relevant variables and upper and lower bound
  bridgesampling_posterior <- JAGS_bridgesampling_posterior(posterior, prior_list, add_parameters, add_bounds)
  if(ncol(bridgesampling_posterior) == 0)
    stop("Bridge sampling cannot proceed without any estimated parameter")


  ### the marglik function
  full_log_posterior <- function(samples.row, data, prior_list, add_parameters, ...){

    parameters      <- JAGS_marglik_parameters(samples.row, prior_list)
    if(!is.null(add_parameters)){
      parameters    <- c(parameters, samples.row[add_parameters])
    }
    marglik_priors  <- JAGS_marglik_priors(samples.row, prior_list)
    marglik_model   <- log_posterior(parameters, data, ...)

    return(marglik_priors + marglik_model)
  }


  ### perform bridgesampling
  marglik <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
      samples        = bridgesampling_posterior,
      data           = data,
      log_posterior  = full_log_posterior,
      prior_list     = prior_list,
      lb             = attr(bridgesampling_posterior, "lb"),
      ub             = attr(bridgesampling_posterior, "ub"),
      silent         = silent,
      maxiter        = maxiter,
      add_parameters = add_parameters,
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



#' @title Prepare 'JAGS' posterior for 'bridgesampling'
#'
#' @description prepares posterior distribution for 'bridgesampling'
#' by removing unnecessary parameters and attaching lower and upper
#' bounds of parameters based on a list of prior distributions.
#'
#' @param posterior matrix of mcmc samples from the posterior
#' distribution
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#'
#' @inheritParams JAGS_bridgesampling
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
  if(!is.null(prior_list)){
    if(!is.list(prior_list))
      stop("'prior_list' must be a list.")
    if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
      stop("'prior_list' must be a list of priors.")
  }
  if(!is.null(add_parameters)){
    if(!is.character(add_parameters))
      stop("'add_parameters' must be a character vector.")
    if(!is.list(add_bounds))
      stop("'add_bounds' must be a list.")
    if(length(add_bounds) != 2 | !all(names(add_bounds) %in% c("lb", "ub")))
      stop("'add_bounds' must contain lower and upper bounds ('lb' and 'ub').")
    if(length(add_bounds[["lb"]]) != length(add_parameters) | length(add_bounds[["ub"]]) != length(add_parameters))
      stop("lb' and 'ub' must have the same lenght as the 'add_parameters'.")
    if(!is.numeric(add_bounds[["lb"]]) | !is.numeric(add_bounds[["ub"]]))
      stop("lb' and 'ub' must be numeric vectors.")
  }


  # get information about the specified parameters
  parameters_names <- .JAGS_bridgesampling_posterior_info(prior_list)

  # add the user defined parameters
  if(!is.null(add_parameters)){
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

.JAGS_bridgesampling_posterior_info                <- function(prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # add the resulting parameters
  parameters    <- NULL
  parameters_lb <- NULL
  parameters_ub <- NULL
  for(i in seq_along(prior_list)){

    add_parameter <- NULL

    if(is.prior.weightfunction(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.weightfunction(prior_list[[i]])

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      add_parameter <- .JAGS_bridgesampling_posterior_info.PP(prior_list[[i]])

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
  if(!is.character(parameter_name) | length(parameter_name) != 1)
    stop("'parameter_name' must be a character vector of length 1.")


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


  if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){

    parameter <- paste0("eta[",1:length(prior$parameters[["alpha"]]),"]")
    attr(parameter, "lb") <- rep(0,   length(parameter))
    attr(parameter, "ub") <- rep(Inf, length(parameter))

  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){

    parameter <- c(paste0("eta1[",1:length(prior$parameters[["alpha1"]]),"]"), paste0("eta2[",1:length(prior$parameters[["alpha2"]]),"]"))
    attr(parameter, "lb") <- rep(0,   length(parameter))
    attr(parameter, "ub") <- rep(Inf, length(parameter))

  }else if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    parameter <- NULL

  }

  names(attr(parameter, "lb")) <- parameter
  names(attr(parameter, "ub")) <- parameter

  return(parameter)
}

#' @title Compute marginal likelihood for 'JAGS' priors
#'
#' @description Computes marginal likelihood for the
#' prior part of a 'JAGS' model within 'bridgesampling'
#' function
#'
#' @param samples samples provided by bridgesampling
#' function
#'
#' @inheritParams JAGS_add_priors
#'
#' @return \code{JAGS_marglik_priors} returns a numeric value
#' of likelihood evaluated at the current posterior sample.
#'
#' @export
JAGS_marglik_priors                <- function(samples, prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # add the resulting parameters
  marglik <- 0
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.weightfunction(samples, prior_list[[i]])

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      marglik <- marglik + .JAGS_marglik_priors.PP(samples, prior_list[[i]])

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
  if(!is.character(parameter_name) | length(parameter_name) != 1)
    stop("'parameter_name' must be a character vector of length 1.")

  if(prior[["distribution"]] == "invgamma"){

    sampling_prior <- prior(
      "distribution" = "gamma",
      "parameters"   = list("shape" = prior$parameters[["shape"]], "rate" = prior$parameters[["scale"]]),
      "truncation"   = list("lower" = prior$truncation[["upper"]]^-1, "upper" = prior$truncation[["lower"]]^-1))
    marglik <- lpdf(sampling_prior, samples[[ paste0("inv_", parameter_name) ]])

  }else if(prior[["distribution"]] == "point"){

    marglik <- 0

  }else{

    marglik <- lpdf(prior, samples[[ parameter_name ]])

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

  if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    marglik <- 0

  }else if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){

    marglik <- sum(stats::dgamma(samples[ paste0("eta[",1:length(prior$parameters[["alpha"]]),"]") ], shape = prior$parameters[["alpha"]], rate = 1, log = TRUE))

  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){

    marglik <-
      sum(stats::dgamma(samples[ paste0("eta1[",1:length(prior$parameters[["alpha1"]]),"]") ], shape = prior$parameters[["alpha1"]], rate = 1, log = TRUE)) +
      sum(stats::dgamma(samples[ paste0("eta2[",1:length(prior$parameters[["alpha2"]]),"]") ], shape = prior$parameters[["alpha2"]], rate = 1, log = TRUE))

  }

  return(marglik)
}



#' @title Extract parameters for 'JAGS' priors
#'
#' @description Extracts transformed parameters from the
#' prior part of a 'JAGS' model inside of a 'bridgesampling'
#' function (returns them as a named list)
#'
#' @param samples samples provided by bridgesampling
#' function
#'
#' @return \code{JAGS_marglik_parameters} returns a named list
#' of (transformed) posterior samples.
#'
#' @inheritParams JAGS_add_priors
#' @export
JAGS_marglik_parameters                <- function(samples, prior_list){

  # return empty string in case that no prior was specified
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

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      parameters <- c(parameters, .JAGS_marglik_parameters.PP(samples, prior_list[[i]]))

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
  if(!is.character(parameter_name) | length(parameter_name) != 1)
    stop("'parameter_name' must be a character vector of length 1.")


  parameter <- list()
  if(prior[["distribution"]] == "invgamma"){
    parameter[[parameter_name]] <- samples[[ paste0("inv_", parameter_name) ]]^-1
  }else if(prior[["distribution"]] == "point"){
    parameter[[parameter_name]] <- prior$parameters[["location"]]
  }else{
    parameter[[parameter_name]] <- samples[[ parameter_name ]]
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
  if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){

    eta     <- samples[ paste0("eta[",1:length(prior$parameters[["alpha"]]),"]") ]
    std_eta <- eta / sum(eta)
    parameter[["omega"]] <- cumsum(std_eta)

  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){

    J1 <- length(prior$parameters[["alpha1"]])
    J2 <- length(prior$parameters[["alpha2"]])
    omega    <- rep(NA, J1 + J2 - 1)
    eta1     <- samples[ paste0("eta1[",1:J1,"]") ]
    std_eta1 <- eta1 / sum(eta1)
    omega[J2:length(omega)] <- cumsum(std_eta1)
    eta2     <- samples[ paste0("eta2[",1:J2,"]") ]
    std_eta2 <- (eta2 / sum(eta2)) * (1 - std_eta1[1])
    omega[1:(J2-1)] <- rev(cumsum(std_eta2[J2:2])) + std_eta1[1]
    parameter[["omega"]] <- omega

  }else if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    parameter[["omega"]] <- prior$parameters[["omega"]]

  }

  return(parameter)
}
