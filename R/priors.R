#' @title Creates a prior distribution
#'
#' @description \code{prior} creates a prior distribution.
#' The prior can be visualized by the \code{plot} function.
#'
#' @param distribution name of the prior distribution. The
#' possible options are
#' \describe{
#'   \item{\code{"point"}}{for a point density characterized by a
#'   \code{location} parameter.}
#'   \item{\code{"normal"}}{for a normal distribution characterized
#'   by a \code{mean} and \code{sd} parameters.}
#'   \item{\code{"lognormal"}}{for a lognormal distribution characterized
#'   by a \code{meanlog} and \code{sdlog} parameters.}
#'   \item{\code{"cauchy"}}{for a Cauchy distribution characterized
#'   by a \code{location} and \code{scale} parameters. Internally
#'   converted into a generalized t-distribution with \code{df = 1}.}
#'   \item{\code{"t"}}{for a generalized t-distribution characterized
#'   by a \code{location}, \code{scale}, and \code{df} parameters.}
#'   \item{\code{"gamma"}}{for a gamma distribution characterized
#'   by either \code{shape} and \code{rate}, or \code{shape} and
#'   \code{scale} parameters. The later is internally converted to
#'   the \code{shape} and \code{rate} parametrization}
#'   \item{\code{"invgamma"}}{for an inverse-gamma distribution
#'   characterized by a \code{shape} and \code{scale} parameters. The
#'   JAGS part uses a 1/gamma distribution with a shape and rate
#'   parameter.}
#'   \item{\code{"beta"}}{for a beta distribution
#'   characterized by an \code{alpha} and \code{beta} parameters.}
#'   \item{\code{"exp"}}{for an exponential distribution
#'   characterized by either \code{rate} or \code{scale}
#'   parameter. The later is internally converted to
#'   \code{rate}.}
#'   \item{\code{"uniform"}}{for a uniform distribution defined on a
#'   range from \code{a} to \code{b}}
#' }
#' @param parameters list of appropriate parameters for a given
#' \code{distribution}.
#' @param truncation list with two elements, \code{lower} and
#' \code{upper}, that define the lower and upper truncation of the
#' distribution. Defaults to \code{list(lower = -Inf, upper = Inf)}.
#' The truncation is automatically set to the bounds of the support.
#' @param prior_weights prior odds associated with a given distribution.
#' The value is passed into the model fitting function, which creates models
#' corresponding to all combinations of prior distributions for each of
#' the model parameters and sets the model priors odds to the product
#' of its prior distributions.
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # create a half-normal standard normal prior distribution
#' p2 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1),
#' truncation = list(lower = 0, upper = Inf))
#'
#' # the prior distribution can be visualized using the plot function
#' # (see ?plot.prior for all options)
#' plot(p1)
#'
#' @return \code{prior} and \code{prior_none} return an object of class 'prior'.
#' A named list containing the distribution name, parameters, and prior weights.
#'
#' @name prior
#' @export prior
#' @export prior_none
#' @seealso [plot.prior()], \link[stats]{Normal}, \link[stats]{Lognormal}, \link[stats]{Cauchy},
#' \link[stats]{Beta}, \link[stats]{Exponential},
#' \link[extraDistr]{LocationScaleT}, \link[extraDistr]{InvGamma}.

#' @rdname prior
prior <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_weights = 1){

  # general input check (detailed checks are performed withing the constructors)
  check_char(distribution, "distribution")
  check_list(parameters, "parameters")
  #sapply(seq_along(parameters), function(i)check_real(parameters[[i]], names(parameters[[i]]), check_length = 0))
  check_list(truncation, "truncation")
  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

  # clean the input name
  distribution <- .prior_clean_input_name(distribution)

  if(distribution %in% c("norm", "normal")){
    distribution <- "normal"
  }else if(distribution %in% c("lnorm", "lognormal")){
    distribution <- "lognormal"
  }else if(distribution %in% c("t", "student")){
    distribution <- "t"
  }else if(distribution %in% c("cauchy")){
    distribution <- "cauchy"
  }else if(distribution %in% c("invgamma", "inversegamma")){
    distribution <- "invgamma"
  }else if(distribution %in% c("gamma")){
    distribution <- "gamma"
  }else if(distribution %in% c("beta")){
    distribution <- "beta"
  }else if(distribution %in% c("exp", "exponential")){
    distribution <- "exp"
  }else if(distribution %in% c("uniform", "unif")){
    distribution <- "uniform"
  }else if(distribution %in% c("point", "spike")){
    distribution <- "point"
  }else if(distribution %in% c("multivariatenorm", "multivariatenormal", "mnorm", "mnormal")){
    distribution <- "mnormal"
  }else if(distribution %in% c("multivariatet", "multivariatestudent", "mt", "mstudent")){
    distribution <- "mt"
  }else if(distribution %in% c("multivariatecauchy", "mcauchy")){
    distribution <- "mcauchy"
  }else{
    stop(paste0("The specified distribution name '", distribution,"' is not known. Please, see '?prior' for more information about supported prior distributions."))
  }

  # check the passed settings
  output <- do.call(paste0(".prior_", distribution), list(parameters = parameters, truncation = truncation))

  # add the prior odds
  output$prior_weights <- prior_weights

  return(output)
}

#' @rdname prior
prior_none <- function(prior_weights = 1){

  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

  out <- list()
  out$distribution <- "none"
  out$prior_weights   <- prior_weights
  class(out)       <- c("prior", "prior.none")

  return(out)
}


#' @title Creates a prior distribution for a weight function
#'
#' @description \code{prior_weightfunction} creates a prior distribution for fitting
#' a RoBMA selection model. The prior can be visualized by the \code{plot} function.
#'
#' @param distribution name of the prior distribution. The
#' possible options are
#' \describe{
#'   \item{\code{"two.sided"}}{for a two-sided weight function
#'   characterized by a vector \code{steps} and vector \code{alpha}
#'   parameters. The \code{alpha} parameter determines an alpha
#'   parameter of Dirichlet distribution which cumulative sum
#'   is used for the weights omega.}
#'   \item{\code{"one.sided"}}{for a one-sided weight function
#'   characterized by either a vector \code{steps} and vector
#'   \code{alpha} parameter, leading to a monotonic one-sided
#'   function, or by a vector \code{steps}, vector \code{alpha1},
#'   and vector \code{alpha2} parameters leading non-monotonic
#'   one-sided weight function. The \code{alpha} / \code{alpha1} and
#'   \code{alpha2} parameters determine an alpha parameter of
#'   Dirichlet distribution which cumulative sum is used for
#'   the weights omega.}
#' }
#' @param parameters list of appropriate parameters for a given
#' \code{distribution}.
#' @param prior_weights prior odds associated with a given distribution.
#' The model fitting function usually creates models corresponding to
#' all combinations of prior distributions for each of the model
#' parameters, and sets the model priors odds to the product of
#' its prior distributions.
#'
#' @examples
#' p1 <- prior_weightfunction("one-sided", parameters = list(steps = c(.05, .10), alpha = c(1, 1, 1)))
#'
#' # the prior distribution can be visualized using the plot function
#' # (see ?plot.prior for all options)
#' plot(p1)
#'
#' @return \code{prior_weightfunction} returns an object of class 'prior'.
#'
#' @export  prior_weightfunction
#' @seealso [plot.prior()]
prior_weightfunction <- function(distribution, parameters, prior_weights = 1){

  # general input check
  check_char(distribution, "distribution")
  check_list(parameters, "parameters")
  sapply(seq_along(parameters), function(i)check_real(parameters[[i]], names(parameters[[i]]), check_length = 0))
  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

  # clean the input name
  distribution <- .prior_clean_input_name(distribution)

  if(distribution == "twosided"){
    distribution <- "two.sided"
  }else if(distribution == "onesided"){
    distribution <- "one.sided"
  }else if(distribution == "twosidedfixed"){
    distribution <- "two.sided.fixed"
  }else if(distribution == "onesidedfixed"){
    distribution <- "one.sided.fixed"
  }else{
    stop(paste0("The specified distribution name '", distribution,"' is not known. Please, see '?prior_weightfunction' for more information about supported prior distributions for weight functions."))
  }

  # check the passed settings
  output <- do.call(paste0(".prior_weightfunction_", distribution), list(parameters = parameters))

  # add the prior odds
  output$prior_weights <- prior_weights

  class(output) <- c("prior", "prior.weightfunction")

  return(output)
}

#' @title Creates a prior distribution for PET or PEESE models
#'
#' @description \code{prior} creates a prior distribution for fitting a PET or
#' PEESE style models in RoBMA. The prior distribution can be visualized
#' by the \code{plot} function.
#'
#' @examples
#' # create a half-Cauchy prior distribution
#' # (PET and PEESE specific functions automatically set lower truncation at 0)
#' p1 <- prior_PET(distribution = "Cauchy", parameters = list(location = 0, scale = 1))
#'
#' plot(p1)
#'
#' @return \code{prior_PET} and \code{prior_PEESE} return an object of class 'prior'.
#'
#' @inheritParams prior
#' @export prior_PET
#' @export prior_PEESE
#' @seealso [plot.prior()], [prior()]
#' @name prior_PP
NULL

#' @rdname prior_PP
prior_PET   <- function(distribution, parameters, truncation = list(lower = 0, upper = Inf), prior_weights = 1){

  output <- prior(distribution, parameters, truncation, prior_weights)

  class(output) <- c(class(output), "prior.PET")

  return(output)
}
#' @rdname prior_PP
prior_PEESE <- function(distribution, parameters, truncation = list(lower = 0, upper = Inf), prior_weights = 1){

  output <- prior(distribution, parameters, truncation, prior_weights)

  class(output) <- c(class(output), "prior.PEESE")

  return(output)
}

#' @title Creates a prior distribution for factors
#'
#' @description \code{prior_factor} creates a prior distribution for fitting
#' models with factor predictors. (Note that results across different operating
#' systems might vary due to differences in JAGS numerical precision.)
#'
#' @param contrast type of contrast for the prior distribution. The possible options are
#' \describe{
#'   \item{\code{"orthonormal"}}{for contrast centered around the grand mean
#'   with equal marginal distributions, making the prior distribution exchangeable
#'   across factor levels. Only supports \code{distribution = "mnormal"} and
#'   \code{distribution = "mt"} which generates the corresponding multivariate normal/t
#'   distributions.}
#'   \item{\code{"treatment"}}{for contrasts using the first level as a comparison
#'   group and setting equal prior distribution on differences between the individual
#'   factor levels and the comparison level.}
#' }
#'
#'
#' @examples
#' # create an orthonormal prior distribution
#' p1 <- prior_factor(distribution = "mnormal", contrast = "orthonormal",
#'                    parameters = list(mean = 0, sd = 1))
#'
#' @return return an object of class 'prior'.
#'
#' @inheritParams prior
#' @export  prior_factor
#' @seealso [prior()]
prior_factor <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_weights = 1, contrast = "orthonormal"){

  # general input check (detailed checks are performed withing the constructors)
  check_char(contrast, "contrast", allow_values = c("orthonormal", "treatment", "dummy"))

  # check its compatibility with the contrasts
  if(contrast == "orthonormal"){

    # add the (yet unspecified) dimensions parameter
    if(is.null(names(parameters))){
      parameters <- c(parameters, NA)
    }else{
      parameters[["K"]] <- NA
    }


    # generate the prior object
    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    if(!is.prior.vector(output))
      stop("'orthonormal' contrasts require vector prior distribution.")
    if(!all(sapply(output[["truncation"]], is.infinite)))
      stop("'orthonormal' contrasts do not support truncation.")

    class(output) <- c(class(output), "prior.factor", "prior.orthonormal")

  }else if(contrast %in% c("treatment", "dummy")){

    # generate the prior object
    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    if(!is.prior.simple(output))
      stop("'treatment' contrasts require univariate prior distribution.")

    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    class(output) <- c(class(output), "prior.factor", "prior.dummy")
  }

  return(output)
}

#### functions for constructing prior distributions ####
.prior_normal    <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("mean", "sd"), "normal")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$mean, "mean")
  .check_parameter(parameters$sd,   "sd")
  .check_parameter_positive(parameters$sd, "sd")

  # add the values to the output
  output$distribution <- "normal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_lognormal <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("meanlog", "sdlog"), "lognormal")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$meanlog, "meanlog")
  .check_parameter(parameters$sdlog,   "sdlog")
  .check_parameter_positive(parameters$sdlog, "sdlog")

  # add the values to the output
  output$distribution <- "lognormal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_cauchy    <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale"), "Cauchy")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter_positive(parameters$scale, "scale")

  # deal with as with a t-distribution
  parameters$df <- 1

  output$distribution <- "t"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_t         <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "df"), "student-t")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter(parameters$df,       "df")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_positive(parameters$df,    "df")

  # add the values to the output
  output$distribution <- "t"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_gamma     <- function(parameters, truncation){

  output <- list()

  # deal with possible scale parametrization
  if(is.null(names(parameters))){
    if("scale" %in% names(parameters)){
      parameters <- .check_and_name_parameters(parameters, c("shape", "scale"), "gamma")
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
  }

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("shape", "rate"), "gamma")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$shape, "shape")
  .check_parameter(parameters$rate,  "rate")
  .check_parameter_positive(parameters$shape, "shape")
  .check_parameter_positive(parameters$rate,  "rate")

  # add the values to the output
  output$distribution <- "gamma"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_invgamma  <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("shape", "scale"), "invgamma")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$shape, "shape")
  .check_parameter(parameters$scale, "scale")
  .check_parameter_positive(parameters$shape, "shape")
  .check_parameter_positive(parameters$scale, "scale")

  # add the values to the output
  output$distribution <- "invgamma"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_exp       <- function(parameters, truncation){

  output <- list()

  # deal with possible scale parametrization
  if(is.null(names(parameters))){
    if("scale" %in% names(parameters)){
      parameters <- .check_and_name_parameters(parameters, c("scale"), "exp")
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
  }

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("rate"), "exp")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$rate, "rate")
  .check_parameter_positive(parameters$rate, "rate")

  # add the values to the output
  output$distribution <- "exp"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_beta      <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("alpha", "beta"), "beta")
  truncation <- .check_and_set_truncation(truncation, lower = 0, upper = 1)

  # check individual parameters
  .check_parameter(parameters$alpha, "alpha")
  .check_parameter(parameters$beta,  "beta")
  .check_parameter_positive(parameters$alpha, "alpha")
  .check_parameter_positive(parameters$beta,  "beta")

  # add the values to the output
  output$distribution <- "beta"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_uniform   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("a", "b"), "uniform")

  # check individual parameters
  .check_parameter(parameters$a, "a")
  .check_parameter(parameters$b, "b")

  if(parameters$a >= parameters$b)
    stop("Parameter 'a' must be lower than the parameter 'b'.")

  # add the values to the output
  output$distribution <- "uniform"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$a, upper = parameters$b)

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_point     <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location"), "point")

  # check individual parameters
  .check_parameter(parameters$location, "location")

  # add the values to the output
  output$distribution <- "point"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$location, upper = parameters$location)

  class(output) <- c("prior", "prior.simple", "prior.point")

  return(output)
}
.prior_mnormal   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("mean", "sd", "K"), "multivariate normal")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$mean, "mean")
  .check_parameter(parameters$sd,   "sd")
  .check_parameter_positive(parameters$sd, "sd")
  .check_parameter_dimensions(parameters$K, "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor

  # add the values to the output
  output$distribution <- "mnormal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}
.prior_mcauchy   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "K"), "multivariate Cauchy")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_dimensions(parameters$K,   "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor

  # deal with as with a t-distribution
  parameters$df <- 1

  output$distribution <- "mt"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}
.prior_mt        <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "df", "K"), "multivariate student-t")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter(parameters$df,       "df")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_positive(parameters$df,    "df")
  .check_parameter_dimensions(parameters$K,   "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor


  # add the values to the output
  output$distribution <- "mt"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}

.prior_weightfunction_two.sided <- function(parameters){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("steps", "alpha"), "two-sided weightfunction")

  # check individual parameters
  .check_parameter_weigthfunction(parameters$steps, parameters$alpha)

  # reverse the ordering of alpha and weights - to correspond to ordering on t-statistics
  parameters$steps <- rev(parameters$steps)
  parameters$alpha <- rev(parameters$alpha)

  # add the values to the output
  output$distribution <- "two.sided"
  output$parameters   <- parameters
  output$truncation   <- list(lower = 0, upper = 1)

  return(output)
}
.prior_weightfunction_one.sided <- function(parameters){

  output <- list()

  if(length(parameters) == 2){

    # check overall settings
    parameters <- .check_and_name_parameters(parameters, c("steps", "alpha"), "one-sided weightfunction")

    # check individual parameters
    .check_parameter_weigthfunction(parameters$steps, parameters$alpha)

    # reverse the ordering of alpha and weights - to correspond to ordering on t-statistics
    parameters$steps <- rev(parameters$steps)
    parameters$alpha <- rev(parameters$alpha)

    # add the values to the output
    output$distribution <- "one.sided"
    output$parameters   <- parameters
    output$truncation   <- list(lower = 0, upper = 1)

  }else{

    # check overall settings
    parameters <- .check_and_name_parameters(parameters, c("steps", "alpha1", "alpha2"), "one-sided weightfunction")

    # check individual parameters
    .check_parameter_weigthfunction(parameters$steps, alpha1 = parameters$alpha1, alpha2 = parameters$alpha2)

    # reverse the ordering of alpha and weights - to correspond to ordering on t-statistics
    parameters$steps  <- rev(parameters$steps)
    parameters$alpha1 <- rev(parameters$alpha1)

    # add the values to the output
    output$distribution <- "one.sided"
    output$parameters   <- parameters
    output$truncation   <- list(lower = 0, upper = 1)

  }

  return(output)
}
.prior_weightfunction_two.sided.fixed <- function(parameters){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("steps", "omega"), "two-sided.fixed weightfunction")

  # check individual parameters
  .check_parameter_weigthfunction(parameters$steps, omega = parameters$omega)

  # reverse the ordering of alpha and weights - to correspond to ordering on t-statistics
  parameters$steps <- rev(parameters$steps)
  parameters$omega <- rev(parameters$omega)

  # add the values to the output
  output$distribution <- "two.sided.fixed"
  output$parameters   <- parameters
  output$truncation   <- list(lower = 0, upper = 1)

  return(output)
}
.prior_weightfunction_one.sided.fixed <- function(parameters){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("steps", "omega"), "two-sided.fixed weightfunction")

  # check individual parameters
  .check_parameter_weigthfunction(parameters$steps, omega = parameters$omega)

  # reverse the ordering of alpha and weights - to correspond to ordering on t-statistics
  parameters$steps <- rev(parameters$steps)
  parameters$omega <- rev(parameters$omega)

  # add the values to the output
  output$distribution <- "one.sided.fixed"
  output$parameters   <- parameters
  output$truncation   <- list(lower = 0, upper = 1)

  return(output)
}

.prior_normal.orthonormal <- function(parameters){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, "sd", "normal.orthonormal")

  # check individual parameters
  .check_parameter(parameters$sd,   "sd")
  .check_parameter_positive(parameters$sd, "sd")

  # add the values to the output
  output$distribution <- "normal.orthonormal"
  output$parameters   <- parameters

  return(output)
}

#### elementary prior related functions ####
#' @title Elementary prior related functions
#'
#' @description Density (pdf / lpdf), distribution
#' function (cdf / ccdf), quantile function (quant),
#' random generation (rng), mean, standard deviation (sd),
#' and marginal variants of the functions (mpdf, mlpf, mcdf,
#' mccdf, mquant) for prior distributions.
#'
#' @param x prior distribution
#' @param y vector of observations
#' @param q vector or matrix of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # generate a random sample from the prior
#' rng(p1, 10)
#'
#' # compute cumulative density function
#' cdf(p1, 0)
#'
#' # obtain quantile
#' quant(p1, .5)
#'
#' # compute probability density
#' pdf(p1, c(0, 1, 2))
#'
#' @return \code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
#' the (marginal) density and the log of (marginal) density,
#' \code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
#' (marginal) distribution and the complement of (marginal) distribution function,
#' \code{quant} (\code{mquant}) give the (marginal) quantile function,
#' and \code{rng} generates random deviates for an object of class 'prior'.
#'
#' @exportS3Method rng prior
#' @exportS3Method cdf prior
#' @exportS3Method ccdf prior
#' @exportS3Method quant prior
#' @exportS3Method lpdf prior
#' @exportS3Method pdf prior
#' @exportS3Method mcdf prior
#' @exportS3Method mccdf prior
#' @exportS3Method mquant prior
#' @exportS3Method mlpdf prior
#' @exportS3Method mpdf prior
#' @name prior_functions
NULL


#### joint distribution functions ####
#' @rdname prior_functions
rng.prior   <- function(x, n, ...){

  prior <- x

  .check_n(n)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    x <- NULL
    # guesstimate the number of samples needed before the truncation
    nn <- round(n * 1/.prior_C(prior) * 1.10)

    while(length(x) < n){
      temp_x <- switch(
        prior[["distribution"]],
        "normal"    = stats::rnorm(nn, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]]),
        "lognormal" = stats::rlnorm(nn, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]]),
        "t"         = extraDistr::rlst(nn, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]]),
        "gamma"     = stats::rgamma(nn, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]]),
        "invgamma"  = extraDistr::rinvgamma(nn, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]]),
        "beta"      = stats::rbeta(nn, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]]),
        "exp"       = stats::rexp(nn, rate = prior$parameters[["rate"]]),
        "uniform"   = stats::runif(nn, min = prior$parameters[["a"]], max = prior$parameters[["b"]]),
        "point"     = rpoint(n, location = prior$parameters[["location"]])
      )
      x <- c(x, temp_x[temp_x >= prior$truncation[["lower"]] & temp_x <= prior$truncation[["upper"]]])
    }

    x <- x[1:n]

  }else if(is.prior.orthonormal(prior)){

    if(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["mean"]],
      "mt"      = prior$parameters[["location"]]) != 0){
      stop("the orthonormal prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- attr(prior, "levels") - 1
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par1 <- rep(0, prior$parameter[["K"]])
    par2 <- diag(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["sd"]]^2,
      "mt"      = prior$parameter[["scale"]]^2
    ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
    )

    x <- x %*% t(contr.orthonormal(1:(prior$parameters[["K"]] + 1)))

  }else if(is.prior.vector(prior)){

    # TODO: generalize this to priors with covariances
    par1 <- rep(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]]
    ), length = prior$parameter[["K"]])
    par2 <- diag(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["sd"]]^2,
      "mt"      = prior$parameter[["scale"]]^2
    ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
    )

  }else if(is.prior.weightfunction(prior)){

    x <- switch(
      prior[["distribution"]],
      "one.sided" = rone.sided(n, alpha = if(!is.null(prior$parameters[["alpha"]])) prior$parameters[["alpha"]], alpha1 = if(!is.null(prior$parameters[["alpha1"]])) prior$parameters[["alpha1"]], alpha2 = if(!is.null(prior$parameters[["alpha2"]])) prior$parameters[["alpha2"]]),
      "two.sided" = rtwo.sided(n, alpha = prior$parameters[["alpha"]]),
      "one.sided.fixed" = rone.sided_fixed(n, omega = prior$parameters[["omega"]]),
      "two.sided.fixed" = rtwo.sided_fixed(n, omega = prior$parameters[["omega"]])
    )

  }

  return(x)
}
#' @rdname prior_functions
cdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    p <- rep(NA, length(q))

    # deal with values below the truncation point
    q_lower    <- q < prior$truncation[["lower"]]
    p[q_lower] <- 0

    # compute for the values above truncation point
    p[!q_lower] <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(q[!q_lower], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(q[!q_lower], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(q[!q_lower], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(q[!q_lower], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(q[!q_lower], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(q[!q_lower], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(q[!q_lower], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(q[!q_lower], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "point"     = ppoint(q[!q_lower], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

    # deal with truncation
    if(prior[["distribution"]] != "point"){
      p[!q_lower] <- (p[!q_lower] - .prior_C1(prior)) /.prior_C(prior)
    }

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal cdfs are implemented for prior weightfunctions.")

  }

  return(p)
}
#' @rdname prior_functions
ccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    p <- rep(NA, length(q))

    # deal with values above the truncation point
    q_higher    <- q > prior$truncation[["upper"]]
    p[q_higher] <- 0

    # compute for the values belove truncation point
    p[!q_higher] <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(q[!q_higher], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = FALSE, log.p = FALSE),
      "lognormal" = stats::plnorm(q[!q_higher], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = FALSE, log.p = FALSE),
      "t"         = extraDistr::plst(q[!q_higher], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = FALSE, log.p = FALSE),
      "gamma"     = stats::pgamma(q[!q_higher], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = FALSE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(q[!q_higher], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = FALSE, log.p = FALSE),
      "beta"      = stats::pbeta(q[!q_higher], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = FALSE, log.p = FALSE),
      "exp"       = stats::pexp(q[!q_higher], rate = prior$parameters[["rate"]], lower.tail = FALSE, log.p = FALSE),
      "uniform"   = stats::punif(q[!q_higher], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = FALSE, log.p = FALSE),
      "point"     = ppoint(q[!q_higher], location = prior$parameters[["location"]], lower.tail = FALSE, log.p = FALSE)
    )

    # deal with truncation
    if(prior[["distribution"]] != "point"){
      p[!q_higher] <- (p[!q_higher] - (1 - .prior_C2(prior))) /.prior_C(prior)
    }

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal ccdf functions are implemented for prior weightfunctions.")

  }

  return(p)
}
#' @rdname prior_functions
lpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    log_lik <- switch(
      prior[["distribution"]],
      "normal"    = stats::dnorm(x, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], log = TRUE),
      "lognormal" = stats::dlnorm(x, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], log = TRUE),
      "t"         = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], log = TRUE),
      "gamma"     = stats::dgamma(x, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], log = TRUE),
      "invgamma"  = extraDistr::dinvgamma(x, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], log = TRUE),
      "beta"      = stats::dbeta(x, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], log = TRUE),
      "exp"       = stats::dexp(x, rate = prior$parameters[["rate"]], log = TRUE),
      "uniform"   = stats::dunif(x, min = prior$parameters[["a"]], max = prior$parameters[["b"]], log = TRUE),
      "point"     = dpoint(x, location = prior$parameters[["location"]], log = TRUE)
    )

    log_lik[x < prior$truncation[["lower"]] | x > prior$truncation[["upper"]]] <- -Inf
    if(prior[["distribution"]] != "point"){
      log_lik <- log_lik - log(.prior_C(prior))
    }

  }else if(is.prior.vector(prior)){

    # TODO: generalize this to priors with covariances
    par1 <- rep(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]]
    ), length = prior$parameter[["K"]])
    par2 <- diag(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["sd"]]^2,
      "mt"      = prior$parameter[["scale"]]^2
    ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::dmvnorm(x, mean = par1, sigma = par2, log = TRUE),
      "mt"         = mvtnorm::dmvt(x, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted", log = TRUE),
    )

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal lpdf are implemented for prior weightfunctions.")

  }

  return(log_lik)
}
#' @rdname prior_functions
pdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  log_lik <- lpdf(prior, x)
  lik     <- exp(log_lik)

  return(lik)
}
#' @rdname prior_functions
quant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    if(.is_prior_default_range(prior)){

      q <- switch(
        prior[["distribution"]],
        "normal"    = stats::qnorm(p, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
        "lognormal" = stats::qlnorm(p, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
        "t"         = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
        "gamma"     = stats::qgamma(p, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
        "invgamma"  = extraDistr::qinvgamma(p, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
        "beta"      = stats::qbeta(p, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
        "exp"       = stats::qexp(p, rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
        "uniform"   = stats::qunif(p, min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
        "point"     = qpoint(p, location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
      )

    }else{

      if(!is.infinite(prior$truncation[["upper"]]) & !is.infinite(prior$truncation[["upper"]])){
        start_value <- (prior$truncation[["upper"]]  - prior$truncation[["lower"]]) / 2
      }else if(!is.infinite(prior$truncation[["upper"]])){
        start_value <- prior$truncation[["upper"]] - 1
      }else if(!is.infinite(prior$truncation[["lower"]])){
        start_value <- prior$truncation[["lower"]] + 1
      }else{
        start_value <- 0
      }

      q <- sapply(p, function(p_i){
        stats::optim(
          par     = start_value,
          fn      = function(x, prior, p_i)(cdf(prior, x) - p_i)^2,
          lower   = prior$truncation[["lower"]],
          upper   = prior$truncation[["upper"]],
          prior   = prior,
          p_i     = p_i,
          method  = "L-BFGS-B",
          control = list(
            factr = 1e3
          )
        )$par
      })

    }


  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal quantile functions are implemented for prior weightfunctions.")

  }

  return(q)
}
.prior_C1 <- function(prior){

  if(is.prior.simple(prior)){

    C1 <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(prior$truncation[["lower"]], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(prior$truncation[["lower"]], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(prior$truncation[["lower"]], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(prior$truncation[["lower"]], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(prior$truncation[["lower"]], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["lower"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["lower"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["lower"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "point"     = ppoint(prior$truncation[["lower"]], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

  }

  return(C1)
}
.prior_C2 <- function(prior){

  if(is.prior.simple(prior)){

    C2 <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(prior$truncation[["upper"]], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(prior$truncation[["upper"]], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(prior$truncation[["upper"]], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(prior$truncation[["upper"]], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(prior$truncation[["upper"]], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["upper"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["upper"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["upper"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "point"     = ppoint(prior$truncation[["upper"]], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

  }

  return(C2)
}
.prior_C  <- function(prior){

  if(is.prior.simple(prior)){

    C <- .prior_C2(prior) - .prior_C1(prior)

  }

  return(C)
}

# tools
.is_prior_default_range   <- function(prior){
  default_range <- switch(
    prior[["distribution"]],
    "normal"    = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "lognormal" = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "t"         = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "gamma"     = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "invgamma"  = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "beta"      = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & isTRUE(all.equal(prior$truncation[["upper"]], 1)),
    "exp"       = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "uniform"   = TRUE,
    "point"     = TRUE,
    "one.sided" = TRUE,
    "two.sided" = TRUE,
    "one.sided.fixed" = TRUE,
    "two.sided.fixed" = TRUE,
    "none"            = TRUE
  )
}

#### marginal distribution functions ####
# weightfunctions require just marginals for plotting and etc
# (also, the joint are not implemented in general, since they differ between JAGS/Stan implementations)
#' @rdname prior_functions
mcdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    p <- cdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    p <- switch(
      prior[["distribution"]],
      "one.sided" = mpone.sided(q, alpha = if(!is.null(prior$parameters[["alpha"]])) prior$parameters[["alpha"]], alpha1 = if(!is.null(prior$parameters[["alpha1"]])) prior$parameters[["alpha1"]], alpha2 = if(!is.null(prior$parameters[["alpha2"]])) prior$parameters[["alpha2"]]),
      "two.sided" = mptwo.sided(q, alpha = prior$parameters[["alpha"]]),
      "one.sided.fixed" = mpone.sided_fixed(q, omega = prior$parameters[["omega"]]),
      "two.sided.fixed" = mptwo.sided_fixed(q, omega = prior$parameters[["omega"]])
    )

  }else if(is.prior.orthonormal(prior)){

    if(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["mean"]],
      "mt"      = prior$parameters[["location"]]) != 0){
      stop("the orthonormal prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- attr(prior, "levels") - 1
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["sd"]],
      "mt"      = prior$parameters[["scale"]]) )^2 ))

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2)
    )

  }

  return(p)
}
#' @rdname prior_functions
mccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    p <- ccdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    p <- switch(
      prior[["distribution"]],
      "one.sided" = mpone.sided(q, alpha = if(!is.null(prior$parameters[["alpha"]])) prior$parameters[["alpha"]], alpha1 = if(!is.null(prior$parameters[["alpha1"]])) prior$parameters[["alpha1"]], alpha2 = if(!is.null(prior$parameters[["alpha2"]])) prior$parameters[["alpha2"]], lower.tail = FALSE),
      "two.sided" = mptwo.sided(q, alpha = prior$parameters[["alpha"]], lower.tail = FALSE),
      "one.sided.fixed" = mpone.sided_fixed(q, omega = prior$parameters[["omega"]], lower.tail = FALSE),
      "two.sided.fixed" = mptwo.sided_fixed(q, omega = prior$parameters[["omega"]], lower.tail = FALSE)
    )

  }else if(is.prior.orthonormal(prior)){

    if(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["mean"]],
      "mt"      = prior$parameters[["location"]]) != 0){
      stop("the orthonormal prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- attr(prior, "levels") - 1
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["sd"]],
      "mt"      = prior$parameters[["scale"]]) )^2 ))

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2, lower.tail = FALSE),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2, lower.tail = FALSE)
    )
  }

  return(p)
}
#' @rdname prior_functions
mlpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    log_lik <- lpdf(prior, x)

  }else if(is.prior.weightfunction(prior)){

    log_lik <- switch(
      prior[["distribution"]],
      "one.sided" = mdone.sided(x, alpha = if(!is.null(prior$parameters[["alpha"]])) prior$parameters[["alpha"]], alpha1 = if(!is.null(prior$parameters[["alpha1"]])) prior$parameters[["alpha1"]], alpha2 = if(!is.null(prior$parameters[["alpha2"]])) prior$parameters[["alpha2"]], log = TRUE),
      "two.sided" = mdtwo.sided(x, alpha = prior$parameters[["alpha"]], log = TRUE),
      "one.sided.fixed" = mdone.sided_fixed(x, omega = prior$parameters[["omega"]], log = TRUE),
      "two.sided.fixed" = mdtwo.sided_fixed(x, omega = prior$parameters[["omega"]], log = TRUE),
    )

  }else if(is.prior.orthonormal(prior)){

    if(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["mean"]],
      "mt"      = prior$parameters[["location"]]) != 0){
      stop("the orthonormal prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- attr(prior, "levels") - 1
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["sd"]],
      "mt"      = prior$parameters[["scale"]]) )^2 ))

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal" = stats::dnorm(x, mean = 0, sd = par2, log = TRUE),
      "mt"      = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = 0, sigma = par2, log = TRUE)
    )
  }

  return(log_lik)
}
#' @rdname prior_functions
mpdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  log_lik <- mlpdf(prior, x)
  lik     <- exp(log_lik)

  return(lik)
}
#' @rdname prior_functions
mquant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.simple(prior)){

    q <- quant(prior, p)

  }else if(is.prior.weightfunction(prior)){

    q <- switch(
      prior[["distribution"]],
      "one.sided" = mqone.sided(p, alpha = if(!is.null(prior$parameters[["alpha"]])) prior$parameters[["alpha"]], alpha1 = if(!is.null(prior$parameters[["alpha1"]])) prior$parameters[["alpha1"]], alpha2 = if(!is.null(prior$parameters[["alpha2"]])) prior$parameters[["alpha2"]]),
      "two.sided" = mqtwo.sided(p, alpha = prior$parameters[["alpha"]]),
      "one.sided.fixed" = mqone.sided_fixed(p, omega = prior$parameters[["omega"]]),
      "two.sided.fixed" = mqtwo.sided_fixed(p, omega = prior$parameters[["omega"]])
    )

  }else if(is.prior.orthonormal(prior)){

    if(switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["mean"]],
      "mt"      = prior$parameters[["location"]]) != 0){
      stop("the orthonormal prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- attr(prior, "levels") - 1
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
      prior[["distribution"]],
      "mnormal" = prior$parameters[["sd"]],
      "mt"      = prior$parameters[["scale"]]) )^2 ))

    q <- switch(
      prior[["distribution"]],
      "mnormal" = stats::qnorm(p, mean = 0, sd = par2),
      "mt"      = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = 0, sigma = par2)
    )

  }

  return(q)
}

### create generic method
#' @title Creates generics for common statistical functions
#'
#' @description Density (pdf / lpdf), distribution
#' function (cdf / ccdf), quantile function (quant),
#' random generation (rng), mean, standard deviation (sd),
#' and marginal variants of the functions (mpdf, mlpf, mcdf,
#' mccdf, mquant).
#'
#' @param x main argument
#' @param ... unused arguments
#'
#' @return \code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
#' the (marginal) density and the log of (marginal) density,
#' \code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
#' (marginal) distribution and the complement of (marginal) distribution function,
#' \code{quant} (\code{mquant}) give the (marginal) quantile function,
#' and \code{rng} generates random deviates for an object of class 'prior'.
#'
#' The \code{pdf} function proceeds to PDF graphics device if \code{x} is a character.
#'
#' @export rng
#' @export cdf
#' @export ccdf
#' @export quant
#' @export lpdf
#' @export pdf
#' @export mcdf
#' @export mccdf
#' @export mquant
#' @export mlpdf
#' @export mpdf
#' @importFrom grDevices pdf
#' @name prior_functions_methods
NULL

#' @rdname prior_functions_methods
rng    <- function(x, ...){
  UseMethod("rng")
}
#' @rdname prior_functions_methods
cdf    <- function(x, ...){
  UseMethod("cdf")
}
#' @rdname prior_functions_methods
ccdf   <- function(x, ...){
  UseMethod("ccdf")
}
#' @rdname prior_functions_methods
quant  <- function(x, ...){
  UseMethod("quant")
}
#' @rdname prior_functions_methods
lpdf   <- function(x, ...){
  UseMethod("lpdf")
}
#' @rdname prior_functions_methods
pdf    <- function(x, ...){
  UseMethod("pdf")
}
#' @rdname prior_functions_methods
mcdf   <- function(x, ...){
  UseMethod("mcdf")
}
#' @rdname prior_functions_methods
mccdf  <- function(x, ...){
  UseMethod("mccdf")
}
#' @rdname prior_functions_methods
mquant <- function(x, ...){
  UseMethod("mquant")
}
#' @rdname prior_functions_methods
mlpdf  <- function(x, ...){
  UseMethod("mlpdf")
}
#' @rdname prior_functions_methods
mpdf   <- function(x, ...){
  UseMethod("mpdf")
}

#' @export
pdf.default  <- function(x, ...){
  grDevices::pdf(x, ...)
}


#### additional prior related function ####
#' @title Prior mean
#'
#' @description Computes mean of a prior
#' distribution.
#'
#' @param x a prior
#' @param ... unused
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute mean of the prior distribution
#' mean(p1)
#'
#' @return a mean of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @rdname mean.prior
#' @export
mean.prior   <- function(x, ...){

  .check_prior(x, "x")

  if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      m <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["mean"]],
        "lognormal" = exp(x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2/2),
        "t"         = ifelse(x$parameters[["df"]] > 1, x$parameters[["location"]], NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]],
        "invgamma"  = ifelse(x$parameters[["shape"]] > 1, x$parameters[["scale"]]/(x$parameters[["shape"]] - 1), NaN),
        "beta"      = x$parameters[["alpha"]] / (x$parameters[["alpha"]] + x$parameters[["beta"]]),
        "exp"       = 1 / x$parameters[["rate"]],
        "uniform"   = (x$parameters[["a"]] + x$parameters[["b"]]) / 2,
        "point"     = x$parameters[["location"]]
      )

    }else{

      # check for undefined values
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 1){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 1){
          return(NaN)
        }
      }

      m <- stats::integrate(
        f       = function(x, prior) x * pdf(prior, x),
        lower   = x$truncation[["lower"]],
        upper   = x$truncation[["upper"]],
        prior   = x
      )$value

    }


  }else if(is.prior.weightfunction(x)){

    m <- .mean.weightfunction(x)

  }

  return(m)
}


### create generic methods from sd and var
#' @title Creates generic for sd function
#'
#' @param x main argument
#' @param ... additional arguments
#'
#' @return \code{sd} returns a standard deviation
#' of the supplied object (if it is either a numeric vector
#' or an object of class 'prior').
#'
#' @seealso \link[stats]{sd}
#' @export
sd  <- function(x, ...){
  UseMethod("sd")
}

#' @title Creates generic for var function
#'
#' @param x main argument
#' @param ... additional arguments
#'
#' @return \code{var} returns a variance
#' of the supplied object (if it is either a numeric vector
#' or an object of class 'prior').
#'
#' @seealso \link[stats]{cor}
#' @export
var <- function(x, ...){
  UseMethod("var")
}

#' @export
sd.default  <- function(x, ...){
  stats::sd(x, ...)
}
#' @export
var.default <- function(x, ...){
  stats::var(x, ...)
}


#' @title Prior var
#'
#' @description Computes variance
#' of a prior distribution.
#'
#' @param x a prior
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute variance of the prior distribution
#' var(p1)
#'
#' @return a variance of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @importFrom stats var
#' @rdname var.prior
#' @export
var.prior   <- function(x, ...){

  .check_prior(x, "x")

  if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      var <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["sd"]]^2,
        "lognormal" = (exp(x$parameters[["sdlog"]]^2) - 1) * exp(2 * x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2),
        "t"         = ifelse(x$parameters[["df"]] > 2, x$parameters[["scale"]]^2 * x$parameters[["df"]] / (x$parameters[["df"]] - 2), NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]]^2,
        "invgamma"  = ifelse(x$parameters[["shape"]] > 2, x$parameters[["scale"]]^2 / (x$parameters[["shape"]] - 1)^2 * (x$parameters[["shape"]] - 2), NaN),
        "beta"      = (x$parameters[["alpha"]] * x$parameters[["beta"]]) / ((x$parameters[["alpha"]] + x$parameters[["beta"]])^2 * (x$parameters[["alpha"]] + x$parameters[["beta"]] + 1)),
        "exp"       = 1 / x$parameters[["rate"]]^2,
        "uniform"   = (x$parameters[["b"]] - x$parameters[["a"]])^2 / 12,
        "point"     = 0
      )

    }else{

      # check for undefined values
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 2){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 2){
          return(NaN)
        }
      }

      E2 <- stats::integrate(
        f       = function(x, prior) x^2 * pdf(prior, x),
        lower   = x$truncation[["lower"]],
        upper   = x$truncation[["upper"]],
        prior   = x
      )$value

      var <- E2 - mean(x)^2
    }


  }else if(is.prior.weightfunction(x)){

    var <- .var.weightfunction(x)

  }

  return(var)
}

#' @title Prior sd
#'
#' @description Computes standard deviation
#' of a prior distribution.
#'
#' @param x a prior
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute sd of the prior distribution
#' sd(p1)
#'
#' @return a standard deviation of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @importFrom stats sd
#' @rdname sd.prior
#' @export
sd.prior     <- function(x, ...){

  .check_prior(x, "x")

  sd <- sqrt(var(x))

  return(sd)
}

.mean.weightfunction <- function(prior){

  if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    m <- prior$parameters[["omega"]]

  }else if(all(names(prior[["parameters"]]) %in% c("steps", "alpha"))){

    alpha       <- prior$parameters[["alpha"]]
    alpha_alpha <- cumsum(alpha)[-length(alpha)]
    alpha_beta  <- cumsum(alpha[length(alpha):1])[(length(alpha)-1):1]

    m <- c(alpha_alpha/(alpha_alpha + alpha_beta), 1)

  }else{

    stop("Analytical solutions for the mean of one-sided non-monotonic weights are not implemented.")

  }

  return(m)
}
.var.weightfunction  <- function(prior){

  if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    m <- rep(0, length(prior$parameters[["omega"]]))

  }else if(all(names(prior[["parameters"]]) %in% c("steps", "alpha"))){

    alpha       <- prior$parameters[["alpha"]]
    alpha_alpha <- cumsum(alpha)[-length(alpha)]
    alpha_beta  <- cumsum(alpha[length(alpha):1])[(length(alpha)-1):1]

    m <- c((alpha_alpha * alpha_beta) / ((alpha_alpha + alpha_beta)^2 * (alpha_alpha + alpha_beta + 1)), 0)

  }else{

    stop("Analytical solutions for the var/sd of one-sided non-monotonic weights are not implemented.")

  }

  return(m)
}
