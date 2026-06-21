#### additional prior related function ####
#' @title Prior mean
#'
#' @description Computes mean of a prior
#' distribution. (In case of orthonormal prior distributions
#' for factors, the mean of for the deviations from intercept
#' is returned.)
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

  if(is.prior.spike_and_slab(x)){

    m <- mean(.get_spike_and_slab_variable(x)) * mean(.get_spike_and_slab_inclusion(x))

  }else if(is.prior.mixture(x)){

    stop("No mean is implemented for prior mixtures.")

  }else if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      m <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["mean"]],
        "lognormal" = exp(x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2/2),
        "t"         = ifelse(x$parameters[["df"]] > 1, x$parameters[["location"]], NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]],
        "invgamma"  = ifelse(x$parameters[["shape"]] > 1, x$parameters[["scale"]]/(x$parameters[["shape"]] - 1), NaN),
        "moment"    = x$parameters[["location"]],
        "invmoment" = ifelse(x$parameters[["df"]] > 1, x$parameters[["location"]], NaN),
        "beta"      = x$parameters[["alpha"]] / (x$parameters[["alpha"]] + x$parameters[["beta"]]),
        "bernoulli" = x$parameters[["probability"]],
        "exp"       = 1 / x$parameters[["rate"]],
        "uniform"   = (x$parameters[["a"]] + x$parameters[["b"]]) / 2,
        "point"     = x$parameters[["location"]]
      )

    }else{

      if(is.prior.discrete(x)){
        discrete <- .prior_simple_truncated_discrete(x)
        return(sum(discrete[["support"]] * discrete[["prob"]]))
      }

      # check for undefined values when the problematic tail is still unbounded
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 1 &&
           (is.infinite(x$truncation[["lower"]]) || is.infinite(x$truncation[["upper"]]))){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 1 && is.infinite(x$truncation[["upper"]])){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invmoment"){
        if(x$parameters[["df"]] <= 1 &&
           (is.infinite(x$truncation[["lower"]]) || is.infinite(x$truncation[["upper"]]))){
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

  }else if(is_prior_phacking(x) || is_prior_bias(x)){

    .selection_prior_stop_unsupported_generic("mean", x)

  }else if(is.prior.vector(x) && !is.prior.factor(x)){

    m <- .prior_vector_mean(x)

  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){

    par1 <- switch(
      x[["distribution"]],
      "mnormal" = x$parameter[["mean"]],
      "mt"      = x$parameter[["location"]],
      "mpoint"  = x$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(x[["distribution"]] == "mt" && x$parameters[["df"]] <= 1){
      m <- NaN
    }else{
      # orthonormal prior distributions must be centered
      m <- 0
    }

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

  if(is.prior.spike_and_slab(x)){

    inclusion_mean <- mean(.get_spike_and_slab_inclusion(x))
    variable_mean  <- mean(.get_spike_and_slab_variable(x))
    variable_E2    <- var(.get_spike_and_slab_variable(x)) + variable_mean^2

    var <- inclusion_mean * variable_E2 - (inclusion_mean * variable_mean)^2

  }else if(is.prior.mixture(x)){

    stop("No var is implemented for prior mixtures.")

  }else if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      var <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["sd"]]^2,
        "lognormal" = (exp(x$parameters[["sdlog"]]^2) - 1) * exp(2 * x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2),
        "t"         = ifelse(x$parameters[["df"]] > 2, x$parameters[["scale"]]^2 * x$parameters[["df"]] / (x$parameters[["df"]] - 2), NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]]^2,
        "invgamma"  = ifelse(x$parameters[["shape"]] > 2, x$parameters[["scale"]]^2 / ((x$parameters[["shape"]] - 1)^2 * (x$parameters[["shape"]] - 2)), NaN),
        "moment"    = x$parameters[["tau"]] * (2 * x$parameters[["order"]] + 1),
        "invmoment" = ifelse(
          x$parameters[["df"]] > 2,
          x$parameters[["tau"]] * exp(lgamma(x$parameters[["df"]] / (2 * x$parameters[["order"]]) - 1 / x$parameters[["order"]]) -
            lgamma(x$parameters[["df"]] / (2 * x$parameters[["order"]]))),
          NaN
        ),
        "beta"      = (x$parameters[["alpha"]] * x$parameters[["beta"]]) / ((x$parameters[["alpha"]] + x$parameters[["beta"]])^2 * (x$parameters[["alpha"]] + x$parameters[["beta"]] + 1)),
        "bernoulli" = (x$parameters[["probability"]] * (1 - x$parameters[["probability"]]) ),
        "exp"       = 1 / x$parameters[["rate"]]^2,
        "uniform"   = (x$parameters[["b"]] - x$parameters[["a"]])^2 / 12,
        "point"     = 0
      )

    }else{

      if(is.prior.discrete(x)){
        discrete <- .prior_simple_truncated_discrete(x)
        E2       <- sum(discrete[["support"]]^2 * discrete[["prob"]])
        return(E2 - mean(x)^2)
      }

      # check for undefined values when the problematic tail is still unbounded
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 2 &&
           (is.infinite(x$truncation[["lower"]]) || is.infinite(x$truncation[["upper"]]))){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 2 && is.infinite(x$truncation[["upper"]])){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invmoment"){
        if(x$parameters[["df"]] <= 2 &&
           (is.infinite(x$truncation[["lower"]]) || is.infinite(x$truncation[["upper"]]))){
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

  }else if(is_prior_phacking(x) || is_prior_bias(x)){

    .selection_prior_stop_unsupported_generic("variance", x)

  }else if(is.prior.vector(x) && !is.prior.factor(x)){

    var <- .prior_vector_var(x)

  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){

    par1 <- switch(
      x[["distribution"]],
      "mnormal" = x$parameter[["mean"]],
      "mt"      = x$parameter[["location"]],
      "mpoint"  = x$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(x[["distribution"]] == "mpoint"){
      var <- 0
    }else if(x[["distribution"]] == "mt" && x$parameters[["df"]] <= 2){
      var <- NaN
    }else{


      if(is.na(x$parameters[["K"]]) && !is.null(attr(x, "levels"))){
        x$parameters[["K"]] <- .get_prior_factor_levels(x)
      }else if(is.na(x$parameters[["K"]])){
        x$parameters[["K"]] <- 1
        warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
      }


      if(is.prior.orthonormal(x)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(x$parameters[["K"]] + 1))[1,] * switch(
          x[["distribution"]],
          "mnormal" = x$parameters[["sd"]],
          "mt"      = x$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(x)){
        par2 <- sqrt(sum( (contr.meandif(1:(x$parameters[["K"]] + 1))[1,] * switch(
          x[["distribution"]],
          "mnormal" = x$parameters[["sd"]],
          "mt"      = x$parameters[["scale"]]) )^2 ))
      }

      # use the univariate functions
      var <- switch(
        x[["distribution"]],
        "mnormal" = var.prior(prior("normal", parameters = list(mean = 0, sd = par2))),
        "mt"      = var.prior(prior("t",      parameters = list(location = 0, scale = par2, df = x$parameters[["df"]]))))
    }

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

  components <- .weightfunction_marginal_components(prior)
  m <- vapply(components, .prior_weightfunction_component_mean, numeric(1))

  return(m)
}
.var.weightfunction  <- function(prior){

  components <- .weightfunction_marginal_components(prior)
  m <- vapply(components, .prior_weightfunction_component_var, numeric(1))

  return(m)
}
