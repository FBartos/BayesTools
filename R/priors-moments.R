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

      if(x[["distribution"]] == "normal"){
        return(.prior_normal_truncated_moments(x, "mean"))
      }

      prior_lpdf <- if(x$distribution %in% c("moment", "invmoment")) .prior_simple_lpdf_evaluator(x) else NULL
      integrand <- if(is.null(prior_lpdf)){
        function(x, prior) x * pdf(prior, x)
      }else{
        function(x, prior) x * exp(prior_lpdf(x))
      }
      m <- stats::integrate(
        f       = integrand,
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

  }else{

    .prior_stop_unsupported_method("mean", x)

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
#' @details The variance of a truncated normal prior is computed in closed
#' form. Its terms cancel far in a tail and for narrow truncation intervals,
#' so the variance is returned only when its estimated relative rounding error
#' in double precision is at most 1e-6; otherwise the function stops with an
#' error. For a one-sided truncation, the error is raised beyond about 33
#' prior standard deviations from the mean. The same applies to \code{sd()}
#' and to spike-and-slab priors with such a slab. The variance can then be
#' estimated from \code{rng()} draws of the prior, which sample the truncated
#' distribution exactly, subject to Monte Carlo error.
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

      if(x[["distribution"]] == "normal"){
        return(.prior_normal_truncated_moments(x, "var"))
      }

      prior_lpdf <- if(x$distribution %in% c("moment", "invmoment")) .prior_simple_lpdf_evaluator(x) else NULL
      integrand <- if(is.null(prior_lpdf)){
        function(x, prior) x^2 * pdf(prior, x)
      }else{
        function(x, prior) x^2 * exp(prior_lpdf(x))
      }
      E2 <- stats::integrate(
        f       = integrand,
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

  }else{

    .prior_stop_unsupported_method("variance", x)

  }

  return(var)
}

#' @title Prior sd
#'
#' @description Computes standard deviation
#' of a prior distribution.
#'
#' @details The standard deviation is the square root of \code{var()}; see
#' \code{\link[=var.prior]{var()}} for the precision check of truncated normal
#' priors.
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

.prior_normal_truncated_moments <- function(prior, moment = c("mean", "var")){

  moment <- match.arg(moment)

  mu    <- prior$parameters[["mean"]]
  sigma <- prior$parameters[["sd"]]
  if(identical(moment, "mean")){
    return(mu + sigma * .prior_normal_truncated_terms(prior)$shift)
  }

  # The variance is returned only when its estimated relative rounding error
  # is at most 1e-6; far tails and narrow intervals exhaust double precision.
  variance <- .prior_normal_truncated_variance(prior)
  if(!isTRUE(variance$relative_error <= 1e-6)){
    bounds <- vapply(signif(variance$bounds, 4), format, character(1))
    error  <- if(is.finite(variance$relative_error)){
      paste0(
        "its estimated relative rounding error is ",
        format(signif(variance$relative_error, 2), trim = TRUE)
      )
    }else{
      "its estimated rounding error exceeds the computed variance"
    }
    stop(
      "The variance of the truncated normal prior is unavailable in double-precision ",
      "arithmetic: ", error, " for a truncation from ", bounds[1], " to ", bounds[2],
      " prior standard deviations from the mean. Estimate the variance from ",
      "rng() draws of the prior, which sample the truncated distribution exactly, ",
      "or use a truncation closer to the prior mean or a wider truncation interval.",
      call. = FALSE
    )
  }
  sigma^2 * variance$value
}

.prior_normal_truncated_terms <- function(prior){

  # Closed-form moments of N(mu, sigma^2) truncated to [lower, upper] with
  # alpha, beta the standardized bounds and Z = Phi(beta) - Phi(alpha):
  # mean = mu + sigma * (phi(alpha) - phi(beta)) / Z and
  # var  = sigma^2 * (1 + (alpha phi(alpha) - beta phi(beta)) / Z
  #                   - ((phi(alpha) - phi(beta)) / Z)^2).
  # The ratios phi / Z use the log-space normalizing constant.
  mu     <- prior$parameters[["mean"]]
  sigma  <- prior$parameters[["sd"]]
  bounds <- (c(prior$truncation[["lower"]], prior$truncation[["upper"]]) - mu) / sigma
  log_Z  <- .prior_normal_log_C(prior)

  log_density  <- stats::dnorm(bounds, log = TRUE)
  ratio        <- ifelse(is.infinite(bounds), 0, exp(log_density - log_Z))
  scaled_ratio <- ifelse(is.infinite(bounds), 0, bounds * ratio)

  list(
    bounds       = bounds,
    log_density  = log_density,
    log_Z        = log_Z,
    ratio        = ratio,
    scaled_ratio = scaled_ratio,
    shift        = ratio[1] - ratio[2]
  )
}

# Standardized variance V = 1 + a r_a - b r_b - m^2 of a truncated normal, with
# r_a = phi(a) / Z (zero for an infinite bound) and m = r_a - r_b, and a
# conservative bound on its relative rounding error. Far in a tail the terms
# are of order a^2 while V is about 1 / a^2; for an interval of width w
# (standardized) they are of order max(1, |a|) / w while V is about w^2 / 12.
# Error model, with eps = .Machine$double.eps:
# - log phi(k) carries an absolute rounding error of at most
#   eps * (2 |log phi(k)| + 2); forming log r_k = log phi(k) - log Z and its
#   exponential adds at most eps * (|log r_k| + 2). Together these give the
#   bound-specific log error z_k.
# - Each log tail mass L behind log Z carries at most eps * (4 |L| + 4), and
#   log Z = L_1 + log1p(-exp(L_2 - L_1)) multiplies the error of L_1 by 1 + q
#   and that of L_2 by q, where q = exp(L_2) / Z. With the final rounding this
#   gives the log error e_Z shared by both ratios.
# The ratios are exact up to r_k exp(eta_k), |eta_k| <= h_k = z_k + e_Z, with
# the same e_Z in both. With g_a = a r_a - 2 m r_a and g_b = 2 m r_b - b r_b
# the sensitivities of V to log r_a and log r_b, perturbing the ratios changes
# V by at most
#   sum_k z_k |g_k| + e_Z |g_a + g_b| + sum_k |g_k| (expm1(h_k) - h_k)
#     + (sum_k r_k expm1(h_k))^2;
# the shared e_Z enters only through g_a + g_b = V - 1 - m^2, which avoids the
# spurious 1 / width cancellation of narrow intervals. Evaluating V rounds by
# at most 2 eps (1 + |a r_a| + |b r_b| + m^2), and rounding the standardized
# bounds (relative error <= eps) changes V by at most
# 2 eps (|a| |dV/da| + |b| |dV/db|), with dV/da = r_a (V - (a - m)^2) and
# dV/db = r_b ((b - m)^2 - V). With delta the sum of these terms, the relative
# error of the computed V is at most delta / (V - delta).
.prior_normal_truncated_variance <- function(prior){

  eps   <- .Machine$double.eps
  terms <- .prior_normal_truncated_terms(prior)
  finite <- is.finite(terms$bounds)
  ratio  <- terms$ratio
  shift  <- terms$shift
  value  <- 1 + terms$scaled_ratio[1] - terms$scaled_ratio[2] - shift^2

  # log tail masses combined by .prior_normal_log_interval_mass()
  mu    <- prior$parameters[["mean"]]
  sigma <- prior$parameters[["sd"]]
  lower <- prior$truncation[["lower"]]
  upper <- prior$truncation[["upper"]]
  if(lower >= mu){
    log_masses <- stats::pnorm(c(lower, upper), mu, sigma, lower.tail = FALSE, log.p = TRUE)
  }else{
    log_masses <- stats::pnorm(c(upper, lower), mu, sigma, lower.tail = TRUE, log.p = TRUE)
  }
  log_mass_error <- ifelse(is.finite(log_masses), eps * (4 * abs(log_masses) + 4), 0)
  q <- exp(log_masses[2] - terms$log_Z)
  log_Z_error <- log_mass_error[1] * (1 + q) + log_mass_error[2] * q +
    eps * (abs(terms$log_Z) + 2 * (1 + q))

  own_error <- ifelse(
    finite,
    eps * (2 * abs(terms$log_density) + 2) +
      eps * (abs(terms$log_density - terms$log_Z) + 2),
    0
  )
  log_ratio_error <- ifelse(finite, own_error + log_Z_error, 0)
  gradient <- c(
    terms$scaled_ratio[1] - 2 * shift * ratio[1],
    2 * shift * ratio[2] - terms$scaled_ratio[2]
  )

  propagated <- sum(own_error * abs(gradient)) + log_Z_error * abs(sum(gradient)) +
    sum(abs(gradient) * (expm1(log_ratio_error) - log_ratio_error)) +
    sum(ratio * expm1(log_ratio_error))^2
  evaluation <- 2 * eps * (1 + sum(abs(terms$scaled_ratio)) + shift^2)
  mean_distance <- ifelse(finite, (terms$bounds - shift)^2, 0)
  standardization <- 2 * eps * sum(ifelse(
    finite,
    abs(terms$bounds) * ratio * abs(value - mean_distance),
    0
  ))
  delta <- propagated + evaluation + standardization

  relative_error <- if(is.finite(value) && is.finite(delta) && value > delta){
    delta / (value - delta)
  }else{
    Inf
  }

  list(
    value          = value,
    relative_error = relative_error,
    bounds         = terms$bounds
  )
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
