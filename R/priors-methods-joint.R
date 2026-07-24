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

  dots  <- list(...)
  if(!is.null(dots[["transform_factor_samples"]])){
    check_bool(dots[["transform_factor_samples"]], "transform_factor_samples")
    transform_factor_samples <- dots[["transform_factor_samples"]]
  }else{
    transform_factor_samples <- TRUE
  }
  if(!is.null(dots[["sample_components"]])){
    check_bool(dots[["sample_components"]], "sample_components")
    sample_components <- dots[["sample_components"]]
  }else{
    sample_components <- FALSE
  }

  if(is.prior.spike_and_slab(prior)){

    inclusion_prob <- rng(.get_spike_and_slab_inclusion(prior), n)
    if(!is.numeric(inclusion_prob) || !is.null(dim(inclusion_prob)) || length(inclusion_prob) != n ||
       anyNA(inclusion_prob) || any(!is.finite(inclusion_prob)) ||
       any(inclusion_prob < 0 | inclusion_prob > 1)){
      stop("'prior_inclusion' must generate scalar probabilities within 0 and 1.", call. = FALSE)
    }
    inclusion <- stats::rbinom(n, size = 1, prob = inclusion_prob)

    if(sample_components)
      return(inclusion)

    x         <- rng(.get_spike_and_slab_variable(prior), n) * inclusion
    attr(x, "inclusion") <- inclusion

  }else if(is.prior.mixture(prior)){

    component_probabilities <- attr(prior, "prior_weights")
    components              <- sample(seq_along(component_probabilities), size = n, replace = TRUE, prob = component_probabilities)

    if(sample_components)
      return(components)

    if(inherits(prior, "prior.bias_mixture")){

      branch_info <- lapply(prior, .selection_branch_info)
      spec        <- selection_backend_spec(prior)

      is_PET        <- sapply(prior, is.prior.PET)
      is_PEESE      <- sapply(prior, is.prior.PEESE)
      has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
      has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

      out_names <- character()
      if(any(has_selection)){
        omega_names <- paste0("omega[", seq_len(spec$step$n_bins), "]")
        out_names   <- c(out_names, omega_names)

        selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
        omega_mapping    <- weightfunctions_mapping(selection_priors, one_sided = TRUE)
        selection_index  <- which(has_selection)
      }
      if(any(has_phacking)){
        out_names <- c(out_names, "alpha", "pi_null")
      }
      if(any(is_PET)){
        out_names <- c(out_names, "PET")
      }
      if(any(is_PEESE)){
        out_names <- c(out_names, "PEESE")
      }

      x <- matrix(0, nrow = n, ncol = length(out_names))
      colnames(x) <- out_names
      if(any(has_selection)){
        x[, omega_names] <- 1
      }

      for(component in unique(components)){
        component_rows <- component == components
        n_component    <- sum(component_rows)

        if(has_selection[component]){
          selection_i <- match(component, selection_index)
          selection_samples <- rng(branch_info[[component]]$selection, n_component)
          x[component_rows, omega_names] <- selection_samples[, paste0("omega[", omega_mapping[[selection_i]], "]"), drop = FALSE]
        }

        if(has_phacking[component]){
          phacking_samples <- rng(branch_info[[component]]$phacking, n_component)
          x[component_rows, c("alpha", "pi_null")] <- phacking_samples[, c("alpha", "pi_null"), drop = FALSE]
        }

        if(is_PET[component]){
          x[component_rows, "PET"] <- rng(prior[[component]], n_component, transform_factor_samples = FALSE)
        }
        if(is_PEESE[component]){
          x[component_rows, "PEESE"] <- rng(prior[[component]], n_component, transform_factor_samples = FALSE)
        }
      }

    }else if(inherits(prior, "prior.factor_mixture")){

      prior_type <- .get_prior_factor_list_type(prior)

      if(transform_factor_samples){
        x <- matrix(NA, nrow = n, ncol = prior_type[["K"]] + 1)
      }else{
        x <- matrix(NA, nrow = n, ncol = prior_type[["K"]])
      }

      for(component in unique(components)){
        x[component == components,] <- rng(prior[[component]], sum(component == components), transform_factor_samples = transform_factor_samples)
      }

    }else if(inherits(prior, "prior.simple_mixture")){

      x <- rep(NA, n)
      for(component in unique(components)){
        x[component == components] <- rng(prior[[component]], sum(component == components))
      }

    }else{
      stop("unsupported prior mixture type")
    }


    attr(x, "components") <- components

  }else if(is.prior.ordered(prior)){

    quantity <- if(is.null(dots[["quantity"]])) "level" else dots[["quantity"]]
    x <- .prior_ordered_rng(
      prior = prior,
      n = n,
      transform_factor_samples = transform_factor_samples,
      quantity = quantity
    )

  }else if(is.prior.simple(prior)){

    x <- .prior_simple_rng(prior, n)

  }else if(transform_factor_samples && (is.prior.orthonormal(prior) | is.prior.meandif(prior))){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par1 <- rep(0, prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
      "mpoint"     = rmpoint(n, location = par1)
    )


    if(is.prior.orthonormal(prior)){
      x <- x %*% t(contr.orthonormal(1:(prior$parameters[["K"]] + 1)))
    }else if(is.prior.meandif(prior)){
      x <- x %*% t(contr.meandif(1:(prior$parameters[["K"]] + 1)))
    }

  }else if(is.prior.vector(prior)){

    if(prior[["distribution"]] != "mpoint")
      .check_vector_truncation_unsupported(prior$truncation)

    if(prior[["distribution"]] == "dirichlet"){
      return(.prior_dirichlet_rng(prior, n))
    }

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    # TODO: generalize this to priors with covariances
    par1 <- rep(par1, length = prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
      "mpoint"     = rmpoint(n, location = par1)
    )

  }else if(is.prior.weightfunction(prior)){

    x <- .weightfunction_rng(prior, n)

  }else if(is_prior_phacking(prior)){

    alpha <- rng(prior$alpha, n)
    x <- cbind(
      alpha   = alpha,
      pi_null = phack_pi_null(alpha, prior$form, prior$source, prior$destination, target = prior$target)
    )

  }else if(is_prior_bias(prior)){

    selection_backend_spec(prior)

    x <- NULL
    if(!is.null(prior$selection)){
      x <- cbind(x, rng(prior$selection, n))
    }
    if(!is.null(prior$phacking)){
      x <- cbind(x, rng(prior$phacking, n))
    }

  }

  return(x)
}
#' @rdname prior_functions
cdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(as.vector(q))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No cdfs are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No cdfs are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_cdf(prior, q)

  }else if(is.prior.vector(prior)){

    stop("No cdfs are implemented for vector priors.")

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal cdfs are implemented for prior weightfunctions.")

  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){

    .selection_prior_stop_unsupported_generic("cdf", prior)

  }

  return(p)
}
#' @rdname prior_functions
ccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(as.vector(q))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No ccdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No ccdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_ccdf(prior, q)

  }else if(is.prior.vector(prior)){

    stop("No cdfs are implemented for vector priors.")

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal ccdf functions are implemented for prior weightfunctions.")

  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){

    .selection_prior_stop_unsupported_generic("ccdf", prior)

  }

  return(p)
}
#' @rdname prior_functions
lpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(as.vector(x))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No lpdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No lpdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    log_lik <- .prior_simple_lpdf(prior, x)

  }else if(is.prior.vector(prior)){

    if(prior[["distribution"]] != "mpoint")
      .check_vector_truncation_unsupported(prior$truncation)

    if(prior[["distribution"]] == "dirichlet"){
      return(.prior_dirichlet_lpdf(prior, x))
    }

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    # TODO: generalize this to priors with covariances
    par1 <- rep(par1, length = prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::dmvnorm(x, mean = par1, sigma = par2, log = TRUE),
      "mt"         = mvtnorm::dmvt(x, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted", log = TRUE),
      "mpoint"     = dmpoint(x, location = par1, log = TRUE)
    )

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal lpdf are implemented for prior weightfunctions.")

  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){

    .selection_prior_stop_unsupported_generic("lpdf", prior)

  }

  return(log_lik)
}
#' @rdname prior_functions
pdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(as.vector(x))
  .check_prior(prior)

  if(is.prior.simple(prior)){
    lik <- .prior_simple_pdf(prior, x)
  }else{
    log_lik <- lpdf.prior(prior, x)
    lik     <- exp(log_lik)
  }

  return(lik)
}
#' @rdname prior_functions
quant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No quantile functions are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No quantile functions are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    q <- .prior_simple_quant(prior, p)

  }else if(is.prior.vector(prior) && !is.prior.factor(prior)){

    q <- mquant(prior, p)

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal quantile functions are implemented for prior weightfunctions.")

  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){

    .selection_prior_stop_unsupported_generic("quantile functions", prior)

  }

  return(q)
}

.prior_simple_base_d <- function(prior, x, log = FALSE){

  switch(
    prior[["distribution"]],
    "normal"    = stats::dnorm(x, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], log = log),
    "lognormal" = stats::dlnorm(x, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], log = log),
    "t"         = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], log = log),
    "gamma"     = stats::dgamma(x, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], log = log),
    "invgamma"  = .dinvgamma_prior(x, shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]], log = log),
    "beta"      = stats::dbeta(x, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], log = log),
    "bernoulli" = stats::dbinom(x, size = 1, prob = prior$parameters[["probability"]], log = log),
    "exp"       = stats::dexp(x, rate = prior$parameters[["rate"]], log = log),
    "uniform"   = stats::dunif(x, min = prior$parameters[["a"]], max = prior$parameters[["b"]], log = log),
    "moment"    = .dmoment_prior(x, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], log = log),
    "invmoment" = .dinvmoment_prior(x, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]], log = log),
    "point"     = dpoint(x, location = prior$parameters[["location"]], log = log)
  )
}

.prior_simple_base_p <- function(prior, q, lower.tail = TRUE){

  switch(
    prior[["distribution"]],
    "normal"    = stats::pnorm(q, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = lower.tail, log.p = FALSE),
    "lognormal" = stats::plnorm(q, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = lower.tail, log.p = FALSE),
    "t"         = extraDistr::plst(q, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "gamma"     = stats::pgamma(q, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "invgamma"  = .pinvgamma_prior(q, shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "beta"      = stats::pbeta(q, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = lower.tail, log.p = FALSE),
    "bernoulli" = stats::pbinom(q, size = 1, prob = prior$parameters[["probability"]], lower.tail = lower.tail, log.p = FALSE),
    "exp"       = stats::pexp(q, rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "uniform"   = stats::punif(q, min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = lower.tail, log.p = FALSE),
    "moment"    = .pmoment_prior(q, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], lower.tail = lower.tail),
    "invmoment" = .pinvmoment_prior(q, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]], lower.tail = lower.tail),
    "point"     = ppoint(q, location = prior$parameters[["location"]], lower.tail = lower.tail, log.p = FALSE)
  )
}

.prior_simple_base_q <- function(prior, p, lower.tail = TRUE){

  switch(
    prior[["distribution"]],
    "normal"    = stats::qnorm(p, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = lower.tail, log.p = FALSE),
    "lognormal" = stats::qlnorm(p, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = lower.tail, log.p = FALSE),
    "t"         = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "gamma"     = stats::qgamma(p, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "invgamma"  = .qinvgamma_prior(p, shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "beta"      = stats::qbeta(p, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = lower.tail, log.p = FALSE),
    "bernoulli" = stats::qbinom(p, size = 1, prob = prior$parameters[["probability"]], lower.tail = lower.tail, log.p = FALSE),
    "exp"       = stats::qexp(p, rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "uniform"   = stats::qunif(p, min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = lower.tail, log.p = FALSE),
    "moment"    = .qmoment_prior(p, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], lower.tail = lower.tail),
    "invmoment" = .qinvmoment_prior(p, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]], lower.tail = lower.tail),
    "point"     = qpoint(p, location = prior$parameters[["location"]], lower.tail = lower.tail, log.p = FALSE)
  )
}

.prior_simple_base_r <- function(prior, n){

  switch(
    prior[["distribution"]],
    "normal"    = stats::rnorm(n, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]]),
    "lognormal" = stats::rlnorm(n, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]]),
    "t"         = extraDistr::rlst(n, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]]),
    "gamma"     = stats::rgamma(n, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]]),
    "invgamma"  = .rinvgamma_prior(n, shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]]),
    "beta"      = stats::rbeta(n, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]]),
    "bernoulli" = stats::rbinom(n, size = 1, prob = prior$parameters[["probability"]]),
    "exp"       = stats::rexp(n, rate = prior$parameters[["rate"]]),
    "uniform"   = stats::runif(n, min = prior$parameters[["a"]], max = prior$parameters[["b"]]),
    "moment"    = .rmoment_prior(n, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]]),
    "invmoment" = .rinvmoment_prior(n, location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]]),
    "point"     = rpoint(n, location = prior$parameters[["location"]])
  )
}

.prior_simple_discrete_support <- function(prior){

  switch(
    prior[["distribution"]],
    "bernoulli" = c(0, 1)
  )
}

.prior_simple_discrete_prob <- function(prior){

  switch(
    prior[["distribution"]],
    "bernoulli" = c(1 - prior$parameters[["probability"]], prior$parameters[["probability"]])
  )
}

.prior_simple_truncated_discrete <- function(prior){

  support <- .prior_simple_discrete_support(prior)
  prob    <- .prior_simple_discrete_prob(prior)
  keep    <- support >= prior$truncation[["lower"]] & support <= prior$truncation[["upper"]]

  support <- support[keep]
  prob    <- prob[keep]

  if(length(support) == 0 || sum(prob) <= 0){
    stop("The truncated discrete prior has no probability mass in its support.", call. = FALSE)
  }

  list(
    support = support,
    prob    = prob / sum(prob)
  )
}

.prior_simple_use_survival_truncation <- function(prior){

  if(prior[["distribution"]] == "point" || is.prior.discrete(prior) ||
     is.infinite(prior$truncation[["lower"]])){
    return(FALSE)
  }

  C1 <- .prior_C1(prior)
  is.finite(C1) && C1 > .5
}

.prior_simple_survival_C <- function(prior){

  .prior_simple_base_p(prior, prior$truncation[["lower"]], lower.tail = FALSE) -
    .prior_simple_base_p(prior, prior$truncation[["upper"]], lower.tail = FALSE)
}

.prior_normal_logdiffexp <- function(log_x, log_y){

  if(length(log_x) == 0 || length(log_y) == 0){
    return(numeric(0))
  }

  n     <- max(length(log_x), length(log_y))
  log_x <- rep(log_x, length.out = n)
  log_y <- rep(log_y, length.out = n)

  log_difference <- rep(NA_real_, n)
  both_zero      <- is.infinite(log_x) & log_x < 0 &
    is.infinite(log_y) & log_y < 0
  log_difference[both_zero] <- -Inf

  known <- !is.na(log_x) & !is.na(log_y) & !both_zero
  if(any(known)){
    delta      <- pmin(log_y[known] - log_x[known], 0)
    correction <- numeric(length(delta))
    far_apart  <- delta < log(.5)

    correction[far_apart]  <- log1p(-exp(delta[far_apart]))
    correction[!far_apart] <- log(-expm1(delta[!far_apart]))
    log_difference[known]   <- log_x[known] + correction
  }

  log_difference
}

.prior_normal_logaddexp <- function(log_x, log_y){

  if(length(log_x) == 0 || length(log_y) == 0){
    return(numeric(0))
  }

  n     <- max(length(log_x), length(log_y))
  log_x <- rep(log_x, length.out = n)
  log_y <- rep(log_y, length.out = n)

  log_sum       <- rep(NA_real_, n)
  both_zero     <- is.infinite(log_x) & log_x < 0 &
    is.infinite(log_y) & log_y < 0
  log_sum[both_zero] <- -Inf

  known <- !is.na(log_x) & !is.na(log_y) & !both_zero
  if(any(known)){
    larger         <- pmax(log_x[known], log_y[known])
    smaller        <- pmin(log_x[known], log_y[known])
    log_sum[known] <- larger + log1p(exp(smaller - larger))
  }

  log_sum
}

.prior_normal_log_interval_mass <- function(prior, lower, upper){

  if(length(lower) == 0 || length(upper) == 0){
    return(numeric(0))
  }

  n     <- max(length(lower), length(upper))
  lower <- rep(lower, length.out = n)
  upper <- rep(upper, length.out = n)

  log_mass    <- rep(NA_real_, n)
  known       <- !is.na(lower) & !is.na(upper)
  use_survival <- known & lower >= prior$parameters[["mean"]]

  if(any(use_survival)){
    log_lower <- stats::pnorm(
      lower[use_survival],
      mean       = prior$parameters[["mean"]],
      sd         = prior$parameters[["sd"]],
      lower.tail = FALSE,
      log.p      = TRUE
    )
    log_upper <- stats::pnorm(
      upper[use_survival],
      mean       = prior$parameters[["mean"]],
      sd         = prior$parameters[["sd"]],
      lower.tail = FALSE,
      log.p      = TRUE
    )
    log_mass[use_survival] <- .prior_normal_logdiffexp(log_lower, log_upper)
  }

  use_cdf <- known & !use_survival
  if(any(use_cdf)){
    log_upper <- stats::pnorm(
      upper[use_cdf],
      mean       = prior$parameters[["mean"]],
      sd         = prior$parameters[["sd"]],
      lower.tail = TRUE,
      log.p      = TRUE
    )
    log_lower <- stats::pnorm(
      lower[use_cdf],
      mean       = prior$parameters[["mean"]],
      sd         = prior$parameters[["sd"]],
      lower.tail = TRUE,
      log.p      = TRUE
    )
    log_mass[use_cdf] <- .prior_normal_logdiffexp(log_upper, log_lower)
  }

  log_mass
}

.prior_normal_log_C <- function(prior){

  .prior_normal_log_interval_mass(
    prior,
    prior$truncation[["lower"]],
    prior$truncation[["upper"]]
  )
}

.prior_simple_cdf <- function(prior, q){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_p(prior, q, lower.tail = TRUE))
  }

  p        <- rep(NA_real_, length(q))
  q_known  <- !is.na(q)
  q_lower  <- q_known & q < prior$truncation[["lower"]]
  q_higher <- q_known & q > prior$truncation[["upper"]]
  q_inside <- q_known & !q_lower & !q_higher

  p[q_lower]  <- 0
  p[q_higher] <- 1

  if(any(q_inside)){
    if(prior[["distribution"]] == "normal"){
      p[q_inside] <- exp(
        .prior_normal_log_interval_mass(
          prior,
          prior$truncation[["lower"]],
          q[q_inside]
        ) -
          .prior_normal_log_C(prior)
      )
      p[q_inside] <- pmin(1, pmax(0, p[q_inside]))
      return(p)
    }

    p[q_inside] <- .prior_simple_base_p(prior, q[q_inside], lower.tail = TRUE)

    if(prior[["distribution"]] != "point"){
      if(.prior_simple_use_survival_truncation(prior)){
        S1          <- .prior_simple_base_p(prior, prior$truncation[["lower"]], lower.tail = FALSE)
        S_q         <- .prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE)
        p[q_inside] <- (S1 - S_q) / .prior_C(prior)
      }else{
        C1          <- .prior_C1(prior)
        p[q_inside] <- (p[q_inside] - C1) / .prior_C(prior)
      }
    }
  }

  return(p)
}

.prior_simple_ccdf <- function(prior, q){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_p(prior, q, lower.tail = FALSE))
  }

  p        <- rep(NA_real_, length(q))
  q_known  <- !is.na(q)
  q_lower  <- q_known & q < prior$truncation[["lower"]]
  q_higher <- q_known & q > prior$truncation[["upper"]]
  q_inside <- q_known & !q_lower & !q_higher

  p[q_lower]  <- 1
  p[q_higher] <- 0

  if(any(q_inside)){
    if(prior[["distribution"]] == "normal"){
      p[q_inside] <- exp(
        .prior_normal_log_interval_mass(
          prior,
          q[q_inside],
          prior$truncation[["upper"]]
        ) -
          .prior_normal_log_C(prior)
      )
      p[q_inside] <- pmin(1, pmax(0, p[q_inside]))
      return(p)
    }

    p[q_inside] <- .prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE)

    if(prior[["distribution"]] != "point"){
      if(.prior_simple_use_survival_truncation(prior)){
        S2          <- .prior_simple_base_p(prior, prior$truncation[["upper"]], lower.tail = FALSE)
        p[q_inside] <- (p[q_inside] - S2) / .prior_C(prior)
      }else{
        C2          <- .prior_C2(prior)
        p[q_inside] <- (p[q_inside] - (1 - C2)) / .prior_C(prior)
      }
    }
  }

  return(p)
}

.prior_simple_lpdf <- function(prior, x){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_d(prior, x, log = TRUE))
  }

  log_lik <- .prior_simple_base_d(prior, x, log = TRUE)
  log_lik[x < prior$truncation[["lower"]] | x > prior$truncation[["upper"]]] <- -Inf

  if(prior[["distribution"]] != "point"){
    if(prior[["distribution"]] == "normal"){
      log_lik <- log_lik - .prior_normal_log_C(prior)
    }else{
      log_lik <- log_lik - log(.prior_C(prior))
    }
  }

  return(log_lik)
}

.prior_simple_pdf <- function(prior, x){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_d(prior, x, log = FALSE))
  }

  if(prior[["distribution"]] == "normal"){
    return(exp(.prior_simple_lpdf(prior, x)))
  }

  lik <- .prior_simple_base_d(prior, x, log = FALSE)
  lik[x < prior$truncation[["lower"]] | x > prior$truncation[["upper"]]] <- 0

  if(prior[["distribution"]] != "point"){
    lik <- lik / .prior_C(prior)
  }

  return(lik)
}

.prior_simple_quant <- function(prior, p){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_q(prior, p))
  }

  if(is.prior.discrete(prior)){
    discrete <- .prior_simple_truncated_discrete(prior)
    cdf      <- cumsum(discrete[["prob"]])
    return(vapply(p, function(p_i) discrete[["support"]][which(cdf >= p_i)[1]], numeric(1)))
  }

  if(prior[["distribution"]] == "normal"){
    q            <- rep(NA_real_, length(p))
    p_known      <- !is.na(p)
    p_lower      <- p_known & p == 0
    p_upper      <- p_known & p == 1
    p_inside     <- p_known & p > 0 & p < 1
    lower        <- prior$truncation[["lower"]]
    upper        <- prior$truncation[["upper"]]
    log_C        <- .prior_normal_log_C(prior)

    q[p_lower] <- lower
    q[p_upper] <- upper

    if(any(p_inside)){
      if(lower >= prior$parameters[["mean"]]){
        log_S2 <- stats::pnorm(
          upper,
          mean       = prior$parameters[["mean"]],
          sd         = prior$parameters[["sd"]],
          lower.tail = FALSE,
          log.p      = TRUE
        )
        log_target <- .prior_normal_logaddexp(
          log_S2,
          log1p(-p[p_inside]) + log_C
        )
        q[p_inside] <- stats::qnorm(
          p          = pmin(log_target, 0),
          mean       = prior$parameters[["mean"]],
          sd         = prior$parameters[["sd"]],
          lower.tail = FALSE,
          log.p      = TRUE
        )
      }else{
        log_F1 <- stats::pnorm(
          lower,
          mean       = prior$parameters[["mean"]],
          sd         = prior$parameters[["sd"]],
          lower.tail = TRUE,
          log.p      = TRUE
        )
        log_target <- .prior_normal_logaddexp(
          log_F1,
          log(p[p_inside]) + log_C
        )
        q[p_inside] <- stats::qnorm(
          p          = pmin(log_target, 0),
          mean       = prior$parameters[["mean"]],
          sd         = prior$parameters[["sd"]],
          lower.tail = TRUE,
          log.p      = TRUE
        )
      }
    }

    return(pmax(lower, pmin(upper, q)))
  }

  C1 <- .prior_C1(prior)
  if(.prior_simple_use_survival_truncation(prior)){
    S1 <- .prior_simple_base_p(prior, prior$truncation[["lower"]], lower.tail = FALSE)
    S2 <- .prior_simple_base_p(prior, prior$truncation[["upper"]], lower.tail = FALSE)
    return(.prior_simple_base_q(prior, S1 - p * (S1 - S2), lower.tail = FALSE))
  }

  .prior_simple_base_q(prior, C1 + p * .prior_C(prior))
}

.prior_simple_quant_optim <- function(prior, p){

  if(!is.infinite(prior$truncation[["lower"]]) & !is.infinite(prior$truncation[["upper"]])){
    start_value <- prior$truncation[["lower"]] + (prior$truncation[["upper"]]  - prior$truncation[["lower"]]) / 2
  }else if(!is.infinite(prior$truncation[["upper"]])){
    start_value <- prior$truncation[["upper"]] - 1
  }else if(!is.infinite(prior$truncation[["lower"]])){
    start_value <- prior$truncation[["lower"]] + 1
  }else{
    start_value <- 0
  }

  sapply(p, function(p_i){
    stats::optim(
      par     = start_value,
      fn      = function(x, prior, p_i)(.prior_simple_cdf(prior, x) - p_i)^2,
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

.prior_simple_rng <- function(prior, n){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_r(prior, n))
  }

  if(is.prior.discrete(prior)){
    discrete <- .prior_simple_truncated_discrete(prior)
    return(sample(discrete[["support"]], size = n, replace = TRUE, prob = discrete[["prob"]]))
  }

  if(.prior_simple_use_survival_truncation(prior)){
    return(.prior_simple_quant(prior, stats::runif(n)))
  }

  C1 <- .prior_C1(prior)
  .prior_simple_base_q(prior, stats::runif(n, min = C1, max = C1 + .prior_C(prior)))
}

.prior_simple_rng_rejection <- function(prior, n){

  x  <- NULL
  nn <- round(n * 1 / .prior_C(prior) * 1.10)

  while(length(x) < n){
    temp_x <- .prior_simple_base_r(prior, nn)
    x      <- c(x, temp_x[temp_x >= prior$truncation[["lower"]] & temp_x <= prior$truncation[["upper"]]])
  }

  x[1:n]
}

.prior_C1 <- function(prior){

  if(is.prior.simple(prior)){

    C1 <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(prior$truncation[["lower"]], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(prior$truncation[["lower"]], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(prior$truncation[["lower"]], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(prior$truncation[["lower"]], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = .pinvgamma_prior(prior$truncation[["lower"]], shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["lower"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "bernoulli" = stats::pbinom(ceiling(prior$truncation[["lower"]]) - 1, size = 1, prob = prior$parameters[["probability"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["lower"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["lower"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "moment"    = .pmoment_prior(prior$truncation[["lower"]], location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], lower.tail = TRUE),
      "invmoment" = .pinvmoment_prior(prior$truncation[["lower"]], location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]], lower.tail = TRUE),
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
      "invgamma"  = .pinvgamma_prior(prior$truncation[["upper"]], shape = prior$parameters[["shape"]], scale = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["upper"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "bernoulli" = stats::pbinom(prior$truncation[["upper"]], size = 1, prob = prior$parameters[["probability"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["upper"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["upper"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "moment"    = .pmoment_prior(prior$truncation[["upper"]], location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], lower.tail = TRUE),
      "invmoment" = .pinvmoment_prior(prior$truncation[["upper"]], location = prior$parameters[["location"]], tau = prior$parameters[["tau"]], order = prior$parameters[["order"]], df = prior$parameters[["df"]], lower.tail = TRUE),
      "point"     = ppoint(prior$truncation[["upper"]], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

  }

  return(C2)
}
.prior_C  <- function(prior){

  if(is.prior.simple(prior)){

    C <- .prior_C2(prior) - .prior_C1(prior)
    if(prior[["distribution"]] != "point" && !is.prior.discrete(prior) &&
       (.prior_simple_use_survival_truncation(prior) || !is.finite(C) || C <= 0)){
      C <- .prior_simple_survival_C(prior)
    }

  }

  return(C)
}

.prior_vector_dimension <- function(prior){

  if(identical(prior$distribution, "dirichlet")){
    return(length(prior$parameters[["alpha"]]))
  }

  K <- prior$parameters[["K"]]
  .check_parameter_dimensions(K, "K", allow_NA = FALSE)
  as.integer(K)
}
.prior_vector_primary_location <- function(prior){

  location <- switch(
    prior[["distribution"]],
    "mnormal" = prior$parameters[["mean"]],
    "mt"      = prior$parameters[["location"]],
    "mpoint"  = prior$parameters[["location"]]
  )

  if(length(location) != 1){
    stop("Only exchangeable vector prior distributions with scalar location parameters are supported.", call. = FALSE)
  }

  location
}
.prior_vector_x_matrix <- function(x, K, name = "x"){

  if(is.matrix(x)){
    if(ncol(x) != K){
      stop(paste0("The '", name, "' argument must have ", K, " columns."), call. = FALSE)
    }
    return(x)
  }

  matrix(rep(x, times = K), ncol = K)
}
.prior_dirichlet_rng <- function(prior, n){

  alpha <- prior$parameters[["alpha"]]
  out <- extraDistr::rdirichlet(n, alpha = alpha)
  colnames(out) <- paste0("V", seq_along(alpha))
  out
}
.prior_dirichlet_x_matrix <- function(x, K, name = "x"){

  if(is.matrix(x)){
    if(ncol(x) != K){
      stop(paste0("The '", name, "' argument must have ", K, " columns."), call. = FALSE)
    }
    return(x)
  }
  if(is.data.frame(x)){
    x <- as.matrix(x)
    if(ncol(x) != K){
      stop(paste0("The '", name, "' argument must have ", K, " columns."), call. = FALSE)
    }
    return(x)
  }
  if(length(x) == K){
    return(matrix(x, nrow = 1L))
  }
  if(length(x) %% K == 0L){
    return(matrix(x, ncol = K, byrow = TRUE))
  }

  stop(paste0("The '", name, "' argument length must be a multiple of the Dirichlet dimension."), call. = FALSE)
}
.prior_dirichlet_lpdf <- function(prior, x){

  alpha <- prior$parameters[["alpha"]]
  K <- length(alpha)
  x_mat <- .prior_dirichlet_x_matrix(x, K, "x")
  log_const <- lgamma(sum(alpha)) - sum(lgamma(alpha))

  apply(x_mat, 1L, function(row){
    if(anyNA(row)){
      return(NA_real_)
    }
    if(any(row < 0 | row > 1) || !isTRUE(all.equal(sum(row), 1, tolerance = 1e-8))){
      return(-Inf)
    }
    zero <- row == 0
    if(any(zero & alpha < 1)){
      return(Inf)
    }
    if(any(zero & alpha > 1)){
      return(-Inf)
    }
    positive <- row > 0
    log_const + sum((alpha[positive] - 1) * log(row[positive]))
  })
}
.prior_vector_marginal_cdf <- function(prior, q, lower.tail = TRUE){

  K     <- .prior_vector_dimension(prior)
  q_mat <- .prior_vector_x_matrix(q, K, "q")

  if(identical(prior$distribution, "dirichlet")){
    alpha <- prior$parameters[["alpha"]]
    alpha0 <- sum(alpha)
    out <- vapply(seq_len(K), function(i){
      stats::pbeta(q_mat[, i], shape1 = alpha[i], shape2 = alpha0 - alpha[i], lower.tail = lower.tail)
    }, numeric(nrow(q_mat)))
    return(matrix(out, nrow = nrow(q_mat), ncol = K))
  }

  loc   <- .prior_vector_primary_location(prior)

  out <- switch(
    prior[["distribution"]],
    "mnormal" = apply(q_mat, 2, stats::pnorm, mean = loc, sd = prior$parameters[["sd"]], lower.tail = lower.tail),
    "mt"      = apply(q_mat, 2, extraDistr::plst, df = prior$parameters[["df"]], mu = loc, sigma = prior$parameters[["scale"]], lower.tail = lower.tail),
    "mpoint"  = apply(q_mat, 2, ppoint, location = loc, lower.tail = lower.tail)
  )

  matrix(out, nrow = nrow(q_mat), ncol = K)
}
.prior_vector_marginal_lpdf <- function(prior, x){

  K     <- .prior_vector_dimension(prior)
  x_mat <- .prior_vector_x_matrix(x, K, "x")

  if(identical(prior$distribution, "dirichlet")){
    alpha <- prior$parameters[["alpha"]]
    alpha0 <- sum(alpha)
    out <- vapply(seq_len(K), function(i){
      stats::dbeta(x_mat[, i], shape1 = alpha[i], shape2 = alpha0 - alpha[i], log = TRUE)
    }, numeric(nrow(x_mat)))
    return(matrix(out, nrow = nrow(x_mat), ncol = K))
  }

  loc   <- .prior_vector_primary_location(prior)

  out <- switch(
    prior[["distribution"]],
    "mnormal" = apply(x_mat, 2, stats::dnorm, mean = loc, sd = prior$parameters[["sd"]], log = TRUE),
    "mt"      = apply(x_mat, 2, extraDistr::dlst, df = prior$parameters[["df"]], mu = loc, sigma = prior$parameters[["scale"]], log = TRUE),
    "mpoint"  = apply(x_mat, 2, dpoint, location = loc, log = TRUE)
  )

  matrix(out, nrow = nrow(x_mat), ncol = K)
}
.prior_vector_marginal_quant <- function(prior, p){

  K   <- .prior_vector_dimension(prior)

  if(identical(prior$distribution, "dirichlet")){
    alpha <- prior$parameters[["alpha"]]
    alpha0 <- sum(alpha)
    q <- vapply(seq_len(K), function(i){
      stats::qbeta(p, shape1 = alpha[i], shape2 = alpha0 - alpha[i])
    }, numeric(length(p)))
    return(matrix(q, nrow = length(p), ncol = K))
  }

  loc <- .prior_vector_primary_location(prior)

  q <- switch(
    prior[["distribution"]],
    "mnormal" = stats::qnorm(p, mean = loc, sd = prior$parameters[["sd"]]),
    "mt"      = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = loc, sigma = prior$parameters[["scale"]]),
    "mpoint"  = qpoint(p, location = loc)
  )

  matrix(rep(q, times = K), ncol = K)
}
.prior_vector_mean <- function(prior){

  K   <- .prior_vector_dimension(prior)

  if(identical(prior$distribution, "dirichlet")){
    alpha <- prior$parameters[["alpha"]]
    return(alpha / sum(alpha))
  }

  loc <- .prior_vector_primary_location(prior)

  m <- switch(
    prior[["distribution"]],
    "mnormal" = loc,
    "mt"      = ifelse(prior$parameters[["df"]] > 1, loc, NaN),
    "mpoint"  = loc
  )

  rep(m, K)
}
.prior_vector_var <- function(prior){

  K <- .prior_vector_dimension(prior)

  if(identical(prior$distribution, "dirichlet")){
    alpha <- prior$parameters[["alpha"]]
    alpha0 <- sum(alpha)
    return(alpha * (alpha0 - alpha) / (alpha0^2 * (alpha0 + 1)))
  }

  v <- switch(
    prior[["distribution"]],
    "mnormal" = prior$parameters[["sd"]]^2,
    "mt"      = ifelse(prior$parameters[["df"]] > 2, prior$parameters[["scale"]]^2 * prior$parameters[["df"]] / (prior$parameters[["df"]] - 2), NaN),
    "mpoint"  = 0
  )

  rep(v, K)
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
    "moment"    = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "invmoment" = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "beta"      = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & isTRUE(all.equal(prior$truncation[["upper"]], 1)),
    "bernoulli" = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & isTRUE(all.equal(prior$truncation[["upper"]], 1)),
    "exp"       = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "uniform"   = isTRUE(all.equal(prior$truncation[["lower"]], prior$parameters[["a"]])) &
      isTRUE(all.equal(prior$truncation[["upper"]], prior$parameters[["b"]])),
    "point"     = TRUE,
    "mpoint"          = TRUE,
    "dirichlet"       = TRUE,
    "weightfunction"  = TRUE,
    "none"            = TRUE
  )
}
