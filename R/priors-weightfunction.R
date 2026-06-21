#' @title Creates a prior distribution for a weight function
#'
#' @description \code{prior_weightfunction} creates a prior distribution for
#' fitting a RoBMA selection model. The \code{side} and \code{steps} arguments
#' define the p-value bins, and the \code{weights} argument defines the prior on
#' the publication weights in those bins.
#'
#' @param side side geometry. Either \code{"one-sided"} or \code{"two-sided"}.
#' @param steps increasing p-value cut points between 0 and 1.
#' @param weights a weight-prior object created by \code{wf_cumulative()},
#' \code{wf_fixed()}, or \code{wf_independent()}.
#' @param reference reference bin. Currently only \code{"most_significant"} is
#' supported and fixes the most significant bin to \code{omega = 1}.
#' @param prior_weights prior odds associated with a given distribution.
#'
#' @examples
#' p1 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05, .10),
#'   weights = wf_cumulative(alpha = c(1, 1, 1))
#' )
#'
#' p2 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05),
#'   weights = wf_independent(prior("beta", list(1, 1)))
#' )
#'
#' @return \code{prior_weightfunction} returns an object of class 'prior'.
#'
#' @export prior_weightfunction
#' @seealso [plot.prior()]
prior_weightfunction <- function(side = "one-sided", steps = c(.025, .05),
                                 weights = wf_cumulative(),
                                 reference = "most_significant",
                                 prior_weights = 1){

  check_char(side, "side")
  check_real(steps, "steps", check_length = 0)
  check_char(reference, "reference", allow_values = "most_significant")
  .check_prior_weight(prior_weights)

  side <- .weightfunction_normalize_side(side)
  steps <- .weightfunction_validate_steps(steps)
  bins <- .weightfunction_bins(steps)
  weights <- .weightfunction_validate_weights(weights, n_bins = nrow(bins), reference = reference)
  truncation <- .weightfunction_weights_truncation(weights)

  output <- list(
    distribution  = "weightfunction",
    side          = side,
    steps         = steps,
    bins          = bins,
    reference     = reference,
    weights       = weights,
    parameters    = list(steps = steps),
    truncation    = truncation,
    prior_weights = prior_weights
  )

  class(output) <- c("prior", "prior.weightfunction")

  return(output)
}

#' @rdname prior_weightfunction
#' @param alpha positive cumulative-Dirichlet concentration parameters. If
#' omitted, a flat Dirichlet prior is used with one concentration parameter per
#' bin.
#' @export
wf_cumulative <- function(alpha = NULL){

  if(!is.null(alpha)){
    check_real(alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = 0, allow_NA = FALSE)
    if(any(!is.finite(alpha))){
      stop("The 'alpha' argument must be finite.", call. = FALSE)
    }
  }

  out <- list(type = "cumulative", alpha = alpha)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.cumulative")
  return(out)
}

#' @rdname prior_weightfunction
#' @param omega fixed non-negative relative publication weights, one per bin.
#' The reference-bin weight must be exactly 1.
#' @export
wf_fixed <- function(omega){

  check_real(omega, "omega", lower = 0, check_length = 0, allow_NA = FALSE)
  if(any(!is.finite(omega))){
    stop("The 'omega' argument must be finite.", call. = FALSE)
  }

  out <- list(type = "fixed", omega = omega)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.fixed")
  return(out)
}

#' @rdname prior_weightfunction
#' @param prior prior distribution for each non-reference weight.
#' @param scale latent scale for independent weights. \code{"omega"} places the
#' prior directly on the non-negative publication weight. \code{"log_omega"}
#' places the prior on \code{log(omega)} and transforms with
#' \code{omega = exp(log_omega)}, allowing weights above one whenever the log
#' prior assigns mass above zero.
#' @export
wf_independent <- function(prior, scale = "omega"){

  .check_prior(prior)
  if(!is.prior.simple(prior) || is.prior.point(prior) || is.prior.discrete(prior)){
    stop("'prior' must be a continuous simple prior distribution.", call. = FALSE)
  }

  check_char(scale, "scale", allow_values = c("omega", "log_omega", "log"))
  if(scale == "log"){
    scale <- "log_omega"
  }

  if(scale == "omega"){
    if(prior$truncation[["lower"]] < 0){
      stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
    }
  }

  out <- list(type = "independent", scale = scale, prior = prior)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.independent")
  return(out)
}

.weightfunction_normalize_side <- function(side){

  side_clean <- .prior_clean_input_name(side)
  if(side_clean %in% c("onesided", "one")){
    return("one-sided")
  }
  if(side_clean %in% c("twosided", "two")){
    return("two-sided")
  }

  stop("'side' must be either 'one-sided' or 'two-sided'.", call. = FALSE)
}
.weightfunction_validate_steps <- function(steps){

  check_real(steps, "steps", check_length = 0, allow_NA = FALSE)

  if(length(steps) == 0){
    stop("'steps' must contain at least one p-value cut point.", call. = FALSE)
  }
  if(any(steps >= 1) || any(steps <= 0)){
    stop("'steps' must be higher than 0 and lower than 1.", call. = FALSE)
  }
  if(anyDuplicated(steps)){
    stop("'steps' must not contain duplicate cut points.", call. = FALSE)
  }
  if(!all(steps == cummax(steps))){
    stop("'steps' must be monotonically increasing.", call. = FALSE)
  }

  steps
}
.weightfunction_bins <- function(steps){

  cuts <- c(0, steps, 1)
  data.frame(
    lower     = cuts[-length(cuts)],
    upper     = cuts[-1],
    reference = c(TRUE, rep(FALSE, length(cuts) - 2L))
  )
}
.weightfunction_validate_weights <- function(weights, n_bins, reference){

  if(!inherits(weights, "weightfunction_weights")){
    stop("'weights' must be created by wf_cumulative(), wf_fixed(), or wf_independent().", call. = FALSE)
  }

  if(weights$type == "cumulative"){
    if(is.null(weights$alpha)){
      weights$alpha <- rep(1, n_bins)
    }
    check_real(weights$alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = n_bins, allow_NA = FALSE)
    if(any(!is.finite(weights$alpha))){
      stop("The 'alpha' argument must be finite.", call. = FALSE)
    }

  }else if(weights$type == "fixed"){
    check_real(weights$omega, "omega", lower = 0, check_length = n_bins, allow_NA = FALSE)
    if(any(!is.finite(weights$omega))){
      stop("The 'omega' argument must be finite.", call. = FALSE)
    }
    if(reference == "most_significant" && !isTRUE(all.equal(weights$omega[1], 1))){
      stop("The reference-bin fixed weight must be exactly 1.", call. = FALSE)
    }

  }else if(weights$type == "independent"){
    .check_prior(weights$prior)
    if(!is.prior.simple(weights$prior) || is.prior.point(weights$prior) || is.prior.discrete(weights$prior)){
      stop("'weights$prior' must be a continuous simple prior distribution.", call. = FALSE)
    }
    if(weights$scale == "omega"){
      if(weights$prior$truncation[["lower"]] < 0){
        stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
      }
    }else if(weights$scale != "log_omega"){
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

  }else{
    stop("Unsupported weightfunction weight prior type.", call. = FALSE)
  }

  weights
}
.weightfunction_weights_truncation <- function(weights){

  if(weights$type == "cumulative"){
    return(list(lower = 0, upper = 1))
  }

  if(weights$type == "fixed"){
    return(list(lower = 0, upper = max(1, weights$omega, na.rm = TRUE)))
  }

  if(weights$type == "independent"){
    if(weights$scale == "omega"){
      upper <- weights$prior$truncation[["upper"]]
    }else if(weights$scale == "log_omega"){
      upper <- weights$prior$truncation[["upper"]]
      upper <- if(is.infinite(upper)) Inf else exp(upper)
    }else{
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

    return(list(lower = 0, upper = max(1, upper)))
  }

  stop("Unsupported weightfunction weight prior type.", call. = FALSE)
}
.weightfunction_n_bins <- function(prior){

  if(!is.prior.weightfunction(prior)){
    stop("'prior' must be a weightfunction prior.", call. = FALSE)
  }

  nrow(prior$bins)
}
.weightfunction_is_fixed <- function(prior){
  prior$weights$type == "fixed"
}
.weightfunction_local_cuts <- function(prior){

  c(prior$bins$lower[1], prior$bins$upper)
}
.weightfunction_reference_index <- function(prior){
  which(prior$bins$reference)
}
.weightfunction_alpha_marginal <- function(alpha, index){

  if(index <= 1L){
    return(list(type = "point", location = 1))
  }

  list(
    type  = "beta",
    alpha = sum(alpha[index:length(alpha)]),
    beta  = sum(alpha[seq_len(index - 1L)])
  )
}
.weightfunction_rng <- function(prior, n){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    out <- matrix(rep(prior$weights$omega, each = n), nrow = n)

  }else if(prior$weights$type == "cumulative"){
    theta <- extraDistr::rdirichlet(n, alpha = prior$weights$alpha)
    out   <- t(apply(theta[,J:1, drop = FALSE], 1, cumsum))[,J:1, drop = FALSE]
    # Keep the reference bin exact; row-wise cumulative sums can leave 1 +/- eps.
    out[,1L] <- 1

  }else if(prior$weights$type == "independent"){
    out <- matrix(1, nrow = n, ncol = J)
    if(J > 1L){
      draws <- matrix(rng(prior$weights$prior, n * (J - 1L)), nrow = n, ncol = J - 1L)
      if(prior$weights$scale == "log_omega"){
        draws <- exp(draws)
      }
      out[,2:J] <- draws
    }
  }

  colnames(out) <- paste0("omega[", seq_len(J), "]")
  out
}
.weightfunction_marginal_components <- function(prior){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    return(lapply(prior$weights$omega, function(x){
      list(type = "point", location = x)
    }))
  }

  if(prior$weights$type == "cumulative"){
    return(lapply(seq_len(J), function(j){
      .weightfunction_alpha_marginal(prior$weights$alpha, j)
    }))
  }

  if(prior$weights$type == "independent"){
    return(lapply(seq_len(J), function(j){
      if(j == 1L){
        list(type = "point", location = 1)
      }else{
        list(type = "prior", prior = prior$weights$prior, scale = prior$weights$scale)
      }
    }))
  }
}
.prior_weightfunction_component_range <- function(component, quantiles = .005){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      if(component$scale == "omega"){
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) mquant(component$prior, quantiles) else lower
        upper <- if(is.infinite(upper)) mquant(component$prior, 1 - quantiles) else upper
        c(lower, upper)
      }else{
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) 0 else exp(lower)
        upper <- if(is.infinite(upper)) exp(mquant(component$prior, 1 - quantiles)) else exp(upper)
        c(lower, upper)
      }
    }
  )
}
.weightfunction_range <- function(prior, quantiles = .005){

  ranges <- do.call(rbind, lapply(
    .weightfunction_marginal_components(prior),
    .prior_weightfunction_component_range,
    quantiles = quantiles
  ))

  x_range <- range(c(0, 1, as.vector(ranges)), finite = TRUE)
  if(x_range[1] == x_range[2]){
    x_range <- range(c(0, 1, x_range), finite = TRUE)
  }

  x_range
}
.prior_weightfunction_component_cdf <- function(component, q){

  switch(
    component$type,
    "point" = ppoint(q, location = component$location),
    "beta"  = stats::pbeta(q, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mcdf(component$prior, q)
      }else{
        p <- rep(NA_real_, length(q))
        q_known <- !is.na(q)
        p[q_known & q <= 0] <- 0
        inside <- q_known & q > 0
        if(any(inside)){
          p[inside] <- mcdf(component$prior, log(q[inside]))
        }
        p
      }
    }
  )
}
.prior_weightfunction_component_ccdf <- function(component, q){
  if(component$type == "point"){
    return(ppoint(q, location = component$location, lower.tail = FALSE))
  }
  1 - .prior_weightfunction_component_cdf(component, q)
}
.prior_weightfunction_component_lpdf <- function(component, x){

  switch(
    component$type,
    "point" = dpoint(x, location = component$location, log = TRUE),
    "beta"  = stats::dbeta(x, shape1 = component$alpha, shape2 = component$beta, log = TRUE),
    "prior" = {
      if(component$scale == "omega"){
        mlpdf(component$prior, x)
      }else{
        out <- rep(-Inf, length(x))
        inside <- x > 0
        out[inside] <- mlpdf(component$prior, log(x[inside])) - log(x[inside])
        out
      }
    }
  )
}
.prior_weightfunction_component_quant <- function(component, p){

  switch(
    component$type,
    "point" = ifelse(is.na(p), NA_real_, component$location),
    "beta"  = stats::qbeta(p, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mquant(component$prior, p)
      }else{
        exp(mquant(component$prior, p))
      }
    }
  )
}
.prior_weightfunction_component_mean <- function(component){

  switch(
    component$type,
    "point" = component$location,
    "beta"  = component$alpha / (component$alpha + component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mean(component$prior)
      }else{
        stats::integrate(
          f     = function(x, prior) {
            y <- exp(x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
      }
    }
  )
}
.prior_weightfunction_component_var <- function(component){

  switch(
    component$type,
    "point" = 0,
    "beta"  = (component$alpha * component$beta) /
      ((component$alpha + component$beta)^2 * (component$alpha + component$beta + 1)),
    "prior" = {
      if(component$scale == "omega"){
        var(component$prior)
      }else{
        m1 <- .prior_weightfunction_component_mean(component)
        m2 <- stats::integrate(
          f     = function(x, prior) {
            y <- exp(2 * x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
        m2 - m1^2
      }
    }
  )
}
