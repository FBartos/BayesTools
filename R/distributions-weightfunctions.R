#' @title Weight functions
#'
#' @description Marginal density, marginal distribution function,
#' marginal quantile function and random generation for weight functions.
#'
#' @param x,q numeric vector of quantiles; the parameter arguments (\code{alpha},
#' \code{omega}) may be matrices with one row per quantile.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param alpha vector or matrix with concentration parameters for
#' the Dirichlet distribution for a monotonic one.sided or a two.sided
#' weight function.
#' @param alpha1 vector or matrix with concentration parameters for
#' the Dirichlet distribution for the expected direction of non-monotonic
#' one.sided of weight function.
#' @param alpha2 vector or matrix with concentration parameters for
#' the Dirichlet distribution for the unexpected direction of non-monotonic
#' one.sided of weight function.
#' @param omega vector or matrix of fixed non-negative relative weights for a
#' one.sided or a two.sided weight function. The first/reference weight must be
#' exactly 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @details The marginal density, distribution function, and quantile function
#' are implemented for the monotonic parameterization (\code{alpha}) and for
#' fixed weights (\code{omega}); the non-monotonic parameterization
#' (\code{alpha1}, \code{alpha2}) supports random generation only. The same
#' marginals of a [prior_weightfunction()] prior are available through
#' [mpdf()], [mcdf()], [mquant()], and [rng()].
#'
#' The marginals of the monotonic parameterization are Beta distributions with
#' the cumulative sums of \code{alpha} as shapes. Finite \code{alpha} whose
#' cumulative sums exceed the largest double signal an error of class
#' \code{BayesTools_weightfunction_beta_shape_unavailable} (parent
#' \code{BayesTools_weightfunction_unavailable}) from the density,
#' distribution, and quantile functions. Degenerate Dirichlet draws of the
#' random generation signal \code{BayesTools_dirichlet_rng_unavailable}.
#'
#' @examples
#' # draw samples from a two-sided weight function
#' rtwo.sided(10, alpha = c(1, 1))
#'
#' # draw samples from a monotone one-sided weight function
#' rone.sided(10, alpha = c(1, 1, 1))
#'
#' # draw samples from a non-monotone one-sided weight function
#' rone.sided(10, alpha1 = c(1, 1), alpha2 = c(1, 1))
#'
#' @return \code{mdone.sided}, \code{mdtwo.sided}, \code{mdone.sided_fixed},
#' and \code{mdtwo.sided_fixed} give the marginal density,
#' \code{mpone.sided}, \code{mptwo.sided}, \code{mpone.sided_fixed},
#' and \code{mptwo.sided_fixed} give the marginal distribution function,
#' \code{mqone.sided}, \code{mqtwo.sided}, \code{mqone.sided_fixed},
#' and \code{mqtwo.sided_fixed} give the marginal quantile function,
#' and \code{rone.sided}, \code{rtwo.sided}, \code{rone.sided_fixed},
#' and \code{rtwo.sided_fixed} generate random deviates.
#'
#' @export mdone.sided
#' @export mdtwo.sided
#' @export mdone.sided_fixed
#' @export mdtwo.sided_fixed
#' @export rone.sided
#' @export rtwo.sided
#' @export rone.sided_fixed
#' @export rtwo.sided_fixed
#' @export mpone.sided
#' @export mptwo.sided
#' @export mpone.sided_fixed
#' @export mptwo.sided_fixed
#' @export mqone.sided
#' @export mqtwo.sided
#' @export mqone.sided_fixed
#' @export mqtwo.sided_fixed
#' @name weightfunctions
NULL

#### wrappers ####
#' @rdname weightfunctions
mdone.sided <- function(x, alpha = NULL, alpha1 = NULL, alpha2 = NULL, log = FALSE){

  parameterization <- .weightfunctions_one_sided_parameterization(alpha, alpha1, alpha2)

  # common input check
  .check_log(log)
  .check_x(x, lower = 0, upper = 1)

  if(parameterization == "monotonic"){
    lik <- .mdone.sided_monotonic(x, alpha, log)
  }else{
    lik <- .mdone.sided_general(x, alpha1, alpha2, log)
  }

  return(lik)
}
#' @rdname weightfunctions
mdtwo.sided <- function(x, alpha, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(x, lower = 0, upper = 1)

  lik <- mdone.sided(x, alpha = alpha, log = log)

  return(lik)
}
#' @rdname weightfunctions
mdone.sided_fixed <- function(x, omega, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(x, lower = 0)

  lik <- .mdone.sided_fixed(x, omega, log)

  return(lik)
}
#' @rdname weightfunctions
mdtwo.sided_fixed <- function(x, omega, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(x, lower = 0)

  lik <- mdone.sided_fixed(x, omega = omega, log = log)

  return(lik)
}

#' @rdname weightfunctions
rone.sided <- function(n, alpha = NULL, alpha1 = NULL, alpha2 = NULL){

  parameterization <- .weightfunctions_one_sided_parameterization(alpha, alpha1, alpha2)

  # common input check
  .check_n(n)

  if(parameterization == "monotonic"){
    x <- .rone.sided_monotonic(n, alpha)
  }else{
    x <- .rone.sided_general(n, alpha1, alpha2)
  }

  return(x)
}
#' @rdname weightfunctions
rtwo.sided <- function(n, alpha){

  # common input check
  .check_n(n)

  x <- rone.sided(n, alpha = alpha)

  return(x)
}
#' @rdname weightfunctions
rone.sided_fixed <- function(n, omega){

  # common input check
  .check_n(n)

  x <- .rone.sided_fixed(n, omega)

  return(x)
}
#' @rdname weightfunctions
rtwo.sided_fixed <- function(n, omega){

  # common input check
  .check_n(n)

  x <- rone.sided_fixed(n, omega = omega)

  return(x)
}

#' @rdname weightfunctions
mpone.sided <- function(q, alpha = NULL, alpha1 = NULL, alpha2 = NULL, lower.tail = TRUE, log.p = FALSE){

  parameterization <- .weightfunctions_one_sided_parameterization(alpha, alpha1, alpha2)

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0, upper = 1)

  if(parameterization == "monotonic"){
    p <- .mpone.sided_monotonic(q, alpha, lower.tail, log.p)
  }else{
    p <- .mpone.sided_general(q, alpha1, alpha2, lower.tail, log.p)
  }

  return(p)
}
#' @rdname weightfunctions
mptwo.sided <- function(q, alpha, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0, upper = 1)

  p <- mpone.sided(q, alpha = alpha, lower.tail = lower.tail, log.p = log.p)

  return(p)
}
#' @rdname weightfunctions
mpone.sided_fixed <- function(q, omega, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0)

  p <- .mpone.sided_fixed(q, omega, lower.tail, log.p)

  return(p)
}
#' @rdname weightfunctions
mptwo.sided_fixed <- function(q, omega, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0)

  p <- mpone.sided_fixed(q, omega = omega, lower.tail = lower.tail, log.p = log.p)

  return(p)
}

#' @rdname weightfunctions
mqone.sided <- function(p, alpha = NULL, alpha1 = NULL, alpha2 = NULL, lower.tail = TRUE, log.p = FALSE){

  parameterization <- .weightfunctions_one_sided_parameterization(alpha, alpha1, alpha2)

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  if(parameterization == "monotonic"){
    q <- .mqone.sided_monotonic(p, alpha, lower.tail, log.p)
  }else{
    q <- .mqone.sided_general(p, alpha1, alpha2, lower.tail, log.p)
  }

  return(q)
}
#' @rdname weightfunctions
mqtwo.sided <- function(p, alpha, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  q <- mqone.sided(p, alpha = alpha, lower.tail = lower.tail, log.p = log.p)

  return(q)
}
#' @rdname weightfunctions
mqone.sided_fixed <- function(p, omega, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  q <- .mqone.sided_fixed(p, omega, lower.tail, log.p)

  return(q)
}
#' @rdname weightfunctions
mqtwo.sided_fixed <- function(p, omega, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  q <- mqone.sided_fixed(p, omega = omega, lower.tail = lower.tail, log.p = log.p)

  return(q)
}


###### helper functions
#### density functions ####
.mdone.sided_general   <- function(x, alpha1, alpha2, log){

  # would require the product of two beta-distributed variables, since the weights in the opposite direction
  # start at the first cutoff modeled by the first alpha parameter
  stop("Not implemented")
}
.mdone.sided_monotonic <- function(x, alpha, log){

  # input check
  .weightfunctions_check_alpha(alpha, "alpha")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  if(length(x) == 1){
    x <- rep(x, nrow(alpha))
  }

  if(nrow(alpha) != length(x) & nrow(alpha) != 1)
    stop("Non matching dimensions of 'alpha' and 'x'.")

  if(nrow(alpha) != length(x) & nrow(alpha) == 1){
    alpha <- do.call(rbind, lapply(1:length(x), function(i)alpha))
  }

  # marginals of reversed cumulative sums of Dirichlet variables are beta
  # distributed; the first/reference publication weight is fixed to one.
  beta_shapes <- .weightfunctions_beta_shapes(alpha)
  alpha_alpha <- beta_shapes[["shape1"]]
  alpha_beta  <- beta_shapes[["shape2"]]

  lik <- cbind(
    dpoint(x, location = 1, log = log),
    do.call(cbind, lapply(2:ncol(alpha), function(i){
      stats::dbeta(x, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], log = log)
    }))
  )

  return(lik)
}
.mdone.sided_fixed     <- function(x, omega, log){

  .weightfunctions_check_omega(omega, "omega")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(omega)){
    omega <- matrix(omega, nrow = 1)
  }
  if(length(x) == 1){
    x <- rep(x, nrow(omega))
  }

  if(nrow(omega) != length(x) & nrow(omega) != 1)
    stop("Non matching dimensions of 'omega' and 'x'.")

  if(nrow(omega) != length(x) & nrow(omega) == 1){
    omega <- do.call(rbind, lapply(1:length(x), function(i)omega))
  }

  lik <- do.call(cbind, lapply(1:ncol(omega), function(i)dpoint(x, location = omega[,i], log = log)))
}

#### random number generators ####
.rone.sided_general   <- function(n, alpha1, alpha2){

  # input check
  .weightfunctions_check_alpha(alpha1, "alpha1")
  .weightfunctions_check_alpha(alpha2, "alpha2")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha1)){
    alpha1 <- matrix(alpha1, nrow = 1)
  }
  if(!is.matrix(alpha2)){
    alpha2 <- matrix(alpha2, nrow = 1)
  }

  if(nrow(alpha1) != nrow(alpha2))
    stop("Non matching dimensions of 'alpha1' and 'alpha2'.")
  if(nrow(alpha1) != n & nrow(alpha1) != 1)
    stop("Incompatible dimensions of requested number of samples and 'alpha'.")

  if(nrow(alpha1) != n & nrow(alpha1) == 1){
    alpha1 <- do.call(rbind, lapply(1:n, function(i)alpha1))
    alpha2 <- do.call(rbind, lapply(1:n, function(i)alpha2))
  }

  x1 <- .prior_dirichlet_draws(alpha1, n)
  x2 <- .prior_dirichlet_draws(alpha2, n) * (1 - x1[,1])

  x <- matrix(ncol = ncol(alpha1) + ncol(alpha2) - 1, nrow = n)

  x[,ncol(alpha2):ncol(x)] <- t(apply(x1, 1, cumsum))
  for(i in 2:ncol(alpha2)){
    x[,i-1] = apply(x2[,i:ncol(alpha2), drop = FALSE], 1, sum) + x1[,1]
  }

  return(x[,ncol(x):1,drop = FALSE])
}
.rone.sided_monotonic <- function(n, alpha){

  # input check
  .weightfunctions_check_alpha(alpha, "alpha")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }

  if(nrow(alpha) != n & nrow(alpha) != 1)
    stop("Incompatible dimensions of requested number of samples and 'alpha'.")

  if(nrow(alpha) != n & nrow(alpha) == 1){
    alpha <- do.call(rbind, lapply(1:n, function(i)alpha))
  }

  x <- .prior_dirichlet_draws(alpha, n)
  x <- t(apply(x[,ncol(x):1, drop = FALSE], 1, cumsum))[,ncol(x):1, drop = FALSE]
  x[,1] <- 1

  return(x)
}
.rone.sided_fixed     <- function(n, omega){

  .weightfunctions_check_omega(omega, "omega")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(omega)){
    omega <- matrix(omega, nrow = 1)
  }

  if(nrow(omega) != n & nrow(omega) != 1)
    stop("Incompatible dimensions of requested number of samples and 'omega'.")

  if(nrow(omega) != n & nrow(omega) == 1){
    omega <- do.call(rbind, lapply(1:n, function(i)omega))
  }

  x <- omega

  return(x)
}

#### marginal distribution functions ####
.mpone.sided_general   <- function(q, alpha1, alpha2, lower.tail, log.p){

  stop("Not implemented")
}
.mpone.sided_monotonic <- function(q, alpha, lower.tail, log.p){

  # input check
  .weightfunctions_check_alpha(alpha, "alpha")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  if(length(q) == 1){
    q <- rep(q, nrow(alpha))
  }

  if(nrow(alpha) != length(q) & nrow(alpha) != 1)
    stop("Non matching dimensions of 'alpha' and 'q'.")

  if(nrow(alpha) != length(q) & nrow(alpha) == 1){
    alpha <- do.call(rbind, lapply(1:length(q), function(i)alpha))
  }

  # marginals of reversed cumulative sums of Dirichlet variables are beta
  # distributed; the first/reference publication weight is fixed to one.
  beta_shapes <- .weightfunctions_beta_shapes(alpha)
  alpha_alpha <- beta_shapes[["shape1"]]
  alpha_beta  <- beta_shapes[["shape2"]]

  p <- cbind(
    ppoint(q, location = 1, lower.tail = lower.tail, log.p = log.p),
    do.call(cbind, lapply(2:ncol(alpha), function(i){
      stats::pbeta(q, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], log.p = log.p, lower.tail = lower.tail)
    }))
  )

  return(p)
}
.mpone.sided_fixed     <- function(q, omega, lower.tail, log.p){

  .weightfunctions_check_omega(omega, "omega")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(omega)){
    omega <- matrix(omega, nrow = 1)
  }
  if(length(q) == 1){
    q <- rep(q, nrow(omega))
  }

  if(nrow(omega) != length(q) & nrow(omega) != 1)
    stop("Non matching dimensions of 'omega' and 'q'.")

  if(nrow(omega) != length(q) & nrow(omega) == 1){
    omega <- do.call(rbind, lapply(1:length(q), function(i)omega))
  }

  p <- do.call(cbind, lapply(1:ncol(omega), function(i)ppoint(q, location = omega[,i], lower.tail = lower.tail, log.p = log.p)))

  return(p)
}

#### marginal quantile functions ####
.mqone.sided_general   <- function(p, alpha1, alpha2, lower.tail, log.p){

  stop("Not implemented")
}
.mqone.sided_monotonic <- function(p, alpha, lower.tail, log.p){

  # input check
  .weightfunctions_check_alpha(alpha, "alpha")


  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  if(length(p) == 1){
    p <- rep(p, nrow(alpha))
  }

  if(nrow(alpha) != length(p) & nrow(alpha) != 1)
    stop("Non matching dimensions of 'alpha' and 'p'.")

  if(nrow(alpha) != length(p) & nrow(alpha) == 1){
    alpha <- do.call(rbind, lapply(1:length(p), function(i)alpha))
  }

  # marginals of reversed cumulative sums of Dirichlet variables are beta
  # distributed; the first/reference publication weight is fixed to one.
  beta_shapes <- .weightfunctions_beta_shapes(alpha)
  alpha_alpha <- beta_shapes[["shape1"]]
  alpha_beta  <- beta_shapes[["shape2"]]

  q <- cbind(
    qpoint(p, location = 1, lower.tail = lower.tail, log.p = log.p),
    do.call(cbind, lapply(2:ncol(alpha), function(i){
      stats::qbeta(p, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], lower.tail = lower.tail, log.p = log.p)
    }))
  )

  # make sure that the range is correct qpoint operates on unrestricted range
  q <- ifelse(q == -Inf, 0, q)
  q <- ifelse(q ==  Inf, 1, q)
  return(q)
}
.mqone.sided_fixed     <- function(p, omega, lower.tail, log.p){

  .weightfunctions_check_omega(omega, "omega")

  # transform to matrices for easier manipulation and checks
  if(!is.matrix(omega)){
    omega <- matrix(omega, nrow = 1)
  }
  if(length(p) == 1){
    p <- rep(p, nrow(omega))
  }

  if(nrow(omega) != length(p) & nrow(omega) != 1)
    stop("Non matching dimensions of 'omega' and 'p'.")

  if(nrow(omega) != length(p) & nrow(omega) == 1){
    omega <- do.call(rbind, lapply(1:length(p), function(i)omega))
  }

  q <- do.call(cbind, lapply(1:ncol(omega), function(i)qpoint(p, location = omega[,i], lower.tail = lower.tail, log.p = log.p)))

  # make sure that the range is correct
  q <- ifelse(q == -Inf, 0, q)
  q <- ifelse(q ==  Inf, omega, q)
  return(q)
}

### helper functions
# The Beta shapes of the marginals of the monotone weight function: the
# reference weight is fixed to one and the weight j has the shapes
# sum(alpha[j:J]) and sum(alpha[1:(j - 1)]). Finite concentrations can have
# cumulative sums beyond the largest double; dbeta(), pbeta() and qbeta() then
# silently evaluate a different distribution, so such shapes are refused.
.weightfunctions_beta_shapes <- function(alpha){

  shape1 <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[, ncol(alpha):1, drop = FALSE]
  shape2 <- cbind(0, t(apply(alpha, 1, cumsum))[, -ncol(alpha), drop = FALSE])

  if(any(!is.finite(shape1[, -1L])) || any(!is.finite(shape2[, -1L]))){
    stop(errorCondition(
      message = paste0(
        "The Beta shapes of the weight-function marginals are not representable: ",
        "the cumulative sums of 'alpha' exceed the largest double. ",
        "Use less extreme 'alpha' values."
      ),
      call = NULL,
      class = c("BayesTools_weightfunction_beta_shape_unavailable", "BayesTools_weightfunction_unavailable")
    ))
  }

  list(shape1 = shape1, shape2 = shape2)
}
.weightfunctions_one_sided_parameterization <- function(alpha, alpha1, alpha2){
  has_alpha  <- !is.null(alpha)
  has_alpha1 <- !is.null(alpha1)
  has_alpha2 <- !is.null(alpha2)

  if(has_alpha && (has_alpha1 || has_alpha2))
    stop("Use exactly one one-sided weightfunction parameterization: either 'alpha' or both 'alpha1' and 'alpha2'.")
  if(has_alpha)
    return("monotonic")
  if(has_alpha1 && has_alpha2)
    return("general")
  if(has_alpha1 || has_alpha2)
    stop("The general one-sided weightfunction parameterization requires both 'alpha1' and 'alpha2'.")

  stop("A one-sided weightfunction parameterization is required: use either 'alpha' or both 'alpha1' and 'alpha2'.")
}
.weightfunctions_check_alpha <- function(alpha, name = "alpha"){
  if(!is.numeric(alpha) | !(is.vector(alpha) | is.matrix(alpha)))
    stop(paste0("'", name, "' must be a numeric vector or a matrix."))
  if(is.vector(alpha))if(length(alpha) < 2)
    stop(paste0("'", name, "' must be a vector of length at least 2."))
  if(is.matrix(alpha))if(ncol(alpha) < 2)
    stop(paste0("'", name, "' must be a matrix with at least 2 columns."))
  if(any(!is.finite(alpha)))
    stop(paste0("'", name, "' must be finite."))
  if(!all(alpha > 0))
    stop(paste0("'", name, "' must be positive."))
}
.weightfunctions_check_omega <- function(omega, name = "omega"){
  if(!is.numeric(omega) | !(is.vector(omega) | is.matrix(omega)))
    stop(paste0("'", name, "' must be a numeric vector or a matrix."))
  if(is.vector(omega))if(length(omega) < 2)
    stop(paste0("'", name, "' must be a vector of length at least 2."))
  if(is.matrix(omega))if(ncol(omega) < 2)
    stop(paste0("'", name, "' must be a matrix with at least 2 columns."))
  if(any(!is.finite(omega)))
    stop(paste0("'", name, "' must be finite."))
  if(!all(omega >= 0))
    stop(paste0("'", name, "' must be non-negative."))
  if(is.vector(omega) && omega[1] != 1)
    stop(paste0("The reference-bin '", name, "' weight must be exactly 1."))
  if(is.matrix(omega) && any(omega[,1] != 1))
    stop(paste0("The reference-bin '", name, "' weight must be exactly 1."))
}

# Step-plot vertex indices of a weightfunction with 'n_cuts' cut points
# (including the 0 and 1 end points): x indexes the cut points and y the bins.
.weightfunction_step_indices <- function(n_cuts){

  if(!is.numeric(n_cuts) || length(n_cuts) != 1L || !is.finite(n_cuts) ||
     n_cuts < 2 || n_cuts != as.integer(n_cuts)){
    stop("'n_cuts' must be an integer of at least 2.", call. = FALSE)
  }
  n_cuts <- as.integer(n_cuts)
  n_bins <- n_cuts - 1L
  if(n_cuts == 2L){
    return(list(x = c(1L, 2L), y = c(1L, 1L)))
  }

  list(
    x = c(1L, sort(rep(2L:(n_cuts - 1L), 2L)), n_cuts),
    y = sort(rep(seq_len(n_bins), 2L))
  )
}
# Plotted mean of a weight that has none (a heavy right tail on the log scale
# makes E[omega] infinite) or none in double precision (it exceeds the largest
# double, .weightfunction_plot_component_mean()): the mean line is omitted for
# that bin, so NA breaks the line, and the y range comes from the finite values.
.weightfunction_plot_mean <- function(x_mean){

  x_mean[!is.finite(x_mean)] <- NA_real_
  x_mean
}
# The mean of a weight component for plotting, from 'mean_function': a mean
# beyond the largest double (the range limit of the moment calculation, which
# mean() reports) is plotted like an infinite one.
.weightfunction_plot_component_mean <- function(mean_function, component){

  tryCatch(
    mean_function(component),
    BayesTools_numerical_range_limit = function(condition) Inf
  )
}
.weightfunction_plot_y_range <- function(x_mean, x_lCI, x_uCI){

  c(0, max(1, .weightfunction_plot_mean(x_mean), x_lCI, x_uCI, na.rm = TRUE))
}
