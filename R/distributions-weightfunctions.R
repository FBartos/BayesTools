#' @title Weight functions
#'
#' @description Marginal density, marginal distribution function,
#' marginal quantile function and random generation for weight functions.
#'
#' @param x,q vector or matrix of quantiles.
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
#' @param omega vector or matrix of fixed probabilities for a one.sided or
#' a two.sided weight function.
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
#'
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

  # common input check
  .check_log(log)
  .check_x(x, lower = 0, upper = 1)

  if(!is.null(alpha)){
    lik <- .mdone.sided_monotonic(x, alpha, log)
  }else if(!is.null(alpha1) & !is.null(alpha2)){
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
  .check_x(x, lower = 0, upper = 1)

  lik <- .mdone.sided_fixed(x, omega, log)

  return(lik)
}
#' @rdname weightfunctions
mdtwo.sided_fixed <- function(x, omega, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(x, lower = 0, upper = 1)

  lik <- mdone.sided_fixed(x, omega = omega, log = log)

  return(lik)
}

#' @rdname weightfunctions
rone.sided <- function(n, alpha = NULL, alpha1 = NULL, alpha2 = NULL){

  # common input check
  .check_n(n)

  if(!is.null(alpha)){
    x <- .rone.sided_monotonic(n, alpha)
  }else if(!is.null(alpha1) & !is.null(alpha2)){
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

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0, upper = 1)

  if(!is.null(alpha)){
    p <- .mpone.sided_monotonic(q, alpha, lower.tail, log.p)
  }else if(!is.null(alpha1) & !is.null(alpha2)){
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
  .check_q(q, lower = 0, upper = 1)

  p <- .mpone.sided_fixed(q, omega, lower.tail, log.p)

  return(p)
}
#' @rdname weightfunctions
mptwo.sided_fixed <- function(q, omega, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q, lower = 0, upper = 1)

  p <- mpone.sided_fixed(q, omega = omega, lower.tail = lower.tail, log.p = log.p)

  return(p)
}

#' @rdname weightfunctions
mqone.sided <- function(p, alpha = NULL, alpha1 = NULL, alpha2 = NULL, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  if(!is.null(alpha)){
    q <- .mqone.sided_monotonic(p, alpha, lower.tail, log.p)
  }else if(!is.null(alpha1) & !is.null(alpha2)){
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

  stop("Not implemented")
  # would require product of two beta-distributed variables, since the weights in the opposite direction
  # start at the first cutoff modeled by the first alpha parameter
  # https://math.stackexchange.com/questions/1073364/product-of-two-beta-distributed-random-variables
  # https://www.dm.fct.unl.pt/sites/www.dm.fct.unl.pt/files/preprints/2012/7_12.pdf
  # probably solvable with Meijer g-function

  # input check
  .weightfunctions_check_alpha(alpha1, "alpha1")
  .weightfunctions_check_alpha(alpha2, "alpha2")

  # transform to matrices for easier manipulation and checks
  if(length(x) == 1){
    x <- rep(x, nrow(alpha1))
  }
  if(!is.matrix(alpha1)){
    alpha1 <- matrix(alpha1, nrow = 1)
  }
  if(!is.matrix(alpha2)){
    alpha2 <- matrix(alpha2, nrow = 1)
  }

  if(nrow(alpha1) != nrow(alpha2))
    stop("Non matching dimensions of 'alpha1' and 'alpha2'.")
  if(nrow(alpha1) != length(x) & nrow(alpha1) != 1)
    stop("Non matching dimensions of 'alpha' and 'x'.")

  if(nrow(alpha1) != length(x) & nrow(alpha1) == 1){
    alpha1 <- do.call(rbind, lapply(1:length(x), function(i)alpha1))
    alpha2 <- do.call(rbind, lapply(1:length(x), function(i)alpha2))
  }


  lik1 <- .mdone.sided_monotonic(x, alpha = alpha1, log = log)
  # the side in the unexpected direction starts with the first step of alpha1, e.g., something like this:
  # lik2 <- .dmone.sided_monotonic(x, alpha = alpha2, log = log) * beta(alpha[,1], t(apply(alpha[,-1], 1, sum)))
  #
  # return(cbind(lik2, lik1))
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

  # marginals of sums of Dirichlet distributed variables are beta distributed
  alpha_alpha <- t(apply(alpha, 1, cumsum))
  alpha_beta  <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[,(ncol(alpha)-1):1, drop = FALSE]

  lik <- do.call(cbind, lapply(1:(ncol(alpha) - 1), function(i)stats::dbeta(x, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], log = log)))
  lik <- cbind(lik, dpoint(x, location = 1))

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

  x1 <- extraDistr::rdirichlet(n, alpha = alpha1)
  x2 <- extraDistr::rdirichlet(n, alpha = alpha2) * (1 - x1[,1])

  x <- matrix(ncol = ncol(alpha1) + ncol(alpha2) - 1, nrow = n)

  x[,ncol(alpha2):ncol(x)] <- t(apply(x1, 1, cumsum))
  for(i in 2:ncol(alpha2)){
    x[,i-1] = apply(x2[,i:ncol(alpha2), drop = FALSE], 1, sum) + x1[,1]
  }

  return(x)
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

  x <- extraDistr::rdirichlet(n, alpha = alpha)
  x <- t(apply(x, 1, cumsum))

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
  if(length(q) == 1){
    q <- rep(q, nrow(alpha1))
  }

  if(nrow(alpha1) != nrow(alpha2))
    stop("Non matching dimensions of 'alpha1' and 'alpha2'.")
  if(nrow(alpha1) != length(q) & nrow(alpha1) != 1)
    stop("Non matching dimensions of 'alpha' and 'q'.")

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

  # marginals of sums of Dirichlet distributed variables are beta distributed
  alpha_alpha <- t(apply(alpha, 1, cumsum))
  alpha_beta  <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[,(ncol(alpha)-1):1, drop = FALSE]

  p <- do.call(cbind, lapply(1:(ncol(alpha) - 1), function(i)stats::pbeta(q, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], log.p = log.p, lower.tail = lower.tail)))
  p <- cbind(p, ppoint(q, location = 1, lower.tail = lower.tail, log.p = log.p))

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
  if(length(p) == 1){
    p <- rep(p, nrow(alpha1))
  }

  if(nrow(alpha1) != nrow(alpha2))
    stop("Non matching dimensions of 'alpha1' and 'alpha2'.")
  if(nrow(alpha1) != length(p) & nrow(alpha1) != 1)
    stop("Non matching dimensions of 'alpha' and 'p'.")

  if(nrow(alpha1) != length(p) & nrow(alpha1) == 1){
    alpha1 <- do.call(rbind, lapply(1:length(p), function(i)alpha1))
    alpha2 <- do.call(rbind, lapply(1:length(p), function(i)alpha2))
  }

}
.mqone.sided_monotonic <- function(p, alpha, lower.tail, log.p){

  # input check
  .weightfunctions_check_alpha(alpha, "alpha")


  # transform to matrices for easier manipulation and checks
  if(!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  if(length(q) == 1){
    p <- rep(p, nrow(alpha))
  }

  if(nrow(alpha) != length(p) & nrow(alpha) != 1)
    stop("Non matching dimensions of 'alpha' and 'p'.")

  if(nrow(alpha) != length(p) & nrow(alpha) == 1){
    alpha <- do.call(rbind, lapply(1:length(p), function(i)alpha))
  }

  # marginals of sums of Dirichlet distributed variables are beta distributed
  alpha_alpha <- t(apply(alpha, 1, cumsum))
  alpha_beta  <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[,(ncol(alpha)-1):1, drop = FALSE]

  q <- do.call(cbind, lapply(1:(ncol(alpha) - 1), function(i)stats::qbeta(p, shape1 = alpha_alpha[,i], shape2 = alpha_beta[,i], lower.tail = lower.tail, log.p = log.p)))
  q <- cbind(q, qpoint(p, location = 1, lower.tail = lower.tail, log.p = log.p))

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
  if(length(q) == 1){
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
  q <- ifelse(q ==  Inf, 1, q)
  return(q)
}

### helper functions
.weightfunctions_check_alpha <- function(alpha, name = "alpha"){
  if(!is.numeric(alpha) | !(is.vector(alpha) | is.matrix(alpha)))
    stop(paste0("'", name, "' must be a numeric vector or a matrix."))
  if(is.vector(alpha))if(length(alpha) < 2)
    stop(paste0("'", name, "' must be a vector of length at least 2."))
  if(is.matrix(alpha))if(ncol(alpha) < 2)
    stop(paste0("'", name, "' must be a matrix with at least 2 columns."))
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
  if(!all(omega >= 0 & omega <= 1))
    stop(paste0("'", name, "' must be between 0 and 1."))
}
