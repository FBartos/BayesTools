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
#' @param omega vector or matrix of fixed non-negative relative weights for a
#' one.sided or a two.sided weight function. The first/reference weight must be
#' exactly 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @details General one-sided marginals use analytic beta distributions for
#' the expected-direction weights and adaptive quadrature for the
#' product-of-beta representation of the unexpected-direction weights.
#' Quantiles numerically invert the corresponding tail probability and verify
#' the achieved probability against a fixed error budget. Results from the
#' general parameterization carry a \code{numerical_provenance} attribute with
#' the method, tolerances, and component-level diagnostics.
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

  inputs <- .weightfunctions_general_inputs(
    x,
    alpha1,
    alpha2,
    value_name = "x"
  )
  x <- inputs$value
  alpha1 <- inputs$alpha1
  alpha2 <- inputs$alpha2

  expected <- .weightfunctions_general_expected(
    x,
    alpha1,
    operation = "density"
  )
  unexpected <- matrix(
    NA_real_,
    nrow = length(x),
    ncol = ncol(alpha2) - 1L
  )
  quadrature <- vector("list", length(x) * ncol(unexpected))
  dim(quadrature) <- dim(unexpected)

  for(row in seq_along(x)){
    shapes <- .weightfunctions_general_shapes(alpha1[row, ], alpha2[row, ])
    for(component in seq_along(shapes)){
      result <- .weightfunctions_general_density(
        x[[row]],
        shapes[[component]]
      )
      unexpected[row, component] <- result$value
      quadrature[[row, component]] <- result$provenance
    }
  }

  out <- cbind(expected, unexpected)
  if(log){
    out <- log(out)
  }
  attr(out, "numerical_provenance") <-
    .weightfunctions_general_provenance("density", quadrature)
  out
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
  alpha_alpha <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[, ncol(alpha):1, drop = FALSE]
  alpha_beta  <- cbind(0, t(apply(alpha, 1, cumsum))[, -ncol(alpha), drop = FALSE])

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

  do.call(cbind, lapply(1:ncol(omega), function(i)dpoint(x, location = omega[,i], log = log)))
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

  x <- extraDistr::rdirichlet(n, alpha = alpha)
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

  inputs <- .weightfunctions_general_inputs(
    q,
    alpha1,
    alpha2,
    value_name = "q"
  )
  q <- inputs$value
  alpha1 <- inputs$alpha1
  alpha2 <- inputs$alpha2

  expected <- .weightfunctions_general_expected(
    q,
    alpha1,
    operation = "distribution",
    lower.tail = lower.tail
  )
  unexpected <- matrix(
    NA_real_,
    nrow = length(q),
    ncol = ncol(alpha2) - 1L
  )
  quadrature <- vector("list", length(q) * ncol(unexpected))
  dim(quadrature) <- dim(unexpected)

  for(row in seq_along(q)){
    shapes <- .weightfunctions_general_shapes(alpha1[row, ], alpha2[row, ])
    for(component in seq_along(shapes)){
      result <- .weightfunctions_general_probability(
        q[[row]],
        shapes[[component]],
        lower.tail = lower.tail
      )
      unexpected[row, component] <- result$value
      quadrature[[row, component]] <- result$provenance
    }
  }

  out <- cbind(expected, unexpected)
  if(log.p){
    out <- log(out)
  }
  attr(out, "numerical_provenance") <-
    .weightfunctions_general_provenance("distribution", quadrature)
  out
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
  alpha_alpha <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[, ncol(alpha):1, drop = FALSE]
  alpha_beta  <- cbind(0, t(apply(alpha, 1, cumsum))[, -ncol(alpha), drop = FALSE])

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

  inputs <- .weightfunctions_general_inputs(
    p,
    alpha1,
    alpha2,
    value_name = "p"
  )
  p <- inputs$value
  alpha1 <- inputs$alpha1
  alpha2 <- inputs$alpha2

  expected <- .weightfunctions_general_expected(
    p,
    alpha1,
    operation = "quantile",
    lower.tail = lower.tail,
    log.p = log.p
  )
  probabilities <- if(log.p) exp(p) else p
  if(log.p && any(is.finite(p) & probabilities == 0)){
    stop(
      "A finite log probability underflowed on the probability scale and ",
      "cannot be inverted faithfully.",
      call. = FALSE
    )
  }
  unexpected <- matrix(
    NA_real_,
    nrow = length(p),
    ncol = ncol(alpha2) - 1L
  )
  roots <- vector("list", length(p) * ncol(unexpected))
  dim(roots) <- dim(unexpected)

  for(row in seq_along(p)){
    shapes <- .weightfunctions_general_shapes(alpha1[row, ], alpha2[row, ])
    for(component in seq_along(shapes)){
      result <- .weightfunctions_general_quantile(
        probabilities[[row]],
        shapes[[component]],
        lower.tail = lower.tail
      )
      unexpected[row, component] <- result$value
      roots[[row, component]] <- result$provenance
    }
  }

  out <- cbind(expected, unexpected)
  attr(out, "numerical_provenance") <-
    .weightfunctions_general_provenance("quantile", roots)
  out
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
  alpha_alpha <- t(apply(alpha[,ncol(alpha):1, drop = FALSE], 1, cumsum))[, ncol(alpha):1, drop = FALSE]
  alpha_beta  <- cbind(0, t(apply(alpha, 1, cumsum))[, -ncol(alpha), drop = FALSE])

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
.weightfunctions_general_control <- function(){
  list(
    relative_tolerance = 1e-8,
    absolute_tolerance = 1e-10,
    subdivisions = 500L,
    quantile_tolerance = 1e-9,
    probability_tolerance = 5e-8
  )
}

.weightfunctions_general_inputs <- function(
    value,
    alpha1,
    alpha2,
    value_name){

  .weightfunctions_check_alpha(alpha1, "alpha1")
  .weightfunctions_check_alpha(alpha2, "alpha2")
  if(!is.matrix(alpha1)){
    alpha1 <- matrix(alpha1, nrow = 1L)
  }
  if(!is.matrix(alpha2)){
    alpha2 <- matrix(alpha2, nrow = 1L)
  }
  if(nrow(alpha1) != nrow(alpha2)){
    stop(
      "Non matching dimensions of 'alpha1' and 'alpha2'.",
      call. = FALSE
    )
  }
  if(length(value) == 1L){
    value <- rep(value, nrow(alpha1))
  }
  if(nrow(alpha1) != length(value) && nrow(alpha1) != 1L){
    stop(
      "Non matching dimensions of 'alpha' and '", value_name, "'.",
      call. = FALSE
    )
  }
  if(nrow(alpha1) == 1L && length(value) > 1L){
    alpha1 <- alpha1[rep.int(1L, length(value)), , drop = FALSE]
    alpha2 <- alpha2[rep.int(1L, length(value)), , drop = FALSE]
  }
  list(value = value, alpha1 = alpha1, alpha2 = alpha2)
}

.weightfunctions_general_shapes <- function(alpha1, alpha2){

  shape_A1 <- alpha1[[1L]]
  shape_A2 <- sum(alpha1[-1L])
  indices <- rev(seq.int(2L, length(alpha2)))
  lapply(indices, function(index){
    list(
      A1 = shape_A1,
      A2 = shape_A2,
      B1 = sum(alpha2[index:length(alpha2)]),
      B2 = sum(alpha2[seq_len(index - 1L)])
    )
  })
}

.weightfunctions_general_expected <- function(
    value,
    alpha,
    operation,
    lower.tail = TRUE,
    log.p = FALSE){

  output <- matrix(
    NA_real_,
    nrow = length(value),
    ncol = ncol(alpha)
  )
  output[, 1L] <- switch(
    operation,
    density = dpoint(value, location = 1),
    distribution = ppoint(
      value,
      location = 1,
      lower.tail = lower.tail
    ),
    quantile = qpoint(
      value,
      location = 1,
      lower.tail = lower.tail,
      log.p = log.p
    )
  )
  for(component in seq.int(2L, ncol(alpha))){
    cumulative_index <- ncol(alpha) - component + 1L
    shape1 <- rowSums(alpha[, seq_len(cumulative_index), drop = FALSE])
    shape2 <- rowSums(alpha[, seq.int(
      cumulative_index + 1L,
      ncol(alpha)
    ), drop = FALSE])
    output[, component] <- switch(
      operation,
      density = stats::dbeta(value, shape1, shape2),
      distribution = stats::pbeta(
        value,
        shape1,
        shape2,
        lower.tail = lower.tail
      ),
      quantile = stats::qbeta(
        value,
        shape1,
        shape2,
        lower.tail = lower.tail,
        log.p = log.p
      )
    )
  }
  if(operation == "quantile"){
    output[output == -Inf] <- 0
    output[output == Inf] <- 1
  }
  output
}

.weightfunctions_general_integrate <- function(fun, context){

  control <- .weightfunctions_general_control()
  integration <- tryCatch(
    stats::integrate(
      fun,
      lower = 0,
      upper = 1,
      subdivisions = control$subdivisions,
      rel.tol = control$relative_tolerance,
      abs.tol = control$absolute_tolerance,
      stop.on.error = FALSE
    ),
    error = function(e) e
  )
  valid <- !inherits(integration, "error") &&
    identical(integration$message, "OK") &&
    is.finite(integration$value) &&
    is.finite(integration$abs.error) &&
    integration$abs.error <= max(
      control$absolute_tolerance,
      control$relative_tolerance * abs(integration$value)
    )
  if(!valid){
    detail <- if(inherits(integration, "error")){
      conditionMessage(integration)
    }else{
      paste0(
        integration$message,
        "; absolute error ",
        format(integration$abs.error, digits = 6)
      )
    }
    stop(
      "General one-sided weight-function ", context,
      " quadrature failed: ", detail, ".",
      call. = FALSE
    )
  }
  list(
    value = integration$value,
    provenance = list(
      method = "adaptive quadrature",
      absolute_error = integration$abs.error,
      message = integration$message
    )
  )
}

.weightfunctions_general_density <- function(x, shapes){

  if(x == 0){
    exponent <- shapes$A1 + shapes$B1 - 1
    value <- if(exponent < 0){
      Inf
    }else if(exponent > 0){
      0
    }else{
      exp(
        lbeta(shapes$A1, shapes$B1) -
          lbeta(shapes$A1, shapes$A2) -
          lbeta(shapes$B1, shapes$B2)
      )
    }
    return(list(
      value = value,
      provenance = list(method = "analytic lower-boundary limit")
    ))
  }
  if(x == 1){
    U1 <- shapes$A2
    U2 <- shapes$A1
    V1 <- shapes$B2
    V2 <- shapes$B1
    minimum_shape <- min(U1, V1)
    value <- if(minimum_shape < 1){
      Inf
    }else if(U1 == V1){
      if(U1 == 1) Inf else 0
    }else if(minimum_shape > 1){
      0
    }else if(U1 == 1){
      exp(-lbeta(U1, U2)) * (V1 + V2 - 1) / (V1 - 1)
    }else{
      exp(-lbeta(V1, V2)) * (U1 + U2 - 1) / (U1 - 1)
    }
    return(list(
      value = value,
      provenance = list(method = "analytic upper-boundary limit")
    ))
  }

  result <- .weightfunctions_general_integrate(
    function(t){
      out <- numeric(length(t))
      interior <- t > 0 & t < 1
      if(any(interior)){
        A <- x * t[interior]
        B <- (x - A) / (1 - A)
        log_integrand <- stats::dbeta(
          A,
          shapes$A1,
          shapes$A2,
          log = TRUE
        ) +
          stats::dbeta(B, shapes$B1, shapes$B2, log = TRUE) -
          log1p(-A) +
          log(x)
        out[interior] <- exp(log_integrand)
      }
      out
    },
    context = paste0("density at x = ", format(x, digits = 17))
  )
  result
}

.weightfunctions_general_probability <- function(q, shapes, lower.tail){

  if(q == 0){
    return(list(
      value = if(lower.tail) 0 else 1,
      provenance = list(method = "analytic support boundary")
    ))
  }
  if(q == 1){
    return(list(
      value = if(lower.tail) 1 else 0,
      provenance = list(method = "analytic support boundary")
    ))
  }

  integral <- .weightfunctions_general_integrate(
    function(t){
      out <- numeric(length(t))
      interior <- t > 0 & t < 1
      if(any(interior)){
        A <- q * t[interior]
        B <- (q - A) / (1 - A)
        conditional <- stats::pbeta(
          B,
          shapes$B1,
          shapes$B2,
          lower.tail = lower.tail
        )
        out[interior] <- q *
          stats::dbeta(A, shapes$A1, shapes$A2) *
          conditional
      }
      out
    },
    context = paste0(
      if(lower.tail) "lower-tail" else "upper-tail",
      " probability at q = ",
      format(q, digits = 17)
    )
  )
  if(!lower.tail){
    tail_A <- stats::pbeta(
      q,
      shapes$A1,
      shapes$A2,
      lower.tail = FALSE
    )
    integral$value <- tail_A + integral$value
    integral$provenance$analytic_tail_A <- tail_A
  }
  integral$value <- min(1, max(0, integral$value))
  integral
}

.weightfunctions_general_quantile <- function(p, shapes, lower.tail){

  boundary <- if(p == 0){
    if(lower.tail) 0 else 1
  }else if(p == 1){
    if(lower.tail) 1 else 0
  }else{
    NULL
  }
  if(!is.null(boundary)){
    return(list(
      value = boundary,
      provenance = list(method = "analytic probability boundary")
    ))
  }

  control <- .weightfunctions_general_control()
  probability_evaluations <- 0L
  objective <- function(q){
    probability_evaluations <<- probability_evaluations + 1L
    .weightfunctions_general_probability(
      q,
      shapes,
      lower.tail = lower.tail
    )$value - p
  }
  root <- tryCatch(
    stats::uniroot(
      objective,
      interval = c(0, 1),
      tol = control$quantile_tolerance
    ),
    error = function(e) e
  )
  if(inherits(root, "error") || !is.finite(root$root)){
    stop(
      "General one-sided weight-function quantile inversion failed for p = ",
      format(p, digits = 17), ": ",
      if(inherits(root, "error")) conditionMessage(root) else
        "the root was not finite",
      ".",
      call. = FALSE
    )
  }
  achieved <- .weightfunctions_general_probability(
    root$root,
    shapes,
    lower.tail = lower.tail
  )$value
  probability_error <- abs(achieved - p)
  if(!is.finite(probability_error) ||
     probability_error > control$probability_tolerance){
    stop(
      "General one-sided weight-function quantile inversion did not meet ",
      "the probability error budget for p = ",
      format(p, digits = 17), "; absolute error ",
      format(probability_error, digits = 6), ".",
      call. = FALSE
    )
  }
  list(
    value = root$root,
    provenance = list(
      method = "CDF root inversion",
      iterations = root$iter,
      probability_evaluations = probability_evaluations + 1L,
      achieved_probability = achieved,
      absolute_probability_error = probability_error
    )
  )
}

.weightfunctions_general_provenance <- function(operation, details){
  list(
    method = paste("general one-sided", operation),
    control = .weightfunctions_general_control(),
    details = details
  )
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
