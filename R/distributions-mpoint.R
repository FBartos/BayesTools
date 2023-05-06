#' @title Multivariate point mass distribution
#'
#' @description Density, distribution function,
#' quantile function and random generation for multivariate point distribution.
#'
#' @param x,q vector or matrix of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location vector of locations corresponding to the location of individual points.
#' Alternatively, a matrix with rows corresponding to the location of individual samples
#' and columns correspond to the location of individual points.
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
#'
#' @examples
#' # draw samples from a multivariate point distribution
#' rmpoint(10, location = c(0, 1))
#'
#' @return \code{dpoint} gives the density, \code{ppoint} gives the
#' distribution function, \code{qpoint} gives the quantile function,
#' and \code{rpoint} generates random deviates.
#'
#' @export dmpoint
#' @export rmpoint
#' @export pmpoint
#' @export qmpoint
#' @name mpoint
NULL

#### wrappers ####
#' @rdname mpoint
dmpoint <- function(x, location, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(as.vector(x))
  check_real(as.vector(location), "location", check_length = 0)

  # transform arguments to the matrix form
  if(!is.matrix(x)){
    x <- matrix(x, nrow = 1)
  }
  if(!is.matrix(location)){
    if(length(location) == 1){
      location <- rep(location, ncol(x))
    }
    location <- matrix(location, nrow = 1)
  }

  if((nrow(location) != nrow(x) & nrow(location) != 1) || ncol(location) != ncol(x))
    stop("Non matching dimensions of 'location' and 'x'.")

  if(nrow(location) != nrow(x) & nrow(location) == 1){
    location <- matrix(location, nrow = nrow(x), ncol = ncol(location), byrow = TRUE)
  }

  lik <- sapply(1:nrow(x), function(i){
    if(isTRUE(all.equal(location[i,], x[i,]))){
      return(Inf)
    }else{
      return(0)
    }
  })

  if(log){
    lik <- log(lik)
  }

  return(lik)
}
#' @rdname mpoint
rmpoint <- function(n, location){

  # common input check
  .check_n(n)
  check_real(as.vector(location), "location", check_length = 0)

  # transform arguments to the matrix form
  if(!is.matrix(location)){
    location <- matrix(location, nrow = 1)
  }

  if(nrow(location) != n & nrow(location) != 1)
    stop("Incompatible dimensions of requested number of samples and 'location'.")

  # 'sample' distribution
  if(nrow(location) == n){
    return(location)
  }else{
    x <- matrix(location, byrow = TRUE, ncol = ncol(location), nrow = n)
  }

  return(x)
}
#' @rdname mpoint
pmpoint <- function(q, location, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(as.vector(q))
  check_real(as.vector(location), "location", check_length = 0)

  # transform arguments to the matrix form
  if(!is.matrix(q)){
    q <- matrix(q, nrow = 1)
  }

  if(!is.matrix(location)){
    if(length(location) == 1){
      location <- rep(location, ncol(q))
    }
    location <- matrix(location, nrow = 1)
  }

  if(nrow(location) != nrow(q) & nrow(location) != 1 || ncol(location) != ncol(q))
    stop("Non matching dimensions of 'location' and 'q'.")

  if(nrow(location) != nrow(q) & nrow(location) == 1){
    location <- matrix(location, nrow = nrow(q), ncol = ncol(location), byrow = TRUE)
  }

  p <- as.numeric(apply(location <= q, 1, all))

  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }

  return(p)
}
#' @rdname mpoint
qmpoint <- function(p, location, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)
  check_real(as.vector(location), "location", check_length = 0)

  # transform arguments to the matrix form
  if(!is.matrix(location)){
    location <- matrix(location, nrow = 1)
  }

  if(nrow(location) != length(p) & nrow(location) != 1)
    stop("Non matching dimensions of 'location' and 'p'.")

  if(nrow(location) != length(p) & nrow(location) == 1){
    location <- matrix(location, nrow = length(p), ncol = ncol(location), byrow = TRUE)
  }

  if(log.p){
    p <- exp(p)
  }

  q <- do.call(rbind, lapply(1:length(p), function(i){
    if(all(p[i] > 0)){
      return(location[i,])
    }else{
      return(matrix(ifelse(lower.tail, -Inf, Inf), nrow = 1, ncol = ncol(location)))
    }
  }))

  return(q)
}
