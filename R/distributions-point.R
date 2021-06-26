#' @title Point mass distribution
#'
#' @description Density, distribution function,
#' quantile function and random generation for point distribution.
#'
#' @param x,q vector or matrix of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location vector of locations.
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
#'
#' @examples
#' # draw samples from a point distribution
#' rpoint(10, location = 1)
#'
#' @export dpoint
#' @export rpoint
#' @export ppoint
#' @export qpoint
#' @name point
NULL

#### wrappers ####
#' @rdname point
dpoint <- function(x, location, log = FALSE){

  # common input check
  .check_log(log)
  .check_x(x)

  if(length(location) != length(x) & length(location) != 1)
    stop("Non matching dimensions of 'location' and 'x'.")

  if(length(location) != length(x) & length(location) == 1){
    location <- rep(location, length(x))
  }

  lik <- sapply(1:length(x), function(i){
    if(isTRUE(all.equal(location[i], x[i]))){
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
#' @rdname point
rpoint <- function(n, location){

  # common input check
  .check_n(n)

  if(length(location) != n & length(location) != 1)
    stop("Incompatible dimensions of requested number of samples and 'location'.")

  if(length(location) != n & length(location) == 1){
    x <- rep(location, n)
  }else{
    x <- location
  }

  return(x)
}
#' @rdname point
ppoint <- function(q, location, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_q(q)

  if(length(location) != length(q) & length(location) == 1){
    location <- rep(location, length(q))
  }

  p <- ifelse(location <= q, 1, 0)

  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }

  return(p)
}
#' @rdname point
qpoint <- function(p, location, lower.tail = TRUE, log.p = FALSE){

  # common input check
  .check_log.p(log.p)
  .check_lower.tail(lower.tail)
  .check_p(p, log.p)

  if(length(location) != length(p) & length(location) == 1){
    location <- rep(location, length(p))
  }

  if(log.p){
    p <- exp(p)
  }

  q <- ifelse(p > 0, location, ifelse(lower.tail, -Inf, Inf))

  return(q)
}
