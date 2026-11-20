# The log of the smaller Laplace tail, shared by the density and the CDF.
# Keep a finite standardized distance if the original subtraction overflows.
.laplace_log_tail <- function(x, location, scale){

  distance <- x - location
  overflow <- is.finite(x) & is.infinite(distance)
  standardized <- distance / scale
  standardized[overflow] <- x[overflow] / scale - location / scale
  -abs(standardized) - log(2)
}

.dlaplace_prior <- function(x, location, scale, log = FALSE){

  log_density <- .laplace_log_tail(x, location, scale) - log(scale)
  if(log) log_density else exp(log_density)
}

.plaplace_prior <- function(q, location, scale, lower.tail = TRUE, log.p = FALSE){

  log_tail <- .laplace_log_tail(q, location, scale)
  complement <- !is.na(q) & if(lower.tail) q > location else q < location
  log_tail[complement] <- .prior_log1mexp(log_tail[complement])
  if(log.p) log_tail else exp(log_tail)
}

.qlaplace_prior <- function(p, location, scale, lower.tail = TRUE, log.p = FALSE){

  if(log.p){
    small <- p <= -log(2)
    log_tail <- p
    complement <- !is.na(p) & !small
    log_tail[complement] <- .prior_log1mexp(p[complement])
    standardized <- log(2) + log_tail
  }else{
    small <- p <= .5
    standardized <- log(2 * pmin(p, 1 - p))
  }
  positive <- !is.na(p) & if(lower.tail) !small else small
  standardized[positive] <- -standardized[positive]
  offset <- scale * standardized
  value <- location + offset
  overflow <- is.finite(standardized) & is.infinite(offset)
  value[overflow] <- (location / scale + standardized[overflow]) * scale
  value
}
