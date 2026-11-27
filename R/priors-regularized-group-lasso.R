# A direction of a d-dimensional group lasso has law sqrt(W) Z, with
# W ~ chi-square(d + 1) and Z standard Normal. Its density at zero is
# phi(0) E[1 / sqrt(W)] = B(d/2,1/2) / (2*pi).
.regularized_group_lasso_log_peak <- function(d){

  lbeta(d / 2, .5) - log(2 * pi)
}

# Log E[exp(-z^2/(2 W))], W ~ chi-square(df). The Bessel identity follows
# from integrating the gamma density. At the origin the expectation is 1.
# For df>2 its departure is bounded by z^2/[2*(df-2)]. For df>=2,
# splitting at W=z also bounds it by z since P(chi-square(2)<z)<=z/2.
.regularized_group_lasso_radial <- function(z, df, scaled = FALSE){

  value <- rep(-Inf, length(z))
  near <- z <= 1e-12
  if(df > 2) near <- near | z <= sqrt(2e-12 * (df - 2))
  value[near] <- if(scaled) z[near] else 0
  inside <- !near & is.finite(z)
  if(any(inside)){
    u <- z[inside]
    order <- df / 2
    # In this far tail the first omitted relative correction to the
    # scaled Bessel kernel is O((order^2/u)^2), below double precision.
    far <- u > 1e9 * max(1, order^2)
    log_k <- numeric(length(u))
    log_k[far] <- (log(pi / 2) - log(u[far])) / 2 + log1p((order^2 - .25) / (2 * u[far]))
    if(any(!far)) log_k[!far] <- log(besselK(u[!far], order, expon.scaled = TRUE))
    log_value <- log(2) + order * (log(u) - log(2)) + log_k - lgamma(order)
    # For high orders K_nu(u) can overflow even with exponential scaling.
    # Normalize the order recurrence itself. With F_nu = E exp(-z^2/(2W))
    # for W~chi-square(2nu), F_(nu+1)=F_nu+z^2 F_(nu-1)/[4nu(nu-1)].
    missing <- !is.finite(log_value)
    if(any(missing) && order > 2){
      start <- if(df %% 2 == 0) 1 else .5
      xx <- u[missing]
      prev <- if(start == .5) rep(0, length(xx)) else log(xx * besselK(xx, 1, expon.scaled = TRUE))
      curr <- if(start == .5) log1p(xx) else 2 * log(xx) + log(besselK(xx, 2, expon.scaled = TRUE)) - log(2)
      nu <- start + 1
      while(nu < order){
        next_value <- .prior_normal_logaddexp(curr, prev + 2 * log(xx) - log(4 * nu * (nu - 1)))
        prev <- curr
        curr <- next_value
        nu <- nu + 1
      }
      log_value[missing] <- curr
    }
    value[inside] <- if(scaled) log_value else log_value - u
    if(any(!is.finite(value[inside])) || (!scaled && any(value[inside] > 1e-10))){
      .regularized_numerical_unavailable("group-lasso scale integral")
    }
  }
  value
}

.regularized_group_lasso_log_kernel <- function(z, d, density){

  if(d %% 2L == 1L && d <= 63L){
    # X=G1-G2 for independent Gamma(k,1) variables, k=(d+1)/2.
    # At integer k, expanding the Gamma density/tail at z+G2 gives a
    # finite positive sum of Poisson probabilities with negative-binomial
    # weights. For larger k the Bessel representation avoids a long sum.
    k <- (d + 1L) / 2L
    l <- seq_len(k) - 1L
    coefficient <- if(density) stats::dnbinom(k - 1L - l, size = k, prob = .5, log = TRUE) else
      stats::pnbinom(k - 1L - l, size = k, prob = .5, log.p = TRUE)
    log_polynomial <- rep(coefficient[1L], length(z))
    positive <- is.finite(z) & z > 0
    for(i in l[l > 0L]){
      log_polynomial[positive] <- .prior_normal_logaddexp(log_polynomial[positive],
        coefficient[i + 1L] + i * log(z[positive]) - lfactorial(i))
    }
    return(log_polynomial - z)
  }
  if(density){
    return(.regularized_group_lasso_log_peak(d) + .regularized_group_lasso_radial(z, d))
  }
  # S(z)/f(z) integrates the monotone ratio f(z+t)/f(z) from zero to
  # infinity. Factor out exp(-z) and integrate on the scale of the tail
  # width. The integrand starts at 1 even if S(z) itself is tiny.
  out <- rep(-Inf, length(z))
  out[z == 0] <- -log(2)
  finite <- is.finite(z) & z > 0
  out[finite] <- vapply(z[finite], function(u){
    width <- max(1, min(sqrt(d + 1), (d + 1) / u))
    base <- .regularized_group_lasso_radial(u, d, scaled = TRUE)
    integral <- stats::integrate(function(v){
      t <- width * v
      exp(.regularized_group_lasso_radial(u + t, d, scaled = TRUE) - base - t)
    }, 0, Inf, rel.tol = 1e-9, abs.tol = 0, stop.on.error = FALSE)
    if(!identical(integral$message, "OK") || !is.finite(integral$value) || integral$value <= 0 ||
       !is.finite(integral$abs.error) || integral$abs.error > 1e-8 * integral$value){
      .regularized_numerical_unavailable("group-lasso tail probability")
    }
    .regularized_group_lasso_log_peak(d) + base - u + log(width) + log(integral$value)
  }, numeric(1))
  out
}
