# Integrate on the positive half of a centered shrinkage prior. A
# difference of two separately approximated tails is ill conditioned for
# small intervals. The density is evaluated on the interval itself.
.regularized_horseshoe_center <- function(prior, b){

  # Craig's representation, as in .regularized_horseshoe_rule(): the
  # integrand of P(0 < beta < b) is 1 - F(r) H(a), where F(r) is the
  # scaled complementary error function and H is the slab Laplace
  # transform. Take that complement *before* integrating the scales.
  spec <- prior$regularization
  previous <- NA_real_
  for(panels in c(1L, 2L, 4L, 8L, 16L, 32L)){
    rule <- .regularized_horseshoe_rule(prior, panels, near_pole = TRUE)
    log_r <- log(b) + rule$log_r
    log_complement <- numeric(length(log_r))
    small <- log_r < log(1e-5)
    # 1-F(r) = 2r/sqrt(pi) - r² + 4r³/(3sqrt(pi)) + O(r^4).
    r <- exp(log_r[small])
    log_complement[small] <- log(2 / sqrt(pi)) + log_r[small] +
      log1p(r * (-sqrt(pi) / 2 + 2 * r / 3))
    r <- exp(log_r[!small])
    large <- r > 100
    log_f <- numeric(length(r))
    z <- r[!large]
    log_f[!large] <- log(2) + z^2 + stats::pnorm(-sqrt(2) * z, log.p = TRUE)
    # F(r) = (1 - 1/(2r²) + 3/(4r^4) + O(r^-6))/(sqrt(pi)*r).
    inv2 <- 0.5 / r[large]^2
    log_f[large] <- -log_r[!small][large] - log(pi) / 2 + log1p(inv2 * (3 * inv2 - 1))
    log_complement[!small] <- .prior_log1mexp(log_f)
    if(!is.null(spec$slab_scale)){
      log_u <- log(2 / spec$slab_df) + 2 * (log(b) + rule$log_a - log(spec$slab_scale))
      log_h <- -spec$slab_df / 2 * (pmax(log_u, 0) + log1p(exp(-abs(log_u))))
      # 1 - F H = (1 - H) + H (1 - F).
      log_complement <- .prior_normal_logaddexp(.prior_log1mexp(log_h), log_h + log_complement)
    }
    estimate <- .regularized_log_sum(rule$log_w + log_complement)
    if(is.finite(estimate) && is.finite(previous) && abs(estimate - previous) <= 1e-6){
      return(estimate)
    }
    previous <- estimate
  }
  .regularized_numerical_unavailable("central interval probability")
}

.regularized_positive_region <- function(prior, a, b){

  if(a == b) return(-Inf)
  if(is.infinite(b)){
    return(.regularized_marginal(prior, a, density = FALSE, log = TRUE, lower.tail = FALSE))
  }
  width <- b - a
  if(a > 0 && width <= a){
    # Narrow intervals, including adjacent representable endpoints. The
    # conditional densities are centered and decreasing on this half-line.
    bounds <- .regularized_marginal(prior, c(a, b), density = TRUE, log = TRUE)
    if(all(is.finite(bounds)) && abs(diff(bounds)) < 1e-8){
      return(log(width) + mean(bounds))
    }
  }else{
    # For separated tails the difference has condition number <= 5/3.
    # For a near-whole interval this avoids searching for a narrow peak
    # in an arbitrarily broad domain.
    tails <- tryCatch(.regularized_marginal(prior, c(a, b), density = FALSE, log = TRUE,
      lower.tail = FALSE), BayesTools_regularized_marginal_unavailable = function(e) NULL)
    if(!is.null(tails) && all(is.finite(tails)) && tails[2] <= tails[1] - log(4)){
      return(.prior_logdiffexp(tails[1], tails[2]))
    }
  }
  if(a == 0 && prior$regularization$family == "horseshoe"){
    return(.regularized_horseshoe_center(prior, b))
  }

  rule <- .prior_density_quadrature_rules()$finite
  # At the pole, evaluate the positive half-interval with an open rule
  # and a quadratic substitution to resolve the integrable singularity.
  # For positive endpoints many decades apart, use log-coordinates.
  logarithmic <- a > 0 && log(b) - log(a) > log(4)
  for(panels in c(1L, 2L, 4L, 8L, 16L)){
    t <- as.vector(outer((1 + rule$nodes) / 2, seq_len(panels) - 1L, "+")) / panels
    if(a == 0){
      nodes <- b * t^2
      log_jacobian <- log(b) + log(2) + log(t)
    }else if(logarithmic){
      span <- log(b) - log(a)
      log_nodes <- log(a) + span * t
      nodes <- exp(log_nodes)
      log_jacobian <- log(span) + log_nodes
    }else{
      nodes <- a + width * t
      log_jacobian <- log(width)
    }
    if(any(nodes == 0)) .regularized_numerical_unavailable("interval probability")
    kernel <- .regularized_marginal(prior, nodes, density = TRUE, log = TRUE) + log_jacobian
    kronrod <- .regularized_log_sum(kernel + log(rep(rule$kronrod, panels))) - log(2 * panels)
    gauss <- .regularized_log_sum(kernel + log(rep(rule$gauss, panels))) - log(2 * panels)
    if(is.finite(kronrod) && is.finite(gauss) && abs(kronrod - gauss) <= 1e-6) return(kronrod)
  }
  .regularized_numerical_unavailable("interval probability")
}

.regularized_region_log_mass <- function(prior, lower, upper){

  if(is.infinite(lower) || is.infinite(upper)){
    if(lower == -Inf && upper == Inf) return(0)
    if(lower == -Inf) return(.regularized_marginal(prior, upper, density = FALSE, log = TRUE))
    return(.regularized_marginal(prior, lower, density = FALSE, log = TRUE, lower.tail = FALSE))
  }
  if(upper <= 0) return(.regularized_positive_region(prior, -upper, -lower))
  if(lower >= 0) return(.regularized_positive_region(prior, lower, upper))
  .model_probability_log_sum(c(.regularized_positive_region(prior, 0, -lower),
                               .regularized_positive_region(prior, 0, upper)))
}
