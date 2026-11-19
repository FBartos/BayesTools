# Marginal distribution of one fitted coefficient. Conditional on the scales
# it is Normal (ridge/horseshoe) or Laplace. Integrate over the declared
# distributions of the scales, independently of posterior or prior draws.
.regularized_zero_log_density <- function(prior){

  spec <- prior$regularization
  if(spec$family == "horseshoe") return(Inf)
  kernel <- if(spec$family == "ridge") stats::dnorm(0, log = TRUE) else -log(2)
  if(!is.null(spec$group_lasso_dimension)) kernel <- .regularized_group_lasso_log_peak(spec$group_lasso_dimension)
  if(!is.null(spec$contrast_ratio)) kernel <- kernel - log1p(spec$contrast_ratio)
  if(is.numeric(spec$scale)) return(kernel - log(spec$multiplier) - log(spec$scale))
  at_zero <- .prior_density_ordinate_primitive(spec$scale, 0)
  if(at_zero$behavior %in% c("regular", "infinite")) return(Inf)
  if(at_zero$behavior != "zero") return(NA_real_)
  moment <- .prior_density_inverse_moment(spec$scale, 512L)
  if(!isTRUE(moment$available)) return(NA_real_)
  kernel - log(spec$multiplier) + moment$log_value
}

.regularized_numerical_unavailable <- function(operation){

  stop(structure(list(message = paste0("The regularized coefficient's ", operation,
    " is unavailable at the requested value and precision."), call = NULL),
    class = c("BayesTools_regularized_marginal_unavailable", "BayesTools_numerical_unavailable",
              "error", "condition")))
}

# Composite 21-point Gauss-Kronrod rule on prior probability coordinates.
# sin² maps the open nodes towards both tails and the Jacobian accounts for
# the change of variables. The same scale nodes serve a whole plot or table.
.regularized_probability_rule <- function(panels){
  rule <- .prior_density_quadrature_rules()$finite
  t <- as.vector(outer((1 + rule$nodes) / 2, seq_len(panels) - 1, "+")) / panels
  list(prob = sinpi(t / 2)^2, weight = rep(rule$kronrod, panels) * pi * sinpi(t) / (4 * panels))
}

.regularized_scale_rule <- function(prior, panels, variance = FALSE){

  rule <- .regularized_probability_rule(panels)
  prob <- rule$prob
  weight <- rule$weight
  spec <- prior$regularization
  global <- .regularized_global_scale_rule(spec$scale, panels, rule, 1)
  log_scale <- log(global$tau)
  log_weight <- global$log_weight
  if(variance){
    # c² = nu * slab_scale² / chi-square(nu)
    # Size-bias the slab by c^2: its finite second moment is the
    # normalizer, and its chi-square degrees of freedom decrease by two.
    df <- spec$slab_df - 2
    chi <- stats::qchisq(prob, df = df)
    if(any(!is.finite(chi) | chi <= 0)) .regularized_numerical_unavailable("variance")
    log_c <- log(spec$slab_scale) + (log(spec$slab_df) - log(chi)) / 2
    a <- rep(log_scale, length(log_c))
    b <- rep(log_c, each = length(log_scale))
    # E_lambda[sd² | tau,c] = c² * tau/(tau+c) under half-Cauchy(1).
    log_scale <- log(spec$slab_scale) + (log(spec$slab_df / df) +
      stats::plogis(a - b, log.p = TRUE)) / 2
    log_weight <- as.vector(outer(log_weight, log(weight), "+"))
  }
  list(log_scale = log_scale + log(spec$multiplier), log_weight = log_weight)
}

# For a truncated law, adding a very small interior probability to the
# base CDF at its truncation bound can round back to the bound. Integrate
# its normalized density directly on positive scale coordinates. The
# ordinary quartiles locate its mass; verify the integral on each piece
# against its known probability before using the nodes in a scale mixture.
.regularized_truncated_scale_rule <- function(prior, panels, power){

  lower <- prior$truncation$lower
  upper <- prior$truncation$upper
  cuts <- c(lower, quant(prior, c(.25, .5, .75)), upper)
  if(any(!is.finite(cuts[1:4])) || any(diff(cuts) <= 0)) .regularized_numerical_unavailable("scale integral")
  density <- .prior_simple_lpdf_evaluator(prior, purpose = "natural_integral")
  for(refinement in c(1L, 2L, 4L, 8L)){
    rule <- .regularized_probability_rule(panels * refinement)
    p <- rule$prob
    pieces <- lapply(seq_len(4L), function(i){
      a <- cuts[i]
      b <- cuts[i + 1L]
      if(is.infinite(b)){
        # The unbounded quartile has stable *upper* base probabilities.
        # Invert those directly to cover even very heavy polynomial tails.
        log_tail <- .prior_simple_base_p(prior, a, lower.tail = FALSE, log.p = TRUE) +
          power * log1p(-p)
        tau <- .prior_simple_base_q(prior, log_tail, lower.tail = FALSE, log.p = TRUE)
        if(any(!is.finite(tau) | tau <= a)) .regularized_numerical_unavailable("scale integral")
        return(list(tau = tau, log_weight = log(rule$weight * power / 4) + (power - 1) * log1p(-p)))
      }else if(a == 0){
        tau <- b * p^power
        log_jacobian <- log(b) + log(power) + (power - 1) * log(p)
      }else{
        tau <- a + (b - a) * p
        log_jacobian <- log(b - a)
      }
      if(any(!is.finite(tau) | tau <= a | tau >= b)) .regularized_numerical_unavailable("scale integral")
      list(tau = tau, log_weight = density(tau) + log(rule$weight) + log_jacobian)
    })
    log_prob <- vapply(pieces, function(piece) .regularized_log_sum(piece$log_weight), numeric(1))
    if(all(is.finite(log_prob)) && max(abs(log_prob + log(4))) <= 1e-8){
      return(list(tau = unlist(lapply(pieces, `[[`, "tau"), use.names = FALSE),
                  log_weight = unlist(lapply(pieces, `[[`, "log_weight"), use.names = FALSE)))
    }
  }
  .regularized_numerical_unavailable("scale normalization")
}

.regularized_global_scale_rule <- function(tau, panels, prob, global_power){

  p <- prob$prob^global_power
  if(!is.prior(tau)) return(list(tau = tau, log_weight = 0))
  global_log_weight <- log(global_power * prob$weight) + (global_power - 1) * log(prob$prob)
  # Inverting a symmetric untruncated CDF at (1+p)/2 rounds tiny
  # probabilities onto zero. The positive-half quantiles avoid that loss.
  half <- tau$truncation$lower == 0 && tau$truncation$upper == Inf
  if(half && identical(tau$distribution, "t") && tau$parameters$location == 0){
    df <- tau$parameters$df
    tau <- tau$parameters$scale * if(df == 1){
      slope <- tanpi(pmin(p, 1 - p) / 2)
      ifelse(p <= .5, slope, 1 / slope)
    }else{
      # t^2/(nu+t^2) ~ Beta(1/2, nu/2). Invert the smaller probability
      # and use the complementary Beta for the upper half. qf() may map
      # the lower half back through 1-p and erase positive quantiles.
      # Grade both tails: a large coefficient can depend on the
      # upper end of the global-scale distribution. Keep its small
      # upper-tail probability directly.
      left <- prob$prob <= .5
      interior <- 2 * pmin(prob$prob, 1 - prob$prob)
      tail <- interior^global_power / 2
      global_log_weight <- log(global_power * prob$weight) + (global_power - 1) * log(interior)
      u <- numeric(length(p))
      u[left] <- stats::qbeta(tail[left], .5, df / 2)
      u[!left] <- stats::qbeta(tail[!left], df / 2, .5)
      log_odds <- log(u) - log1p(-u)
      log_odds[!left] <- -log_odds[!left]
      exp((log(df) + log_odds) / 2)
    }
  }else if(half && identical(tau$distribution, "normal") && tau$parameters$mean == 0){
    tau <- tau$parameters$sd * sqrt(stats::qchisq(p, 1))
  }else if(!.is_prior_default_range(tau)){
    return(.regularized_truncated_scale_rule(tau, panels, global_power))
  }else tau <- quant(tau, p)
  if(any(!is.finite(tau) | tau <= 0)) .regularized_numerical_unavailable("scale integral")
  list(tau = tau, log_weight = global_log_weight)
}

# Craig's angular representation of the Normal upper tail, integrated
# analytically over the half-Cauchy local precision and the Gamma slab
# precision. At positive x put a=x/(sqrt(2)*sin(theta)*multiplier).
# The integrand is F(a/tau)*H(a), where F(r)=exp(r^2)*erfc(r) and
# H(a)=(1+2*a^2/(nu*slab_scale^2))^(-nu/2), or 1 for a classical
# horseshoe. Differentiating this tail gives the density with the same
# one-dimensional angular rule (two dimensions for stochastic tau).
.regularized_horseshoe_rule <- function(prior, panels, near_pole = FALSE){

  prob <- .regularized_probability_rule(panels)
  spec <- prior$regularization
  global <- .regularized_global_scale_rule(spec$scale, panels, prob, if(near_pole) 6 else 2)
  angular_power <- if(near_pole) 6 else 3
  # The derivative is concentrated near theta=0 for coefficients close
  # to the horseshoe pole. Cluster angular nodes there with theta=(pi/2)*p^3.
  log_angle <- -log(sinpi(prob$prob^angular_power / 2)) - log(2) / 2 - log(spec$multiplier)
  angle_log_weight <- log(angular_power / 2 * prob$weight) + (angular_power - 1) * log(prob$prob)
  list(log_r = as.vector(outer(log_angle, log(global$tau), "-")),
       log_a = log_angle,
       log_w = as.vector(outer(angle_log_weight, global$log_weight, "+")))
}

.regularized_horseshoe_log_value <- function(log_x, spec, rule, density){

  log_r <- log_x + rule$log_r
  small <- log_r <= log(8)
  log_f <- numeric(length(log_r))
  if(density) log_g <- log_f
  r <- exp(log_r[small])
  log_f[small] <- log(2) + r^2 + stats::pnorm(-sqrt(2) * r, log.p = TRUE)
  if(density){
    # G(r)=-r F'(r)=2r/sqrt(pi)-2r^2 F(r).
    log_g[small] <- log(2 / sqrt(pi)) + log_r[small] +
      log1p(-sqrt(pi) * r * exp(log_f[small]))
  }
  if(any(!small)){
    # Asymptotic expansion of F and its derivative; successive terms
    # decrease for r>8. The first omitted term is <1e-10 of F (or G).
    q <- -0.5 * exp(-2 * log_r[!small])
    term <- f <- g <- rep(1, length(q))
    for(k in seq_len(10L)){
      term <- term * (2 * k - 1) * q
      f <- f + term
      if(density) g <- g + (2 * k + 1) * term
    }
    first <- -log_r[!small] - log(pi) / 2
    log_f[!small] <- first + log(f)
    if(density) log_g[!small] <- first + log(g)
  }
  if(is.null(spec$slab_scale)){
    integrand <- if(density) log_g - log_x else log_f
  }else{
    log_u <- log(2) + 2 * (log_x + rule$log_a - log(spec$slab_scale)) - log(spec$slab_df)
    log_h <- -spec$slab_df / 2 * (pmax(log_u, 0) + log1p(exp(-abs(log_u))))
    if(density){
      log_slope <- log(spec$slab_df) + stats::plogis(log_u, log.p = TRUE)
      integrand <- log_h + .prior_normal_logaddexp(log_g, log_f + log_slope) - log_x
    }else{
      integrand <- log_h + log_f
    }
  }
  .regularized_log_sum(rule$log_w + integrand)
}

.regularized_log_sum <- function(x){

  peak <- max(x)
  if(!is.finite(peak)) return(peak)
  peak + log(sum(exp(x - peak)))
}

.regularized_lasso_log_kernel <- function(z, spec, density){

  if(!is.null(spec$group_lasso_dimension)){
    return(.regularized_group_lasso_log_kernel(z, spec$group_lasso_dimension, density))
  }
  ratio <- spec$contrast_ratio
  log_kernel <- -log(2) - z
  if(is.null(ratio)) return(log_kernel)
  # The density of L1 + r*L2 divides exp(-z) - r*exp(-z/r)
  # by 2*(1-r^2). Factor the larger exponential and use expm1 for
  # the difference; the limit of the divided difference at r = 1 is z.
  ok <- is.finite(z)
  u <- z[ok]
  change <- if(ratio == 1) u else
    -expm1(-u * (1 - ratio) / ratio) * (ratio / (1 - ratio))
  correction <- if(density) log1p(change) - log1p(ratio) else
    log1p(change * (ratio / (1 + ratio)))
  log_kernel[ok] <- log_kernel[ok] + correction
  log_kernel
}

# A learned Normal/Laplace width is a single positive scale integral. Put
# breakpoints at the transition |beta|/tau as well as at the scale prior's
# quantiles. The adaptive log-coordinate integrator can then resolve the
# decades between a near-zero coefficient and a typical global scale.
.regularized_global_integral <- function(prior, x, density){

  spec <- prior$regularization
  scale <- spec$scale
  width <- spec$multiplier
  log_x <- log(x)
  scale_lpdf <- .prior_simple_lpdf_evaluator(scale, purpose = "natural_integral")
  setup <- .prior_conditional_normal_breakpoint_setup(scale, unlist(scale$truncation, use.names = FALSE))
  normal <- list(multiplier = scale, bounds = unlist(scale$truncation, use.names = FALSE),
    additive_mean = 0, additive_sd = 0, product_mean = 0, product_sd = width)
  breaks <- lapply(x, function(value) .prior_conditional_normal_breakpoints(normal, value, setup = setup))
  integrand <- list(
    shared = function(tau) list(log_scale = log(tau) + log(width), log_prior = scale_lpdf(tau)),
    log_value = function(s, index){
      z <- exp(log_x[index] - s$log_scale)
      log_kernel <- if(spec$family == "lasso"){
        .regularized_lasso_log_kernel(z, spec, density)
      }else if(density) stats::dnorm(z, log = TRUE) else stats::pnorm(-z, log.p = TRUE)
      if(density) log_kernel <- log_kernel - s$log_scale
      log_kernel + s$log_prior
    })
  log(.prior_density_quadrature_batch(integrand, breaks))
}

# With classical half-Cauchy global and local scales their standardized
# product has log-scale density 2*v/(pi^2*sinh(v)). It follows by
# integrating the product of their two half-Cauchy densities. Splitting
# in this single coordinate at the conditional kernel's transition
# resolves even far-tail probabilities used to bracket a plot quantile.
.regularized_classical_integral <- function(x, log_scale, density){

  log_x <- log(x)
  breaks <- lapply(log_x - log_scale, function(transition){
    sort(unique(c(-Inf, 0, transition, Inf)))
  })
  kernel <- function(v, index){
    magnitude <- abs(v)
    log_mass <- rep(log(2 / pi^2), length(v))
    away <- magnitude > 0
    u <- magnitude[away]
    log_mass[away] <- log(4 / pi^2) + log(u) - u - log(-expm1(-2 * u))
    log_sd <- log_scale + v
    z <- exp(log_x[index] - log_sd)
    log_kernel <- if(density) stats::dnorm(z, log = TRUE) - log_sd else stats::pnorm(-z, log.p = TRUE)
    exp(log_mass + log_kernel)
  }
  log(.prior_density_quadrature_batch(kernel, breaks))
}

# Compare refinements on the log scale, so even small tail probabilities
# must agree to relative precision. A failed computation never supplies an
# apparent zero density or a finite ordinate at the pole.
.regularized_log_integral <- function(prior, x, operation){

  spec <- prior$regularization
  if(!is.null(spec$group_lasso_dimension) && is.numeric(spec$scale) && operation != "variance"){
    log_scale <- log(spec$multiplier) + log(spec$scale)
    result <- .regularized_group_lasso_log_kernel(exp(log(x) - log_scale), spec$group_lasso_dimension,
      density = operation == "density")
    if(any(!is.finite(result))) .regularized_numerical_unavailable(operation)
    return(if(operation == "density") result - log_scale else result)
  }
  result <- rep(NA_real_, length(x))
  pending <- seq_along(x)
  previous <- result
  horseshoe <- prior$regularization$family == "horseshoe" && operation != "variance"
  near_pole <- FALSE
  if(horseshoe){
    s <- prior$regularization$scale
    if(is.null(prior$regularization$slab_scale) && is.prior(s) && identical(s$distribution, "t") &&
       s$parameters$df == 1 && s$parameters$location == 0 &&
       s$truncation$lower == 0 && s$truncation$upper == Inf){
      log_scale <- log(s$parameters$scale) + log(prior$regularization$multiplier)
      result <- .regularized_classical_integral(x, log_scale, density = operation == "density")
      pending <- which(!is.finite(result))
      if(!length(pending)) return(result)
    }
    typical <- if(is.prior(s)) quant(s, .5) else s
    near_pole <- min(log(x)) - log(typical) - log(prior$regularization$multiplier) < -9
  }
  if(!horseshoe && operation != "variance" && is.prior(prior$regularization$scale)){
    result <- .regularized_global_integral(prior, x, density = operation == "density")
    pending <- which(!is.finite(result))
    if(!length(pending)) return(result)
  }
  panels_to_try <- if(near_pole) c(1L, 2L, 4L, 8L, 16L, 32L) else
    if(horseshoe) c(1L, 2L, 4L, 8L, 16L) else c(1L, 2L, 4L, 8L)
  for(panels in panels_to_try){
    rule <- if(horseshoe) .regularized_horseshoe_rule(prior, panels, near_pole) else
      .regularized_scale_rule(prior, panels, variance = operation == "variance")
    estimate <- vapply(x[pending], function(value){
      if(horseshoe) return(.regularized_horseshoe_log_value(log(value), prior$regularization,
        rule, density = operation == "density"))
      if(operation == "variance") return(.regularized_log_sum(rule$log_weight + 2 * rule$log_scale))
      z <- exp(log(value) - rule$log_scale)
      kernel <- if(operation == "density"){
        if(prior$regularization$family == "lasso") .regularized_lasso_log_kernel(
          z, prior$regularization, TRUE) - rule$log_scale else
          stats::dnorm(z, log = TRUE) - rule$log_scale
      }else{
        if(prior$regularization$family == "lasso") .regularized_lasso_log_kernel(
          z, prior$regularization, FALSE) else
          stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
      }
      .regularized_log_sum(kernel + rule$log_weight)
    }, numeric(1))
    converged <- is.finite(estimate) & is.finite(previous[pending]) &
      abs(estimate - previous[pending]) <= 1e-6
    result[pending[converged]] <- estimate[converged]
    previous[pending] <- estimate
    pending <- pending[!converged]
    if(!length(pending)) return(result)
  }
  .regularized_numerical_unavailable(operation)
}

.regularized_marginal <- function(prior, x, density, log, lower.tail = TRUE){

  out <- x + 0.0
  finite <- is.finite(x)
  out[is.infinite(x)] <- -Inf
  zero <- finite & x == 0
  if(any(zero)){
    out[zero] <- if(density) .regularized_zero_log_density(prior) else -log(2)
    if(anyNA(out[zero])) .regularized_numerical_unavailable("density at zero")
  }
  inside <- finite & x != 0
  if(any(inside)){
    values <- unique(abs(x[inside]))
    logs <- .regularized_log_integral(prior, values, if(density) "density" else "tail probability")
    out[inside] <- logs[match(abs(x[inside]), values)]
  }
  if(!density){
    complement <- !is.na(x) & if(lower.tail) x > 0 else x < 0
    out[complement] <- log1p(-exp(out[complement]))
  }
  if(log) out else exp(out)
}

.regularized_quantile <- function(prior, p){

  vapply(p, function(prob){
    if(is.na(prob)) return(prob)
    if(prob == 0 || prob == 1) return(if(prob == 0) -Inf else Inf)
    if(prob == 0.5) return(0)
    log_tail <- log(min(prob, 1 - prob))
    difference <- function(log_x){
      .regularized_marginal(prior, exp(log_x), FALSE, TRUE, lower.tail = FALSE) - log_tail
    }
    scale <- prior$regularization$scale
    center <- log(if(is.numeric(scale)) scale else quant(scale, 0.5)) + log(prior$regularization$multiplier)
    root <- stats::uniroot(difference, center + c(-1, 1), extendInt = "downX", tol = 1e-8)
    sign(prob - 0.5) * exp(root$root)
  }, numeric(1), USE.NAMES = FALSE)
}

.regularized_moment <- function(prior, variance){

  spec <- prior$regularization
  if(spec$family == "horseshoe"){
    # The half-Cauchy local scale has no first moment. As for a Cauchy
    # coefficient, the mean and its centered variance are undefined.
    if(is.null(spec$slab_scale)) return(NaN)
    if(spec$slab_df <= if(variance) 2 else 1) .regularized_numerical_unavailable("moment")
    if(!variance) return(0)
    return(exp(.regularized_log_integral(prior, 1, "variance")))
  }
  scale_mean <- if(is.numeric(spec$scale)) spec$scale else mean(spec$scale)
  if(!variance) return(if(is.finite(scale_mean)) 0 else NaN)
  second <- if(is.numeric(spec$scale)) spec$scale^2 else var(spec$scale) + scale_mean^2
  variance <- spec$multiplier^2 * second * if(spec$family == "lasso") 2 else 1
  if(!is.null(spec$group_lasso_dimension)) variance <- spec$multiplier^2 * second * (spec$group_lasso_dimension + 1)
  if(is.null(spec$contrast_ratio)) variance else variance * (1 + spec$contrast_ratio^2)
}
