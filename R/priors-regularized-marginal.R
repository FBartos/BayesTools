# Marginal distribution of one fitted coefficient. Conditional on the scales
# it is Normal (ridge/horseshoe) or Laplace. Integrate over the declared
# distributions of the scales, independently of posterior or prior draws.
.regularized_zero_log_density <- function(prior){

  spec <- prior$regularization
  if(spec$family == "horseshoe") return(Inf)
  kernel <- if(spec$family == "ridge") stats::dnorm(0, log = TRUE) else -log(2)
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
  if(is.prior(spec$scale)){
    global <- quant(spec$scale, prob)
    w <- weight
  }else{
    global <- spec$scale
    w <- 1
  }
  if(any(!is.finite(global) | global <= 0)) .regularized_numerical_unavailable("scale integral")
  log_scale <- log(global)
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
    w <- as.vector(outer(w, weight))
  }
  list(log_scale = log_scale + log(spec$multiplier), log_weight = log(w))
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
  global_power <- if(near_pole) 6 else 2
  angular_power <- if(near_pole) 6 else 3
  p <- prob$prob^global_power
  tau <- spec$scale
  if(is.prior(tau)){
    # Inverting a symmetric untruncated CDF at (1+p)/2 rounds tiny
    # probabilities onto zero. The positive-half quantiles avoid that loss.
    half <- tau$truncation$lower == 0 && tau$truncation$upper == Inf
    if(half && identical(tau$distribution, "t") && tau$parameters$location == 0){
      df <- tau$parameters$df
      tau <- tau$parameters$scale * if(df == 1){
        slope <- tanpi(pmin(p, 1 - p) / 2)
        ifelse(p <= .5, slope, 1 / slope)
      }else sqrt(stats::qf(p, 1, df))
    }else if(half && identical(tau$distribution, "normal") && tau$parameters$mean == 0){
      tau <- tau$parameters$sd * sqrt(stats::qchisq(p, 1))
    }else tau <- quant(tau, p)
  }
  if(any(!is.finite(tau) | tau <= 0)) .regularized_numerical_unavailable("scale integral")
  # The derivative is concentrated near theta=0 for coefficients close
  # to the horseshoe pole. Cluster angular nodes there with theta=(pi/2)*p^3.
  log_angle <- -log(sinpi(prob$prob^angular_power / 2)) - log(2) / 2 - log(spec$multiplier)
  global_weight <- if(is.prior(spec$scale)) global_power * prob$prob^(global_power - 1) * prob$weight else 1
  list(log_r = as.vector(outer(log_angle, log(tau), "-")),
       log_a = log_angle,
       log_w = log(as.vector(outer(angular_power / 2 * prob$weight * prob$prob^(angular_power - 1), global_weight))))
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
        -log(2) - z
      }else if(density) stats::dnorm(z, log = TRUE) else stats::pnorm(-z, log.p = TRUE)
      if(density) log_kernel <- log_kernel - s$log_scale
      log_kernel + s$log_prior
    })
  log(.prior_density_quadrature_batch(integrand, breaks))
}

# Compare refinements on the log scale, so even small tail probabilities
# must agree to relative precision. A failed computation never supplies an
# apparent zero density or a finite ordinate at the pole.
.regularized_log_integral <- function(prior, x, operation){

  result <- rep(NA_real_, length(x))
  pending <- seq_along(x)
  previous <- result
  horseshoe <- prior$regularization$family == "horseshoe" && operation != "variance"
  near_pole <- FALSE
  if(horseshoe){
    s <- prior$regularization$scale
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
        if(prior$regularization$family == "lasso") -log(2) - z - rule$log_scale else
          stats::dnorm(z, log = TRUE) - rule$log_scale
      }else{
        if(prior$regularization$family == "lasso") -log(2) - z else
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
    if(is.null(spec$slab_scale)) return(if(variance) Inf else NaN)
    if(spec$slab_df <= if(variance) 2 else 1) .regularized_numerical_unavailable("moment")
    if(!variance) return(0)
    return(exp(.regularized_log_integral(prior, 1, "variance")))
  }
  scale_mean <- if(is.numeric(spec$scale)) spec$scale else mean(spec$scale)
  if(!variance) return(if(is.finite(scale_mean)) 0 else NaN)
  second <- if(is.numeric(spec$scale)) spec$scale^2 else var(spec$scale) + scale_mean^2
  spec$multiplier^2 * second * if(spec$family == "lasso") 2 else 1
}
