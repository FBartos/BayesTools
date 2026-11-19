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
.regularized_scale_rule <- function(prior, panels, variance = FALSE){

  rule <- .prior_density_quadrature_rules()$finite
  t <- as.vector(outer((1 + rule$nodes) / 2, seq_len(panels) - 1, "+")) / panels
  prob <- sin(pi * t / 2)^2
  weight <- rep(rule$kronrod, panels) * pi * sinpi(t) / (4 * panels)
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
  if(spec$family == "horseshoe" && !variance){
    log_local <- log(stats::qcauchy((1 - prob) / 2, lower.tail = FALSE))
    log_scale <- as.vector(outer(log_scale, log_local, "+"))
    w <- as.vector(outer(w, weight))
  }
  if(!is.null(spec$slab_scale)){
    # c² = nu * slab_scale² / chi-square(nu)
    log_c <- log(spec$slab_scale) + (log(spec$slab_df) -
      log(stats::qchisq(prob, df = spec$slab_df))) / 2
    a <- rep(log_scale, length(log_c))
    b <- rep(log_c, each = length(log_scale))
    if(variance){
      # E_lambda[sd² | tau,c] = c² * tau/(tau+c) under half-Cauchy(1).
      log_scale <- b + stats::plogis(a - b, log.p = TRUE) / 2
    }else{
      log_scale <- pmin(a, b) - log1p(exp(-2 * abs(a - b))) / 2
    }
    w <- as.vector(outer(w, weight))
  }
  list(log_scale = log_scale + log(spec$multiplier), log_weight = log(w))
}

.regularized_log_sum <- function(x){

  peak <- max(x)
  if(!is.finite(peak)) return(peak)
  peak + log(sum(exp(x - peak)))
}

# Compare refinements on the log scale, so even small tail probabilities
# must agree to relative precision. A failed computation never supplies an
# apparent zero density or a finite ordinate at the pole.
.regularized_log_integral <- function(prior, x, operation){

  result <- rep(NA_real_, length(x))
  pending <- seq_along(x)
  previous <- result
  for(panels in c(1L, 2L, 4L, 8L)){
    rule <- .regularized_scale_rule(prior, panels, variance = operation == "variance")
    estimate <- vapply(x[pending], function(value){
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
