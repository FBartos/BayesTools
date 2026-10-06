.prior_numerical_condition <- function(operation, family, requested_scale, indices,
                                       reason, range = FALSE, error = FALSE,
                                       rng = FALSE){

  structure(
    list(
      message = paste0("The ", family, " prior ", operation, " is numerically ",
                       if(range) "outside the requested representable range" else "unavailable",
                       ". ", reason, "."),
      call = NULL,
      operation = operation,
      family = family,
      requested_scale = requested_scale,
      indices = indices,
      reason = reason
    ),
    class = c(if(rng) "BayesTools_prior_rng_unavailable",
              if(range) "BayesTools_numerical_range_limit" else "BayesTools_numerical_unavailable",
              "BayesTools_numerical_condition", if(error) "error" else "warning", "condition")
  )
}

.prior_numerical_signal <- function(operation, family, requested_scale, indices,
                                    reason, range = FALSE, error = FALSE,
                                    rng = FALSE){

  condition <- .prior_numerical_condition(operation, family, requested_scale,
                                         indices, reason, range, error, rng)
  if(error) stop(condition)
  warning(condition)
  invisible(NULL)
}

.prior_numerical_result <- function(output, input, domain, operation, family,
                                    requested_scale, interior = rep(TRUE, length(output)),
                                    rng = FALSE, bounded = FALSE, warn_range = TRUE){

  if(!isTRUE(domain)) return(output)
  known <- !is.na(input) & interior
  unavailable <- which(known & is.nan(output))
  if(length(unavailable)){
    .prior_numerical_signal(operation, family, requested_scale, unavailable,
      "The declared interior calculation could not be resolved at supported precision",
      error = rng || bounded, rng = rng)
    return(output)
  }
  if(!warn_range) return(output)
  range <- which(known & (is.infinite(output) |
    (output == 0 & requested_scale == "natural" &
     (operation == "density" |
      (family == "invgamma" & operation %in% c("quantile", "sampling"))))))
  if(length(range)){
    .prior_numerical_signal(operation, family, requested_scale, range,
      "Use an available logarithmic result or inspect the declared numerical limit",
      range = TRUE, error = bounded, rng = rng && bounded)
  }
  output
}

.prior_numerical_nonlocal_domain <- function(location, tau, order, df = NULL){

  is.finite(location) && is.finite(tau) && tau > 0 &&
    is.finite(order) && order >= 1 && order == floor(order) &&
    (is.null(df) || (is.finite(df) && df > 0))
}

.prior_numerical_finite <- function(values, family, operation = "sampling"){

  if(!is.character(family) || length(family) != 1L) family <- "declared"
  failed <- which(!is.finite(values) | (family == "invgamma" & values <= 0))
  if(length(failed)){
    .prior_numerical_signal(operation, family, "finite", failed,
      "This consumer requires finite representable values",
      error = TRUE, rng = identical(operation, "sampling"))
  }
  values
}

.prior_numerical_without_warnings <- function(expression){

  withCallingHandlers(expression,
    BayesTools_numerical_condition = function(condition){
      if(inherits(condition, "warning")) invokeRestart("muffleWarning")
    })
}

.prior_nonlocal_log_interval_mass <- function(prior, lower, upper){

  .BayesTools_require_native_nonlocal()
  parameters <- prior$parameters
  .Call("BayesTools_nonlocal_log_interval_mass", lower, upper,
        parameters$location, parameters$tau, as.numeric(parameters$order),
        if(prior$distribution == "invmoment") parameters$df else 0,
        prior$distribution == "invmoment", PACKAGE = "BayesTools")
}

.prior_nonlocal_truncated_quantile <- function(prior, p){

  .BayesTools_require_native_nonlocal()
  parameters <- prior$parameters
  out <- .Call("BayesTools_nonlocal_truncated_quantile", p,
        prior$truncation$lower, prior$truncation$upper,
        parameters$location, parameters$tau, as.numeric(parameters$order),
        if(prior$distribution == "invmoment") parameters$df else 0,
        prior$distribution == "invmoment", PACKAGE = "BayesTools")
  log_mass <- attr(out, "log_normalizer", exact = TRUE)
  attr(out, "log_normalizer") <- NULL
  if(any(p > 0 & p < 1, na.rm = TRUE) && !is.finite(log_mass)){
    .prior_numerical_signal("normalization", prior$distribution, "log", 1L,
      "The positive truncation mass could not be resolved at supported precision", error = TRUE)
  }
  .prior_numerical_result(out, p, TRUE, "quantile", prior$distribution, "natural",
                          interior = p > 0 & p < 1)
}

.prior_log1mexp <- function(x){

  out <- rep(NaN, length(x))
  valid <- !is.na(x) & x <= 0
  far <- valid & x < log(.5)
  out[far] <- log1p(-exp(x[far]))
  out[valid & !far] <- log(-expm1(x[valid & !far]))
  out[is.na(x)] <- x[is.na(x)]
  out
}

.prior_logdiffexp <- function(x, y){

  if(length(x) != length(y)) stop("Log probability pairs must have equal lengths.", call. = FALSE)
  out <- rep(NaN, length(x))
  valid <- !is.na(x) & !is.na(y) & x >= y
  zero <- valid & x == y
  out[zero] <- -Inf
  out[valid & !zero] <- x[valid & !zero] + .prior_log1mexp(y[valid & !zero] - x[valid & !zero])
  out[is.na(y)] <- y[is.na(y)]
  out[is.na(x)] <- x[is.na(x)]
  out
}

.prior_simple_log_interval_mass <- function(prior, lower, upper){

  if(length(lower) == 1L) lower <- rep(lower, length(upper))
  if(length(upper) == 1L) upper <- rep(upper, length(lower))
  if(length(lower) != length(upper)) stop("Interval bounds must have compatible lengths.", call. = FALSE)
  if(prior$distribution == "normal") return(.prior_normal_log_interval_mass(prior, lower, upper))
  if(prior$distribution %in% c("moment", "invmoment")) return(.prior_nonlocal_log_interval_mass(prior, lower, upper))
  tails <- .prior_numerical_without_warnings(list(
    .prior_simple_base_p(prior, lower, log.p = TRUE),
    .prior_simple_base_p(prior, upper, log.p = TRUE),
    .prior_simple_base_p(prior, lower, lower.tail = FALSE, log.p = TRUE),
    .prior_simple_base_p(prior, upper, lower.tail = FALSE, log.p = TRUE)))
  log_pl <- tails[[1L]]
  log_pu <- tails[[2L]]
  log_ql <- tails[[3L]]
  log_qu <- tails[[4L]]
  from_p <- .prior_logdiffexp(log_pu, log_pl)
  from_q <- .prior_logdiffexp(log_ql, log_qu)
  out <- .prior_log1mexp(.prior_normal_logaddexp(log_pl, log_qu))
  prefer_p <- !is.na(log_pu) & log_pu <= log(.5)
  prefer_q <- !prefer_p & !is.na(log_ql) & log_ql <= log(.5)
  out[prefer_p] <- from_p[prefer_p]
  out[prefer_q] <- from_q[prefer_q]
  failed <- !is.finite(out)
  out[failed & is.finite(from_p)] <- from_p[failed & is.finite(from_p)]
  failed <- !is.finite(out)
  out[failed & is.finite(from_q)] <- from_q[failed & is.finite(from_q)]
  out[!is.na(lower) & !is.na(upper) & lower == upper] <- -Inf
  unresolved <- !is.na(lower) & !is.na(upper) & lower < upper & !is.finite(out)
  out[unresolved] <- NaN
  out[is.na(upper)] <- upper[is.na(upper)]
  out[is.na(lower)] <- lower[is.na(lower)]
  out
}

.prior_simple_log_C <- function(prior){

  shape <- prior$parameters$shape
  if(prior$distribution == "gamma" && is.numeric(shape) && length(shape) == 1L &&
     is.finite(shape) && shape > 0 && shape < .Machine$double.xmin){
    if(prior$truncation$lower <= 0 && prior$truncation$upper == Inf) return(0)
    .prior_numerical_signal("normalization", "gamma", "log", 1L,
      "Ordinary Gamma truncation tails are unavailable in the subnormal-shape floating-point regime",
      error = TRUE)
  }
  log_mass <- .prior_simple_log_interval_mass(prior, prior$truncation$lower, prior$truncation$upper)
  if(length(log_mass) != 1L || !is.finite(log_mass)){
    .prior_numerical_signal("normalization", prior$distribution, "log", 1L,
      "The positive truncation mass could not be resolved at supported precision", error = TRUE)
  }
  log_mass
}

.prior_simple_truncated_probability <- function(prior, q, lower.tail){

  lower <- prior$truncation$lower
  upper <- prior$truncation$upper
  out <- rep(NaN, length(q))
  out[is.na(q)] <- q[is.na(q)]
  known <- !is.na(q)
  out[known & q <= lower] <- if(lower.tail) 0 else 1
  out[known & q >= upper] <- if(lower.tail) 1 else 0
  interior <- known & q > lower & q < upper
  if(any(interior)){
    log_mass <- .prior_simple_log_C(prior)
    interval <- if(lower.tail) .prior_simple_log_interval_mass(prior, lower, q[interior]) else
      .prior_simple_log_interval_mass(prior, q[interior], upper)
    out[interior] <- exp(interval - log_mass)
    out[interior & !is.na(out) & (out < 0 | out > 1)] <- NaN
  }
  .prior_numerical_result(out, q, TRUE, "distribution", prior$distribution, "natural", interior)
}

.prior_simple_truncated_quantile <- function(prior, p){

  if(prior$distribution %in% c("moment", "invmoment")) return(.prior_nonlocal_truncated_quantile(prior, p))
  lower <- prior$truncation$lower
  upper <- prior$truncation$upper
  out <- rep(NaN, length(p))
  out[is.na(p)] <- p[is.na(p)]
  known <- !is.na(p)
  out[known & p == 0] <- lower
  out[known & p == 1] <- upper
  interior <- known & p > 0 & p < 1
  if(any(interior)){
    log_mass <- .prior_simple_log_C(prior)
    log_p <- .prior_normal_logaddexp(.prior_numerical_without_warnings(.prior_simple_base_p(prior, lower, log.p = TRUE)),
                                    log(p[interior]) + log_mass)
    log_q <- .prior_normal_logaddexp(.prior_numerical_without_warnings(.prior_simple_base_p(prior, upper, lower.tail = FALSE, log.p = TRUE)),
                                    log1p(-p[interior]) + log_mass)
    use_p <- !is.na(log_p) & (is.na(log_q) | log_p <= log_q)
    values <- rep(NaN, sum(interior))
    if(any(use_p)) values[use_p] <- .prior_numerical_without_warnings(.prior_simple_base_q(prior, log_p[use_p], log.p = TRUE))
    if(any(!use_p)) values[!use_p] <- .prior_numerical_without_warnings(.prior_simple_base_q(prior, log_q[!use_p], lower.tail = FALSE, log.p = TRUE))
    values[!is.finite(values) | values <= lower | values >= upper] <- NaN
    out[interior] <- values
  }
  .prior_numerical_result(out, p, TRUE, "quantile", prior$distribution, "natural", interior)
}

.prior_numerical_ordinate <- function(expression, value, method, provenance){

  tryCatch(expression,
    BayesTools_numerical_unavailable = function(condition){
      .prior_density_ordinate_imprecise(value,
        paste0("The declared numerical calculation (", condition$reason, ")"),
        method, c(provenance, list(numerical_condition = unclass(condition))))
    })
}
