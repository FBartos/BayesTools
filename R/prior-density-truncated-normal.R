# Closed-form density of a Gaussian convolution with one truncated normal term.
#
# A Gaussian convolution X = a_m + w T + G (Gaussian terms and points with
# mean a_m, G ~ N(0, a_s), a_s > 0, plus one other continuous scalar term T
# with weight w; .prior_density_ordinate_gaussian_convolution_spec()) whose
# other term is a normal N(m, s) truncated to [l, u] (at least one finite
# bound) has a closed-form density. With T' = w T, a normal N(m_1, s_1),
# m_1 = w m, s_1 = |w| s, truncated to the image [l', u'] of [l, u], the
# product of the Gaussian densities
#   phi(t; m_1, s_1) phi(x - a_m - t; 0, a_s) = phi(x; a_m + m_1, s) phi(t; mu_x, tau),
# s^2 = s_1^2 + a_s^2, mu_x = m_1 + s_1^2 / s^2 (x - a_m - m_1),
# tau = s_1 a_s / s (T' given X = x is normal before truncation), integrates
# over [l', u'] to
#   f(x) = phi(x; a_m + m_1, s) [Phi(beta(x)) - Phi(alpha(x))] / [Phi(B) - Phi(A)],
# alpha(x) = (l' - mu_x) / tau, beta(x) = (u' - mu_x) / tau,
# A = (l - m) / s, B = (u - m) / s (the normalizer of T is invariant to w).
# Both interval masses are evaluated in log space from upper- or lower-tail
# probabilities (.prior_normal_log_interval_mass()), so tail values keep their
# relative precision. The density is positive and finite everywhere (a
# positive-variance Gaussian convolution): the ordinate is regular and exact at
# every value. Region probabilities keep the conditional-normal quadrature of
# the Gaussian convolution (they need the bivariate normal distribution
# function).

# Whether a Gaussian convolution specification has a truncated normal other
# term (the closed form applies).
.prior_truncated_normal_convolution_eligible <- function(spec){

  multiplier <- spec$multiplier
  isTRUE(spec$product_sd == 0) && isTRUE(spec$additive_sd > 0) &&
    is.finite(spec$additive_sd) && is.finite(spec$additive_mean) &&
    is.numeric(spec$product_mean) && length(spec$product_mean) == 1L &&
    is.finite(spec$product_mean) && spec$product_mean != 0 &&
    is.prior(multiplier) && is.prior.simple(multiplier) &&
    identical(multiplier$distribution, "normal") &&
    .prior_density_ordinate_parameters_numeric(multiplier) &&
    is.numeric(spec$bounds) && length(spec$bounds) == 2L && !anyNA(spec$bounds) &&
    spec$bounds[1L] < spec$bounds[2L] && any(is.finite(spec$bounds))
}

# Value-independent parameters of the closed form.
.prior_truncated_normal_convolution_parameters <- function(spec){

  weight <- spec$product_mean
  term_mean <- weight * spec$multiplier$parameters$mean
  term_sd <- abs(weight) * spec$multiplier$parameters$sd
  sd <- .prior_density_ordinate_stable_norm(c(term_sd, spec$additive_sd))
  list(
    term_mean      = term_mean,
    term_sd        = term_sd,
    term_bounds    = sort(weight * spec$bounds),
    mean           = spec$additive_mean + term_mean,
    sd             = sd,
    shrinkage      = (term_sd / sd)^2,
    conditional_sd = term_sd * (spec$additive_sd / sd),
    log_normalizer = .prior_normal_log_C(spec$multiplier)
  )
}

# Log density at the values 'x' (vectorized).
.prior_truncated_normal_convolution_log_density <- function(spec, x,
                                                            parameters = .prior_truncated_normal_convolution_parameters(spec)){

  conditional_mean <- parameters$term_mean +
    parameters$shrinkage * (x - spec$additive_mean - parameters$term_mean)
  log_mass <- .prior_normal_log_interval_mass(
    prior("normal", list(mean = 0, sd = 1)),
    (parameters$term_bounds[1L] - conditional_mean) / parameters$conditional_sd,
    (parameters$term_bounds[2L] - conditional_mean) / parameters$conditional_sd
  )
  stats::dnorm(x, parameters$mean, parameters$sd, log = TRUE) + log_mass -
    parameters$log_normalizer
}

# Structural provenance of the closed form (value independent).
.prior_truncated_normal_convolution_provenance <- function(spec){

  parameters <- .prior_truncated_normal_convolution_parameters(spec)
  list(
    kind                  = "truncated_normal_convolution",
    additive              = c(mean = spec$additive_mean, sd = spec$additive_sd),
    weight                = unname(spec$product_mean),
    truncated_normal      = .prior_density_ordinate_prior_provenance(spec$multiplier),
    convolution           = c(mean = parameters$mean, sd = parameters$sd,
                              conditional_sd = parameters$conditional_sd),
    independent_sources   = spec$sources,
    structural_regularity = "positive_variance_gaussian_convolution"
  )
}

.prior_truncated_normal_convolution_ordinate <- function(spec, value){

  log_density <- .prior_truncated_normal_convolution_log_density(spec, value)
  .prior_density_ordinate_result(
    value       = value,
    behavior    = "regular",
    log_density = log_density,
    exact       = TRUE,
    method      = "truncated_normal_convolution",
    reason      = if(identical(log_density, -Inf)){
      paste0(
        "The continuous prior density is structurally regular, but its ",
        "log-density underflows in ordinary floating-point arithmetic."
      )
    },
    provenance  = .prior_truncated_normal_convolution_provenance(spec)
  )
}

# Continuous density at the values 'x' (plotted densities).
.prior_truncated_normal_convolution_density <- function(spec, x){

  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if(any(finite)){
    out[finite] <- exp(.prior_truncated_normal_convolution_log_density(spec, x[finite]))
  }
  out
}
