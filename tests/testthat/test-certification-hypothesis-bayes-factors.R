skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

# hypothesis_BF() on a conjugate model. With a known residual SD of 1 and
# N(0, 1) priors on the intercept, the slope, and the two mean-difference
# coefficients of a three-level factor, the coefficients have the posterior
# N(m, V) with V = (X'X + I)^-1 and m = V X'y, so every point and region Bayes
# factor of the slope or of a linear combination of factor levels is exact.
# Tolerances are four Monte Carlo SDs of each estimator at the fit's
# effective sample size (ESS): 1 / sqrt(ESS p (1 - p)) for the log posterior
# odds of a region; sqrt(z^2 + (z^2 - 1)^2 / 2) / sqrt(ESS) for the log
# normal-approximation ordinate at z posterior SDs from the null; and for the
# kernel ordinate in addition sqrt(R(K) / (ESS h f)) with the Gaussian kernel
# (R(K) = 1 / (2 sqrt(pi))), its normal-reference bandwidth h, and its
# smoothing bias h^2 (z^2 - 1) / (2 s^2). The reported error%(BF) of each
# region test bounds the observed relative error within four times.

hypothesis_data <- local({
  set.seed(1)
  data <- data.frame(x1 = stats::rnorm(45), f = factor(rep(c("a", "b", "c"), length.out = 45)))
  data$y <- 0.2 + 0.25 * data$x1 + c(-0.3, 0.1, 0.2)[as.integer(data$f)] + stats::rnorm(45)
  data
})
hypothesis_contrast <- contr.meandif(levels(hypothesis_data$f))
hypothesis_posterior <- local({
  X <- cbind(1, hypothesis_data$x1, hypothesis_contrast[as.integer(hypothesis_data$f), ])
  V <- solve(crossprod(X) + diag(ncol(X)))
  list(mean = drop(V %*% crossprod(X, hypothesis_data$y)), cov = V)
})

# the coefficient weights of the slope and of the level effects
slope_weights <- c(0, 1, 0, 0)
level_weights <- function(level) c(0, 0, hypothesis_contrast[level, ])

# the posterior and prior of a linear combination of the coefficients
linear_combination <- function(w) {
  c(mean     = sum(w * hypothesis_posterior$mean),
    sd       = sqrt(drop(t(w) %*% hypothesis_posterior$cov %*% w)),
    prior_sd = sqrt(sum(w[-1L]^2)))
}
exact_point <- function(w) {
  l <- linear_combination(w)
  stats::dnorm(0, 0, l[["prior_sd"]], log = TRUE) - stats::dnorm(0, l[["mean"]], l[["sd"]], log = TRUE)
}
exact_region <- function(w, lower, upper) {
  l <- linear_combination(w)
  p <- stats::pnorm(upper, l[["mean"]], l[["sd"]]) - stats::pnorm(lower, l[["mean"]], l[["sd"]])
  q <- stats::pnorm(upper, 0, l[["prior_sd"]]) - stats::pnorm(lower, 0, l[["prior_sd"]])
  c(log_BF = stats::qlogis(p) - stats::qlogis(q), p = p)
}

test_that("hypothesis_BF() matches the exact Bayes factors of a conjugate model", {

  unit <- prior("normal", list(mean = 0, sd = 1))
  fit  <- JAGS_fit(
    "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) } }",
    list(y = hypothesis_data$y, N = nrow(hypothesis_data)),
    formula_list = list(mu = ~ x1 + f), formula_data_list = list(mu = hypothesis_data),
    formula_prior_list = list(mu = list(intercept = unit, x1 = unit,
                                        f = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif"))),
    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
  )
  n_draws <- 40000
  ess <- min(as.data.frame(JAGS_estimates_table(fit))[["ESS"]])

  regions <- list(
    "x1 > 0"              = exact_region(slope_weights, 0, Inf),
    "x1 > 0.1 & x1 < 0.4" = exact_region(slope_weights, 0.1, 0.4),
    "f[a] < f[c]"         = exact_region(level_weights(1) - level_weights(3), -Inf, 0),
    "f[a] + f[b] > 0"     = exact_region(level_weights(1) + level_weights(2), 0, Inf)
  )
  points <- list(
    "x1 = 0"      = slope_weights,
    "f[b] = f[c]" = level_weights(2) - level_weights(3),
    "f[b] = 0"    = level_weights(2)
  )

  region_table <- as.data.frame(hypothesis_BF(fit, names(regions), seed = 1))
  region_exact <- vapply(regions, `[[`, numeric(1), "log_BF")
  region_p     <- vapply(regions, `[[`, numeric(1), "p")
  region_error <- log(region_table[["BF10"]]) - region_exact
  expect_true(all(abs(region_error) < 4 / sqrt(ess * region_p * (1 - region_p))))
  expect_true(all(abs(expm1(region_error)) < 4 * region_table[["BF_error"]] / 100))

  ordinate_sd <- vapply(points, function(w) {
    l <- linear_combination(w)
    z <- l[["mean"]] / l[["sd"]]
    sqrt(z^2 + (z^2 - 1)^2 / 2) / sqrt(ess)
  }, numeric(1))
  kernel_allowance <- vapply(points, function(w) {
    l <- linear_combination(w)
    z <- l[["mean"]] / l[["sd"]]
    h <- 0.9 * l[["sd"]] * n_draws^(-1 / 5)
    f <- stats::dnorm(0, l[["mean"]], l[["sd"]])
    4 * sqrt(1 / (2 * sqrt(pi)) / (ess * h * f)) + abs(h^2 * (z^2 - 1) / (2 * l[["sd"]]^2))
  }, numeric(1))
  point_exact <- vapply(points, exact_point, numeric(1))

  normal_table <- as.data.frame(hypothesis_BF(fit, names(points), seed = 1, density_method = "normal"))
  expect_true(all(abs(log(normal_table[["BF10"]]) - point_exact) < 4 * ordinate_sd))

  kernel_table <- as.data.frame(hypothesis_BF(fit, names(points), seed = 1))
  expect_true(all(abs(log(kernel_table[["BF10"]]) - point_exact) < 4 * ordinate_sd + kernel_allowance))
})
