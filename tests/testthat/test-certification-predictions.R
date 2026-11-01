skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

# JAGS_predict_formula() against exact Gaussian predictions. With a known
# residual SD of 1, N(0, 1) priors on the intercept and the slope, and random
# group intercepts with a known SD (a point prior at 0.6), the coefficients and
# the group intercepts theta = (beta, u) have the posterior N(m, S) with
# S = (A'A + P^-1)^-1, m = S A'y, A = [X Z], and P = diag(1, 1, 0.6^2, ...).
# A prediction row x (with group indicator z) then has the exact posterior
# N(x'm, x'S x) for the fixed target, N((x, z)'m, (x, z)'S (x, z)) for the
# conditional target of a fitted group, and N(x'm_beta, x'S_beta x + 0.6^2)
# for a new group under new_levels = "sample"; the marginal target with the
# covariance method has the fixed mean and the random-effect covariance
# 0.6^2 z z'. Sampled means and SDs are held to four Monte Carlo SDs at the
# effective sample size of each prediction's draws: SD / sqrt(ESS) for a mean
# and SD / sqrt(2 ESS) for an SD (the predictions are normal). Every ESS
# exceeds 1000, so no bound exceeds 0.13 of a prediction SD.

prediction_data <- local({
  set.seed(4)
  data <- data.frame(x = stats::rnorm(30), g = factor(rep(letters[1:6], each = 5)))
  data$y <- 0.2 + 0.5 * data$x + stats::rnorm(6, 0, 0.6)[as.integer(data$g)] + stats::rnorm(30)
  data
})
prediction_sd <- 0.6
prediction_exact <- local({
  A <- cbind(1, prediction_data$x, stats::model.matrix(~ g - 1, prediction_data))
  P <- diag(c(1, 1, rep(prediction_sd^2, nlevels(prediction_data$g))))
  S <- solve(crossprod(A) + solve(P))
  list(mean = drop(S %*% crossprod(A, prediction_data$y)), cov = S)
})
# prediction rows: two fitted groups and one new group "z"
prediction_new <- data.frame(x = c(-1, 0.5, 2), g = factor(c("a", "d", "z"), levels = c(letters[1:6], "z")))

prediction_moments <- function(weights, extra_variance = 0) {
  c(mean = sum(weights * prediction_exact$mean),
    sd   = sqrt(drop(t(weights) %*% prediction_exact$cov %*% weights) + extra_variance))
}
prediction_weights <- function(x, group = NULL) {
  u <- numeric(nlevels(prediction_data$g))
  if (!is.null(group)) u[match(group, levels(prediction_data$g))] <- 1
  c(1, x, u)
}

test_that("formula predictions match the exact Gaussian posterior predictions", {

  unit <- prior("normal", list(mean = 0, sd = 1))
  fit  <- JAGS_fit(
    "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) } }",
    list(y = prediction_data$y, N = nrow(prediction_data)),
    formula_list = list(mu = ~ x + (1 | g)), formula_data_list = list(mu = prediction_data),
    formula_prior_list = list(mu = list(intercept = unit, x = unit)),
    formula_random_prior_list = list(mu = prior_random(g = random_block(sd = prior("point", list(location = prediction_sd))))),
    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
  )
  expect_moments <- function(draws, exact) {
    ess <- apply(draws, 1L, function(row) coda::effectiveSize(coda::mcmc(row)))
    expect_true(all(is.finite(ess) & ess > 1000))
    observed <- cbind(mean = rowMeans(draws), sd = apply(draws, 1L, stats::sd))
    expect_true(all(abs(observed[, "mean"] - exact[, "mean"]) < 4 * exact[, "sd"] / sqrt(ess)))
    expect_true(all(abs(observed[, "sd"] - exact[, "sd"]) < 4 * exact[, "sd"] / sqrt(2 * ess)))
  }

  fixed <- JAGS_predict_formula(fit, "mu", data = prediction_new, formula_target = "fixed")
  expect_moments(fixed$value, t(vapply(prediction_new$x, function(x) prediction_moments(prediction_weights(x)), numeric(2))))

  conditional <- JAGS_predict_formula(fit, "mu", data = prediction_new, formula_target = "conditional",
                                      new_levels = "sample", seed = 1)
  expect_moments(conditional$value, rbind(
    prediction_moments(prediction_weights(-1, "a")),
    prediction_moments(prediction_weights(0.5, "d")),
    prediction_moments(prediction_weights(2), extra_variance = prediction_sd^2)
  ))

  # the marginal target: fixed means and the random-effect covariance
  marginal <- JAGS_predict_formula(fit, "mu", data = prediction_new, formula_target = "marginal",
                                   marginal_method = "covariance", new_levels = "sample")
  expect_moments(marginal$mean, t(vapply(prediction_new$x, function(x) prediction_moments(prediction_weights(x)), numeric(2))))
  same_group <- outer(prediction_new$g, prediction_new$g, "==") * prediction_sd^2
  # one covariance per draw (draw x row x row), each the known one
  expect_lt(max(abs(sweep(marginal$vcov$samples, c(2L, 3L), same_group))), 1e-10)
})
