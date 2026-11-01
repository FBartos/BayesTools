skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Model averaging over all eight subsets of three predictors against exact
# enumeration. With a known residual SD of 1 and N(0, 1) priors on the
# intercept and every included slope, a model's marginal likelihood is
# N(y; 0, X X' + I) and its slopes have the posterior N(V X'y, V) with
# V = (X'X + I)^-1, so posterior model probabilities, inclusion Bayes
# factors, and model-averaged posterior means are exact. Bridge sampling of
# these four- to five-dimensional targets from 4 x 5000 draws has errors
# below 0.001 (checked across seeds), so the log marginal likelihoods are
# held to 0.02, the posterior probabilities, which change by at most
# 0.25 x the log marginal likelihood error, to 0.005, and the inclusion
# Bayes factors to 0.05 on the log scale. Sampled quantities are held to
# four Monte Carlo SDs at their effective sample size.

averaging_data <- local({
  set.seed(2)
  data <- data.frame(x1 = stats::rnorm(40), x2 = stats::rnorm(40), x3 = stats::rnorm(40))
  data$y <- 0.1 + 0.45 * data$x1 + 0.2 * data$x2 + stats::rnorm(40)
  data
})
averaging_terms   <- c("x1", "x2", "x3")
averaging_subsets <- as.matrix(expand.grid(x1 = c(FALSE, TRUE), x2 = c(FALSE, TRUE), x3 = c(FALSE, TRUE)))
averaging_exact <- local({
  gaussian_log_density <- function(y, covariance) {
    factor <- chol(covariance)
    z <- backsolve(factor, y, transpose = TRUE)
    -0.5 * (length(y) * log(2 * pi) + sum(z^2)) - sum(log(diag(factor)))
  }
  models <- lapply(seq_len(nrow(averaging_subsets)), function(i) {
    included <- averaging_terms[averaging_subsets[i, ]]
    X <- cbind(1, as.matrix(averaging_data[, included, drop = FALSE]))
    V <- solve(crossprod(X) + diag(ncol(X)))
    mean <- stats::setNames(numeric(3), averaging_terms)
    second_moment <- mean
    mean[included] <- drop(V %*% crossprod(X, averaging_data$y))[-1L]
    second_moment[included] <- diag(V)[-1L] + mean[included]^2
    list(logml = gaussian_log_density(averaging_data$y, tcrossprod(X) + diag(nrow(X))),
         mean = mean, second_moment = second_moment)
  })
  logml <- vapply(models, `[[`, numeric(1), "logml")
  post  <- exp(logml - max(logml)) / sum(exp(logml - max(logml)))
  inclusion <- colSums(post * averaging_subsets)
  mean <- colSums(post * t(vapply(models, `[[`, numeric(3), "mean")))
  list(
    logml       = logml,
    post_prob   = post,
    log_BF      = stats::qlogis(inclusion),
    mean        = mean,
    sd          = sqrt(colSums(post * t(vapply(models, `[[`, numeric(3), "second_moment"))) - mean^2),
    inclusion   = inclusion
  )
})

averaging_likelihood <- "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) } }"
averaging_data_list  <- list(y = averaging_data$y, N = nrow(averaging_data))
averaging_unit       <- prior("normal", list(mean = 0, sd = 1))

test_that("bridge-sampled model averaging reproduces exact enumeration", {

  log_posterior <- function(parameters, data) sum(stats::dnorm(data$y, parameters$mu, 1, log = TRUE))
  models <- lapply(seq_len(nrow(averaging_subsets)), function(i) {
    included <- averaging_terms[averaging_subsets[i, ]]
    formula  <- if (length(included) > 0L) stats::reformulate(included) else ~ 1
    priors   <- c(list(intercept = averaging_unit), stats::setNames(rep(list(averaging_unit), length(included)), included))
    fit <- JAGS_fit(averaging_likelihood, averaging_data_list,
                    formula_list = list(mu = formula), formula_data_list = list(mu = averaging_data),
                    formula_prior_list = list(mu = priors), chains = 4, adapt = 500, burnin = 500, sample = 5000, seed = 1)
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = averaging_data_list, seed = 1),
         prior_weights = 1)
  })
  expect_equal_each(vapply(models, function(model) model$marglik$logml, numeric(1)), averaging_exact$logml, tolerance = 0.02)

  models <- models_inference(models)
  expect_equal_each(vapply(models, function(model) model$inference$post_prob, numeric(1)), averaging_exact$post_prob, tolerance = 0.005)

  parameters   <- paste0("mu_", averaging_terms)
  is_null_list <- stats::setNames(lapply(averaging_terms, function(term) !averaging_subsets[, term]), parameters)
  inference    <- ensemble_inference(models, parameters = parameters, is_null_list = is_null_list, conditional = FALSE)
  expect_equal_each(log(vapply(inference, `[[`, numeric(1), "BF")), unname(averaging_exact$log_BF), tolerance = 0.05)

  mixed <- mix_posteriors(models, parameters = parameters, is_null_list = is_null_list, seed = 1)
  draws <- vapply(mixed, length, numeric(1))
  expect_true(all(abs(vapply(mixed, mean, numeric(1)) - averaging_exact$mean) < 4 * averaging_exact$sd / sqrt(draws)))
})

test_that("a spike-and-slab fit samples the exact posterior model probabilities", {

  spike <- prior_spike_and_slab(averaging_unit)
  fit <- JAGS_fit(averaging_likelihood, averaging_data_list,
                  formula_list = list(mu = ~ x1 + x2 + x3), formula_data_list = list(mu = averaging_data),
                  formula_prior_list = list(mu = list(intercept = averaging_unit, x1 = spike, x2 = spike, x3 = spike)),
                  chains = 4, adapt = 1000, burnin = 1000, sample = 20000, seed = 1)
  estimates  <- as.data.frame(JAGS_estimates_table(fit))
  indicators <- paste0("(mu) ", averaging_terms, " (inclusion)")
  ess <- min(estimates[indicators, "ESS"])

  samples  <- JAGS_estimates_table(fit, return_samples = TRUE)[, indicators] > 0.5
  sampled  <- vapply(seq_len(nrow(averaging_subsets)), function(i) {
    mean(colSums(t(samples) == averaging_subsets[i, ]) == 3L)
  }, numeric(1))
  p <- averaging_exact$post_prob
  expect_true(all(abs(sampled - p) < 4 * sqrt(p * (1 - p) / ess)))

  inclusion <- averaging_exact$inclusion
  expect_true(all(abs(colMeans(samples) - inclusion) < 4 * sqrt(inclusion * (1 - inclusion) / ess)))
})
