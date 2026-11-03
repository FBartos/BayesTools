skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Coefficient priors scaled by the residual SD ('multiply_by' = "sigma", the
# default priors of the BayesFactor package): y_i = mu + sigma theta x_i +
# e_i with e_i ~ N(0, sigma^2), theta ~ Cauchy(0, r), mu ~ uniform(-100, 100),
# and sigma ~ Gamma(0.001, 0.001) truncated to [0.01, 100]. The intercept
# integrates analytically (its posterior lies far inside the uniform range),
# leaving p(y | theta, sigma) = (2 pi sigma^2)^(-(n - 1) / 2) n^(-1 / 2) / 200
# exp(-S / (2 sigma^2)) with S = Syy - 2 sigma theta Sxy + sigma^2 theta^2
# Sxx of the centred data. theta = r tan(phi) with phi uniform under the
# Cauchy prior, and log sigma, are integrated by trapezoid sums on fine grids
# (halving both steps changes the log Bayes factor by about 1e-12; the
# truncated gamma prior reproduces BayesFactor's Jeffreys-prior values within
# 1.1e-4: 0.235804 vs 0.235803 and 12.94298 vs 12.94309). Two cases: a
# two-level orthonormal factor (the sleep groups, r = 1 / 2) and a slope of a
# predictor standardized by its population SD (attitude complaints,
# r = sqrt(2) / 4). The bridge log Bayes factors of 4 x 10000 draws are held to
# 0.05; for the factor, the Savage-Dickey ratio, whose prior ordinate is the
# Cauchy density of the unscaled theta at 0, is held to four Monte Carlo SDs
# of the kernel ordinate at the exact posterior density of theta at 0 (as in
# the hypothesis-bayes-factors case) plus 0.05 for its smoothing bias (the
# slope's posterior density at 0 is about 2e-6, too far in the tail for a
# kernel ordinate).

scaled_exact <- function(y, x, r, step_sigma = 0.002, n_phi = 8000L) {
  n   <- length(y)
  yc  <- y - mean(y)
  xc  <- x - mean(x)
  Syy <- sum(yc^2); Sxy <- sum(xc * yc); Sxx <- sum(xc^2)
  log_sigma  <- seq(log(0.01), log(100), by = step_sigma)
  sigma      <- exp(log_sigma)
  log_prior  <- stats::dgamma(sigma, 0.001, 0.001, log = TRUE) + log_sigma + log(step_sigma) -
    log(stats::pgamma(100, 0.001, 0.001) - stats::pgamma(0.01, 0.001, 0.001))
  phi   <- seq(-pi / 2, pi / 2, length.out = n_phi + 1L)[2:n_phi]
  theta <- r * tan(phi)
  log_lik <- function(sigma, theta) {
    -(n - 1) / 2 * log(2 * pi * sigma^2) - 0.5 * log(n) - log(200) -
      (Syy - 2 * sigma * theta * Sxy + sigma^2 * theta^2 * Sxx) / (2 * sigma^2)
  }
  logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))
  # per sigma: the integral over theta, and the posterior moments of theta
  rows <- t(vapply(sigma, function(s) {
    l <- log_lik(s, theta) + log(phi[2L] - phi[1L]) - log(pi)
    w <- exp(l - max(l))
    c(log = logsumexp(l), mean = sum(w * theta) / sum(w), second = sum(w * theta^2) / sum(w))
  }, numeric(3)))
  log_joint <- rows[, "log"] + log_prior
  log_ml_1  <- logsumexp(log_joint)
  weights   <- exp(log_joint - log_ml_1)
  log_ml_0  <- logsumexp(log_lik(sigma, 0) + log_prior)
  mean_theta <- sum(weights * rows[, "mean"])
  list(
    log_BF  = log_ml_1 - log_ml_0,
    # the posterior density of theta at 0: p(y | theta = 0) p(theta = 0) / p(y)
    density_at_zero = exp(log_ml_0 - log_ml_1) * stats::dcauchy(0, 0, r),
    sd_theta = sqrt(sum(weights * rows[, "second"]) - mean_theta^2)
  )
}

scaled_cases <- local({
  sl <- sleep
  sl$y <- as.numeric(scale(sl$extra))
  at <- attitude
  at$y <- as.numeric(scale(at$rating))
  at$complaints_s <- (at$complaints - mean(at$complaints)) / sqrt(mean((at$complaints - mean(at$complaints))^2))
  list(
    factor = list(data = sl, formula = ~ group, term = "group", r = 1 / 2, statement = "group[1] = 0",
                  x = contr.orthonormal(levels(sl$group))[as.integer(sl$group), 1L],
                  coefficient = function(r) prior_factor("mcauchy", list(location = 0, scale = r), contrast = "orthonormal")),
    slope  = list(data = at, formula = ~ complaints_s, term = "complaints_s", r = sqrt(2) / 4, statement = NULL,
                  x = at$complaints_s,
                  coefficient = function(r) prior("cauchy", list(location = 0, scale = r)))
  )
})

test_that("bridge sampling and Savage-Dickey match the exact Bayes factors of sigma-scaled priors", {

  sigma_prior <- prior("gamma", list(shape = 0.001, rate = 0.001), truncation = list(lower = 0.01, upper = 100))
  flat_prior  <- prior("uniform", list(a = -100, b = 100))
  likelihood  <- "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], pow(sigma, -2)) } }"
  log_posterior <- function(parameters, data) sum(stats::dnorm(data$y, parameters$mu, parameters$sigma, log = TRUE))

  for (case in scaled_cases) {
    data_list <- list(y = case$data$y, N = nrow(case$data))
    coefficient_prior <- case$coefficient(case$r)
    attr(coefficient_prior, "multiply_by") <- "sigma"
    fit_1 <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = sigma_prior),
                      formula_list = list(mu = case$formula), formula_data_list = list(mu = case$data),
                      formula_prior_list = list(mu = stats::setNames(list(flat_prior, coefficient_prior), c("intercept", case$term))),
                      chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    fit_0 <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = sigma_prior),
                      formula_list = list(mu = ~ 1), formula_data_list = list(mu = case$data),
                      formula_prior_list = list(mu = list(intercept = flat_prior)),
                      chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    marglik <- function(fit) JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list,
                                                 prior_list = list(sigma = sigma_prior), seed = 1)
    exact <- scaled_exact(case$data$y, case$x, case$r)

    expect_lt(abs(marglik(fit_1)$logml - marglik(fit_0)$logml - exact$log_BF), 0.05)

    # Savage-Dickey on the unscaled coefficient
    if (is.null(case$statement)) next
    ess  <- min(as.data.frame(JAGS_estimates_table(fit_1))[["ESS"]], na.rm = TRUE)
    h    <- 0.9 * exact$sd_theta * 40000^(-1 / 5)
    kernel_sd <- sqrt(1 / (2 * sqrt(pi)) / (ess * h * exact$density_at_zero))
    savage_dickey <- log(as.data.frame(hypothesis_BF(fit_1, case$statement, seed = 1))[["BF10"]])
    expect_lt(abs(savage_dickey - exact$log_BF), 4 * kernel_sd + 0.05)
  }
})
