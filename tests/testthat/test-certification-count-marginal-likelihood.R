skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Bridge sampling of count GLMMs against a reference marginal likelihood:
# Poisson and negative-binomial (Poisson-gamma) models with random group
# intercepts. Given the coefficients, the group SD tau, and the size phi, each
# group's intercept u integrates by 20-node Gauss-Hermite quadrature centred
# at its conditional mode (damped Newton steps with the analytic derivatives
# of the concave log integrand) and scaled by the curvature there (in
# development, 40 nodes changed the reference by at most 3e-12, and the rule
# agreed with statmod::gauss.quad() to 2e-13). The coefficients, log tau, and
# log phi are importance-sampled from a multivariate t (5 df) fitted to the
# posterior draws, whose relative standard error is the standard error of the
# reference on the log scale. The bridge
# estimates of 4 x 5000 draws are held to 0.05 plus four reference standard
# errors.

count_data <- local({
  set.seed(6)
  data <- data.frame(x = stats::rnorm(60), g = factor(rep(seq_len(15), each = 4)))
  mu <- exp(0.5 + 0.4 * data$x + stats::rnorm(15, 0, 0.5)[as.integer(data$g)])
  data$y <- stats::rnbinom(60, size = 3, mu = mu)
  data
})

# the n-node Gauss-Hermite rule (Golub-Welsch): the nodes are the eigenvalues
# of the Jacobi matrix with off-diagonal sqrt(k / 2), k = 1, ..., n - 1, and
# the weights sqrt(pi) times the squared first components of its eigenvectors
gauss_hermite <- function(n) {
  jacobi <- matrix(0, n, n)
  off    <- sqrt(seq_len(n - 1L) / 2)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = sqrt(pi) * decomposition$vectors[1L, ]^2)
}

count_reference <- function(fit, nb, size = 20000L, n_nodes = 20L) {
  group <- as.integer(count_data$g)
  y     <- count_data$y
  X     <- cbind(1, count_data$x)
  nodes <- gauss_hermite(n_nodes)
  log_groups <- function(beta, tau, phi) {
    eta <- drop(X %*% beta)
    log_f <- function(u) {
      lambda <- exp(eta + u[group, , drop = FALSE])
      log_p  <- if (nb) stats::dnbinom(y, size = phi, mu = lambda, log = TRUE) else stats::dpois(y, lambda, log = TRUE)
      rowsum(log_p, group, reorder = TRUE) + stats::dnorm(u, 0, tau, log = TRUE)
    }
    derivatives <- function(u) {
      lambda <- exp(eta + u[group])
      d1 <- if (nb) y - (y + phi) * lambda / (phi + lambda) else y - lambda
      d2 <- if (nb) -(y + phi) * phi * lambda / (phi + lambda)^2 else -lambda
      cbind(rowsum(d1, group, reorder = TRUE) - u / tau^2, rowsum(d2, group, reorder = TRUE) - 1 / tau^2)
    }
    u <- numeric(max(group))
    for (iteration in 1:50) {
      d    <- derivatives(u)
      step <- pmax(pmin(-d[, 1L] / d[, 2L], 1), -1)
      u    <- u + step
      if (max(abs(step)) < 1e-10) break
    }
    s <- sqrt(-1 / derivatives(u)[, 2L])
    values <- log_f(u + sqrt(2) * s %*% t(nodes$nodes)) +
      matrix(log(nodes$weights) + nodes$nodes^2, max(group), n_nodes, byrow = TRUE)
    top <- apply(values, 1L, max)
    sum(top + log(rowSums(exp(values - top))) + log(sqrt(2) * s))
  }
  draws <- JAGS_estimates_table(fit, return_samples = TRUE)
  U <- cbind(draws[, c("(mu) intercept", "(mu) x")], log(draws[, "(mu) sd(intercept)"]), if (nb) log(draws[, "phi"]))
  set.seed(1)
  center <- colMeans(U)
  scale  <- 1.5 * stats::cov(U)
  proposals <- mvtnorm::rmvt(size, sigma = scale, df = 5, delta = center, type = "shifted")
  log_q     <- mvtnorm::dmvt(proposals, delta = center, sigma = scale, df = 5, log = TRUE, type = "shifted")
  log_w <- vapply(seq_len(size), function(i) {
    p   <- proposals[i, ]
    tau <- exp(p[[3L]])
    phi <- if (nb) exp(p[[4L]]) else NULL
    log_groups(p[1:2], tau, phi) + stats::dnorm(p[[1L]], 0, 2.5, log = TRUE) + stats::dnorm(p[[2L]], 0, 1, log = TRUE) +
      log(2) + stats::dnorm(tau, 0, 1, log = TRUE) + p[[3L]] +
      if (nb) stats::dgamma(phi, 2, 0.5, log = TRUE) + p[[4L]] else 0
  }, numeric(1)) - log_q
  weights <- exp(log_w - max(log_w))
  c(logml = max(log_w) + log(mean(weights)), se = stats::sd(weights) / sqrt(size) / mean(weights))
}

test_that("bridge sampling of Poisson and negative-binomial GLMMs matches the reference marginal likelihoods", {

  intercept_prior <- prior("normal", list(mean = 0, sd = 2.5))
  unit_prior      <- prior("normal", list(mean = 0, sd = 1))
  sd_prior        <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  phi_prior       <- prior("gamma", list(shape = 2, rate = 0.5))
  data_list <- list(y = count_data$y, N = nrow(count_data))
  likelihoods <- list(
    poisson = "model{ for(i in 1:N){ y[i] ~ dpois(exp(mu[i])) } }",
    nb      = "model{ for(i in 1:N){ y[i] ~ dnegbin(phi / (phi + exp(mu[i])), phi) } }"
  )
  log_posteriors <- list(
    poisson = function(parameters, data) sum(stats::dpois(data$y, exp(parameters$mu), log = TRUE)),
    nb      = function(parameters, data) sum(stats::dnbinom(data$y, size = parameters$phi, mu = exp(parameters$mu), log = TRUE))
  )
  logml <- list()
  for (model in names(likelihoods)) {
    prior_list <- if (identical(model, "nb")) list(phi = phi_prior) else NULL
    fit <- JAGS_fit(likelihoods[[model]], data_list, prior_list = prior_list,
                    formula_list = list(mu = ~ x + (1 | g)), formula_data_list = list(mu = count_data),
                    formula_prior_list = list(mu = list(intercept = intercept_prior, x = unit_prior)),
                    formula_random_prior_list = list(mu = prior_random(g = random_block(sd = sd_prior))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 5000, seed = 1)
    bridge    <- JAGS_bridgesampling(fit, log_posterior = log_posteriors[[model]], data = data_list,
                                     prior_list = prior_list, seed = 1)$logml
    reference <- count_reference(fit, nb = identical(model, "nb"))
    expect_lt(abs(bridge - reference[["logml"]]), 0.05 + 4 * reference[["se"]])
    logml[[model]] <- c(bridge = bridge, reference = reference[["logml"]])
  }
  # the overdispersion Bayes factor
  expect_lt(abs((logml$nb - logml$poisson)[["bridge"]] - (logml$nb - logml$poisson)[["reference"]]), 0.1)
})
