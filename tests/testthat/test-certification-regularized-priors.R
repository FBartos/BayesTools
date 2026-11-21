skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("rstan")
skip_if_not_installed("posterior")
skip_if_not_installed("bridgesampling")

regularized_prior <- function(family, contrast = NULL, units = 1) {

  global <- prior("cauchy", list(0, 0.2 * units), truncation = list(0, Inf))
  if (family == "finite-slab-horseshoe") {
    return(prior_regularized("horseshoe", global, contrast = contrast, slab_scale = 2 * units, slab_df = 4))
  }
  prior_regularized(family, global, contrast = contrast)
}

# Integrate the independent N(0, 1) intercept first. Its marginal covariance
# A = I + 11' has determinant n+1 and inverse I - 11'/(n+1), including for
# the unbalanced month design. Then integrate the N(0, tau^2 I) slopes.
ridge_reference <- function(X, y) {

  n <- length(y)
  weighted_y <- y - sum(y) / (n + 1)
  weighted_X <- sweep(X, 2L, colSums(X) / (n + 1), "-")
  gram <- eigen(crossprod(X, weighted_X), symmetric = TRUE)
  xy <- drop(crossprod(gram$vectors, crossprod(X, weighted_y)))
  # tau = 0.2 tan(angle) makes the half-Cauchy density 2/pi in angle.
  integral <- function(component = NULL) {
    integrand <- function(angle) {
      tau2 <- (0.2 * tan(angle))^2
      variance <- tau2 / (1 + tau2 * gram$values)
      weight <- exp(sum(variance * xy^2 - log1p(tau2 * gram$values)) / 2)
      if (is.null(component)) return(weight)
      weight * sum(gram$vectors[component, ] * variance * xy)
    }
    integrate(Vectorize(integrand), 0, pi / 2, rel.tol = 1e-10)$value
  }
  normalizer <- integral()
  intercept_only <- -0.5 * (n * log(2 * pi) + log1p(n) + sum(y * weighted_y))
  list(mean = vapply(seq_len(ncol(X)), integral, numeric(1)) / normalizer,
       logml = intercept_only + log((2 / pi) * normalizer))
}

# Check the exact target used by RStan's bridge sampler, independently of its
# posterior draws. The transform Jacobian is excluded on both sides because
# these are densities in the constrained sampling parameters.
check_stan_log_density <- function(stan, program, X, y, unit) {

  alpha <- -0.15
  tau <- 0.35
  par <- list(alpha = alpha, tau = tau)
  prior <- dnorm(alpha, log = TRUE) + log(2) + dcauchy(tau, scale = 0.2, log = TRUE)
  coefficient <- seq(-0.4, 0.6, length.out = ncol(X))
  if (program %in% c("horseshoe", "finite-slab-horseshoe")) {
    par$lambda <- seq(0.4, 1.3, length.out = max(unit))
    prior <- prior + sum(log(2) + dcauchy(par$lambda, log = TRUE))
    variance <- (tau * par$lambda[unit])^2
    if (program == "finite-slab-horseshoe") {
      par$slab_aux <- 0.8
      prior <- prior + dgamma(1 / par$slab_aux, 2, rate = 2, log = TRUE) - 2 * log(par$slab_aux)
      variance <- 1 / (1 / variance + 1 / (4 * par$slab_aux))
    }
    covariance <- diag(length(y)) + tcrossprod(sweep(X, 2L, sqrt(variance), "*"))
    likelihood <- mvtnorm::dmvnorm(y - alpha, sigma = covariance, log = TRUE)
  } else {
    if (program == "grouped-lasso") {
      par$beta <- coefficient
      dimension <- ncol(X) - 1L
      radius <- sqrt(sum(head(coefficient, dimension)^2))
      # R ~ Gamma(d, scale=tau), with a uniform direction on the unit sphere.
      sphere_area <- 2 * pi^(dimension / 2) / gamma(dimension / 2)
      prior <- prior + dgamma(radius, dimension, scale = tau, log = TRUE) -
        log(sphere_area) - (dimension - 1) * log(radius) -
        log(2 * tau) - abs(tail(coefficient, 1L)) / tau
    } else {
      par$z <- coefficient
      prior <- prior + if (program == "ridge") sum(dnorm(par$z, log = TRUE)) else sum(-log(2) - abs(par$z))
      coefficient <- tau * par$z
    }
    likelihood <- sum(dnorm(y, alpha + drop(X %*% coefficient), 1, log = TRUE))
  }
  unconstrained <- rstan::unconstrain_pars(stan, par)
  expect_equal(rstan::log_prob(stan, unconstrained, adjust_transform = FALSE),
               prior + likelihood, tolerance = 1e-10, info = program)
}

regularized_draw_summary <- function(draws, columns, basis) {

  values <- lapply(columns, function(column) draws[, , column])
  if (is.null(basis)) {
    values[["P(both within 0.2)"]] <-
      1 * (abs(values$Agriculture) < 0.2 & abs(values$Examination) < 0.2)
  } else {
    coefficients <- values[colnames(basis)]
    for (month in rownames(basis)) {
      values[[month]] <- Reduce(`+`, Map(`*`, coefficients, basis[month, ]))
    }
    radius <- sqrt(Reduce(`+`, lapply(coefficients, `^`, 2)))
    values[["P(month radius < 0.75)"]] <- 1 * (radius < 0.75)
  }
  # These scale priors are heavy-tailed. Compare their log-scale means.
  scales <- intersect(c("global", "month local", "wind local"), names(values))
  values[scales] <- lapply(values[scales], log)
  t(vapply(values, function(x) {
    c(mean = mean(x), mcse = posterior::mcse_mean(x), rhat = posterior::rhat(x))
  }, numeric(3)))
}

regularized_families <- c("ridge", "lasso", "horseshoe", "finite-slab-horseshoe")

regularized_bridge_checks <- function(jags, stan, y, seed, case) {

  log_likelihood <- function(parameters, data) {
    sum(dnorm(data$y, parameters$mu, 1, log = TRUE))
  }
  jags_bridge <- JAGS_bridgesampling(jags, log_posterior = log_likelihood,
    data = list(y = y, N = length(y)), prior_list = list(), seed = seed)
  stan_bridge <- withr::with_seed(seed + 1L,
    bridgesampling::bridge_sampler(stan, method = "normal", cores = 1L, silent = TRUE, maxiter = 10000L))
  expect_true(isTRUE(jags_bridge$success) && isTRUE(jags_bridge$repetitions$within_maxiter), info = case)
  expect_length(jags_bridge$diagnostics$upstream_warnings, 0L)
  expect_lt(stan_bridge$niter, 10000L, label = paste(case, "RStan bridge iterations"))
  for (target in c("q11", "q21")) {
    values <- stan_bridge[[target]]
    expect_true(length(values) > 0L && all(is.finite(values)), info = paste(case, "RStan", target))
  }
  logml <- c(BayesTools = jags_bridge$logml, RStan = stan_bridge$logml)
  mcse <- c(BayesTools = jags_bridge$repetitions$mcse, RStan = stan_bridge$mcse_logml)
  expect_true(length(mcse) == 2L && all(is.finite(mcse) & mcse > 0), info = case)
  # The precision target is an estimated log-evidence MCSE below 0.05
  # (about 5% relative Monte Carlo standard error in marginal likelihood).
  expect_lt(max(mcse), 0.05, label = paste(case, "maximum log evidence MCSE"))
  # Cross-engine fits and proposal draws use separate RNG streams. Use the
  # same simultaneous 4.5-SE threshold as the posterior comparisons below.
  expect_lt(abs(diff(logml)) / sqrt(sum(mcse^2)), 4.5, label = paste(case, "log evidence difference / SE"))
  list(logml = logml, mcse = mcse)
}

test_that("both centered codings retain their coded-coefficient prior scale", {

  months <- data.frame(Month = factor(month.abb[5:9]))
  conversion <- sqrt(1 - 1 / nlevels(months$Month))
  for (family in regularized_families) {
    level_draws <- function(contrast, units) {
      bound <- JAGS_formula(~ Month, "mu", months, list(
        intercept = prior("normal", list(0, 1)), sh = regularized_prior(family, contrast, units)
      ))
      set.seed(50931)
      unname(rng(bound$prior_list$mu_Month, 256L))
    }
    orthonormal <- level_draws("orthonormal", 1)
    expect_equal(level_draws("meandif", 1) * conversion, orthonormal, tolerance = 1e-10)
    expect_equal(level_draws("meandif", conversion), orthonormal, tolerance = 1e-10)
  }
})

test_that("regression shrinkage matches Stan in coefficient, level, joint, and scale summaries", {

  swiss <- as.data.frame(scale(datasets::swiss))
  slopes <- setdiff(names(swiss), "Fertility")
  swiss_X <- as.matrix(swiss[slopes])

  ozone <- na.omit(datasets::airquality[c("Ozone", "Month", "Wind")])
  ozone$Ozone <- as.numeric(scale(ozone$Ozone))
  ozone$Month <- factor(ozone$Month, labels = month.abb[5:9])
  ozone$Wind <- as.numeric(scale(ozone$Wind))
  # Independent reference: build the two bases with base R.
  helmert <- contr.helmert(levels(ozone$Month))
  orthonormal <- sweep(helmert, 2, sqrt(colSums(helmert^2)), "/")
  colnames(orthonormal) <- paste0("b", seq_len(ncol(orthonormal)))
  bases <- list(orthonormal = orthonormal,
                meandif = orthonormal / sqrt(1 - 1 / nrow(orthonormal)))

  # The numerical verification uses the very same Stan programs as the two
  # readable scenarios. Fit fresh so a change to BayesTools is exercised.
  programs <- c(regularized_families, "grouped-lasso")
  models <- lapply(setNames(programs, programs), function(program) {
    file <- testthat::test_path("..", "scenarios", "stan", paste0(program, ".stan"))
    rstan::stan_model(file, auto_write = FALSE)
  })
  intercept <- prior("normal", list(0, 1))
  for (design in c("swiss", names(bases))) {
    basis <- bases[[design]]
    if (is.null(basis)) {
      dat <- swiss
      y <- dat$Fertility
      X <- swiss_X
      formula <- reformulate(slopes)
      contrast <- NULL
      beta_names <- jags_terms <- slopes
      unit <- seq_len(ncol(X))
      group <- "shrinkage"
      jags_seed <- 89211
      stan_seed <- 26011
    } else {
      dat <- ozone
      y <- dat$Ozone
      X <- cbind(basis[as.integer(dat$Month), ], dat$Wind)
      formula <- ~ Month + Wind
      contrast <- design
      beta_names <- c(colnames(basis), "Wind")
      jags_terms <- c(sprintf("Month[%d]", seq_len(ncol(basis))), "Wind")
      unit <- c(rep(1L, ncol(basis)), 2L)
      group <- "sh"
      jags_seed <- 73429
      stan_seed <- 6283
    }
    for (family in regularized_families) {
      priors <- list(intercept = intercept)
      priors[[group]] <- regularized_prior(family, contrast)
      jags <- JAGS_fit("model{ for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) } }",
        data = list(y = y, N = length(y)), prior_list = list(),
        formula_list = list(mu = formula), formula_data_list = list(mu = dat),
        formula_prior_list = list(mu = priors),
        chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = jags_seed)
      jags_draws <- aperm(as.array(coda::as.mcmc.list(jags)), c(1, 3, 2))

      program <- if (family == "lasso" && !is.null(basis)) "grouped-lasso" else family
      stan_data <- list(N = nrow(X), K = ncol(X), X = X, y = y)
      horseshoe <- family %in% c("horseshoe", "finite-slab-horseshoe")
      if (horseshoe) stan_data <- c(stan_data, list(J = max(unit), unit = unit))
      stan <- rstan::sampling(models[[program]], data = stan_data,
        chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = stan_seed, refresh = 0,
        control = list(adapt_delta = 0.995, max_treedepth = 14))
      stan_draws <- as.array(stan, pars = c("alpha", "beta", "tau", if (horseshoe) "lambda"))
      check_stan_log_density(stan, program, X, y, unit)

      labels <- c("intercept", beta_names, "global")
      jags_columns <- setNames(paste0("mu_", c("intercept", jags_terms, paste0(group, "_global"))), labels)
      stan_columns <- setNames(c("alpha", sprintf("beta[%d]", seq_len(ncol(X))), "tau"), labels)
      if (horseshoe && !is.null(basis)) {
        jags_columns <- c(jags_columns, "month local" = "mu_Month_local_1", "wind local" = "mu_Wind_local_1")
        stan_columns <- c(stan_columns, "month local" = "lambda[1]", "wind local" = "lambda[2]")
      }
      a <- regularized_draw_summary(jags_draws, jags_columns, basis)
      b <- regularized_draw_summary(stan_draws, stan_columns, basis)
      z <- (a[, "mean"] - b[, "mean"]) / sqrt(a[, "mcse"]^2 + b[, "mcse"]^2)
      diagnostics <- rstan::get_sampler_params(stan, inc_warmup = FALSE)
      divergences <- sum(vapply(diagnostics, function(x) sum(x[, "divergent__"]), numeric(1)))
      case <- paste(design, family)
      expect_equal(divergences, 0, info = case)
      expect_lt(max(a[, "rhat"], b[, "rhat"]), 1.02, label = paste(case, "maximum Rhat"))
      expect_lt(max(abs(z)), 4.5, label = paste(case, "maximum |z| at", names(which.max(abs(z)))))
      seed <- 61270L + 10L * match(family, regularized_families)
      bridge <- regularized_bridge_checks(jags, stan, y, seed, case)
      if (family == "ridge") {
        ridge_bridge <- bridge
        exact <- ridge_reference(X, y)
        expect_lt(max(abs(bridge$logml - exact$logml) / bridge$mcse), 4.5,
                  label = paste(case, "ridge versus exact log evidence"))
        for (estimate in list(a, b)) {
          error <- (estimate[beta_names, "mean"] - exact$mean) / estimate[beta_names, "mcse"]
          expect_lt(max(abs(error)), 4.5, label = paste(case, "ridge versus exact posterior"))
        }
      } else {
        log_bf <- bridge$logml - ridge_bridge$logml
        # Families reuse fitting seeds within each engine. Adding their marginal
        # MCSE estimates gives a conservative estimate of the SE of the log BF.
        bf_mcse_bound <- bridge$mcse + ridge_bridge$mcse
        expect_lt(abs(diff(log_bf)) / sqrt(sum(bf_mcse_bound^2)), 4.5,
                  label = paste(case, "log BF versus ridge, difference / SE bound"))
      }
    }
  }
})
