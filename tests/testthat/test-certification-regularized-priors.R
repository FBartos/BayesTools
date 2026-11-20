skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("rstan")
skip_if_not_installed("posterior")

regularized_prior <- function(family, contrast = NULL, units = 1) {

  global <- prior("cauchy", list(0, 0.2 * units), truncation = list(0, Inf))
  if (family == "finite-slab-horseshoe") {
    return(prior_regularized("horseshoe", global, contrast = contrast, slab_scale = 2 * units, slab_df = 4))
  }
  prior_regularized(family, global, contrast = contrast)
}

# With centered predictors and intercept independent of the slopes, integrate
# the Gaussian posterior mean of each ridge slope over the posterior of tau.
# tau = 0.2 tan(angle) makes the half-Cauchy prior uniform in angle.
ridge_posterior_mean <- function(X, y) {

  gram <- eigen(crossprod(X), symmetric = TRUE)
  xy <- drop(crossprod(gram$vectors, crossprod(X, y)))
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
  vapply(seq_len(ncol(X)), integral, numeric(1)) / integral()
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
  exact_ridge <- ridge_posterior_mean(swiss_X, swiss$Fertility)

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
      if (design == "swiss" && family == "ridge") {
        for (estimate in list(a, b)) {
          error <- (estimate[slopes, "mean"] - exact_ridge) / estimate[slopes, "mcse"]
          expect_lt(max(abs(error)), 4.5, label = "ridge versus exact posterior")
        }
      }
    }
  }
})
