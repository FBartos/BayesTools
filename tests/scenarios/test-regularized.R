if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized")

# Gaussian regression of fertility in 47 Swiss provinces. Every continuous
# variable is standardized once, and both engines use residual SD 1.
testthat::test_that("regularized regression of Swiss fertility", {
  dat <- as.data.frame(scale(datasets::swiss))
  slopes <- setdiff(names(dat), "Fertility")
  design <- as.matrix(dat[, slopes])
  models <- c("ridge", "lasso", "horseshoe", "regularized_horseshoe")
  stan_file <- file.path(.scenario_config()$root, "stan", "regularized.stan")
  stan_data <- list(N = nrow(dat), K = length(slopes), X = design, y = dat$Fertility)
  stan_model <- NULL
  stan <- list()
  for (model in models) {
    stan[[model]] <- scenario_fit(paste0("stan_", model), {
      if (is.null(stan_model)) {
        stan_model <- rstan::stan_model(file = stan_file, auto_write = FALSE)
      }
      fit <- rstan::sampling(stan_model,
        data = c(stan_data, list(model_id = match(model, models))),
        chains = 3, cores = 3, warmup = 1500, iter = 6500,
        control = list(adapt_delta = 0.995, max_treedepth = 14), seed = 26011, refresh = 0
      )
      list(draws = as.array(fit, pars = c("alpha", "beta", "tau")),
           diagnostics = lapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), colSums))
    }, cache_version = strtoi(substr(tools::md5sum(stan_file), 1, 7), base = 16L) + 1L)
  }

  # For ridge, conditioning on tau gives a Gaussian posterior. Integrating its
  # means with the Gaussian marginal likelihood gives an independent reference.
  # tau = 0.2 * tan(angle) makes the half-Cauchy measure uniform on (0, pi/2).
  gram <- eigen(crossprod(design), symmetric = TRUE)
  xy <- drop(crossprod(gram$vectors, crossprod(design, dat$Fertility)))
  ridge_integral <- function(component = NULL) {
    integrand <- function(angle) {
      tau2 <- (0.2 * tan(angle))^2
      variance <- tau2 / (1 + tau2 * gram$values)
      weight <- exp(sum(variance * xy^2 - log1p(tau2 * gram$values)) / 2)
      if (is.null(component)) return(weight)
      weight * sum(gram$vectors[component, ] * variance * xy)
    }
    integrate(Vectorize(integrand), 0, pi / 2, rel.tol = 1e-10)$value
  }
  exact_ridge <- vapply(seq_along(slopes), ridge_integral, numeric(1)) / ridge_integral()

  formula <- reformulate(slopes)
  likelihood <- "model { for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) } }"
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
  for (model in models) {
    jags <- scenario_fit(paste0("jags_", model), {
      shrinkage <- if (model == "regularized_horseshoe") {
        prior_regularized("horseshoe", scale = global, slab_scale = 2, slab_df = 4)
      } else {
        prior_regularized(model, scale = global)
      }
      fit <- JAGS_fit(likelihood, data = list(y = dat$Fertility, N = nrow(dat)),
                      prior_list = list(), formula_list = list(mu = formula),
                      formula_data_list = list(mu = dat),
                      formula_prior_list = list(mu = list(
                        intercept = prior("normal", list(0, 1)), shrinkage = shrinkage
                      )),
                      chains = 3, adapt = 1500, burnin = 3000, sample = 15000,
                      seed = 89211)
      as.array(coda::as.mcmc.list(fit))
    })
    # coda arrays are iterations x variables x chains.
    jags <- aperm(jags, c(1, 3, 2))
    stan_draws <- stan[[model]]$draws
    jags_columns <- paste0("mu_", c("intercept", slopes, "shrinkage_global"))
    stan_columns <- c("alpha", sprintf("beta[%d]", seq_along(slopes)), "tau")
    summarize <- function(draws, columns) {
      values <- lapply(columns, function(nm) draws[, , nm])
      names(values) <- c("intercept", slopes, "log(global)")
      # A log scale permits a meaningful mean even with a half-Cauchy upper tail.
      values[["log(global)"]] <- log(values[["log(global)"]])
      values[["P(both within 0.2)"]] <-
        (abs(values$Agriculture) < 0.2 & abs(values$Examination) < 0.2) * 1
      t(vapply(values, function(x) {
        c(mean = mean(x), mcse = posterior::mcse_mean(x),
          q05 = unname(quantile(x, 0.05)), q95 = unname(quantile(x, 0.95)),
          rhat = posterior::rhat(x))
      }, numeric(5)))
    }
    reference <- summarize(stan_draws, stan_columns)
    estimate <- summarize(jags, jags_columns)
    combined_mcse <- sqrt(reference[, "mcse"]^2 + estimate[, "mcse"]^2)
    z <- (estimate[, "mean"] - reference[, "mean"]) / combined_mcse
    if (model == "ridge") {
      for (summary in list(reference, estimate)) {
        testthat::expect_lt(max(abs(summary[slopes, "mean"] - exact_ridge) /
                                  summary[slopes, "mcse"]), 4.5)
      }
    }
    divergences <- sum(vapply(stan[[model]]$diagnostics, "[[", numeric(1), "divergent__"))
    comparison <- cbind(Stan = reference[, "mean"], JAGS = estimate[, "mean"],
                        MCSE = combined_mcse, z = z,
                        Stan05 = reference[, "q05"], JAGS05 = estimate[, "q05"],
                        Stan95 = reference[, "q95"], JAGS95 = estimate[, "q95"],
                        Rhat = pmax(reference[, "rhat"], estimate[, "rhat"]))
    scenario_text(paste0("comparison-", model), {
      cat("Stan post-warmup divergences:", divergences, "\n")
      print(round(comparison, 4))
    })
    testthat::expect_equal(divergences, 0)
    testthat::expect_true(all(is.finite(z) & abs(z) < 4.5))
    testthat::expect_lt(max(comparison[, "Rhat"]), 1.02)
  }
})
