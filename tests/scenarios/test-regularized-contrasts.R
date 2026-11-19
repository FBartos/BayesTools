if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized-contrasts")

# New York daily ozone measurements, May through September 1973.
# A five-level month effect and a standardized wind slope share a global scale.
testthat::test_that("shrinkage of monthly ozone effects", {
  dat <- na.omit(datasets::airquality[c("Ozone", "Month", "Wind")])
  dat$Month <- factor(dat$Month, labels = month.abb[5:9])
  dat$Wind <- as.numeric(scale(dat$Wind))
  dat$Ozone <- as.numeric(scale(dat$Ozone))
  n_months <- nlevels(dat$Month)
  conversion <- sqrt(1 - 1 / n_months)
  # Construct the reference bases from base R, independently of JAGS_formula().
  helmert <- contr.helmert(n_months)
  orthonormal <- sweep(helmert, 2, sqrt(colSums(helmert^2)), "/")
  bases <- list(orthonormal = orthonormal, meandif = orthonormal / conversion)
  families <- c("ridge", "lasso", "horseshoe", "regularized_horseshoe")
  stan_file <- file.path(.scenario_config()$root, "stan", "regularized-contrasts.stan")
  stan_model <- NULL
  reference <- list()
  for (contrast in names(bases)) {
    X <- cbind(bases[[contrast]][as.integer(dat$Month), ], dat$Wind)
    for (family in families) {
      key <- paste(contrast, family, sep = "_")
      reference[[key]] <- scenario_fit(paste0("stan_", key), {
        if (is.null(stan_model)) {
          stan_model <- rstan::stan_model(file = stan_file, auto_write = FALSE)
        }
        fit <- rstan::sampling(stan_model,
          data = list(N = nrow(X), K = ncol(X), X = X, y = dat$Ozone,
                      family = match(family, families)),
          chains = 3, cores = 3, iter = 6500, warmup = 1500, seed = 6283, refresh = 0,
          control = list(adapt_delta = 0.995, max_treedepth = 14)
        )
        pars <- c("alpha", "beta", "tau", if (family %in% families[3:4]) "lambda")
        diagnostics <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
        list(draws = as.array(fit, pars = pars), divergences = sum(vapply(
          diagnostics, function(chain) sum(chain[, "divergent__"]), numeric(1))))
      }, cache_version = strtoi(substr(tools::md5sum(stan_file), 1, 7), base = 16L) + 1L)
    }
  }

  intercept <- prior("normal", list(0, 1))
  prior_group <- function(family, contrast, units = 1) {
    global <- prior("cauchy", list(0, 0.2 * units), truncation = list(0, Inf))
    if (family == "regularized_horseshoe") {
      prior_regularized("horseshoe", global, contrast = contrast, slab_scale = 2 * units, slab_df = 4)
    } else {
      prior_regularized(family, global, contrast = contrast)
    }
  }
  for (family in families) {
    # Holding primitive draws fixed, the same coded scale induces larger
    # level effects under meandif. Rescale both global and finite slab for
    # that basis to recover the orthonormal prior on the month effects.
    level_draws <- function(contrast, units) {
      term <- JAGS_formula(~ Month, "mu", dat, list(
        intercept = intercept, sh = prior_group(family, contrast, units)
      ))$prior_list$mu_Month
      set.seed(50931)
      unname(rng(term, 256L))
    }
    ortho <- level_draws("orthonormal", 1)
    testthat::expect_equal(level_draws("meandif", 1) * conversion, ortho, tolerance = 1e-10)
    testthat::expect_equal(level_draws("meandif", conversion), ortho, tolerance = 1e-10)
  }

  summarize <- function(draws, columns, C) {
    values <- lapply(columns, function(id) draws[, , id])
    names(values) <- names(columns)
    coefficients <- values[paste0("b", seq_len(ncol(C)))]
    for (i in seq_len(nrow(C))) {
      values[[paste0("month[", levels(dat$Month)[i], "]")]] <-
        Reduce(`+`, Map(`*`, coefficients, C[i, ]))
    }
    radius <- sqrt(Reduce(`+`, lapply(coefficients, `^`, 2)))
    values[["P(month radius < 0.75)"]] <- (radius < 0.75) * 1
    for (id in intersect(c("global", "month local", "wind local"), names(values))) {
      values[[paste0("log(", id, ")")]] <- log(values[[id]])
      values[[id]] <- NULL
    }
    t(vapply(values, function(x) {
      c(mean = mean(x), mcse = posterior::mcse_mean(x), rhat = posterior::rhat(x))
    }, numeric(3)))
  }
  for (contrast in names(bases)) {
    C <- bases[[contrast]]
    for (family in families) {
      key <- paste(contrast, family, sep = "_")
      jags <- scenario_fit(paste0("jags_", key), {
        fit <- JAGS_fit("model { for(i in 1:N) { y[i] ~ dnorm(mu[i], 1) } }",
          data = list(y = dat$Ozone, N = nrow(dat)), prior_list = list(),
          formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = dat),
          formula_prior_list = list(mu = list(intercept = intercept, sh = prior_group(family, contrast))),
          chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429
        )
        aperm(as.array(coda::as.mcmc.list(fit)), c(1, 3, 2))
      })
      idx <- seq_len(ncol(C))
      labels <- c("intercept", paste0("b", idx), "wind", "global")
      jags_columns <- setNames(c("mu_intercept", sprintf("mu_Month[%d]", idx), "mu_Wind", "mu_sh_global"), labels)
      stan_columns <- setNames(c("alpha", sprintf("beta[%d]", seq_len(n_months)), "tau"), labels)
      if (family %in% families[3:4]) {
        jags_columns <- c(jags_columns, "month local" = "mu_Month_local_1", "wind local" = "mu_Wind_local_1")
        stan_columns <- c(stan_columns, "month local" = "lambda[1]", "wind local" = "lambda[2]")
      }

      a <- summarize(reference[[key]]$draws, stan_columns, C)
      b <- summarize(jags, jags_columns, C)
      se <- sqrt(a[, "mcse"]^2 + b[, "mcse"]^2)
      z <- (b[, "mean"] - a[, "mean"]) / se
      table <- cbind(Stan = a[, "mean"], JAGS = b[, "mean"], MCSE = se, z = z,
                     Rhat = pmax(a[, "rhat"], b[, "rhat"]))
      scenario_text(paste0("comparison-", key), print(round(table, 4)))
      testthat::expect_equal(reference[[key]]$divergences, 0)
      testthat::expect_lt(max(abs(z)), 4.5)
      testthat::expect_lt(max(table[, "Rhat"]), 1.02)
    }
  }
})
