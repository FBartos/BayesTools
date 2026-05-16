skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS LKJ-Cholesky Fit Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile MCMC checks for the package-shipped JAGS LKJ-Cholesky module.
#
# TAGS: @fit, @JAGS, @Stan, @LKJ, @Cholesky
# ============================================================================ #

skip_on_cran()
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

.fit_jags_lkj_cholesky_prior <- function(K, eta, sample = 3000, seed = 1L,
                                          backend = "module") {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = K,
    eta = eta,
    include_correlation = TRUE,
    include_primitives = TRUE,
    backend = backend
  )

  user_silent.jags <- runjags::runjags.getOption("silent.jags")
  user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  fit <- suppressWarnings(runjags::run.jags(
    model = paste0("model{\n", module$syntax, "\n}"),
    monitor = module$monitor,
    n.chains = 2,
    adapt = 500,
    burnin = 500,
    sample = sample,
    thin = 1,
    method = "rjags",
    summarise = FALSE,
    plots = FALSE,
    inits = lapply(seq_len(2), function(i){
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  list(
    fit = fit,
    module = module,
    samples = as.matrix(fit$mcmc)
  )
}

.fit_jags_lkj_cholesky_prior_cached <- local({
  cache <- new.env(parent = emptyenv())

  function(K, eta, sample = 3000, seed = 1L, backend = "module") {
    key <- paste(K, eta, sample, seed, backend, sep = "|")
    if(!exists(key, envir = cache, inherits = FALSE)){
      assign(
        key,
        .fit_jags_lkj_cholesky_prior(
          K = K,
          eta = eta,
          sample = sample,
          seed = seed,
          backend = backend
        ),
        envir = cache
      )
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

.eval_jags_lkj_cholesky_transform <- function(u, K) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  K2 <- K * K
  syntax <- c(
    paste0("  Omega_L_flat[1:", K2, "] <- bt_lkj_cholesky(u, ", K, ")"),
    paste0("  Omega_R_flat[1:", K2, "] <- bt_lkj_corr(u, ", K, ")")
  )
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      flat_index <- BayesTools:::.bt_lkj_cholesky_flat_index(i, j, K)
      syntax <- c(
        syntax,
        paste0("  Omega_L[", i, ",", j, "] <- Omega_L_flat[", flat_index, "]"),
        paste0("  Omega_R[", i, ",", j, "] <- Omega_R_flat[", flat_index, "]")
      )
    }
  }

  con <- textConnection(paste0("model{\n", paste(syntax, collapse = "\n"), "\n}\n"))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u),
    n.chains = 1,
    n.adapt = 0,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = c("Omega_L", "Omega_R"),
    n.iter = 1,
    quiet = TRUE,
    progress.bar = "none"
  )

  as.matrix(samples)
}

.jags_lkj_vector_column <- function(samples, base_name, index) {
  indexed_name <- paste0(base_name, "[", index, "]")
  if(indexed_name %in% colnames(samples)){
    return(indexed_name)
  }
  if(index == 1L && base_name %in% colnames(samples)){
    return(base_name)
  }
  stop("Missing monitored JAGS column: ", indexed_name, call. = FALSE)
}

.jags_lkj_matrix_draw <- function(samples, row, prefix, K) {
  out <- matrix(NA_real_, K, K)
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      out[i, j] <- samples[row, paste0(prefix, "[", i, ",", j, "]")]
    }
  }
  out
}

.jags_lkj_offdiag_columns <- function(samples, prefix, K) {
  out <- character(0)
  for(row in 2:K){
    for(column in seq_len(row - 1L)){
      out <- c(out, paste0(prefix, "[", column, ",", row, "]"))
    }
  }
  if(!all(out %in% colnames(samples))){
    stop("Missing monitored JAGS correlation columns.", call. = FALSE)
  }

  out
}

.eval_jags_lkj_cpc_deviance <- function(u, alpha) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  con <- textConnection(paste0(
    "model{\n",
    "  u[1:", length(u), "] ~ dbt_lkj_cpc(alpha)\n",
    "}\n"
  ))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u, alpha = alpha),
    n.chains = 2,
    n.adapt = 0,
    quiet = TRUE
  )
  dic <- rjags::dic.samples(
    model = model,
    n.iter = 1,
    type = "pD",
    quiet = TRUE,
    progress.bar = "none"
  )

  list(
    deviance = as.numeric(dic[["deviance"]]),
    penalty = as.numeric(dic[["penalty"]])
  )
}

.lkj_bivariate_data <- function(N = 28L, rho = 0.55, seed = 4711L) {
  set.seed(seed)
  z_1 <- stats::rnorm(N)
  z_2 <- rho * z_1 + sqrt(1 - rho^2) * stats::rnorm(N)
  cbind(z_1, z_2)
}

.lkj_bivariate_loglik <- function(rho, y) {
  denom <- 1 - rho^2
  quad <- (y[, 1]^2 - 2 * rho * y[, 1] * y[, 2] + y[, 2]^2) / denom
  sum(-0.5 * log(denom) - 0.5 * quad)
}

.lkj_bivariate_grid_oracle <- function(y, eta, grid_size = 20001L) {
  rho <- seq(-.999, .999, length.out = grid_size)
  log_density <- vapply(rho, .lkj_bivariate_loglik, numeric(1), y = y) +
    (eta - 1) * log1p(-rho^2)
  log_density <- log_density - max(log_density)
  weights <- exp(log_density)
  weights <- weights / sum(weights)
  cdf <- cumsum(weights)
  posterior_mean <- sum(weights * rho)

  list(
    mean = posterior_mean,
    sd = sqrt(sum(weights * (rho - posterior_mean)^2)),
    quantiles = stats::approx(cdf, rho, xout = c(.1, .5, .9), ties = "ordered")$y
  )
}

.fit_jags_lkj_bivariate_posterior <- function(y, eta, sample = 5000,
                                              seed = 812L) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  model <- paste0(
    "model{\n",
    "  alpha[1] <- eta\n",
    "  u[1:1] ~ dbt_lkj_cpc(alpha)\n",
    "  R_flat[1:4] <- bt_lkj_corr(u, 2)\n",
    "  rho <- R_flat[2]\n",
    "  for(i in 1:N){\n",
    "    loglik[i] <- -0.5 * log(1 - pow(rho, 2)) - ",
    "(pow(y[i,1], 2) - 2 * rho * y[i,1] * y[i,2] + pow(y[i,2], 2)) / ",
    "(2 * (1 - pow(rho, 2)))\n",
    "    zeros[i] ~ dpois(C - loglik[i])\n",
    "  }\n",
    "}\n"
  )

  user_silent.jags <- runjags::runjags.getOption("silent.jags")
  user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  fit <- suppressWarnings(runjags::run.jags(
    model = model,
    data = list(y = y, eta = eta, N = nrow(y), zeros = rep(0, nrow(y)), C = 100),
    monitor = "rho",
    n.chains = 2,
    adapt = 1000,
    burnin = 1000,
    sample = sample,
    thin = 1,
    method = "rjags",
    summarise = FALSE,
    plots = FALSE,
    inits = lapply(seq_len(2), function(i){
      list(u = c(0.5), .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  as.numeric(as.matrix(fit$mcmc)[, "rho"])
}

test_that("JAGS LKJ-Cholesky module matches primitive and correlation marginals", {
  settings <- data.frame(
    K = c(2L, 3L, 4L, 5L),
    eta = c(0.35, 1, 2.5, 0.75),
    sample = c(3000L, 3000L, 3000L, 4000L)
  )

  for(setting_i in seq_len(nrow(settings))){
    K <- settings$K[setting_i]
    eta <- settings$eta[setting_i]
    result <- .fit_jags_lkj_cholesky_prior_cached(
      K = K,
      eta = eta,
      sample = settings$sample[setting_i],
      seed = 100 + setting_i
    )
    samples <- result$samples
    pairs <- result$module$pairs

    u_columns <- character(nrow(pairs))
    for(p in seq_len(nrow(pairs))){
      u_name <- .jags_lkj_vector_column(samples, "Omega_lkj_u", p)
      u_columns[p] <- u_name
      alpha <- pairs$alpha[p]
      expect_equal(mean(samples[, u_name]), 0.5, tolerance = 0.05)
      expect_true(abs(stats::var(samples[, u_name]) - 1 / (4 * (2 * alpha + 1))) < 0.012)
      expect_equal(
        stats::quantile(samples[, u_name], c(.1, .5, .9), names = FALSE),
        stats::qbeta(c(.1, .5, .9), alpha, alpha),
        tolerance = 0.06
      )
    }

    if(length(u_columns) >= 2L){
      primitive_cor <- stats::cor(samples[, u_columns, drop = FALSE])
      expect_lt(max(abs(primitive_cor[upper.tri(primitive_cor)])), 0.10)
    }

    rho_columns <- .jags_lkj_offdiag_columns(samples, "Omega_R", K)
    alpha_rho <- eta - 1 + K / 2
    rho_quantiles <- 2 * stats::qbeta(c(.1, .5, .9), alpha_rho, alpha_rho) - 1
    rho_means <- vapply(rho_columns, function(name) mean(samples[, name]), numeric(1))
    rho_sds <- vapply(rho_columns, function(name) stats::sd(samples[, name]), numeric(1))
    for(rho_name in rho_columns){
      expect_equal(mean(samples[, rho_name]), 0, tolerance = 0.045)
      expect_true(abs(stats::var(samples[, rho_name]) - 1 / (2 * alpha_rho + 1)) < 0.025)
      expect_equal(
        stats::quantile(samples[, rho_name], c(.1, .5, .9), names = FALSE),
        rho_quantiles,
        tolerance = 0.08
      )
    }
    expect_lt(max(rho_means) - min(rho_means), 0.08)
    expect_lt(max(rho_sds) - min(rho_sds), 0.08)

    check_rows <- unique(round(seq(1, nrow(samples), length.out = min(25, nrow(samples)))))
    for(row in check_rows){
      L <- .jags_lkj_matrix_draw(samples, row, "Omega_L", K)
      R <- .jags_lkj_matrix_draw(samples, row, "Omega_R", K)

      expect_equal(L[upper.tri(L)], rep(0, K * (K - 1) / 2), tolerance = 1e-10)
      expect_true(all(diag(L) > 0))
      expect_equal(rowSums(L^2), rep(1, K), tolerance = 1e-8)
      expect_equal(R, L %*% t(L), tolerance = 1e-8)
      expect_equal(diag(R), rep(1, K), tolerance = 1e-8)
      expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
    }
  }
})

test_that("compiled LKJ-Cholesky functions match hand-coded transform oracles", {
  settings <- list(
    list(
      K = 2L,
      u = c(0.31),
      expected_L = matrix(c(
        1, 0,
        -0.38, sqrt(1 - .38^2)
      ), nrow = 2, byrow = TRUE)
    ),
    list(
      K = 3L,
      u = c(0.61, 0.22, 0.74),
      expected_L = matrix(c(
        1, 0, 0,
        .22, sqrt(1 - .22^2), 0,
        -.56, .48 * sqrt(1 - .56^2), sqrt(1 - .56^2) * sqrt(1 - .48^2)
      ), nrow = 3, byrow = TRUE)
    )
  )

  for(setting in settings){
    samples <- .eval_jags_lkj_cholesky_transform(setting$u, setting$K)
    expected_L <- setting$expected_L
    expected_R <- expected_L %*% t(expected_L)
    observed_L <- .jags_lkj_matrix_draw(samples, 1, "Omega_L", setting$K)
    observed_R <- .jags_lkj_matrix_draw(samples, 1, "Omega_R", setting$K)

    expect_equal(observed_L, expected_L, tolerance = 1e-12)
    expect_equal(observed_R, expected_R, tolerance = 1e-12)
    expect_equal(observed_L[upper.tri(observed_L)], rep(0, setting$K * (setting$K - 1L) / 2L), tolerance = 1e-12)
    expect_true(all(diag(observed_L) > 0))
    expect_equal(rowSums(observed_L^2), rep(1, setting$K), tolerance = 1e-12)
    expect_equal(diag(observed_R), rep(1, setting$K), tolerance = 1e-12)
  }
})

test_that("compiled LKJ CPC distribution logDensity matches beta density for observed nodes", {
  K <- 4L
  eta <- 0.75
  u <- c(0.08, 0.91, 0.47, 0.12, 0.63, 0.88)
  alpha <- BayesTools:::.bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)$alpha
  deviance <- .eval_jags_lkj_cpc_deviance(u = u, alpha = alpha)
  expected_deviance <- -2 * sum(stats::dbeta(u, alpha, alpha, log = TRUE))

  expect_equal(deviance$deviance, expected_deviance, tolerance = 1e-12)
  expect_equal(deviance$penalty, 0, tolerance = 1e-12)
})

test_that("JAGS LKJ posterior for bivariate normal correlation matches grid oracle", {
  eta <- 0.8
  y <- .lkj_bivariate_data()
  oracle <- .lkj_bivariate_grid_oracle(y = y, eta = eta)
  rho <- .fit_jags_lkj_bivariate_posterior(y = y, eta = eta, sample = 5000)

  expect_equal(mean(rho), oracle$mean, tolerance = 0.06)
  expect_equal(stats::sd(rho), oracle$sd, tolerance = 0.06)
  expect_equal(
    stats::quantile(rho, c(.1, .5, .9), names = FALSE),
    oracle$quantiles,
    tolerance = 0.09
  )
})

test_that("JAGS_fit loads BayesTools module from generated LKJ metadata", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = 2,
    eta = 1.1,
    include_correlation = TRUE
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 10
  ))

  expect_s3_class(fit, "BayesTools_fit")
  expect_equal(attr(fit, "jags_modules"), "BayesTools")
  samples <- as.matrix(fit$mcmc)
  expect_true(all(c("Omega_L[1,1]", "Omega_R[1,2]") %in% colnames(samples)))
})

test_that("one-dimensional generated LKJ module runs without native LKJ calls", {
  module <- JAGS_lkj_corr_cholesky(
    name = "One",
    K = 1,
    eta = 3,
    include_correlation = TRUE
  )

  expect_false(grepl("dbt_lkj_cpc", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_cholesky", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_corr", module$syntax, fixed = TRUE))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 11
  ))

  expect_s3_class(fit, "BayesTools_fit")
  samples <- as.matrix(fit$mcmc)
  expect_equal(colnames(samples), c("One_L", "One_R"))
  expect_equal(unname(samples[, "One_L"]), rep(1, nrow(samples)))
  expect_equal(unname(samples[, "One_R"]), rep(1, nrow(samples)))
})
