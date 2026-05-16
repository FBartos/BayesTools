skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS LKJ-Cholesky Module
# ============================================================================ #
#
# PURPOSE:
#   Deterministic tests for the package-shipped JAGS LKJ-Cholesky syntax module.
#
# TAGS: @unit, @JAGS, @LKJ, @Cholesky
# ============================================================================ #

test_that("JAGS_lkj_corr_cholesky validates inputs and exposes metadata", {
  expect_error(JAGS_lkj_corr_cholesky("bad-name", 2), "valid JAGS node prefix")
  expect_error(JAGS_lkj_corr_cholesky("Omega", 0), "K")
  expect_error(JAGS_lkj_corr_cholesky("Omega", 2, eta = 0), "eta")
  expect_error(JAGS_lkj_corr_cholesky("Omega", 2, include_correlation = NA), "include_correlation")
  expect_error(JAGS_lkj_corr_cholesky("Omega", 2, backend = "bad"), "one of")

  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = 3,
    eta = 1.5,
    include_correlation = TRUE,
    include_primitives = TRUE
  )

  expect_s3_class(module, "BayesTools_JAGS_lkj_corr_cholesky")
  expect_equal(module$K, 3)
  expect_equal(module$eta, 1.5)
  expect_equal(module$backend, "module")
  expect_equal(module$jags_module, "BayesTools")
  expect_equal(module$required_packages, "BayesTools")
  expect_equal(module$cholesky_name, "Omega_L")
  expect_equal(module$correlation_name, "Omega_R")
  expect_equal(module$primitive_names, paste0("Omega_lkj_u[", 1:3, "]"))
  expect_equal(module$primitive_bounds$lb, stats::setNames(rep(0, 3), module$primitive_names))
  expect_equal(module$primitive_bounds$ub, stats::setNames(rep(1, 3), module$primitive_names))
  expect_equal(module$monitor, c("Omega_L", "Omega_R", paste0("Omega_lkj_u[", 1:3, "]"), paste0("Omega_lkj_cpc[", 1:3, "]")))

  expect_equal(module$pairs$i, c(1L, 1L, 2L))
  expect_equal(module$pairs$j, c(2L, 3L, 3L))
  expect_equal(module$pairs$alpha, c(2, 2, 1.5))
})

test_that("LKJ primitive helpers enforce open u support before native transforms", {

  near_lower <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(1e-8, K = 2)
  near_upper <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_R(1 - 1e-8, K = 2)
  expect_true(all(is.finite(as.vector(near_lower))))
  expect_true(all(is.finite(as.vector(near_upper))))

  for(value in c(0, 1, NA_real_, Inf, -Inf)){
    expect_error(
      BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(value, K = 2),
      "strictly between 0 and 1",
      fixed = TRUE
    )
    expect_error(
      BayesTools:::.bt_lkj_cholesky_cpc_u_to_R(value, K = 2),
      "strictly between 0 and 1",
      fixed = TRUE
    )
  }

  invalid_matrix <- matrix(c(0.25, 0, 0.75), nrow = 1)
  expect_error(
    BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(invalid_matrix, K = 3),
    "strictly between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_lkj_cholesky_cpc_u_to_R(invalid_matrix, K = 3),
    "strictly between 0 and 1",
    fixed = TRUE
  )

  prior_values <- BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(
    matrix(c(0.25, 0.5, 0.75, 0.25, 0, 0.75), nrow = 2, byrow = TRUE),
    K = 3,
    eta = 1
  )
  expect_true(is.finite(prior_values[1]))
  expect_equal(prior_values[2], -Inf)
})

test_that("LKJ alpha generation is available without native routines", {
  testthat::local_mocked_bindings(
    .BayesTools_require_native_lkj = function(){
      stop("native LKJ routines requested", call. = FALSE)
    },
    .package = "BayesTools"
  )

  expect_equal(
    BayesTools:::.bt_lkj_cholesky_alpha(K = 5, eta = 0.75),
    c(2.25, 2.25, 1.75, 2.25, 1.75, 1.25, 2.25, 1.75, 1.25, 0.75)
  )

  module <- JAGS_lkj_corr_cholesky(
    name = "Sigma",
    K = 3,
    eta = 1,
    include_correlation = TRUE,
    include_primitives = FALSE,
    backend = "syntax"
  )

  expect_equal(module$pairs$alpha, c(1.5, 1.5, 1))
  expect_match(module$syntax, "Sigma_lkj_u[1] ~ dbeta(1.5, 1.5)", fixed = TRUE)
  expect_error(
    BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(c(0.5), K = 2),
    "native LKJ routines requested",
    fixed = TRUE
  )
})

test_that("LKJ CPC construction returns valid lower Cholesky factors", {
  cpc <- c(0.2, -0.4, 0.5)
  L <- BayesTools:::.bt_lkj_cholesky_cpc_to_L(cpc, K = 3)

  expected <- matrix(
    c(
      1, 0, 0,
      0.2, sqrt(1 - 0.2^2), 0,
      -0.4, 0.5 * sqrt(1 - (-0.4)^2), sqrt(1 - (-0.4)^2) * sqrt(1 - 0.5^2)
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(L, expected, tolerance = 1e-12)
  expect_equal(L[upper.tri(L)], rep(0, 3), tolerance = 1e-12)
  expect_true(all(diag(L) > 0))
  expect_equal(rowSums(L^2), rep(1, 3), tolerance = 1e-12)

  R <- BayesTools:::.bt_lkj_cholesky_corr(L)
  expect_equal(R, L %*% t(L), tolerance = 1e-12)
  expect_equal(diag(R), rep(1, 3), tolerance = 1e-12)
  expect_equal(R, t(R), tolerance = 1e-12)
  expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("LKJ CPC construction remains stable near primitive boundaries", {
  settings <- list(
    list(K = 2L, u = c(1e-8)),
    list(K = 2L, u = c(1 - 1e-8)),
    list(K = 4L, u = c(1e-6, .999999, .5, .01, .99, .37)),
    list(K = 5L, u = c(.13, .87, .49, .02, .98, .31, .69, .42, .58, .76))
  )

  for(setting in settings){
    L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(setting$u, K = setting$K)
    R <- BayesTools:::.bt_lkj_cholesky_corr(L)

    expect_true(all(is.finite(L)))
    expect_equal(L[upper.tri(L)], rep(0, setting$K * (setting$K - 1L) / 2L), tolerance = 1e-12)
    expect_true(all(diag(L) > 0))
    expect_equal(rowSums(L^2), rep(1, setting$K), tolerance = 1e-10)
    expect_equal(R, L %*% t(L), tolerance = 1e-10)
    expect_equal(R, t(R), tolerance = 1e-10)
    expect_equal(diag(R), rep(1, setting$K), tolerance = 1e-10)
    expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > -1e-10))
  }
})

test_that("native LKJ helpers handle vectorized R-side transforms", {
  u <- rbind(
    c(0.61, 0.22, 0.74),
    c(0.42, 0.81, 0.35)
  )
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u, K = 3)
  R <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_R(u, K = 3)
  alpha <- BayesTools:::.bt_lkj_cholesky_alpha(K = 3, eta = 1.4)
  log_prior <- BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(u, K = 3, eta = 1.4)

  expected_L_1 <- matrix(c(
    1, 0, 0,
    .22, sqrt(1 - .22^2), 0,
    -.56, .48 * sqrt(1 - .56^2), sqrt(1 - .56^2) * sqrt(1 - .48^2)
  ), nrow = 3, byrow = TRUE)

  expect_equal(dim(L), c(2L, 3L, 3L))
  expect_equal(dim(R), c(2L, 3L, 3L))
  expect_equal(L[1, , ], expected_L_1, tolerance = 1e-12)
  expect_equal(R[1, , ], expected_L_1 %*% t(expected_L_1), tolerance = 1e-12)
  expect_equal(log_prior, rowSums(matrix(
    stats::dbeta(as.vector(u), rep(alpha, each = nrow(u)), rep(alpha, each = nrow(u)), log = TRUE),
    nrow = nrow(u),
    ncol = ncol(u)
  )), tolerance = 1e-12)
})

test_that("LKJ CPC beta density induces Stan Cholesky density kernel", {
  K <- 5L
  eta <- 0.7
  u_1 <- c(.13, .87, .49, .02, .98, .31, .69, .42, .58, .76)
  u_2 <- c(.73, .27, .51, .91, .09, .64, .36, .56, .44, .19)
  cpc_1 <- 2 * u_1 - 1
  cpc_2 <- 2 * u_2 - 1
  pairs <- BayesTools:::.bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)

  log_row_jacobian <- function(cpc) {
    out <- 0
    for(p in seq_len(nrow(pairs))){
      row <- pairs$j[p]
      column <- pairs$i[p]
      exponent <- (row - 1L - column) / 2
      if(exponent > 0){
        out <- out + exponent * log1p(-cpc[p]^2)
      }
    }
    out
  }

  L_1 <- BayesTools:::.bt_lkj_cholesky_cpc_to_L(cpc_1, K = K)
  L_2 <- BayesTools:::.bt_lkj_cholesky_cpc_to_L(cpc_2, K = K)
  induced_log_density_difference <-
    BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(u_1, K = K, eta = eta) -
    log_row_jacobian(cpc_1) -
    BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(u_2, K = K, eta = eta) +
    log_row_jacobian(cpc_2)
  stan_kernel_difference <- sum(
    (K - seq_len(K) + 2 * eta - 2) *
      (log(diag(L_1)) - log(diag(L_2)))
  )

  expect_equal(induced_log_density_difference, stan_kernel_difference, tolerance = 1e-10)
})

test_that("LKJ primitive beta construction has the expected marginal behavior", {
  set.seed(17)
  eta <- 2
  draws <- BayesTools:::.bt_lkj_cholesky_rng(n = 8000, K = 2, eta = eta)
  rho <- draws[, 2, 1]

  expect_equal(mean(rho), 0, tolerance = 0.03)
  expect_equal(stats::var(rho), 1 / (2 * eta + 1), tolerance = 0.03)
  expect_true(all(abs(rho) < 1))
})

test_that("LKJ primitive log prior uses beta stochastic coordinates", {
  u <- c(0.6, 0.25, 0.8, 0.4, 0.7, 0.55)
  pairs <- BayesTools:::.bt_lkj_cholesky_cpc_pairs(K = 4, eta = 0.75)

  expect_equal(
    BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(u, K = 4, eta = 0.75),
    sum(stats::dbeta(u, pairs$alpha, pairs$alpha, log = TRUE)),
    tolerance = 1e-12
  )

  u_matrix <- rbind(u, rev(u))
  expect_equal(
    BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(u_matrix, K = 4, eta = 0.75),
    c(
      sum(stats::dbeta(u, pairs$alpha, pairs$alpha, log = TRUE)),
      sum(stats::dbeta(rev(u), pairs$alpha, pairs$alpha, log = TRUE))
    ),
    tolerance = 1e-12
  )
})

test_that("JAGS LKJ-Cholesky module backend uses compiled distribution and transformations", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Sigma",
    K = 3,
    eta = 1,
    include_correlation = TRUE,
    include_primitives = FALSE
  )

  expect_equal(module$backend, "module")
  expect_match(module$syntax, "Sigma_lkj_alpha[1] <- 1.5", fixed = TRUE)
  expect_match(module$syntax, "Sigma_lkj_alpha[3] <- 1", fixed = TRUE)
  expect_match(module$syntax, "Sigma_lkj_u[1:3] ~ dbt_lkj_cpc(Sigma_lkj_alpha)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_L_flat[1:9] <- bt_lkj_cholesky(Sigma_lkj_u, 3)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R_flat[1:9] <- bt_lkj_corr(Sigma_lkj_u, 3)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_L[3,2] <- Sigma_L_flat[8]", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R[1,3] <- Sigma_R_flat[3]", fixed = TRUE)

  expect_false(grepl("dbeta", module$syntax, fixed = TRUE))
  expect_false(grepl("sqrt", module$syntax, fixed = TRUE))
  expect_false(grepl("pow", module$syntax, fixed = TRUE))
  expect_false(grepl("inprod", module$syntax, fixed = TRUE))
  expect_false(grepl("dmnorm", module$syntax, fixed = TRUE))
  expect_false(grepl("dwish", module$syntax, fixed = TRUE))
})

test_that("compiled backend keeps generated syntax compact for larger dimensions", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Sigma",
    K = 8,
    eta = 1.25,
    include_correlation = TRUE,
    include_primitives = FALSE
  )
  fallback <- JAGS_lkj_corr_cholesky(
    name = "Sigma",
    K = 8,
    eta = 1.25,
    include_correlation = TRUE,
    include_primitives = FALSE,
    backend = "syntax"
  )

  module_lines <- strsplit(module$syntax, "\n", fixed = TRUE)[[1]]
  fallback_lines <- strsplit(fallback$syntax, "\n", fixed = TRUE)[[1]]

  expect_lt(length(module_lines), length(fallback_lines))
  expect_match(module$syntax, "Sigma_lkj_u[1:28] ~ dbt_lkj_cpc(Sigma_lkj_alpha)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_L_flat[1:64] <- bt_lkj_cholesky(Sigma_lkj_u, 8)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R_flat[1:64] <- bt_lkj_corr(Sigma_lkj_u, 8)", fixed = TRUE)
  expect_false(grepl("sqrt", module$syntax, fixed = TRUE))
  expect_false(grepl("inprod", module$syntax, fixed = TRUE))
  expect_true(grepl("sqrt", fallback$syntax, fixed = TRUE))
  expect_true(grepl("inprod", fallback$syntax, fixed = TRUE))
})

test_that("JAGS LKJ-Cholesky syntax fallback is unrolled and deterministic after primitives", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Sigma",
    K = 3,
    eta = 1,
    include_correlation = TRUE,
    include_primitives = FALSE,
    backend = "syntax"
  )

  expect_equal(module$backend, "syntax")
  expect_match(module$syntax, "Sigma_lkj_u[1] ~ dbeta(1.5, 1.5)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_lkj_u[3] ~ dbeta(1, 1)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_L[3,2] <- Sigma_lkj_cpc[3] * sqrt(1 - pow(Sigma_lkj_cpc[2], 2))", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R[3,1] <- inprod(Sigma_L[3,1:3], Sigma_L[1,1:3])", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R[1,3] <- Sigma_R[3,1]", fixed = TRUE)

  expect_false(grepl("dmnorm", module$syntax, fixed = TRUE))
  expect_false(grepl("dwish", module$syntax, fixed = TRUE))
  expect_false(grepl("for(", module$syntax, fixed = TRUE))
  expect_false(grepl("for ", module$syntax, fixed = TRUE))
})

test_that("one-dimensional LKJ-Cholesky module degenerates to identity", {
  module <- JAGS_lkj_corr_cholesky("One", K = 1, eta = 3)

  expect_equal(module$primitive_names, character(0))
  expect_equal(module$pairs$index, integer(0))
  expect_equal(module$backend, "module")
  expect_match(module$syntax, "One_L[1,1] <- 1", fixed = TRUE)
  expect_match(module$syntax, "One_R[1,1] <- 1", fixed = TRUE)
  expect_false(grepl("dbt_lkj_cpc", module$syntax, fixed = TRUE))

  L <- BayesTools:::.bt_lkj_cholesky_cpc_to_L(numeric(0), K = 1)
  expect_equal(L, matrix(1, 1, 1))
})
