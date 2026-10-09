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
  near_upper <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(1 - 1e-8, K = 2)
  expect_true(all(is.finite(as.vector(near_lower))))
  expect_true(all(is.finite(as.vector(near_upper))))

  extreme_u <- rep(
    c(.Machine$double.eps, 1 - .Machine$double.eps),
    length.out = BayesTools:::.bt_lkj_cholesky_n_pairs(8L)
  )
  extreme_L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(extreme_u, K = 8L)
  extreme_R <- tcrossprod(extreme_L)
  expect_true(all(is.finite(extreme_L)))
  expect_true(all(diag(extreme_L) >= 0))
  expect_identical(extreme_R, t(extreme_R))
  expect_equal(diag(extreme_R), rep(1, 8L), tolerance = 32 * .Machine$double.eps)

  for(value in c(0, 1, NA_real_, Inf, -Inf)){
    expect_error(
      BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(value, K = 2),
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
  prior_values <- BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(
    matrix(c(0.25, 0.5, 0.75, 0.25, 0, 0.75), nrow = 2, byrow = TRUE),
    K = 3,
    eta = 1
  )
  expect_true(is.finite(prior_values[1]))
  expect_equal(prior_values[2], -Inf)
})

test_that("native LKJ transforms enforce open u support directly", {

  BayesTools:::.BayesTools_require_native_lkj()
  for(value in c(0, 1, NA_real_, Inf)){
    for(coordinate in c(1L, 3L)){
      invalid_vector <- rep(0.5, 3L)
      invalid_vector[coordinate] <- value
      invalid_matrix <- matrix(0.5, nrow = 3L, ncol = 3L)
      invalid_matrix[coordinate, coordinate] <- value
      for(u in list(invalid_vector, invalid_matrix)){
        expect_error(
          .Call("BayesTools_lkj_cholesky_from_u", u, 3L, PACKAGE = "BayesTools"),
          "'u' values must be finite and strictly between 0 and 1.",
          fixed = TRUE
        )
      }
    }
  }

  expect_equal(
    .Call("BayesTools_lkj_cholesky_from_u", numeric(0), 1L, PACKAGE = "BayesTools"),
    matrix(1, 1, 1)
  )
  expect_equal(
    .Call("BayesTools_lkj_cholesky_from_u", matrix(numeric(0), 3L, 0L), 1L,
          PACKAGE = "BayesTools"),
    array(1, c(3L, 1L, 1L))
  )
})

test_that("native LKJ helpers check geometry before allocating", {

  BayesTools:::.BayesTools_require_native_lkj()
  expect_error(
    .Call("BayesTools_lkj_cholesky_from_u", rep(.5, 32768L), 65537L,
          PACKAGE = "BayesTools"),
    "Native LKJ geometry is unavailable for 'K' = 65537.",
    fixed = TRUE
  )
  for(K in c(65537L, 100000L)){
    expect_error(
      .Call("BayesTools_lkj_alpha", K, 1.25, PACKAGE = "BayesTools"),
      paste0("Native LKJ geometry is unavailable for 'K' = ", K, "."),
      fixed = TRUE
    )
  }
  expect_error(
    .Call("BayesTools_lkj_cholesky_from_u", numeric(0), 100000L,
          PACKAGE = "BayesTools"),
    "Native LKJ geometry is unavailable for 'K' = 100000.",
    fixed = TRUE
  )

  empty_u <- structure(numeric(0), dim = c(0L, 2147450880L))
  empty_L <- .Call("BayesTools_lkj_cholesky_from_u", empty_u, 65536L,
                   PACKAGE = "BayesTools")
  expect_identical(dim(empty_L), c(0L, 65536L, 65536L))
  expect_length(empty_L, 0L)

  for(K in c(1L, 2L, 4L)){
    expected_alpha <- unlist(lapply(seq_len(K - 1L), function(row_offset){
      1.25 + (K - seq_len(row_offset) - 1L) / 2
    }), use.names = FALSE)
    if(K == 1L){
      expected_alpha <- numeric(0)
    }
    expect_identical(
      .Call("BayesTools_lkj_alpha", K, 1.25, PACKAGE = "BayesTools"),
      expected_alpha
    )
  }
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
    include_primitives = FALSE
  )

  expect_equal(module$pairs$alpha, c(1.5, 1.5, 1))
  expect_match(module$syntax, "Sigma_lkj_u[1:3] ~ dbt_lkj_cpc(Sigma_lkj_alpha)", fixed = TRUE)
  expect_error(
    BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(c(0.5), K = 2),
    "native LKJ routines requested",
    fixed = TRUE
  )
})

test_that("LKJ CPC construction returns valid lower Cholesky factors", {
  cpc <- c(0.2, -0.4, 0.5)
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L((cpc + 1) / 2, K = 3)

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

  R <- L %*% t(L)
  expect_equal(R, L %*% t(L), tolerance = 1e-12)
  expect_equal(diag(R), rep(1, 3), tolerance = 1e-12)
  expect_equal(R, t(R), tolerance = 1e-12)
  expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("LKJ primitives round-trip through the Cholesky factor", {

  set.seed(20260923)
  for(K in 2:5){
    n_pairs <- K * (K - 1L) / 2L
    alpha <- BayesTools:::.bt_lkj_cholesky_alpha(K = K, eta = 1.5)
    u <- sapply(alpha, function(a) stats::rbeta(25L, a, a))
    u <- matrix(u, nrow = 25L, ncol = n_pairs)
    L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u, K = K)
    expect_equal(
      BayesTools:::.bt_lkj_cholesky_L_to_cpc_u(L, K = K),
      u,
      tolerance = 1e-12
    )
    expect_equal(
      BayesTools:::.bt_lkj_cholesky_L_to_cpc_u(L[1L, , ], K = K),
      u[1L, ],
      tolerance = 1e-12
    )
    # A factor of an unscaled covariance is not a correlation Cholesky factor.
    expect_true(all(is.na(
      BayesTools:::.bt_lkj_cholesky_L_to_cpc_u(2 * L[1:2, , , drop = FALSE], K = K)
    )))
  }
})

test_that("LKJ primitives round-trip through the Cholesky factor near boundaries", {

  # Partial correlations close to +/-1 make 1 - sum_{k < j} L_ik^2 cancel
  # catastrophically for K >= 4 (round-trip errors up to 0.5 in u and NA rows);
  # the inverse must stay accurate to rounding error.
  set.seed(20260923)
  edges <- c(1e-12, 1e-10, 1e-8, 1e-6, 1 - 1e-6, 1 - 1e-8, 1 - 1e-10, 1 - 1e-12)
  for(K in 2:5){
    n_pairs <- K * (K - 1L) / 2L
    u <- rbind(
      matrix(rep(edges, each = n_pairs), ncol = n_pairs, byrow = TRUE),
      matrix(sample(c(edges, 0.2, 0.5, 0.9), 50L * n_pairs, replace = TRUE),
             ncol = n_pairs)
    )
    L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u, K = K)
    u_back <- BayesTools:::.bt_lkj_cholesky_L_to_cpc_u(L, K = K)
    expect_false(anyNA(u_back))
    expect_lt(max(abs(u_back - u)), 1e-12)
  }
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
    R <- L %*% t(L)

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
  R <- array(NA_real_, dim = dim(L))
  for(draw in seq_len(dim(L)[1L])){
    R[draw, , ] <- tcrossprod(L[draw, , ])
  }
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

.lkj_unit_diagonal_random_term <- function(){
  data <- data.frame(
    x = c(-1, 0.5, 1.2, -0.3, 0.8, -1.4),
    z = c(0.2, -0.7, 1.1, 0.4, -1.5, 0.9),
    id = factor(rep(c("a", "b", "c"), each = 2L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + z + (1 + x + z | id),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      z = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(
      sd = prior("normal", list(0, 1), list(0, Inf)),
      cor = prior_lkj(eta = 1)
    ))
  )
  formula_result$formula_design$random_effects[[1L]]
}

test_that("R-side LKJ correlation reconstructions have an exactly unit diagonal", {

  # The rounded row sums of squares of an LKJ Cholesky factor scatter within
  # 1 +/- a few eps. Correlation matrices rebuilt from the LKJ primitives or
  # their factor carry the exact unit diagonal of the module's bt_lkj_corr();
  # their off-diagonal entries stay the plain L L' products.
  random_term <- .lkj_unit_diagonal_random_term()
  K <- random_term$n_columns
  expect_identical(K, 3L)
  correlation <- random_term$correlation
  u_names <- correlation$primitive_names
  matrix_names <- function(name){
    as.vector(outer(seq_len(K), seq_len(K), Vectorize(function(row, column){
      paste0(name, "[", row, ",", column, "]")
    })))
  }
  L_names <- matrix_names(correlation$cholesky_name)
  R_names <- matrix_names(correlation$correlation_name)
  on_diagonal <- as.vector(diag(K) == 1)

  set.seed(20260924)
  n_draws <- 200L
  u <- matrix(
    stats::rbeta(n_draws * length(u_names), 1.2, 1.2),
    ncol = length(u_names),
    dimnames = list(NULL, u_names)
  )
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u, K = K)
  plain_R <- t(vapply(seq_len(n_draws), function(draw){
    as.vector(tcrossprod(L[draw, , ]))
  }, numeric(K * K)))
  # The draws exercise the rounding: plain row sums of squares miss 1.
  expect_true(any(plain_R[, on_diagonal] != 1))

  # Prior draws of the monitored correlation columns: the registered LKJ
  # deterministic node evaluated on the primitives.
  prior_samples <- cbind(u, BayesTools:::.bt_deterministic_node_evaluate(
    BayesTools:::.bt_dnode_lkj_from_random_term(random_term),
    BayesTools:::.bt_deterministic_lookup(u)
  ))
  expect_true(all(prior_samples[, R_names[on_diagonal]] == 1))
  expect_identical(
    unname(prior_samples[, R_names[!on_diagonal]]),
    unname(plain_R[, !on_diagonal])
  )

  # Source correlations of SD unscaling, from the primitives and from a
  # monitored Cholesky factor.
  group_key <- sub("^mu__xREx__", "", random_term$parameter_stem)
  for(posterior in list(u, prior_samples[, L_names, drop = FALSE])){
    source_R <- BayesTools:::.random_sd_correlation_draws(
      posterior = posterior,
      prefix = "mu",
      group_key = group_key,
      n_terms = K
    )
    source_R <- matrix(source_R, nrow = n_draws, ncol = K * K)
    expect_true(all(source_R[, on_diagonal] == 1))
    expect_identical(source_R[, !on_diagonal], unname(plain_R[, !on_diagonal]))
  }

  # Correlation columns completed for random-effect summaries.
  completed <- BayesTools:::.bt_random_effect_summary_complete_correlation_samples(
    random_term = random_term,
    model_samples = u
  )
  expect_true(all(completed[, R_names[on_diagonal]] == 1))
  expect_identical(
    unname(completed[, R_names[!on_diagonal]]),
    unname(plain_R[, !on_diagonal])
  )

  # Draws without a factor keep missing correlations.
  L[1L, , ] <- NA_real_
  R <- BayesTools:::.bt_lkj_cholesky_L_to_R(L)
  expect_true(all(is.na(R[1L, , ])))
  expect_true(all(R[-1L, 2L, 2L] == 1))
  expect_identical(
    BayesTools:::.bt_lkj_cholesky_L_to_R(L[2L, , ]),
    R[2L, , ]
  )
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

  L_1 <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u_1, K = K)
  L_2 <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u_2, K = K)
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
  alpha <- BayesTools:::.bt_lkj_cholesky_alpha(K = 2, eta = eta)
  u <- matrix(stats::rbeta(8000, alpha, alpha), ncol = 1L)
  draws <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(u, K = 2)
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

  module_lines <- strsplit(module$syntax, "\n", fixed = TRUE)[[1]]

  expect_lt(length(module_lines), 200L)
  expect_match(module$syntax, "Sigma_lkj_u[1:28] ~ dbt_lkj_cpc(Sigma_lkj_alpha)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_L_flat[1:64] <- bt_lkj_cholesky(Sigma_lkj_u, 8)", fixed = TRUE)
  expect_match(module$syntax, "Sigma_R_flat[1:64] <- bt_lkj_corr(Sigma_lkj_u, 8)", fixed = TRUE)
  expect_false(grepl("sqrt", module$syntax, fixed = TRUE))
  expect_false(grepl("inprod", module$syntax, fixed = TRUE))
})

test_that("one-dimensional LKJ-Cholesky module degenerates to identity", {
  module <- JAGS_lkj_corr_cholesky("One", K = 1, eta = 3)

  expect_equal(module$primitive_names, character(0))
  expect_equal(module$pairs$index, integer(0))
  expect_equal(module$backend, "module")
  expect_match(module$syntax, "One_L[1,1] <- 1", fixed = TRUE)
  expect_match(module$syntax, "One_R[1,1] <- 1", fixed = TRUE)
  expect_false(grepl("dbt_lkj_cpc", module$syntax, fixed = TRUE))

  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(numeric(0), K = 1)
  expect_equal(L, matrix(1, 1, 1))
})
