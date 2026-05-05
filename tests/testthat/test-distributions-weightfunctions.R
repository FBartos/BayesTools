skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Distributions - Weight Functions
# ============================================================================ #
#
# PURPOSE:
#   Contract tests for marginal distributions and random generation for
#   reference-first weight functions.
#
# TAGS: @evaluation, @distributions, @weightfunctions
# ============================================================================ #

.cumdirichlet_beta_shapes <- function(alpha){
  J <- length(alpha)
  lapply(seq_len(J), function(j){
    if(j == 1L){
      return(list(type = "point", location = 1))
    }
    list(
      type   = "beta",
      shape1 = sum(alpha[j:J]),
      shape2 = sum(alpha[seq_len(j - 1L)])
    )
  })
}

.expect_monotone_marginals <- function(alpha, q_seq){
  shapes <- .cumdirichlet_beta_shapes(alpha)

  expected_density <- do.call(cbind, lapply(shapes, function(shape){
    if(shape$type == "point"){
      dpoint(q_seq, location = shape$location)
    }else{
      stats::dbeta(q_seq, shape$shape1, shape$shape2)
    }
  }))
  expected_cdf <- do.call(cbind, lapply(shapes, function(shape){
    if(shape$type == "point"){
      ppoint(q_seq, location = shape$location)
    }else{
      stats::pbeta(q_seq, shape$shape1, shape$shape2)
    }
  }))
  expected_quantile <- do.call(cbind, lapply(shapes, function(shape){
    if(shape$type == "point"){
      qpoint(q_seq, location = shape$location)
    }else{
      stats::qbeta(q_seq, shape$shape1, shape$shape2)
    }
  }))
  expected_quantile <- ifelse(expected_quantile == -Inf, 0, expected_quantile)
  expected_quantile <- ifelse(expected_quantile ==  Inf, 1, expected_quantile)

  expect_equal(mdone.sided(q_seq, alpha = alpha), expected_density)
  expect_equal(mdtwo.sided(q_seq, alpha = alpha), expected_density)
  expect_equal(mpone.sided(q_seq, alpha = alpha), expected_cdf)
  expect_equal(mptwo.sided(q_seq, alpha = alpha), expected_cdf)
  expect_equal(mqone.sided(q_seq, alpha = alpha), expected_quantile)
  expect_equal(mqtwo.sided(q_seq, alpha = alpha), expected_quantile)
}

test_that("monotone weightfunction marginal helpers use reference-first cumulative Dirichlet math", {

  q_seq <- seq(0, 1, .1)

  .expect_monotone_marginals(c(1, 1), q_seq)
  .expect_monotone_marginals(c(.3, 1, 2), q_seq)
  .expect_monotone_marginals(c(4, 3, 5, 2), q_seq)

  expect_equal(
    mdone.sided(c(0, 1), alpha = c(1, 1), log = TRUE),
    log(matrix(c(0, Inf, 1, 1), ncol = 2))
  )
  expect_equal(
    mpone.sided(q_seq, alpha = c(1, 1), log.p = TRUE),
    unname(log(cbind(c(rep(0, length(q_seq) - 1), 1), q_seq)))
  )
  expect_equal(
    mpone.sided(q_seq, alpha = c(1, 1), lower.tail = FALSE),
    cbind(1 - c(rep(0, length(q_seq) - 1), 1), 1 - q_seq)
  )
  expect_equal(
    mqone.sided(q_seq, alpha = c(1, 1), lower.tail = FALSE),
    cbind(rep(1, length(q_seq)), stats::qbeta(q_seq, 1, 1, lower.tail = FALSE))
  )
})

test_that("monotone weightfunction RNG fixes the first reference bin", {

  set.seed(1)
  samples <- rone.sided(10000, alpha = c(1, 1))
  expect_true(all(samples[,1] == 1))
  expect_equal(mean(samples[,2]), .5, tolerance = .02)

  set.seed(2)
  samples <- rone.sided(10000, alpha = c(.3, 1, 2))
  expect_true(all(samples[,1] == 1))
  expect_equal(mean(samples[,2]), 3 / 3.3, tolerance = .02)
  expect_equal(mean(samples[,3]), 2 / 3.3, tolerance = .02)

  set.seed(3)
  samples1 <- rone.sided(10, alpha = c(3, 5, .4))
  set.seed(3)
  samples2 <- rtwo.sided(10, alpha = c(3, 5, .4))
  expect_equal(samples1, samples2)

  set.seed(4)
  general <- rone.sided(1000, alpha1 = c(1, 1), alpha2 = c(1, 1))
  expect_equal(general[,1], rep(1, nrow(general)), tolerance = 1e-12)
  expect_true(all(general >= -1e-12 & general <= 1 + 1e-12))
})

test_that("fixed weightfunction helpers support non-negative relative weights", {

  omega <- c(1, .5, 1.5)
  q <- c(0, .5, 1, 1.5)

  expect_equal(
    mdone.sided_fixed(q, omega = omega),
    cbind(
      dpoint(q, location = 1),
      dpoint(q, location = .5),
      dpoint(q, location = 1.5)
    )
  )
  expect_equal(
    mpone.sided_fixed(q, omega = omega),
    cbind(
      ppoint(q, location = 1),
      ppoint(q, location = .5),
      ppoint(q, location = 1.5)
    )
  )
  expect_equal(
    mqone.sided_fixed(c(0, .5, 1), omega = omega),
    matrix(c(0, 1, 1, 0, .5, .5, 0, 1.5, 1.5), ncol = 3)
  )
  expect_equal(
    mqone.sided_fixed(c(0, .5, 1), omega = omega, lower.tail = FALSE),
    matrix(c(1, 1, 1, .5, .5, .5, 1.5, 1.5, 1.5), ncol = 3)
  )
  expect_equal(
    rone.sided_fixed(2, omega = omega),
    matrix(c(1, 1, .5, .5, 1.5, 1.5), nrow = 2)
  )
  expect_equal(mdone.sided_fixed(q, omega = omega), mdtwo.sided_fixed(q, omega = omega))
  expect_equal(mpone.sided_fixed(q, omega = omega), mptwo.sided_fixed(q, omega = omega))
  expect_equal(mqone.sided_fixed(c(0, .5, 1), omega = omega), mqtwo.sided_fixed(c(0, .5, 1), omega = omega))
})

test_that("non-monotonic one-sided density, cdf, and quantile helpers remain explicit", {

  expect_error(mdone.sided(0.5, alpha1 = c(1, 1), alpha2 = c(1, 1)), "Not implemented")
  expect_error(mpone.sided(0.5, alpha1 = c(1, 1), alpha2 = c(1, 1)), "Not implemented")
  expect_error(mqone.sided(0.5, alpha1 = c(1, 1), alpha2 = c(1, 1)), "Not implemented")
})

test_that("one-sided weightfunction wrappers require exactly one public parameterization", {

  expect_error(mdone.sided(0.5), "parameterization is required")
  expect_error(rone.sided(1), "parameterization is required")
  expect_error(mpone.sided(0.5), "parameterization is required")
  expect_error(mqone.sided(0.5), "parameterization is required")

  expect_error(mdone.sided(0.5, alpha = c(1, 1), alpha1 = c(1, 1)), "exactly one")
  expect_error(rone.sided(1, alpha = c(1, 1), alpha1 = c(1, 1), alpha2 = c(1, 1)), "exactly one")
  expect_error(mpone.sided(0.5, alpha = c(1, 1), alpha2 = c(1, 1)), "exactly one")
  expect_error(mqone.sided(0.5, alpha = c(1, 1), alpha1 = c(1, 1), alpha2 = c(1, 1)), "exactly one")

  expect_error(mdone.sided(0.5, alpha1 = c(1, 1)), "requires both 'alpha1' and 'alpha2'")
  expect_error(rone.sided(1, alpha2 = c(1, 1)), "requires both 'alpha1' and 'alpha2'")
  expect_error(mpone.sided(0.5, alpha1 = c(1, 1)), "requires both 'alpha1' and 'alpha2'")
  expect_error(mqone.sided(0.5, alpha2 = c(1, 1)), "requires both 'alpha1' and 'alpha2'")
})

test_that("one-sided weightfunction wrappers accept monotonic and general parameterizations", {

  expect_equal(dim(mdone.sided(0.5, alpha = c(1, 1))), c(1, 2))
  expect_equal(nrow(rone.sided(1, alpha = c(1, 1))), 1)
  expect_equal(dim(mpone.sided(0.5, alpha = c(1, 1))), c(1, 2))
  expect_equal(dim(mqone.sided(0.5, alpha = c(1, 1))), c(1, 2))

  general <- rone.sided(1, alpha1 = c(1, 1), alpha2 = c(1, 1))
  expect_equal(dim(general), c(1, 3))
  expect_equal(general[,1], 1)
})

test_that("weightfunction helper input validation is explicit", {

  expect_error(mdone.sided(0.5, alpha = "not_numeric"), "'alpha' must be a numeric vector or a matrix.")
  expect_error(mpone.sided(0.5, alpha = list(1, 1)), "'alpha' must be a numeric vector or a matrix.")
  expect_error(mqone.sided(0.5, alpha = data.frame(a = 1, b = 1)), "'alpha' must be a numeric vector or a matrix.")
  expect_error(rone.sided(5, alpha = "not_numeric"), "'alpha' must be a numeric vector or a matrix.")
  expect_error(mdone.sided(0.5, alpha = 1), "'alpha' must be a vector of length at least 2.")
  expect_error(mdone.sided(0.5, alpha = matrix(1, nrow = 2, ncol = 1)), "'alpha' must be a matrix with at least 2 columns.")
  expect_error(mdone.sided(0.5, alpha = c(-1, 1)), "'alpha' must be positive.")

  expect_error(mdone.sided_fixed(0.5, omega = "not_numeric"), "'omega' must be a numeric vector or a matrix.")
  expect_error(mpone.sided_fixed(0.5, omega = list(1, .5)), "'omega' must be a numeric vector or a matrix.")
  expect_error(mqone.sided_fixed(0.5, omega = data.frame(a = 1, b = .5)), "'omega' must be a numeric vector or a matrix.")
  expect_error(rone.sided_fixed(5, omega = "not_numeric"), "'omega' must be a numeric vector or a matrix.")
  expect_error(mdone.sided_fixed(0.5, omega = 1), "'omega' must be a vector of length at least 2.")
  expect_error(mdone.sided_fixed(0.5, omega = matrix(1, nrow = 2, ncol = 1)), "'omega' must be a matrix with at least 2 columns.")
  expect_error(mdone.sided_fixed(0.5, omega = c(-.1, 1)), "'omega' must be non-negative.")
  expect_error(mdone.sided_fixed(0.5, omega = c(.5, 1)), "reference-bin")
})

test_that("weightfunction helper dimension checks happen after valid parameter checks", {

  alpha_mat <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  omega_mat <- matrix(c(1, .5, 1, .6), nrow = 2, byrow = TRUE)

  expect_error(mdone.sided(c(.5, .6, .7), alpha = alpha_mat), "Non matching dimensions of 'alpha' and 'x'.")
  expect_error(mdone.sided_fixed(c(.5, .6, .7), omega = omega_mat), "Non matching dimensions of 'omega' and 'x'.")
  expect_error(mpone.sided(c(.5, .6, .7), alpha = alpha_mat), "Non matching dimensions of 'alpha' and 'q'.")
  expect_error(mpone.sided_fixed(c(.5, .6, .7), omega = omega_mat), "Non matching dimensions of 'omega' and 'q'.")
  expect_error(mqone.sided(c(.25, .5, .75), alpha = alpha_mat), "Non matching dimensions of 'alpha' and 'p'.")
  expect_error(mqone.sided_fixed(c(.25, .5, .75), omega = omega_mat), "Non matching dimensions of 'omega' and 'p'.")
  expect_error(rone.sided(5, alpha = alpha_mat), "Incompatible dimensions of requested number of samples and 'alpha'.")
  expect_error(rone.sided_fixed(5, omega = omega_mat), "Incompatible dimensions of requested number of samples and 'omega'.")

  alpha1_mat <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  alpha2_mat <- matrix(c(1, 1, 2, 2, 3, 3), nrow = 3, byrow = TRUE)
  expect_error(rone.sided(5, alpha1 = alpha1_mat, alpha2 = alpha2_mat), "Non matching dimensions of 'alpha1' and 'alpha2'.")
})

test_that("weightfunction helper matrix broadcasting works", {

  expect_equal(dim(mdone.sided(0.5, alpha = matrix(c(1, 1), nrow = 1))), c(1, 2))
  expect_equal(dim(mdone.sided_fixed(0.5, omega = matrix(c(1, .5), nrow = 1))), c(1, 2))
  expect_equal(dim(mdone.sided_fixed(c(.3, .5, .7), omega = c(1, .5))), c(3, 2))
  expect_equal(dim(mpone.sided(0.5, alpha = matrix(c(1, 1), nrow = 1))), c(1, 2))
  expect_equal(dim(mpone.sided_fixed(c(.3, .5, .7), omega = c(1, .5))), c(3, 2))
  expect_equal(dim(mqone.sided(0.5, alpha = matrix(c(1, 1), nrow = 1))), c(1, 2))
  expect_equal(dim(mqone.sided_fixed(c(.25, .5, .75), omega = c(1, .5))), c(3, 2))

  set.seed(1)
  expect_equal(nrow(rone.sided(2, alpha = matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE))), 2)
  expect_equal(nrow(rone.sided_fixed(2, omega = matrix(c(1, .3, 1, .5), nrow = 2, byrow = TRUE))), 2)
  expect_equal(nrow(rone.sided(5, alpha1 = matrix(c(1, 1), nrow = 1), alpha2 = matrix(c(1, 1), nrow = 1))), 5)
})

test_that("monotone weightfunction matrix parameters use row-specific beta shapes", {

  alpha <- matrix(c(
    1, 2, 3,
    2, 3, 5
  ), nrow = 2, byrow = TRUE)
  q <- c(.25, .75)
  p <- c(.2, .8)

  expected_density <- cbind(
    dpoint(q, 1),
    stats::dbeta(q, c(5, 8), c(1, 2)),
    stats::dbeta(q, c(3, 5), c(3, 5))
  )
  expected_cdf <- cbind(
    ppoint(q, 1),
    stats::pbeta(q, c(5, 8), c(1, 2)),
    stats::pbeta(q, c(3, 5), c(3, 5))
  )
  expected_quantile <- cbind(
    1,
    stats::qbeta(p, c(5, 8), c(1, 2)),
    stats::qbeta(p, c(3, 5), c(3, 5))
  )

  expect_equal(mdone.sided(q, alpha = alpha), expected_density, tolerance = 1e-12)
  expect_equal(mdtwo.sided(q, alpha = alpha), expected_density, tolerance = 1e-12)
  expect_equal(mpone.sided(q, alpha = alpha), expected_cdf, tolerance = 1e-12)
  expect_equal(mptwo.sided(q, alpha = alpha), expected_cdf, tolerance = 1e-12)
  expect_equal(mqone.sided(p, alpha = alpha), expected_quantile, tolerance = 1e-12)
  expect_equal(mqtwo.sided(p, alpha = alpha), expected_quantile, tolerance = 1e-12)
})

test_that("fixed weightfunction matrix parameters use row-specific point locations", {

  omega <- matrix(c(
    1, 0, 2,
    1, .5, 1.5
  ), nrow = 2, byrow = TRUE)
  q <- c(0, 1.5)
  p <- c(.25, .75)

  expected_density <- cbind(
    dpoint(q, c(1, 1)),
    dpoint(q, c(0, .5)),
    dpoint(q, c(2, 1.5))
  )
  expected_cdf <- cbind(
    ppoint(q, c(1, 1)),
    ppoint(q, c(0, .5)),
    ppoint(q, c(2, 1.5))
  )
  expected_upper_tail_log <- cbind(
    ppoint(q, c(1, 1), lower.tail = FALSE, log.p = TRUE),
    ppoint(q, c(0, .5), lower.tail = FALSE, log.p = TRUE),
    ppoint(q, c(2, 1.5), lower.tail = FALSE, log.p = TRUE)
  )
  expected_quantile <- cbind(
    qpoint(p, c(1, 1)),
    qpoint(p, c(0, .5)),
    qpoint(p, c(2, 1.5))
  )

  expect_equal(mdone.sided_fixed(q, omega = omega), expected_density)
  expect_equal(mdtwo.sided_fixed(q, omega = omega), expected_density)
  expect_equal(mdone.sided_fixed(q, omega = omega, log = TRUE), log(expected_density))
  expect_equal(mpone.sided_fixed(q, omega = omega), expected_cdf)
  expect_equal(mptwo.sided_fixed(q, omega = omega), expected_cdf)
  expect_equal(
    mpone.sided_fixed(q, omega = omega, lower.tail = FALSE, log.p = TRUE),
    expected_upper_tail_log
  )
  expect_equal(mqone.sided_fixed(p, omega = omega), expected_quantile)
  expect_equal(mqtwo.sided_fixed(p, omega = omega), expected_quantile)
})
