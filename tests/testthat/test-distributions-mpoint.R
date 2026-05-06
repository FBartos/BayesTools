skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Distributions - Multivariate Point
# ============================================================================ #
#
# PURPOSE:
#   Tests for multivariate point distribution functions (dmpoint, rmpoint,
#   pmpoint, qmpoint).
#
# DEPENDENCIES:
#   - None (pure R)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# TAGS: @evaluation, @distributions, @multivariate
# ============================================================================ #

test_that("Density function works", {

  expect_equal(dmpoint(1, 1), Inf)
  expect_equal(dmpoint(0, 1),   0)
  expect_equal(dmpoint(1 + 1e-9, 1), 0)
  expect_equal(dmpoint(NA_real_, 1), NA_real_)
  expect_equal(dmpoint(matrix(c(0, NA_real_), nrow = 1), c(0, 1)), NA_real_)

  expect_equal(dmpoint(matrix(c(0, 1), ncol = 1), 1), c(0,   Inf))
  expect_equal(dmpoint(c(0, 1), c(0, 1)), Inf)
  expect_equal(dmpoint(matrix(c(0, 1), ncol = 2, nrow = 2, byrow = TRUE), c(0, 1)), c(Inf, Inf))
  expect_equal(dmpoint(matrix(c(0, 1), ncol = 2, nrow = 2, byrow = TRUE), matrix(c(0, 1), nrow = 1)), c(Inf, Inf))
  expect_equal(dmpoint(matrix(c(0, 1), ncol = 2, nrow = 2, byrow = TRUE), matrix(c(0, 1, 0, 0), ncol = 2, nrow = 2, byrow = TRUE)), c(Inf, 0))

  expect_equal(dmpoint(matrix(c(0, 1), ncol = 1), 1, log = TRUE),       log(c(0,   Inf)))
  expect_equal(dmpoint(matrix(c(0, 1), ncol = 2, nrow = 2, byrow = TRUE), c(0, 1), log = TRUE), log(c(Inf, Inf)))

  expect_error(dmpoint(c(0, 1), c(0, 1, 2)), "Non matching dimensions of 'location' and 'x'.")
  expect_error(dmpoint(c(0, 1), matrix(c(0, 1, 2), nrow = 1)), "Non matching dimensions of 'location' and 'x'.")
  expect_error(dmpoint(matrix(c(0, 1), nrow = 1), matrix(c(0, 1, 2), nrow = 1)), "Non matching dimensions of 'location' and 'x'.")
})

test_that("Multivariate point distribution accepts named inputs and rejects zero-length inputs", {

  x <- c(first = 0, second = 1)
  q <- matrix(c(0, 1, 2, 2), nrow = 2, byrow = TRUE, dimnames = list(c("below", "above"), c("a", "b")))
  p <- c(zero = 0, half = .5)

  expect_equal(dmpoint(x, c(0, 1)), Inf)
  expect_null(names(dmpoint(x, c(0, 1))))

  expect_equal(pmpoint(q, c(1, 1)), c(0, 1))
  expect_null(names(pmpoint(q, c(1, 1))))

  expect_equal(qmpoint(p, c(0, 1)), matrix(c(-Inf, -Inf, 0, 1), nrow = 2, byrow = TRUE))
  expect_null(dimnames(qmpoint(p, c(0, 1))))

  expect_error(dmpoint(numeric(), 1), "cannot be NULL")
  expect_error(pmpoint(numeric(), 1), "cannot be NULL")
  expect_error(qmpoint(numeric(), 1), "cannot be NULL")
  expect_error(rmpoint(0, 1), "equal or higher than 1")
})

test_that("Multivariate point distribution handles NA, NaN, and infinite values", {

  x <- matrix(c(
    NA_real_, 1,
    NaN,      1,
    Inf,      1,
    -Inf,     1,
    Inf,      1
  ), ncol = 2, byrow = TRUE)
  location <- matrix(c(
    0,   1,
    0,   1,
    0,   1,
    0,   1,
    Inf, 1
  ), ncol = 2, byrow = TRUE)

  expect_equal(dmpoint(x, location), c(NA_real_, NA_real_, 0, 0, Inf))
  expect_equal(pmpoint(x, location), c(NA_real_, NA_real_, 1, 0, 1))
  expect_equal(pmpoint(x, location, lower.tail = FALSE), c(NA_real_, NA_real_, 0, 1, 0))

  expect_equal(qmpoint(c(NA_real_, NaN, 0, 1), c(0, 1)), matrix(c(
    NA_real_, NA_real_,
    NA_real_, NA_real_,
    -Inf,     -Inf,
    0,        1
  ), ncol = 2, byrow = TRUE))
  expect_error(qmpoint(Inf, c(0, 1)), "equal or lower than 1")
  expect_error(qmpoint(-Inf, c(0, 1)), "equal or higher than 0")
})

test_that("Random generator works", {

  expect_equal(rmpoint(10, 1),      matrix(1, nrow = 10, ncol = 1))
  expect_equal(rmpoint(1, c(1, 2)), matrix(c(1, 2), nrow = 1))
  expect_equal(rmpoint(3, c(1, 2)), matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE))
  expect_equal(rmpoint(3, matrix(c(1, 2), nrow = 1)), matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE))
  expect_equal(rmpoint(3, matrix(1:6, nrow = 3, ncol = 2)), matrix(1:6, nrow = 3, ncol = 2))

  expect_error(rmpoint(2, matrix(1:6, nrow = 3, ncol = 2)), "Incompatible dimensions of requested number of samples and 'location'.")
})

test_that("Distribution function works", {

  expect_equal(pmpoint(1, 1), 1)
  expect_equal(pmpoint(c(0, 1, 2), 1), 0)
  expect_equal(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 2), 1), c(0, 1, 1))
  expect_equal(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 2), matrix(1, ncol = 2)), c(0, 1, 1))
  expect_equal(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 2), matrix(1, ncol = 2, nrow = 3)), c(0, 1, 1))

  expect_equal(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 2), matrix(1, ncol = 2, nrow = 3), lower.tail = FALSE), c(1, 0, 0))
  expect_equal(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 2), matrix(1, ncol = 2, nrow = 3), log.p = TRUE), log(c(0, 1, 1)))

  expect_error(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 3), c(1, 2)), "Non matching dimensions of 'location' and 'q'.")
  expect_error(pmpoint(matrix(c(0, 1, 2), nrow = 3, ncol = 3), matrix(1, ncol = 2, nrow = 3)), "Non matching dimensions of 'location' and 'q'.")
})

test_that("Quantile function works", {

  expect_equal(qmpoint(1, 1), matrix(1))
  expect_equal(qmpoint(0, 1), matrix(-Inf))
  expect_equal(qmpoint(c(1, 1, 1), c(1, 1, 1)), matrix(1, ncol = 3, nrow = 3))
  expect_equal(qmpoint(c(1, 1), matrix(1:6, ncol = 3)), matrix(1:6, ncol = 3))
  expect_equal(qmpoint(NA_real_, c(1, 2)), matrix(c(NA_real_, NA_real_), nrow = 1))

  expect_equal(qmpoint(c(0, 1), matrix(1:6, ncol = 3, byrow = TRUE), lower.tail = TRUE), matrix(c(-Inf, -Inf, -Inf, 4:6), byrow = 3, ncol = 3))
  expect_equal(qmpoint(log(c(0.5, 1)), matrix(1:6, ncol = 3), log.p = TRUE), matrix(1:6, ncol = 3))

  expect_error(qmpoint(c(1, 1), matrix(1:9, ncol = 3)), "Non matching dimensions of 'location' and 'p'.")
})

test_that("Multivariate point PMF, CDF, and quantile options are exact", {

  x <- matrix(c(1, 1, 0, 1, 2, 2), ncol = 2, byrow = TRUE)
  q <- matrix(c(0, 1, 1, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  location <- c(1, 1)
  p <- c(0, .5, 1)

  expect_equal(dmpoint(x, location), c(Inf, 0, 0))
  expect_equal(dmpoint(x, location, log = TRUE), c(Inf, -Inf, -Inf))

  expect_equal(pmpoint(q, location), c(0, 0, 1, 1))
  expect_equal(pmpoint(q, location, lower.tail = FALSE), c(1, 1, 0, 0))
  expect_equal(pmpoint(q, location, log.p = TRUE), c(-Inf, -Inf, 0, 0))
  expect_equal(pmpoint(q, location, lower.tail = FALSE, log.p = TRUE), c(0, 0, -Inf, -Inf))

  expect_equal(qmpoint(p, location), matrix(c(-Inf, -Inf, 1, 1, 1, 1), ncol = 2, byrow = TRUE))
  expect_equal(qmpoint(p, location, lower.tail = FALSE), matrix(c(Inf, Inf, 1, 1, 1, 1), ncol = 2, byrow = TRUE))
  expect_equal(qmpoint(log(p), location, log.p = TRUE), matrix(c(-Inf, -Inf, 1, 1, 1, 1), ncol = 2, byrow = TRUE))
  expect_equal(qmpoint(log(p), location, lower.tail = FALSE, log.p = TRUE), matrix(c(Inf, Inf, 1, 1, 1, 1), ncol = 2, byrow = TRUE))
})
