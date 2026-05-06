skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Distributions - Point
# ============================================================================ #
#
# PURPOSE:
#   Tests for point distribution functions (dpoint, rpoint, ppoint, qpoint).
#
# DEPENDENCIES:
#   - None (pure R)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# TAGS: @evaluation, @distributions, @point
# ============================================================================ #

test_that("Density function works", {

  expect_equal(dpoint(1, 1), Inf)
  expect_equal(dpoint(0, 1),   0)
  expect_equal(dpoint(1 + 1e-9, 1), 0)
  expect_equal(dpoint(NA_real_, 1), NA_real_)
  expect_equal(dpoint(1, NA_real_), NA_real_)

  expect_equal(dpoint(c(0, 1), 1),       c(0,   Inf))
  expect_equal(dpoint(c(0, 1), c(0, 1)), c(Inf, Inf))

  expect_equal(dpoint(c(0, 1), 1, log = TRUE),       log(c(0,   Inf)))
  expect_equal(dpoint(c(0, 1), c(0, 1), log = TRUE), log(c(Inf, Inf)))

  expect_error(dpoint(0, c(1, 1)), "Non matching dimensions of 'location' and 'x'.")
})

test_that("Point distribution preserves names and rejects unsupported inputs", {
  x <- c(below = 0, at = 1)
  p <- c(zero = 0, half = .5)

  expect_named(dpoint(x, 1), names(x))
  expect_named(ppoint(x, 1), names(x))
  expect_named(qpoint(p, 1), names(p))

  expect_error(dpoint(numeric(), 1), "cannot be NULL")
  expect_error(ppoint(numeric(), 1), "cannot be NULL")
  expect_error(qpoint(numeric(), 1), "cannot be NULL")
  expect_error(rpoint(0, 1), "equal or higher than 1")

  attributed <- structure(c(0, 1), source = "external")
  expect_error(dpoint(attributed, 1), "numeric vector")
  expect_error(ppoint(attributed, 1), "numeric vector")
  expect_error(qpoint(attributed, 1), "numeric vector")
})

test_that("Random generator works", {

  expect_equal(rpoint(10, 1),      rep(1, 10))
  expect_equal(rpoint(2, c(1, 2)), c(1, 2))

  expect_error(rpoint(10, c(1, 2)), "Incompatible dimensions of requested number of samples and 'location'.")
})

test_that("Distribution function works", {

  expect_equal(ppoint(1, 1), 1)
  expect_equal(ppoint(c(0, 1, 2), 1), c(0, 1, 1))

  expect_equal(ppoint(c(2, 2, 1), c(1, 1, 2)), c(1, 1, 0))
  expect_equal(ppoint(c(2, 2, 1), c(1, 1, 2), lower.tail = FALSE), c(0, 0, 1))
  expect_equal(ppoint(c(2, 2, 1), c(1, 1, 2), log.p = TRUE), log(c(1, 1, 0)))

  expect_error(ppoint(c(1, 2), c(1, 2, 3)), "Non matching dimensions of 'location' and 'q'.")
})

test_that("Quantile function works", {

  expect_equal(qpoint(1, 1), 1)
  expect_equal(qpoint(0, 1), -Inf)
  expect_equal(qpoint(c(0, .5, 1), 1, lower.tail = FALSE), c(Inf, 1, 1))
  expect_equal(qpoint(log(c(1, .5, 0)), c(0, 1, 2), log.p = TRUE), c(0, 1, -Inf))

  expect_error(qpoint(rep(0.5, 3), c(1, 2)))
})
