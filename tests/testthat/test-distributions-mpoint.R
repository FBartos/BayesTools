context("Distributions - Multivariate point")

test_that("Density function works", {

  expect_equal(dmpoint(1, 1), Inf)
  expect_equal(dmpoint(0, 1),   0)

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

  expect_equal(qmpoint(c(0, 1), matrix(1:6, ncol = 3, byrow = TRUE), lower.tail = TRUE), matrix(c(-Inf, -Inf, -Inf, 4:6), byrow = 3, ncol = 3))
  expect_equal(qmpoint(log(c(0.5, 1)), matrix(1:6, ncol = 3), log.p = TRUE), matrix(1:6, ncol = 3))

  expect_error(qmpoint(c(1, 1), matrix(1:9, ncol = 3)), "Non matching dimensions of 'location' and 'p'.")
})

