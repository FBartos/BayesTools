context("Distributions - Point")

test_that("Density function works", {

  expect_equal(dpoint(1, 1), Inf)
  expect_equal(dpoint(0, 1),   0)

  expect_equal(dpoint(c(0, 1), 1),       c(0,   Inf))
  expect_equal(dpoint(c(0, 1), c(0, 1)), c(Inf, Inf))

  expect_equal(dpoint(c(0, 1), 1, log = TRUE),       log(c(0,   Inf)))
  expect_equal(dpoint(c(0, 1), c(0, 1), log = TRUE), log(c(Inf, Inf)))

  expect_error(dpoint(0, c(1, 1)), "Non matching dimensions of 'location' and 'x'.")
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

