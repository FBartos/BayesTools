context("Prior density")

test_that("Prior density function density", {
  set.seed(1)
  expect_doppelganger("prior-density-1-1", function()plot(density(prior("normal", list(0, 1)))))
  expect_doppelganger("prior-density-1-2", function()plot(density(prior("normal", list(0, 1)), x_seq   = seq(0, 1, 0.1))))
  expect_doppelganger("prior-density-1-3", function()plot(density(prior("normal", list(0, 1)), x_range = c(-1, 1))))
  expect_doppelganger("prior-density-1-4", function()plot(density(prior("normal", list(0, 1)), x_range_quant  = .17)))
  expect_doppelganger("prior-density-1-5", function()plot(density(prior("normal", list(0, 1)), force_samples  = TRUE)))
  expect_doppelganger("prior-density-1-6", function()plot(density(prior("normal", list(0, 1)), transformation = "tanh")))
  expect_doppelganger("prior-density-1-7", function()plot(density(prior("normal", list(0, 1)), transformation = "exp")))
  expect_doppelganger("prior-density-1-8", function()plot(density(prior("normal", list(0, 1)), transformation = "lin", transformation_arguments = list(b = 0.5))))

  set.seed(2)
  expect_doppelganger("prior-density-2-1", function()plot(density(prior("normal", list(0, 1), list(0, Inf)))))
  expect_doppelganger("prior-density-2-2", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_seq   = seq(0, 1, 0.1))))
  expect_doppelganger("prior-density-2-3", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_range = c(-1, 1))))
  expect_doppelganger("prior-density-2-4", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_range_quant  = .17)))
  expect_doppelganger("prior-density-2-5", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), force_samples  = TRUE)))
  expect_doppelganger("prior-density-2-6", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "tanh")))
  expect_doppelganger("prior-density-2-7", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "exp")))
  expect_doppelganger("prior-density-2-8", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "lin", transformation_arguments = list(b = 0.5))))

  # weightfunctions
  expect_doppelganger("prior-density-3-1", function()plot(density(prior_weightfunction("one.sided", list(c(.05), c(1, 1))))))
  expect_doppelganger("prior-density-3-2", function()plot(density(prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .25))))))
})
