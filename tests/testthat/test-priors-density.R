context("Prior density")

test_that("Prior density function density", {
  set.seed(1)
  vdiffr::expect_doppelganger("prior-density-1-1", function()plot(density(prior("normal", list(0, 1)))))
  vdiffr::expect_doppelganger("prior-density-1-2", function()plot(density(prior("normal", list(0, 1)), x_seq   = seq(0, 1, 0.1))))
  vdiffr::expect_doppelganger("prior-density-1-3", function()plot(density(prior("normal", list(0, 1)), x_range = c(-1, 1))))
  vdiffr::expect_doppelganger("prior-density-1-4", function()plot(density(prior("normal", list(0, 1)), x_range_quant  = .17)))
  vdiffr::expect_doppelganger("prior-density-1-5", function()plot(density(prior("normal", list(0, 1)), force_samples  = TRUE)))
  vdiffr::expect_doppelganger("prior-density-1-6", function()plot(density(prior("normal", list(0, 1)), transformation = "tanh")))
  vdiffr::expect_doppelganger("prior-density-1-7", function()plot(density(prior("normal", list(0, 1)), transformation = "exp")))
  vdiffr::expect_doppelganger("prior-density-1-8", function()plot(density(prior("normal", list(0, 1)), transformation = "lin", transformation_arguments = list(b = 0.5))))

  set.seed(2)
  vdiffr::expect_doppelganger("prior-density-2-1", function()plot(density(prior("normal", list(0, 1), list(0, Inf)))))
  vdiffr::expect_doppelganger("prior-density-2-2", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_seq   = seq(0, 1, 0.1))))
  vdiffr::expect_doppelganger("prior-density-2-3", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_range = c(-1, 1))))
  vdiffr::expect_doppelganger("prior-density-2-4", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), x_range_quant  = .17)))
  vdiffr::expect_doppelganger("prior-density-2-5", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), force_samples  = TRUE)))
  vdiffr::expect_doppelganger("prior-density-2-6", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "tanh")))
  vdiffr::expect_doppelganger("prior-density-2-7", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "exp")))
  vdiffr::expect_doppelganger("prior-density-2-8", function()plot(density(prior("normal", list(0, 1), list(0, Inf)), transformation = "lin", transformation_arguments = list(b = 0.5))))

  # weightfunctions
  vdiffr::expect_doppelganger("prior-density-3-1", function()plot(density(prior_weightfunction("one.sided", list(c(.05), c(1, 1))))))
  vdiffr::expect_doppelganger("prior-density-3-2", function()plot(density(prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .25))))))

  # factor priors
  vdiffr::expect_doppelganger("prior-density-4-1", function()plot(density(prior_factor("normal", list(0, 1), list(0, Inf), contrast = "independent"))))
  vdiffr::expect_doppelganger("prior-density-4-2", function()plot(density(prior_factor("normal", list(0, 1), list(0, Inf), contrast = "treatment"))))
  vdiffr::expect_doppelganger("prior-density-4-3", function()suppressWarnings(plot(density(prior_factor("mnormal", list(0, 1), contrast = "orthonormal")))))
  vdiffr::expect_doppelganger("prior-density-4-4", function()suppressWarnings(plot(density(prior_factor("mnormal", list(0, 1), contrast = "meandif")))))
})


