# ============================================================================ #
# TEST FILE: Prior Density Function
# ============================================================================ #
#
# PURPOSE:
#   Visual regression tests for the density.prior S3 method including
#   various transformation options.
#
# DEPENDENCIES:
#   - vdiffr: Visual regression testing
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("vdiffr")
#
# TAGS: @evaluation, @visual, @priors, @density
# ============================================================================ #

skip_if_not_installed("vdiffr")

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
  vdiffr::expect_doppelganger("prior-density-3-1", function()plot(density(prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1))))))
  vdiffr::expect_doppelganger("prior-density-3-2", function()plot(density(prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .25))))))

  # factor priors
  vdiffr::expect_doppelganger("prior-density-4-1", function()plot(density(prior_factor("normal", list(0, 1), list(0, Inf), contrast = "independent"))))
  vdiffr::expect_doppelganger("prior-density-4-2", function()plot(density(prior_factor("normal", list(0, 1), list(0, Inf), contrast = "treatment"))))
  vdiffr::expect_doppelganger("prior-density-4-3", function()suppressWarnings(plot(density(prior_factor("mnormal", list(0, 1), contrast = "orthonormal")))))
  vdiffr::expect_doppelganger("prior-density-4-4", function()suppressWarnings(plot(density(prior_factor("mnormal", list(0, 1), contrast = "meandif")))))
  vdiffr::expect_doppelganger("prior-density-4-5", function()suppressWarnings(plot(density(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"), force_samples = TRUE))))

  # PET
  vdiffr::expect_doppelganger("prior-density-5-1", function()plot(density(prior_PET("normal", list(0, 1)))))
  vdiffr::expect_doppelganger("prior-density-5-2", function()plot(density(prior_PET("normal", list(0, 1)), force_samples = TRUE)))
  vdiffr::expect_doppelganger("prior-density-5-3", function()plot(density(prior_PET("normal", list(0, 1)), force_samples = TRUE, transformation = "tanh")))

  # no plotting etc implemented
  xd <- density(prior_spike_and_slab(prior("normal", list(0, 1))))
  expect_s3_class(xd, "density.prior.spike_and_slab")

})


