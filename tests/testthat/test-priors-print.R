skip_if_not_test_profile(c("unit", "visual"))

# ============================================================================ #
# TEST FILE: Prior Print Function
# ============================================================================ #
#
# PURPOSE:
#   Tests for the print.prior S3 method including input validation,
#   formatting options, and output correctness.
#
# DEPENDENCIES:
#   - None (pure R)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# TAGS: @evaluation, @priors, @print
# ============================================================================ #


test_that("Prior print function input validation", {

  p <- prior("normal", list(0, 1))

  # Check invalid inputs
  expect_error(print(p, short_name = "no"), "'short_name'")
  expect_error(print(p, parameter_names = "no"), "'parameter_names'")
  expect_error(print(p, digits_estimates = "two"), "'digits_estimates'")
  expect_error(print(p, plot = "yes"), "'plot'")
  expect_error(print(p, silent = "shh"), "'silent'")
  expect_error(print(p, inline = "no"), "'inline'")

})


test_that("Prior print function works", {

  # check the default options
  p1 <- prior("normal", list(0, 1))

  expect_equal(utils::capture.output(print(p1)), "Normal(0, 1)")
  expect_equal(utils::capture.output(print(p1, short_name = TRUE)), "N(0, 1)")
  expect_equal(utils::capture.output(print(p1, parameter_names = TRUE)), "Normal(mean = 0, sd = 1)")
  expect_equal(utils::capture.output(print(p1, silent = TRUE)), character())

  # check dealing with truncation
  p2 <- prior("Cauchy", list(0, 1), list(0, Inf))
  expect_equal(utils::capture.output(print(p2)), "Cauchy(0, 1)[0, Inf]")

  p3 <- prior("gamma", list(1, 1))
  expect_equal(utils::capture.output(print(p3)), "Gamma(1, 1)")

  p4 <- prior("gamma", list(1, 1), list(0, Inf))
  expect_equal(utils::capture.output(print(p4)), "Gamma(1, 1)")

  # check prefixes
  p5 <- prior_PET("normal", list(1, 1))
  p6 <- prior_PEESE("gamma", list(1, 1))
  expect_equal(utils::capture.output(print(p5)), "PET ~ Normal(1, 1)[0, Inf]")
  expect_equal(utils::capture.output(print(p6)), "PEESE ~ Gamma(1, 1)")

  # check weightfunctions
  p7  <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
  p8  <- prior_weightfunction("one-sided", c(0.05, .95), wf_independent(prior("beta", list(1, 1))))
  p9  <- prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1)))
  p10 <- prior_weightfunction("one-sided", c(0.10), wf_fixed(c(1, .7)))
  expect_equal(utils::capture.output(print(p7)),  "omega[one-sided: .05] ~ CumDirichlet(1, 1)")
  expect_equal(utils::capture.output(print(p8)),  "omega[one-sided: .05, .95] ~ Independent(Beta(1, 1))")
  expect_equal(utils::capture.output(print(p9)),  "omega[two-sided: .05] ~ CumDirichlet(1, 1)")
  expect_equal(utils::capture.output(print(p10)), "omega[one-sided: .1] = (1, 0.7)")
  expect_equal(utils::capture.output(print(p7,  parameter_names = TRUE)), "omega[one-sided: .05] ~ CumDirichlet(alpha = 1, 1)")
  expect_equal(utils::capture.output(print(p8,  parameter_names = TRUE)), "omega[one-sided: .05, .95] ~ Independent(Beta(alpha = 1, beta = 1))")
  expect_equal(utils::capture.output(print(p9,  parameter_names = TRUE)), "omega[two-sided: .05] ~ CumDirichlet(alpha = 1, 1)")
  expect_equal(utils::capture.output(print(p10, parameter_names = TRUE)), "omega[one-sided: .1] = (1, 0.7)")

  # check vector priors
  p11 <- prior(distribution = "mnormal", parameters = list(mean = 0, sd = 1, K = 3))
  p12 <- prior(distribution = "mcauchy", parameters = list(0, 1, 5))
  p13 <- prior(distribution = "mt",      parameters = list(location = 1, scale = .5, df = 3, K = 2))

  expect_equal(utils::capture.output(print(p11)),  "mNormal(0, 1)")
  expect_equal(utils::capture.output(print(p12)),  "mCauchy(0, 1)")
  expect_equal(utils::capture.output(print(p13)),  "mStudent-t(1, 0.5, 3)")

  # check factor priors
  p14 <- prior_factor(distribution = "mnormal", contrast = "orthonormal", parameters = list(0, 1))
  p15 <- prior_factor(distribution = "normal", contrast = "treatment", parameters = list(mean = 0, sd = 1))
  p16 <- prior_factor(distribution = "beta",   contrast = "treatment", parameters = list(alpha = 1, beta = 1))
  p17 <- prior_factor(distribution = "beta",   contrast = "independent", parameters = list(alpha = 1, beta = 1))
  p18 <- prior_factor(distribution = "mnormal", contrast = "meandif", parameters = list(0, 0.5))
  p19 <- prior_factor(distribution = "point", contrast = "orthonormal", parameters = list(location = 0))
  p20 <- prior_factor(distribution = "spike", contrast = "meandif", parameters = list(location = 0))

  expect_equal(utils::capture.output(print(p14)),  "orthonormal contrast: mNormal(0, 1)")
  expect_equal(utils::capture.output(print(p15)),  "treatment contrast: Normal(0, 1)")
  expect_equal(utils::capture.output(print(p16)),  "treatment contrast: Beta(1, 1)")
  expect_equal(utils::capture.output(print(p17)),  "independent contrast: Beta(1, 1)")
  expect_equal(utils::capture.output(print(p18)),  "mean difference contrast: mNormal(0, 0.5)")
  expect_equal(utils::capture.output(print(p19)),  "orthonormal contrast: mSpike(0)")
  expect_equal(utils::capture.output(print(p20)),  "mean difference contrast: mSpike(0)")

  # check plot names
  empty_plot <- function(){
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  }
  if (bayestools_test_profile_includes("visual")) {
    testthat::skip_if_not_installed("vdiffr")

    vdiffr::expect_doppelganger("priors-print-1", function(){
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(mar = oldpar[["mar"]]))
      par(mar = c(0, 0, 0, 0))
      empty_plot()
      text(0.5, 1,   print(p1, plot = TRUE))
      text(0.5, 0.9, print(p1, short_name = TRUE, plot = TRUE))
      text(0.5, 0.8, print(p1, parameter_names = TRUE, plot = TRUE))
      text(0.5, 0.7, print(p1, silent = TRUE, plot = TRUE))
      text(0.5, 0.6, print(p2, plot = TRUE))
      text(0.5, 0.5, print(p2, short_name = TRUE, plot = TRUE))
      text(0.5, 0.4, print(p2, parameter_names = TRUE, plot = TRUE))
      text(0.5, 0.3, print(p2, silent = TRUE, plot = TRUE))
      text(0.5, 0.2, print(p3, silent = TRUE, plot = TRUE))
      text(0.5, 0.1, print(p4, silent = TRUE, plot = TRUE))
    })

    vdiffr::expect_doppelganger("priors-print-2", function(){
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(mar = oldpar[["mar"]]))
      par(mar = c(0, 0, 0, 0))
      empty_plot()
      text(0.5, 1,   print(p5,  plot = TRUE))
      text(0.5, 0.9, print(p6,  plot = TRUE))
      text(0.5, 0.8, print(p7,  plot = TRUE))
      text(0.5, 0.7, print(p8,  plot = TRUE))
      text(0.5, 0.6, print(p9,  plot = TRUE))
      text(0.5, 0.5, print(p10, plot = TRUE))
      text(0.5, 0.4, print(p7,  parameter_names = TRUE, plot = TRUE))
      text(0.5, 0.3, print(p8,  parameter_names = TRUE, plot = TRUE))
      text(0.5, 0.2, print(p9,  parameter_names = TRUE, plot = TRUE))
      text(0.5, 0.1, print(p10, parameter_names = TRUE, plot = TRUE))
    })

    vdiffr::expect_doppelganger("priors-print-3", function(){
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(mar = oldpar[["mar"]]))
      par(mar = c(0, 0, 0, 0))
      empty_plot()
      text(0.5, 1,   print(p11, plot = TRUE))
      text(0.5, 0.9, print(p12, plot = TRUE))
      text(0.5, 0.8, print(p13, plot = TRUE))
      text(0.5, 0.7, print(p14, plot = TRUE))
      text(0.5, 0.6, print(p15, plot = TRUE))
      text(0.5, 0.5, print(p16, plot = TRUE))
      text(0.5, 0.4, print(p17, plot = TRUE))
      text(0.5, 0.3, print(p18, plot = TRUE))
      text(0.5, 0.2, print(p19, plot = TRUE))
      text(0.5, 0.1, print(p20, plot = TRUE))
    })
  }

  p21 <- prior_spike_and_slab(prior("gamma", list(1, 2), list(0, Inf)),
                              prior_inclusion = prior("beta", list(3, 2)))
  p22 <- prior_mixture(
    list(
      prior("normal", list(0,  1)),
      prior("normal", list(-3, 1)),
      prior("gamma",  list(5, 10))
    )
  )
  p23 <- prior_mixture(
    list(
      prior("normal", list(0,  1), prior_weights = 1),
      prior("normal", list(-3, 1), prior_weights = 5),
      prior("gamma",  list(5, 10), prior_weights = 1)
    ),
    is_null = c(T, F, T)
  )
  p24 <- prior_mixture(
    list(
      prior("normal", list(0,  1), prior_weights = 1),
      prior("normal", list(-3, 1), prior_weights = 5)
    ),
    components = c("b", "a")
  )

  expect_equal(utils::capture.output(print(p21)), "Gamma(1, 2) * Beta(3, 2)")
  expect_equal(utils::capture.output(print(p21, short_name = TRUE)), "G(1, 2) * B(3, 2)")
  expect_equal(utils::capture.output(print(p22, parameter_names = TRUE)), c(
    "alternative:", "  (1/3) * Normal(mean = 0, sd = 1)", "  (1/3) * Normal(mean = -3, sd = 1)", "  (1/3) * Gamma(shape = 5, rate = 10)"
  ))
  expect_equal(utils::capture.output(print(p23, short_name = TRUE)), c(
    "alternative:", "  (5/7) * N(-3, 1)", "null:", "  (1/7) * N(0, 1)", "  (1/7) * G(5, 10)"
  ))
  expect_equal(utils::capture.output(print(p24)), c(
    "b:", "  (1/6) * Normal(0, 1)",  "a:", "  (5/6) * Normal(-3, 1)"
  ))
  if (bayestools_test_profile_includes("visual")) {
    testthat::skip_if_not_installed("vdiffr")

    vdiffr::expect_doppelganger("priors-print-4", function(){
      empty_plot()
      text(0.5, 1, print(p21, plot = TRUE))
      text(0.5, 0.9, print(p22, plot = TRUE))
      text(0.5, 0.8, print(p23, plot = TRUE))
      text(0.5, 0.7, print(p24, plot = TRUE))
    })
  }

  # prior expressions
  pe1 <- prior("normal", parameters = list(0, expression(x)))
  expect_equal(utils::capture.output(print(pe1)), "Normal(0, x)")

  if (bayestools_test_profile_includes("visual")) {
    testthat::skip_if_not_installed("vdiffr")

    vdiffr::expect_doppelganger("priors-print-e1", function(){
      empty_plot()
      text(0.5, 1, print(pe1, plot = TRUE))
    })
  }
})


test_that("Prior print for prior_none", {

  p_none <- prior_none()
  output <- utils::capture.output(print(p_none))
  expect_type(output, "character")

  # Silent output
  expect_equal(utils::capture.output(print(p_none, silent = TRUE)), character())

})


test_that("Prior print with inline option for mixtures", {

  p_mix <- prior_mixture(
    list(
      prior("normal", list(0, 1)),
      prior("normal", list(0, 2))
    )
  )

  # Test inline option
  output_inline <- print(p_mix, silent = TRUE, inline = TRUE)
  expect_type(output_inline, "character")

})


test_that("Prior print for additional distributions", {

  # Beta distribution with different parameters
  p_beta <- prior("beta", list(alpha = 2, beta = 5))
  expect_equal(utils::capture.output(print(p_beta)), "Beta(2, 5)")
  expect_equal(utils::capture.output(print(p_beta, short_name = TRUE)), "B(2, 5)")

  # Exponential distribution
  p_exp <- prior("exp", list(rate = 2))
  expect_equal(utils::capture.output(print(p_exp)), "Exponential(2)")
  expect_equal(utils::capture.output(print(p_exp, short_name = TRUE)), "E(2)")

  # Uniform distribution
  p_unif <- prior("uniform", list(a = -1, b = 1))
  expect_equal(utils::capture.output(print(p_unif)), "Uniform(-1, 1)")
  expect_equal(utils::capture.output(print(p_unif, short_name = TRUE)), "U(-1, 1)")

  # Lognormal distribution
  p_ln <- prior("lognormal", list(meanlog = 0, sdlog = 1))
  expect_equal(utils::capture.output(print(p_ln)), "Lognormal(0, 1)")
  expect_equal(utils::capture.output(print(p_ln, short_name = TRUE)), "Ln(0, 1)")

  # Inverse gamma distribution
  p_ig <- prior("invgamma", list(shape = 1, scale = 1))
  expect_equal(utils::capture.output(print(p_ig)), "InvGamma(1, 1)")
  expect_equal(utils::capture.output(print(p_ig, short_name = TRUE)), "Ig(1, 1)")

})


test_that("Prior print digits_estimates parameter", {

  p <- prior("normal", list(mean = 1.2345678, sd = 0.9876543))

  # Default (2 digits)
  expect_match(utils::capture.output(print(p)), "Normal\\(1\\.23, 0\\.99\\)")

  # 4 digits
  expect_match(utils::capture.output(print(p, digits_estimates = 4)), "Normal\\(1\\.2346, 0\\.9877\\)")

  # 0 digits
  expect_match(utils::capture.output(print(p, digits_estimates = 0)), "Normal\\(1, 1\\)")

})

