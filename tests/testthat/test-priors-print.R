context("Prior print function")


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
  p7  <- prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  p8  <- prior_weightfunction("one.sided", list(c(0.05, .95), c(1, 1), c(1, 1)))
  p9  <- prior_weightfunction("two.sided", list(c(0.05), c(1, 1)))
  p10 <- prior_weightfunction("one.sided.fixed", list(c(0.10), c(.7, 1)))
  expect_equal(utils::capture.output(print(p7)),  "omega[one-sided: .05] ~ CumDirichlet(1, 1)")
  expect_equal(utils::capture.output(print(p8)),  "omega[one-sided: .95, .05] ~ CumDirichlet(1, 1), revCumDirichlet(1, 1)")
  expect_equal(utils::capture.output(print(p9)),  "omega[two-sided: .05] ~ CumDirichlet(1, 1)")
  expect_equal(utils::capture.output(print(p10)), "omega[one-sided: .1] = (1, 0.7)")
  expect_equal(utils::capture.output(print(p7,  parameter_names = TRUE)), "omega[one-sided: .05] ~ CumDirichlet(alpha = 1, 1)")
  expect_equal(utils::capture.output(print(p8,  parameter_names = TRUE)), "omega[one-sided: .95, .05] ~ CumDirichlet(alpha1 = 1, 1), revCumDirichlet(alpha2 = 1, 1)")
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
  p15 <- prior_factor(distribution = "normal", contrast = "dummy", parameters = list(mean = 0, sd = 1))
  p16 <- prior_factor(distribution = "beta",   contrast = "dummy", parameters = list(alpha = 1, beta = 1))

  expect_equal(utils::capture.output(print(p14)),  "orthonormal contrast: mNormal(0, 1)")
  expect_equal(utils::capture.output(print(p15)),  "treatment contrast: Normal(0, 1)")
  expect_equal(utils::capture.output(print(p16)),  "treatment contrast: Beta(1, 1)")

  # check plot names
  empty_plot <- function(){
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  }
  expect_doppelganger("priors-print-1", function(){
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

  expect_doppelganger("priors-print-2", function(){
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

  expect_doppelganger("priors-print-3", function(){
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
  })
})
