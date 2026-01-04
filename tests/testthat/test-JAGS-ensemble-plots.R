# ============================================================================ #
# TEST FILE: JAGS Ensemble Plot Functions
# ============================================================================ #
#
# PURPOSE:
#   Visual regression tests for ensemble plot functions (plot_prior_list,
#   plot_posterior, plot_models, etc.).
#
# DEPENDENCIES:
#   - vdiffr: Visual regression testing
#   - rjags, bridgesampling: For tests using pre-fitted models
#   - common-functions.R: temp_fits_dir, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - First section (prior plots): Can run on CRAN (pure R)
#   - Second section (posterior plots): skip_if_no_fits(), skip_if_not_installed()
#   - skip_on_os(): Multivariate sampling differs across OSes
#
# MODELS/FIXTURES:
#   - fit_simple_spike, fit_simple_normal, fit_summary*, fit_marginal_*
#
# TAGS: @evaluation, @visual, @JAGS, @model-averaging, @plots
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-ensemble-plots")
source(testthat::test_path("common-functions.R"))

test_that("helper functions work", {

  # join duplicate
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1)),
    p4  = prior("normal", list(0, 1))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1), prior_weights = 2),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  ))


  # no duplicate
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  ))


  # multiple duplicates
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1)),
    p1  = prior("normal", list(0, 1)),
    p4  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1), prior_weights = 3),
    p2  = prior("lognormal", list(0, .5), prior_weights = 2),
    p3  = prior("point", list(1))
  ))
})


test_that("prior plot functions (simple) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior("normal", list(0, 1))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike
  prior_list       <- list(
    p1 = prior("spike", list(.5))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior (+ check truncation)
  prior_list       <- list(
    p1 = prior("normal", list(0, 1), truncation = list(0, Inf)),
    p2 = prior("normal", list(0, 1.001), truncation = list(0, Inf))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  ### mixtures
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # with additional settings
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-9", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })

  # and more spikes
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(-5))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-11", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # verify aggregation
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(.5))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-simple-13", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, lwd = 2)
    lines_prior_list(prior_list, lty = 2, col = "red", lwd = 2)
  })
})

# ============================================================================ #
# SECTION: Tests requiring pre-fitted models (skip on CRAN)
# ============================================================================ #
skip_on_cran()
skip_if_no_fits()
skip_if_not_installed("vdiffr")

test_that("prior plot functions (PET-PEESE) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_PET("cauchy",   list(0, 1))
  )
  prior_list_mu   <- list(
    m1 = prior("spike", list(0))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = ggplot2::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = ggplot2::alpha("blue", .20))
  })
  prior_list       <- list(
    p1 = prior_PEESE("cauchy",   list(0, 2))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })

  # spike
  prior_list       <- list(
    p1 = prior_PET("point", list(.1))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  prior_list       <- list(
    p1 = prior_PEESE("point", list(.05))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    PET1  = prior_PET("cauchy",   list(0, 1)),
    PET2  = prior_PET("cauchy",   list(0, 1.001))
  )
  prior_list_mu   <- list(
    m1 = prior("spike", list(0)),
    m2 = prior("spike", list(0))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$PET1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  ### the prior joining should give the same prior
  prior_list       <- list(
    PEESE1 = prior_PEESE("cauchy",   list(0, 1)),
    PEESE2 = prior_PEESE("cauchy",   list(0, 1.001))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", col.fill = scales::alpha("red", .20), lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", col.fill = scales::alpha("blue", .20), lwd = 3, lty = 2, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })

  # with additional settings
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-9", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, xlab = "xlab", ylab = "ylab", main = "main", prior_list_mu = prior_list_mu)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-10", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main", prior_list_mu = prior_list_mu)
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5)),
    p3 = prior_none(prior_weights = 4)
  )
  prior_list_mu    <- list(
    m1 = prior("spike",  list(0)),
    m2 = prior("spike",  list(0)),
    m3 = prior("normal", list(.3, .15))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5), prior_list_mu = prior_list_mu)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-PETPEESE-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5), prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
})

test_that("prior plot functions (weightfunctions) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1)))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })

  # spike
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1))),
    p2 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1.0001)))
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.025), c(1, 1))),
    p2 = prior_weightfunction("two.sided", list(c(.05),  c(1, 1))),
    p3 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)), prior_weights = 10)
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

  # with additional settings
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-9", function(){
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", main = "main")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main")
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.5), c(1, 1))),
    p2 = prior_none(),
    p3 = prior_none()
  )
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-wf-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

})

test_that("prior plot functions (orthonormal) work", {

  ### simple cases
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  prior_list$p1$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) +
      geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike & slab mixture
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal"),
    p2 = prior("spike", list(0))
  )
  prior_list$p1$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
  })

  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal"),
    p2 = prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "orthonormal"),
    p3 = prior("spike", list(0))
  )
  prior_list$p1$parameters$K <- 3
  prior_list$p2$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-orthonormal-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 2) +
      geom_prior_list(prior_list, col = "blue", lwd = 1, lty = 2)
  })

})

test_that("prior plot functions (treatment) work", {

  ### simple cases
  prior_list       <- list(
    p1 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "treatment")
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) +
      geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike & slab mixture
  prior_list       <- list(
    p1 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "treatment"),
    p2 = prior("spike", list(0))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
  })

  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 2)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment"),
    p2 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "treatment"),
    p3 = prior("spike", list(0))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-treatment-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 2) +
      geom_prior_list(prior_list, col = "blue", lwd = 1, lty = 2)
  })

})

test_that("prior plot functions (independent) work", {

  ### simple cases
  prior_list       <- list(
    p1 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "independent")
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-independent-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike & slab mixture
  prior_list       <- list(
    p1 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "independent"),
    p2 = prior("spike", list(0))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-independent-2", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
  })

  ### mixtures
  prior_list       <- list(
    p1 = prior_factor("normal", list(mean = 0, sd = 1), contrast = "independent"),
    p2 = prior_factor("beta", list(alpha = 4, beta = 5), contrast = "independent"),
    p3 = prior("spike", list(0))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-prior-independent-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-independent-4", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 2) +
      geom_prior_list(prior_list, col = "blue", lwd = 1, lty = 2)
  })

})

test_that("prior plot functions (meandif) work", {

  ### simple cases
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
  )
  prior_list$p1$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) +
      geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike & slab mixture
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif"),
    p2 = prior("spike", list(0))
  )
  prior_list$p1$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
  })

  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif"),
    p2 = prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "meandif"),
    p3 = prior("spike", list(0))
  )
  prior_list$p1$parameters$K <- 3
  prior_list$p2$parameters$K <- 3

  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-prior-meandif-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 2) +
      geom_prior_list(prior_list, col = "blue", lwd = 1, lty = 2)
  })

})

test_that("posterior plot functions (simple) work", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik0 <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  fit1 <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik1 <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)


  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "m", lwd = 2, col = "red", par_name = expression(mu))
    lines_prior_list(attr(mixed_posteriors$m, "prior_list"), col = "blue")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-2", {
    plot_posterior(mixed_posteriors, "m", plot_type = "ggplot", lwd = 2, col = "red") + geom_prior_list(attr(mixed_posteriors$m, "prior_list"), col = "blue")
  })

  # checks truncation
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-3", function(){
    plot_posterior(mixed_posteriors, "s")
    lines_prior_list(attr(mixed_posteriors$s, "prior_list"), col = "blue")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-4", {
    plot_posterior(mixed_posteriors, "s", plot_type = "ggplot")
  })

  # check transformation
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-5", function(){
    plot_posterior(mixed_posteriors, "s", transformation = "exp")
    lines_prior_list(attr(mixed_posteriors$s, "prior_list"), col = "blue", transformation = "exp")
  })

  # prior and posterior
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-6", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "m", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-7", function(){
    plot_posterior(mixed_posteriors, "m", plot_type = "ggplot", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })

  # check transformation
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-simple-8", function(){
    plot_posterior(mixed_posteriors, "s", transformation = "exp", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })
})

test_that("posterior plot functions (PET-PEESE) work", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_pet.RDS"))
  marglik0 <- readRDS(file.path(temp_marglik_dir, "fit_pet.RDS"))
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_peese.RDS"))
  marglik1 <- readRDS(file.path(temp_marglik_dir, "fit_peese.RDS"))

  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("mu", "PET", "PEESE"), is_null_list = list("mu" = c(T, T), "PET" = c(F,T),  "PEESE" = c(T,F)), seed = 1)

  # Reconstruct priors for plotting (since we don't have them in the test scope directly)
  priors_list0 <- list(
    mu    = prior("spike", list(0)),
    PET   = prior_PET("normal", list(0, .2))
  )
  priors_list1 <- list(
    mu    = prior("spike", list(0)),
    PEESE = prior_PEESE("normal", list(0, .8))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-1", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), par_name = "PET-PEESE", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-2", {
    plot_posterior(mixed_posteriors, "PETPEESE", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), ylim = c(0, .5)) + geom_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })

  # check transformation
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-5", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), main = "PET-PEESE (1/2x)")
    lines_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })

  # prior and posterior
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-6", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-7", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # check transformation
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-8", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), lwd = 2, col = "red", n_points = 50, n_samples = 1000, col.fill = scales::alpha("red", .20), prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # add an overhelming missing model
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_missing.RDS"))
  marglik2 <- readRDS(file.path(temp_marglik_dir, "fit_missing.RDS"))

  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, prior_weights = 4)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("mu" ,"PET", "PEESE"), is_null_list = list("mu" = c(T, T, F),"PET" = c(F,T,F),  "PEESE" = c(T,F,F)), seed = 1)
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-PETPEESE-9", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", ylim = c(0, 3), lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

})

test_that("posterior plot functions (weightfunctions) work", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_wf_onesided.RDS"))
  marglik0 <- readRDS(file.path(temp_marglik_dir, "fit_wf_onesided.RDS"))
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_wf_twosided.RDS"))
  marglik1 <- readRDS(file.path(temp_marglik_dir, "fit_wf_twosided.RDS"))

  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = "omega", is_null_list = list("omega" = c(F,F)), seed = 1)

  # Reconstruct priors
  priors_list0 <- list(
    omega = prior_weightfunction("one.sided", list(c(.025), c(1, 1)))
  )
  priors_list1 <- list(
    omega = prior_weightfunction("two.sided", list(c(.05),  c(1, 1)))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-1", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), par_name = "Selection Models", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$omega, priors_list1$omega), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-2", {
    plot_posterior(mixed_posteriors, "omega", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20)) + geom_prior_list(list(priors_list0$omega, priors_list1$omega), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })

  # rescale-x
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-3", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, rescale_x = TRUE, col = "red", col.fill = scales::alpha("red", .20), par_name = "Selection Models", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$omega, priors_list1$omega), rescale_x = TRUE, n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-4", {
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20)) + geom_prior_list(list(priors_list0$omega, priors_list1$omega), rescale_x = TRUE, n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })

  # prior and posterior
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-6", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-7", function(){
    plot_posterior(mixed_posteriors, "omega", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # rescale-x
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-8", function(){
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-9", function(){
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # add an overhelming missing model
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_wf_missing.RDS"))
  marglik2 <- readRDS(file.path(temp_marglik_dir, "fit_wf_missing.RDS"))

  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, prior_weights = 5)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = "omega", is_null_list = list("omega" = c(F,F,F)), seed = 1)
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-wf-10", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

})

test_that("posterior plot functions (orthonormal) work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik0 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik1 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_intercept", "mu_x_fac3o"),
    is_null_list = list(
      "mu_intercept"  = c(TRUE,  TRUE),
      "mu_x_fac3o"    = c(FALSE, TRUE)
    ),
    seed = 1, n_samples = 10000)
  mixed_posteriors2 <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_fac3o"),
    conditional  = TRUE,
    is_null_list = list(
      "mu_x_fac3o"    = c(TRUE, FALSE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-o-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3o")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-o-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3o", col = c("red", "green", "blue"))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-o-3", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3o", lty = c(2, 3, 4), col = "blue", lwd = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-o-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3o", legend = FALSE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-o-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors2, "mu_x_fac3o")
  })
})

test_that("posterior plot functions (treatment) work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment.RDS"))
  # Create dummy marginal likelihood since this model doesn't have one
  marglik0 <- structure(list(logml = -10), class = "bridge")

  # Create a second model with different prior for comparison
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment.RDS"))
  marglik1 <- structure(list(logml = -12), class = "bridge")

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    is_null_list = list(
      "p1"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)
  mixed_posteriors2 <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    conditional  = TRUE,
    is_null_list = list(
      "p1"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "p1")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors2, "p1")
  })
})

test_that("posterior plot functions (independent) work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_factor_independent.RDS"))
  # Create dummy marginal likelihood since this model doesn't have one
  marglik0 <- structure(list(logml = -15), class = "bridge")

  # Create a second model with different prior for comparison
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_factor_independent.RDS"))
  marglik1 <- structure(list(logml = -17), class = "bridge")

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    is_null_list = list(
      "p1[1]"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)

  mixed_posteriors2 <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    conditional  = TRUE,
    is_null_list = list(
      "p1[1]"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-i-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "p1")
  })

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-i-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors2, "p1")
  })
})

test_that("posterior plot functions (meandif) work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit0 <- readRDS(file.path(temp_fits_dir, "fit_factor_meandif.RDS"))
  # Create dummy marginal likelihood since this model doesn't have one
  marglik0 <- structure(list(logml = -20), class = "bridge")

  # Create a second model with different prior for comparison
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_factor_meandif.RDS"))
  marglik1 <- structure(list(logml = -22), class = "bridge")

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    is_null_list = list(
      "p"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)
  mixed_posteriors2 <- mix_posteriors(
    model_list   = models,
    parameters   = c("p1"),
    conditional  = TRUE,
    is_null_list = list(
      "p"  = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-m-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "p1")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-m-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors2, "p1")
  })

})

test_that("posterior plot model averaging based on complex single JAGS models  (formulas + spike factors + mixture)", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit1 <- readRDS(file.path(temp_fits_dir, "fit_complex_mixed.RDS"))

  mixed_posteriors <- as_mixed_posteriors(
    mode       = fit1,
    parameters = names(attr(fit1, "prior_list"))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-intercept", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_intercept", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_cont1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_cont1", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_fac2t", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_fac3t", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3t", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-sigma", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "sigma", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-bias-PET", function(){
    PET <- mixed_posteriors$bias[,"PET",drop=FALSE]
    attributes(PET) <- c(attributes(PET), attributes(mixed_posteriors$bias)[!names(attributes(mixed_posteriors$bias)) %in% c("dimnames", "dim")])
    attr(PET, "prior_list")[!sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(1:sum(!sapply(attr(PET, "prior_list"), is.prior.PET)), function(i) prior("point", list(0)))
    attr(PET, "prior_list")[ sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(attr(PET, "prior_list")[sapply(attr(PET, "prior_list"), is.prior.PET)], function(p) {
      class(p) <- class(p)[!class(p) %in% "prior.PET"]
      return(p)
    })
    plot_posterior(list(PET = PET), "PET", prior = T, dots_prior = list(col = "grey"))
  })


  mixed_posteriors_conditional1 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "mu_intercept",
    conditional = "mu_intercept",
    force_plots = TRUE
  )

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-intercept-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional1, "mu_intercept", prior = TRUE, dots_prior = list(col = "grey"))
  })

  mixed_posteriors_conditional2 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "mu_x_cont1",
    conditional = "mu_x_cont1",
    force_plots = TRUE
  )

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_cont1-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional2, "mu_x_cont1", prior = TRUE, dots_prior = list(col = "grey"))
  })

  mixed_posteriors_conditional3 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "mu_x_fac2t",
    conditional = "mu_x_fac2t",
    force_plots = TRUE
  )

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_fac2t-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional3, "mu_x_fac2t", prior = TRUE, dots_prior = list(col = "grey"))
  })

  mixed_posteriors_conditional4 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "mu_x_fac3t",
    conditional = "mu_x_fac3t",
    force_plots = TRUE
  )

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-x_fac3t-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional4, "mu_x_fac3t", prior = TRUE, dots_prior = list(col = "grey"))
  })


  mixed_posteriors_conditional5a <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias"
  )

  mixed_posteriors_conditional5b <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "bias",
    force_plots = TRUE
  )

  mixed_posteriors_conditional6a <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "PET",
    force_plots = TRUE
  )

  mixed_posteriors_conditional6b <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "omega",
    force_plots = TRUE
  )


  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-weightfunction", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional5a, parameter = "weightfunction", prior = TRUE, col = "black", col.fill = ggplot2::alpha("grey", 0.2),  dots_prior = list(col = "red", col.fill = ggplot2::alpha("red", 0.5)))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-weightfunction-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors_conditional6b, parameter = "weightfunction", prior = TRUE, col = "black", col.fill = ggplot2::alpha("grey", 0.2),  dots_prior = list(col = "red", col.fill = ggplot2::alpha("red", 0.5)))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-PET-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4), mfrow = c(3, 1))
    hist(mixed_posteriors_conditional5a$bias[,"PET"], breaks = 50, col = "grey", main = "PET", xlab = "PET")
    hist(mixed_posteriors_conditional5b$bias[,"PET"], breaks = 50, col = "grey", main = "PET", xlab = "PET")
    hist(mixed_posteriors_conditional6a$bias[,"PET"], breaks = 50, col = "grey", main = "PET", xlab = "PET")
  })

  vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-omega-con", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4), mfrow = c(3, 3))
    hist(mixed_posteriors_conditional5a$bias[,"omega[0,0.025]"], breaks = 50, col = "grey", main = "omega[0,0.025]", xlab = "omega[0,0.025]")
    hist(mixed_posteriors_conditional5a$bias[,"omega[0.025,0.05]"], breaks = 50, col = "grey", main = "omega[0.025,0.05]", xlab = "omega[0.025,0.05]")
    hist(mixed_posteriors_conditional5a$bias[,"omega[0.05,0.975]"], breaks = 50, col = "grey", main = "omega[0.05,0.975]", xlab = "omega[0.05,0.975]")

    hist(mixed_posteriors_conditional5b$bias[,"omega[0,0.025]"], breaks = 50, col = "grey", main = "omega[0,0.025]", xlab = "omega[0,0.025]")
    hist(mixed_posteriors_conditional5b$bias[,"omega[0.025,0.05]"], breaks = 50, col = "grey", main = "omega[0.025,0.05]", xlab = "omega[0.025,0.05]")
    hist(mixed_posteriors_conditional5b$bias[,"omega[0.05,0.975]"], breaks = 50, col = "grey", main = "omega[0.05,0.975]", xlab = "omega[0.05,0.975]")

    hist(mixed_posteriors_conditional6b$bias[,"omega[0,0.025]"], breaks = 50, col = "grey", main = "omega[0,0.025]", xlab = "omega[0,0.025]")
    hist(mixed_posteriors_conditional6b$bias[,"omega[0.025,0.05]"], breaks = 50, col = "grey", main = "omega[0.025,0.05]", xlab = "omega[0.025,0.05]")
    hist(mixed_posteriors_conditional6b$bias[,"omega[0.05,0.975]"], breaks = 50, col = "grey", main = "omega[0.05,0.975]", xlab = "omega[0.05,0.975]")
  })

})

test_that("posterior plot model averaging based on simple single JAGS models  (formulas)", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit1 <- readRDS(file.path(temp_fits_dir, "fit_simple_formula_mixed.RDS"))

  mixed_posteriors <- as_mixed_posteriors(
    mode       = fit1,
    parameters = names(attr(fit1, "prior_list"))
  )

  vdiffr::expect_doppelganger("model-averaging-simple-plot-ss-posterior-intercept", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_intercept", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-simple-plot-ss-posterior-x_cont1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_cont1", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-simple-plot-ss-posterior-x_fac2t", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-simple-plot-ss-posterior-x_fac3t", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3t", prior = T, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-simple-plot-ss-posterior-sigma", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "sigma", prior = T, dots_prior = list(col = "grey"))
  })
})

test_that("posterior plot model averaging based on complex bias mixture model (PET + PEESE + weightfunction)", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  fit1 <- readRDS(file.path(temp_fits_dir, "fit_complex_bias.RDS"))

  mixed_posteriors <- as_mixed_posteriors(
    mode       = fit1,
    parameters = names(attr(fit1, "prior_list"))
  )

  vdiffr::expect_doppelganger("model-averaging-plot-complex-bias-posterior-mu", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu", prior = TRUE, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-complex-bias-posterior-bias-PET", function(){
    PET <- mixed_posteriors$bias[,"PET",drop=FALSE]
    attributes(PET) <- c(attributes(PET), attributes(mixed_posteriors$bias)[!names(attributes(mixed_posteriors$bias)) %in% c("dimnames", "dim")])
    attr(PET, "prior_list")[!sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(1:sum(!sapply(attr(PET, "prior_list"), is.prior.PET)), function(i) prior("point", list(0)))
    attr(PET, "prior_list")[ sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(attr(PET, "prior_list")[sapply(attr(PET, "prior_list"), is.prior.PET)], function(p) {
      class(p) <- class(p)[!class(p) %in% "prior.PET"]
      return(p)
    })
    plot_posterior(list(PET = PET), "PET", prior = TRUE, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-complex-bias-posterior-bias-PEESE", function(){
    PEESE <- mixed_posteriors$bias[,"PEESE",drop=FALSE]
    attributes(PEESE) <- c(attributes(PEESE), attributes(mixed_posteriors$bias)[!names(attributes(mixed_posteriors$bias)) %in% c("dimnames", "dim")])
    attr(PEESE, "prior_list")[!sapply(attr(PEESE, "prior_list"), is.prior.PEESE)] <- lapply(1:sum(!sapply(attr(PEESE, "prior_list"), is.prior.PEESE)), function(i) prior("point", list(0)))
    attr(PEESE, "prior_list")[ sapply(attr(PEESE, "prior_list"), is.prior.PEESE)] <- lapply(attr(PEESE, "prior_list")[sapply(attr(PEESE, "prior_list"), is.prior.PEESE)], function(p) {
      class(p) <- class(p)[!class(p) %in% "prior.PEESE"]
      return(p)
    })
    plot_posterior(list(PEESE = PEESE), "PEESE", prior = TRUE, dots_prior = list(col = "grey"))
  })


  mixed_posteriors_conditional1 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "PET",
    force_plots = TRUE
  )

  mixed_posteriors_conditional2 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "PEESE",
    force_plots = TRUE
  )

  mixed_posteriors_conditional3 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "PETPEESE",
    force_plots = TRUE
  )

  mixed_posteriors_conditional4 <- as_mixed_posteriors(
    mode       = fit1,
    parameters  = "bias",
    conditional = "omega",
    force_plots = TRUE
  )

  vdiffr::expect_doppelganger("model-averaging-plot-complex-bias-conditional-posterior-PET", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))

    PET <- mixed_posteriors_conditional1$bias[,"PET",drop=FALSE]
    attributes(PET) <- c(attributes(PET), attributes(mixed_posteriors_conditional1$bias)[!names(attributes(mixed_posteriors_conditional1$bias)) %in% c("dimnames", "dim")])
    attr(PET, "prior_list")[!sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(1:sum(!sapply(attr(PET, "prior_list"), is.prior.PET)), function(i) prior("point", list(0)))
    attr(PET, "prior_list")[ sapply(attr(PET, "prior_list"), is.prior.PET)] <- lapply(attr(PET, "prior_list")[sapply(attr(PET, "prior_list"), is.prior.PET)], function(p) {
      class(p) <- class(p)[!class(p) %in% "prior.PET"]
      return(p)
    })

    plot_posterior(list(PET = PET), "PET", prior = TRUE, dots_prior = list(col = "grey"))
  })

  vdiffr::expect_doppelganger("model-averaging-plot-complex-bias-conditional-posterior-PEESE", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))

    PEESE <- mixed_posteriors_conditional2$bias[,"PEESE",drop=FALSE]
    attributes(PEESE) <- c(attributes(PEESE), attributes(mixed_posteriors_conditional1$bias)[!names(attributes(mixed_posteriors_conditional1$bias)) %in% c("dimnames", "dim")])
    attr(PEESE, "prior_list")[!sapply(attr(PEESE, "prior_list"), is.prior.PEESE)] <- lapply(1:sum(!sapply(attr(PEESE, "prior_list"), is.prior.PEESE)), function(i) prior("point", list(0)))
    attr(PEESE, "prior_list")[ sapply(attr(PEESE, "prior_list"), is.prior.PEESE)] <- lapply(attr(PEESE, "prior_list")[sapply(attr(PEESE, "prior_list"), is.prior.PEESE)], function(p) {
      class(p) <- class(p)[!class(p) %in% "prior.PEESE"]
      return(p)
    })

    plot_posterior(list(PEESE = PEESE), "PEESE", prior = TRUE, dots_prior = list(col = "grey"))
  })

})






