context("Model-averaging plot functions")
set.seed(1)

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

# skip the rest as it takes too long
skip_on_cran()

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

  set.seed(1)
  data <- NULL
  priors_list0 <- list(
    m  = prior("spike", list(0)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  priors_list1 <- list(
    m  = prior("normal", list(0, .3)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
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

  set.seed(1)
  data <- NULL
  priors_list0 <- list(
    mu    = prior("spike", list(0)),
    PET   = prior_PET("normal", list(0, .2))
  )
  priors_list1 <- list(
    mu    = prior("spike", list(0)),
    PEESE = prior_PEESE("normal", list(0, .8))
  )
  model_syntax <-
    "model
    {
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("mu", "PET", "PEESE"), is_null_list = list("mu" = c(T, T), "PET" = c(F,T),  "PEESE" = c(T,F)), seed = 1)



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
  priors_list2 <- list(
    mu = prior("normal", list(.2, .2), prior_weights = 4)
  )
  # fit the models
  fit2 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)
  # automatically mix posteriors
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

  set.seed(1)
  data <- NULL
  priors_list0 <- list(
    omega = prior_weightfunction("one.sided", list(c(.025), c(1, 1)))
  )
  priors_list1 <- list(
    omega = prior_weightfunction("two.sided", list(c(.05),  c(1, 1)))
  )
  model_syntax <-
    "model
    {
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = "omega", is_null_list = list("omega" = c(F,F)), seed = 1)



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
  priors_list2 <- list(
    mu = prior("normal", list(0, .8), prior_weights = 4)
  )
  # fit the models
  fit2 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)
  # automatically mix posteriors
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

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac3o = factor(rep(c("A", "B", "C"), 40), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(120, .4  + ifelse(data_formula$x_fac3o == "A", 0.0, ifelse(data_formula$x_fac3o == "B", -0.5, 0.5)), 1),
    N = 120
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ 1)
  formula_list1 <- list(mu = ~ x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 0.5))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

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
    plot_posterior(mixed_posteriors, "mu_x_fac3o", col = c("red", "green", "blue"), dots_prior = list(col = "grey"), prior = TRUE)
  })

  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-o-1", plot_posterior(mixed_posteriors, "mu_x_fac3o", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-o-2", plot_posterior(mixed_posteriors, "mu_x_fac3o", col = c("red", "green", "blue"), plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-o-3", plot_posterior(mixed_posteriors, "mu_x_fac3o", lty = c(2, 3, 4), col = "blue", lwd = 2, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-o-4", plot_posterior(mixed_posteriors, "mu_x_fac3o", legend = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-o-5", plot_posterior(mixed_posteriors, "mu_x_fac3o", col = c("red", "green", "blue"), prior = TRUE, plot_type = "ggplot"))

})


test_that("posterior plot functions (treatment) work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac2t = factor(rep(c("A", "B"), 60), levels = c("A", "B"))
  )
  data <- list(
    y = rnorm(120, .4  + ifelse(data_formula$x_fac2t == "A", 0.0, 0.5), 1),
    N = 120
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ 1)
  formula_list1 <- list(mu = ~ x_fac2t)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac2t"         = prior_factor("normal", contrast = "treatment", list(0, 0.5))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_intercept", "mu_x_fac2t"),
    is_null_list = list(
      "mu_intercept"  = c(TRUE,  TRUE),
      "mu_x_fac2t"    = c(FALSE, TRUE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", col = "red")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-3", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", lty = 2, col = "blue", lwd = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", legend = FALSE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-t-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2t", col = "red", prior = TRUE)
  })

  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-t-1", plot_posterior(mixed_posteriors, "mu_x_fac2t", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-t-2", plot_posterior(mixed_posteriors, "mu_x_fac2t", col = "red", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-t-3", plot_posterior(mixed_posteriors, "mu_x_fac2t", lty = 2, col = "blue", lwd = 2, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-t-4", plot_posterior(mixed_posteriors, "mu_x_fac2t", legend = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-t-5", plot_posterior(mixed_posteriors, "mu_x_fac2t", col = "red", prior = TRUE, plot_type = "ggplot"))

})


test_that("posterior plot functions (independent) work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac2i = factor(rep(c("A", "B"), 150), levels = c("A", "B"))
  )
  data <- list(
    y = rnorm(300, ifelse(data_formula$x_fac2i == "A", 0.0, 0.5), 1),
    N = 300
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_fac2i - 1)
  formula_list1 <- list(mu = ~ x_fac2i - 1)

  formula_prior_list0 <- list(
    mu    = list(
      "x_fac2i" = prior_factor("spike", contrast = "independent", list(0))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "x_fac2i" = prior_factor("normal", contrast = "independent", list(0, 1/4))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_fac2i"),
    is_null_list = list(
      "mu_x_fac2i"    = c(FALSE, TRUE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-i-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2i")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-i-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2i", col = c("green", "yellow"))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-i-3", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac2i", col = "red", prior = TRUE)
  })

  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-i-1", plot_posterior(mixed_posteriors, "mu_x_fac2i", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-i-2", plot_posterior(mixed_posteriors, "mu_x_fac2i", prior = TRUE, plot_type = "ggplot"))

})


test_that("posterior plot functions (meandif) work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac3md = factor(rep(c("A", "B", "C"), 40), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(120, .4  + ifelse(data_formula$x_fac3md == "A", 0.0, ifelse(data_formula$x_fac3md == "B", -0.5, 0.5)), 1),
    N = 120
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ 1)
  formula_list1 <- list(mu = ~ x_fac3md)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3md"         = prior_factor("mnormal", contrast = "meandif", list(0, 0.5))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_intercept", "mu_x_fac3md"),
    is_null_list = list(
      "mu_intercept"  = c(TRUE,  TRUE),
      "mu_x_fac3md"    = c(FALSE, TRUE)
    ),
    seed = 1, n_samples = 10000)

  vdiffr::expect_doppelganger("model-averaging-plot-posterior-md-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3md")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-md-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3md", col = c("red", "green", "blue"))
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-md-3", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3md", lty = c(2, 3, 4), col = "blue", lwd = 2)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-md-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3md", legend = FALSE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-posterior-md-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "mu_x_fac3md", col = c("red", "green", "blue"), dots_prior = list(col = "grey"), prior = TRUE)
  })

  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-md-1", plot_posterior(mixed_posteriors, "mu_x_fac3md", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-md-2", plot_posterior(mixed_posteriors, "mu_x_fac3md", col = c("red", "green", "blue"), plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-md-3", plot_posterior(mixed_posteriors, "mu_x_fac3md", lty = c(2, 3, 4), col = "blue", lwd = 2, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-md-4", plot_posterior(mixed_posteriors, "mu_x_fac3md", legend = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("model-averaging-ggplot-posterior-md-5", plot_posterior(mixed_posteriors, "mu_x_fac3md", col = c("red", "green", "blue"), prior = TRUE, plot_type = "ggplot"))

})


test_that("models plot functions work", {

  set.seed(1)
  ### prior distribution related functions
  p0 <- prior(distribution = "point",  parameters = list(location = 0))
  p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
  p2 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1), truncation = list(0, Inf))

  data <- list(
    x = rnorm(10),
    N = 10
  )

  ## create and fit models
  priors_list0 <- list(mu  = p0)
  priors_list1 <- list(mu  = p1)
  priors_list2 <- list(tau = p2)

  # define likelihood for the data
  model_syntax <-
    "model{
    for(i in 1:N){
      x[i] ~ dnorm(mu, 1)
    }
  }"
  model_syntax2 <-
    "model{
    for(i in 1:N){
      x[i] ~ dnorm(0, pow(tau, -2))
    }
  }"

  # define log posterior for bridge sampling
  log_posterior <- function(parameters, data){
    sum(dnorm(data$x, parameters$mu, 1, log = TRUE))
  }
  log_posterior2 <- function(parameters, data){
    sum(dnorm(data$x, 0, parameters$tau, log = TRUE))
  }
  # fit the models
  fit0 <- JAGS_fit(model_syntax,  data, priors_list0, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 0)
  fit1 <- JAGS_fit(model_syntax,  data, priors_list1, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 1)
  fit2 <- JAGS_fit(model_syntax2, data, priors_list2, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 2)
  # get marginal likelihoods
  marglik0 <- list(
    logml = sum(dnorm(data$x, mean(p0), 1, log = TRUE))
  )
  class(marglik0) <- "bridge"
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior2, data = data, prior_list = priors_list2)
  ## create an object with the models
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
  )
  # compare and summarize the models
  models            <- models_inference(models)
  inference         <- ensemble_inference(model_list = models, parameters = c("mu", "tau"), is_null_list = list("mu" = c(1, 3), "tau" = c(1, 2)))
  mixed_posteriors  <- mix_posteriors(model_list = models, parameters = c("mu", "tau"), is_null_list = list("mu" = c(1, 3), "tau" = c(1, 2)), seed = 1)


  vdiffr::expect_doppelganger("model-averaging-plot-models-1", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-2", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-3", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", prior = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-4", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-5", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", conditional = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-6", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE, conditional = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-7", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-8", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE, plot_type = "ggplot")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-9", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot", y_axis2 = FALSE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-10", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot", show_estimates = FALSE)
  })

})


test_that("models plot functions work (formulas + factors)", {


  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(60, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_fac2t)
  formula_list1 <- list(mu = ~ x_cont1 + x_fac3t)
  formula_list2 <- list(mu = ~ x_fac3o)
  formula_list3 <- list(mu = ~ x_cont1 * x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac2t"         = prior_factor("normal",  contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list2 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )
  formula_prior_list3 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
      "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)
  fit2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2, seed = 3)
  fit3 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3, seed = 4)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)
  marglik2 <- JAGS_bridgesampling(
    fit2, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2)
  marglik3 <- JAGS_bridgesampling(
    fit3, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3)


  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1), prior_weights = 1),
    list(fit = fit2, marglik = marglik2, fit_summary = runjags_estimates_table(fit2), prior_weights = 1),
    list(fit = fit3, marglik = marglik3, fit_summary = runjags_estimates_table(fit3), prior_weights = 1)
  )
  models <- models_inference(models)


  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac3t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac3o"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1__xXx__x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"               = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"               = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac3t"               = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"               = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1__xXx__x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)



  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-1", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_cont1")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-2", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac2t")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-3", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac2t", prior = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(2, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3t")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(2, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3t", prior = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-6", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3o")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3o", prior = TRUE)
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_cont1__xXx__x_fac3o")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-9", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_cont1__xXx__x_fac3o", prior = TRUE)
  })

  # ggplot versions
  p1 <- plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_cont1", plot_type = "ggplot")
  p2 <- plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3o", prior = TRUE, plot_type = "ggplot")

  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-10", p1)
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-11", p2[[1]])
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-12", p2[[2]])
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-13", p2[[3]])

})


test_that("models plot functions work (formulas + spike factors)", {


  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac3md = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, ifelse(data_formula$x_fac3md == "A", 0.0, ifelse(data_formula$x_fac3md == "B", -0.2, 0.4))),
    N = 300
  )


  formula_list <- list(
    mu    = ~ x_fac3md
  )
  formula_data_list <- list(
    mu    = data_formula
  )
  formula_prior_list0 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_fac3md"  = prior_factor("spike", contrast = "meandif", list(0))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_fac3md"  = prior_factor("mnormal", contrast = "meandif", list(0, 0.25))
    )
  )

  prior_list <- list(
    sigma = prior("lognormal", list(0, 1))
  )
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  log_posterior <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = parameters[["sigma"]], log = TRUE)))
  }

  marglik0 <- JAGS_bridgesampling(
    fit                = fit0,
    log_posterior      = log_posterior,
    data               = data,
    prior_list         = prior_list,
    formula_list       = formula_list,
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit                = fit1,
    log_posterior      = log_posterior,
    data               = data,
    prior_list         = prior_list,
    formula_list       = formula_list,
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list1)


  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0, remove_spike_0 = FALSE, transform_factors = TRUE), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1, remove_spike_0 = FALSE, transform_factors = TRUE), prior_weights = 1)
  )
  models <- models_inference(models)

    inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE,  FALSE)
    ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)



  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-s-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3md")
  })
  vdiffr::expect_doppelganger("model-averaging-plot-models-formula-s-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3, 1))
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu_x_fac3md", prior = TRUE)
  })

})


test_that("posterior plot model averaging based on complex single JAGS models  (formulas + spike factors + mixture)", {

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, -0.15 + 0.20 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.2)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )

  # create model with mix of a formula and free parameters ---
  formula_list1 <- list(
    mu    = ~ x_cont1 + x_fac2t + x_fac3t
  )
  formula_data_list1 <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"  = prior_mixture(
        list(
          prior("spike",   list(0),       prior_weights = 2),
          prior("normal",  list(-1, 0.5), prior_weights = 1),
          prior("normal",  list( 1, 0.5), prior_weights = 1)
        ),
        is_null = c(T, F, F)
      ),
      "x_cont1"    = prior_spike_and_slab(prior("normal",  list(0, 1), prior_weights = 1)),
      "x_fac2t"    = prior_mixture(list(
          prior("spike", list(0), prior_weights = 1),
          prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
        ),
        is_null = c(T, F)
      ),
      "x_fac3t"    = prior_mixture(list(
          prior("spike", list(0), prior_weights = 1),
          prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
        ),
        is_null = c(T, F)
      )
    )
  )

  attr(formula_prior_list1$mu$x_cont1, "multiply_by") <- "sigma"
  prior_list1 <- list(
    "sigma" = prior_mixture(
      list(
        prior("normal",    list(0, 1), truncation = list(0, Inf)),
        prior("lognormal", list(0, 1))
      ),
      components = c("normal", "lognormal")
    ),
    "bias"  = prior_mixture(list(
      prior_none(prior_weights = 1),
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/3),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/3),
      prior_PET("normal", list(0, 1), prior_weights = 1/3)
    ), is_null = c(T, F, F, F))
  )
  model_syntax1 <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  if("RoBMA" %in% rownames(installed.packages())){
    require("RoBMA")
  }else{
    skip()
  }

  fit1 <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1)

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

  # TODO: implement
  # vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-bias", function(){
  #   oldpar <- graphics::par(no.readonly = TRUE)
  #   on.exit(graphics::par(mar = oldpar[["mar"]]))
  #   par(mar = c(4, 4, 1, 4), mfrow = c(2, 2))
  #   plot_posterior(mixed_posteriors, "PET", prior = T, dots_prior = list(col = "grey"))
  # })


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

#   # TODO: at some point
#   vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-PETPEESE", function(){
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(mar = oldpar[["mar"]]))
#     par(mar = c(4, 4, 1, 4))
#     plot_posterior(mixed_posteriors_conditional5a, parameter = "PETPEESE", prior = TRUE, col = "black", col.fill = ggplot2::alpha("grey", 0.2),  dots_prior = list(col = "red", col.fill = ggplot2::alpha("red", 0.5)))
#   })
#
#   vdiffr::expect_doppelganger("model-averaging-plot-ss-posterior-PETPEESE-con", function(){
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(mar = oldpar[["mar"]]))
#     par(mar = c(4, 4, 1, 4))
#     plot_posterior(mixed_posteriors_conditional6b, parameter = "PETPEESE", prior = TRUE, col = "black", col.fill = ggplot2::alpha("grey", 0.2),  dots_prior = list(col = "red", col.fill = ggplot2::alpha("red", 0.5)))
#   })

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
