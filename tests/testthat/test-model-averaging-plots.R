context("Model-averaging plot functions")
set.seed(1)

test_that("helper functions", {

  # join duplicate
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1)),
    p4  = prior("normal", list(0, 1))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1), prior_odds = 2),
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
    p1  = prior("normal", list(0, 1), prior_odds = 3),
    p2  = prior("lognormal", list(0, .5), prior_odds = 2),
    p3  = prior("point", list(1))
  ))
})


test_that("prior plot functions (simple)", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior("normal", list(0, 1))
  )

  expect_doppelganger("model-averaging-plot-prior-simple-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike
  prior_list       <- list(
    p1 = prior("spike", list(.5))
  )

  expect_doppelganger("model-averaging-plot-prior-simple-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior (+ check truncation)
  prior_list       <- list(
    p1 = prior("normal", list(0, 1), truncation = list(0, Inf)),
    p2 = prior("normal", list(0, 1.001), truncation = list(0, Inf))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-7", function(){
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-simple-9", function(){
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })
  expect_doppelganger("model-averaging-plot-prior-simple-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })

  # and more spikes
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(-5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-11", function(){
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # verify aggregation
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(.5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-13", function(){
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, lwd = 2)
    lines_prior_list(prior_list, lty = 2, col = "red", lwd = 2)
  })
})


test_that("prior plot functions (PET-PEESE)", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_PET("cauchy",   list(0, 1))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = ggplot2::alpha("red", .20), n_samples = 1000, n_points = 50)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = ggplot2::alpha("blue", .20))
  })
  prior_list       <- list(
    p1 = prior_PEESE("cauchy",   list(0, 2))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50)
  })

  # spike
  prior_list       <- list(
    p1 = prior_PET("point", list(.1))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, n_samples = 1000, n_points = 50)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  prior_list       <- list(
    p1 = prior_PEESE("point", list(.05))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, n_samples = 1000, n_points = 50) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, n_samples = 1000, n_points = 50)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    PET1  = prior_PET("cauchy",   list(0, 1)),
    PET2  = prior_PET("cauchy",   list(0, 1.001))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50)
    lines(prior_list$PET1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  ### the prior joining should give the same prior
  prior_list       <- list(
    PEESE1 = prior_PEESE("cauchy",   list(0, 1)),
    PEESE2 = prior_PEESE("cauchy",   list(0, 1.001))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", col.fill = scales::alpha("red", .20), lwd = 4, n_samples = 1000, n_points = 50) + geom_prior_list(prior_list, col = "blue", col.fill = scales::alpha("blue", .20), lwd = 3, lty = 2, n_samples = 1000, n_points = 50)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50)
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50)
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-9", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, xlab = "xlab", ylab = "ylab", main = "main")
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-10", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main")
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5)),
    p3 = prior_none(prior_odds = 3)
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5))
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50)
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5)) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50)
  })
})


test_that("prior plot functions (weightfunctions)", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1)))
  )

  expect_doppelganger("model-averaging-plot-prior-wf-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-prior-wf-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })

  # spike
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)))
  )

  expect_doppelganger("model-averaging-plot-prior-wf-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1))),
    p2 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1.0001)))
  )
  expect_doppelganger("model-averaging-plot-prior-wf-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-prior-wf-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.025), c(1, 1))),
    p2 = prior_weightfunction("two.sided", list(c(.05),  c(1, 1))),
    p3 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)), prior_odds = 10)
  )
  expect_doppelganger("model-averaging-plot-prior-wf-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-wf-9", function(){
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", main = "main")
  })
  expect_doppelganger("model-averaging-plot-prior-wf-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main")
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.5), c(1, 1))),
    p2 = prior("normal", list(0, 1)),
    p3 = prior_none()
  )
  expect_doppelganger("model-averaging-plot-prior-wf-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

})
