context("Prior plot function")

test_that("Prior plot (simple) function works", {

  # check the default options
  p1 <- prior("normal", list(0, 1))

  vdiffr::expect_doppelganger("priors-plot-1-1", function()plot(p1))
  vdiffr::expect_doppelganger("priors-plot-1-2", function()plot(p1, short_name = TRUE))
  vdiffr::expect_doppelganger("priors-plot-1-3", function()plot(p1, parameter_names = TRUE))
  vdiffr::expect_doppelganger("priors-plot-1-4", function()plot(p1, xlab = "xlab", ylab = "ylab", main = "main"))
  vdiffr::expect_doppelganger("priors-plot-1-5", function()plot(p1, lwd = 3, lty = 3, col = "blue"))
  vdiffr::expect_doppelganger("priors-plot-1-6", function()plot(p1, cex.lab = 2, col.lab = "blue", cex.axis = .5, col.axis = "red", cex.main = 1.75, col.main = "green", main = "main"))
  vdiffr::expect_doppelganger("priors-plot-1-7", function()plot(p1, par_name = "name"))
  vdiffr::expect_doppelganger("priors-plot-1-8", function()plot(p1, par_name = bquote(mu)))

  vdiffr::expect_doppelganger("priors-plot-2-1", plot(p1, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-2", plot(p1, short_name = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-3", plot(p1, parameter_names = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-4", plot(p1, xlab = "xlab", ylab = "ylab", main = "main", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-5", plot(p1, lwd = 3, lty = 3, col = "blue", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-6", plot(p1, cex.lab = 2, col.lab = "blue", cex.axis = .5, col.axis = "red", cex.main = 1.75, col.main = "green", main = "main", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-7", plot(p1, par_name = "name", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-2-8", plot(p1, par_name = bquote(mu), plot_type = "ggplot"))

  # check dealing with truncation and range changes
  set.seed(1)
  p2 <- prior("Cauchy", list(0, 1), list(0, Inf))
  vdiffr::expect_doppelganger("priors-plot-3-1", function()plot(p2))
  vdiffr::expect_doppelganger("priors-plot-3-2", function()plot(p2, xlim = c(0, 3)))
  vdiffr::expect_doppelganger("priors-plot-3-3", function()plot(p2, x_seq = seq(-1, 2, .01)))
  vdiffr::expect_doppelganger("priors-plot-3-4", function()plot(p2, x_range_quant = .10))
  vdiffr::expect_doppelganger("priors-plot-3-5", function()plot(p2, force_samples = TRUE, xlim = c(0, 10)))

  # check transformations
  vdiffr::expect_doppelganger("priors-plot-4-1", function()plot(p1, transformation = "tanh", main = "tanh(x)"))
  vdiffr::expect_doppelganger("priors-plot-4-2", function()plot(p1, transformation = "exp", main = "exp(x)"))
  vdiffr::expect_doppelganger("priors-plot-4-3", function()plot(p1, transformation = "lin", transformation_arguments = list(a = -5, b = 0.5), main = "-5 + 0.5x"))

  # check bounded distributions from both sides
  p3 <- prior("beta", list(3, 2))
  p4 <- prior("uniform", list(.4, 2))
  vdiffr::expect_doppelganger("priors-plot-5-1", function()plot(p3))
  vdiffr::expect_doppelganger("priors-plot-5-2", function()plot(p4))

})

test_that("Prior plot (discrete) function works", {

  # check the default options
  p1 <- prior("bernoulli", list(.33))

  vdiffr::expect_doppelganger("priors-plot-15-1", function()plot(p1))
  vdiffr::expect_doppelganger("priors-plot-15-2", function()plot(p1, short_name = TRUE))
  vdiffr::expect_doppelganger("priors-plot-15-3", function()plot(p1, parameter_names = TRUE))
  vdiffr::expect_doppelganger("priors-plot-15-4", function()plot(p1, lwd = 3, lty = 3, col = "blue"))
  vdiffr::expect_doppelganger("priors-plot-15-5", function()plot(p1, par_name = "name"))

  vdiffr::expect_doppelganger("priors-plot-15-6",  plot(p1, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-7",  plot(p1, short_name = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-8",  plot(p1, parameter_names = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-9",  plot(p1, xlab = "xlab", ylab = "ylab", main = "main", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-10", plot(p1, lwd = 3, lty = 3, col = "blue", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-11", plot(p1, par_name = bquote(mu), plot_type = "ggplot"))

})

test_that("Prior plot (point) function works", {
  # check point distributions
  p5 <- prior("point", list(1.5))
  vdiffr::expect_doppelganger("priors-plot-6-1", function()plot(p5))
  vdiffr::expect_doppelganger("priors-plot-6-2", function()plot(p5, short_name = TRUE))
  vdiffr::expect_doppelganger("priors-plot-6-3", function()plot(p5, parameter_names = TRUE))
  vdiffr::expect_doppelganger("priors-plot-6-4", function()plot(p5, xlab = "xlab", ylab = "ylab", main = "main"))
  vdiffr::expect_doppelganger("priors-plot-6-5", function()plot(p5, lwd = 3, lty = 3, col = "blue"))
  vdiffr::expect_doppelganger("priors-plot-6-6", function()plot(p5, cex.lab = 2, col.lab = "blue", cex.axis = .5, col.axis = "red", cex.main = 1.75, col.main = "green", main = "main"))
  vdiffr::expect_doppelganger("priors-plot-6-7", plot(p5, plot_type = "ggplot"))
})

test_that("Prior plot (spike and slab) function works", {

  # check the default options
  p1 <- prior("bernoulli", list(.33))

  vdiffr::expect_doppelganger("priors-plot-15-1", function()plot(p1))
  vdiffr::expect_doppelganger("priors-plot-15-2", function()plot(p1, short_name = TRUE))
  vdiffr::expect_doppelganger("priors-plot-15-3", function()plot(p1, parameter_names = TRUE))
  vdiffr::expect_doppelganger("priors-plot-15-4", function()plot(p1, lwd = 3, lty = 3, col = "blue"))
  vdiffr::expect_doppelganger("priors-plot-15-5", function()plot(p1, par_name = "name"))

  vdiffr::expect_doppelganger("priors-plot-15-6",  plot(p1, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-7",  plot(p1, short_name = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-8",  plot(p1, parameter_names = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-9",  plot(p1, xlab = "xlab", ylab = "ylab", main = "main", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-10", plot(p1, lwd = 3, lty = 3, col = "blue", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-15-11", plot(p1, par_name = bquote(mu), plot_type = "ggplot"))

})

test_that("Prior plot (weightfunction) function works", {

  set.seed(1)
  # check weightfunctions
  p7  <- prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  p8  <- prior_weightfunction("one.sided", list(c(0.05, .95), c(1, 1), c(1, 1)))
  p9  <- prior_weightfunction("two.sided", list(c(0.05), c(1, 1)))
  p10 <- prior_weightfunction("one.sided.fixed", list(c(0.10), c(1, .7)))

  vdiffr::expect_doppelganger("priors-plot-8-1", function()plot(p7))
  vdiffr::expect_doppelganger("priors-plot-8-2", function()plot(p8))
  vdiffr::expect_doppelganger("priors-plot-8-3", function()plot(p9))
  vdiffr::expect_doppelganger("priors-plot-8-4", function()plot(p10))
  vdiffr::expect_doppelganger("priors-plot-8-5", plot(p7, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-8-6", function()plot(p7, individual =  T))
  vdiffr::expect_doppelganger("priors-plot-8-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    plot(p7, individual =  T, show_figures = NULL)
  })
  vdiffr::expect_doppelganger("priors-plot-8-8", function()plot(p7, rescale_x = TRUE))
})

test_that("Prior plot (PET-PEESE) function works", {

  p6   <- prior_PET("normal", list(1, 1))
  p6.1 <- prior_PEESE("Cauchy", list(0, .1))
  vdiffr::expect_doppelganger("priors-plot-7-1", function()plot(p6))
  vdiffr::expect_doppelganger("priors-plot-7-2", plot(p6, plot_type = "ggplot"))

  vdiffr::expect_doppelganger("priors-plot-7-1-1", function()plot(p6.1))
  vdiffr::expect_doppelganger("priors-plot-7-2-1", plot(p6.1, plot_type = "ggplot"))

})

test_that("Prior plot (orthonormal) function works", {

  p11   <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  p12   <- prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "orthonormal")
  p20   <- prior_factor("point", list(location = 0), contrast = "orthonormal")
  p11.5 <- p11.3 <- p11.2 <- p11
  p12.9 <- p12
  p20.3 <- p20.5 <- p20
  p11.2$parameters$K <- 2
  p11.3$parameters$K <- 3
  p11.5$parameters$K <- 5
  p12.9$parameters$K <- 9
  p20.3$parameters$K <- 3
  p20.5$parameters$K <- 5

  vdiffr::expect_doppelganger("priors-plot-11-2", function()plot(p11.2))
  vdiffr::expect_doppelganger("priors-plot-11-3", plot(p11.3, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-11-5", function()plot(p11.5))
  vdiffr::expect_doppelganger("priors-plot-12-9", function()plot(p12.9))
  vdiffr::expect_doppelganger("priors-plot-20-3", function()plot(p20.3))
  vdiffr::expect_doppelganger("priors-plot-20-5", plot(p20.5, plot_type = "ggplot"))

})

test_that("Prior plot (treatment) function works", {

  p13   <- prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
  p14   <- prior_factor("beta", list(alpha = 2, beta = 3), contrast = "treatment")
  p21   <- prior_factor("spike", list(location = 1), contrast = "treatment")

  vdiffr::expect_doppelganger("priors-plot-13-1", function()plot(p13))
  vdiffr::expect_doppelganger("priors-plot-13-2", plot(p13, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-14",   function()plot(p14))
  vdiffr::expect_doppelganger("priors-plot-21",   function()plot(p21))

})

test_that("Prior plot (independent) function works", {

  p15   <- prior_factor("gamma", list(shape = 2, rate = 3), contrast = "independent")
  p16   <- prior_factor("uniform", list(a = -0.5, b = 1), contrast = "independent")
  p22   <- prior_factor("spike", list(location = 0), contrast = "independent")

  vdiffr::expect_doppelganger("priors-plot-16-1", function()plot(p15))
  vdiffr::expect_doppelganger("priors-plot-16-2", plot(p15, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("priors-plot-17",   function()plot(p16))
  vdiffr::expect_doppelganger("priors-plot-22",   function()plot(p22))

})

test_that("Prior plot (meandif) function works", {

  p18   <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
  p19   <- prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "meandif")
  p23   <- prior_factor("point", list(location = 0), contrast = "meandif")
  p18.5 <- p18.3 <- p18.2 <- p18
  p19.5 <- p19.3 <- p19.2 <- p19
  p23.3 <- p23.5 <- p23
  p18.2$parameters$K <- 2
  p18.3$parameters$K <- 3
  p18.5$parameters$K <- 5
  p19.2$parameters$K <- 2
  p19.3$parameters$K <- 3
  p19.5$parameters$K <- 5

  vdiffr::expect_doppelganger("priors-plot-18", function(){

    plot(p18.2)
    lines(p18.3, col = "blue", lty = 2, lwd = 2)
    lines(p18.5, col = "red", lty = 3, lwd = 2)

  })

  vdiffr::expect_doppelganger("priors-plot-19", function(){

    plot(p19.2)
    lines(p19.3, col = "blue", lty = 2, lwd = 2)
    lines(p19.5, col = "red", lty = 3, lwd = 2)

  })

  vdiffr::expect_doppelganger("priors-plot-23", function(){

    plot(p23.3)
    lines(p23.5, col = "red", lty = 3, lwd = 2)

  })

})
