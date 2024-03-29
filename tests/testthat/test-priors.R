context("Prior distribution functions")

# each test checks that a corresponding prior distribution can be created and the following functions work:
# - random number generator
# - quantile function
# - density function
# - distribution function
# - print function
# - mean and sd functions
test_prior          <- function(prior, skip_moments = FALSE){
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  if(is.prior.discrete(prior)){
    barplot(table(samples)/length(samples), main = print(prior, plot = T), width = 1/(max(samples)+1), space = 0, xlim = c(-0.25, max(samples)+0.25))
  }else if(is.prior.spike_and_slab(prior)){
    xh         <- hist(samples[samples != 0], breaks = 50, plot = FALSE)
    xh$density <- xh$density * mean(samples != 0)
    plot(xh, main = print(prior, plot = T), freq = FALSE)
  }else{
    hist(samples, main = print(prior, plot = T), breaks = 50, freq = FALSE)
  }
  # tests density function
  lines(prior, individual = TRUE)
  # tests quantile function
  if(!is.prior.spike_and_slab(prior)){
    abline(v = quant(prior, 0.5), col = "blue", lwd = 2)
  }
  # tests that pdf(q(x)) == x
  if(!is.prior.point(prior) & !is.prior.discrete(prior) & !is.prior.spike_and_slab(prior)){
    expect_equal(.25, cdf(prior, quant(prior, 0.25)), tolerance = 1e-4)
    expect_equal(.25, ccdf(prior, quant(prior, 0.75)), tolerance = 1e-4)
  }
  # test mean and sd functions
  if(!skip_moments){
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples),   sd(prior),   tolerance = 1e-2)
  }
  return(invisible())
}
test_weightfunction <- function(prior, skip_moments = FALSE){
  set.seed(1)
  # tests rng and print function (for plot)
  samples   <- rng(prior, 10000)
  densities <- density(prior, individual = TRUE)

  if(!all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))){
    quantiles <- mquant(prior, 0.5)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
  par(mfcol = c(1, ncol(samples)-1))

  for(i in 1:(ncol(samples)-1)){
    hist(samples[,i], main = print(prior, plot = T), breaks = 50, freq = FALSE)
    lines(densities[[i]])
    if(!all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))){
      abline(v = quantiles[i], col = "blue", lwd = 2)
    }
    if(!grepl("fixed", prior$distribution) & !all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))){
      expect_equal(.25, mcdf(prior, mquant(prior, 0.25)[,i])[,i], tolerance = 1e-5)
      expect_equal(.25, mccdf.prior(prior, mquant(prior, 0.75)[,i])[,i], tolerance = 1e-5)
    }
    if(!skip_moments){
      expect_equal(apply(samples, 2, mean), mean(prior), tolerance = 1e-2)
      expect_equal(apply(samples, 2, sd),   sd(prior)  , tolerance = 1e-2)
    }
  }
  return(invisible())
}
test_orthonormal    <- function(prior, skip_moments = FALSE){
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  samples <- samples[abs(samples) < 10]
  hist(samples, main = print(prior, plot = T), breaks = 50, freq = FALSE)
  # tests density function
  lines(prior, individual = TRUE)
  # tests quantile function
  abline(v = mquant(prior, 0.5), col = "blue", lwd = 2)
  # tests that pdf(q(x)) == x
  if(!is.prior.point(prior)){
    expect_equal(.25, mcdf(prior, mquant(prior, 0.25)), tolerance = 1e-5)
    expect_equal(.25, mccdf(prior, mquant(prior, 0.75)), tolerance = 1e-5)
  }
  # test mean and sd functions
  if(!skip_moments){
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples),   sd(prior),   tolerance = 1e-2)
  }
  return(invisible())
}
test_meandif        <- function(prior, skip_moments = FALSE){
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  samples <- samples[abs(samples) < 10]
  hist(samples, main = print(prior, plot = T), breaks = 50, freq = FALSE)
  # tests density function
  lines(prior, individual = TRUE)
  # tests quantile function
  abline(v = mquant(prior, 0.5), col = "blue", lwd = 2)
  # tests that pdf(q(x)) == x
  if(!is.prior.point(prior)){
    expect_equal(.25, mcdf(prior, mquant(prior, 0.25)), tolerance = 1e-5)
    expect_equal(.25, mccdf(prior, mquant(prior, 0.75)), tolerance = 1e-5)
  }
  # test mean and sd functions
  if(!skip_moments){
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples),   sd(prior),   tolerance = 1e-2)
  }
  return(invisible())
}

test_that("Normal prior distribution works", {

  vdiffr::expect_doppelganger("prior-normal-1", function()test_prior(prior("normal", list(0, 1))))
  vdiffr::expect_doppelganger("prior-normal-2", function()test_prior(prior("normal", list(1, 1), list(0, Inf))))
  vdiffr::expect_doppelganger("prior-normal-3", function()test_prior(prior("normal", list(2, 1), list(-Inf, 0))))
  vdiffr::expect_doppelganger("prior-normal-4", function()test_prior(prior("normal", list(0, 1), list(1, 2))))

})

test_that("Log-normal prior distribution works", {

  vdiffr::expect_doppelganger("prior-lognormal-1", function()test_prior(prior("lognormal", list(0, .5))))
  vdiffr::expect_doppelganger("prior-lognormal-2", function()test_prior(prior("lognormal", list(1, .3), list(1, Inf))))

})

test_that("Student-t prior distribution works", {

  vdiffr::expect_doppelganger("prior-t-1", function()test_prior(prior("t", list(0, .5, 5))))
  vdiffr::expect_doppelganger("prior-t-2", function()test_prior(prior("t", list(1, .3, 8), list(1, Inf))))

})

test_that("Cauchy prior distribution works", {

  vdiffr::expect_doppelganger("prior-cauchy-1", function()test_prior(prior("Cauchy", list(0, 1)), skip_moments = TRUE))
  vdiffr::expect_doppelganger("prior-cauchy-2", function()test_prior(prior("Cauchy", list(1, 0.1), list(-Inf, 0)), skip_moments = TRUE))

})

test_that("Gamma prior distribution works", {

  vdiffr::expect_doppelganger("prior-gamma-1", function()test_prior(prior("gamma", list(1, 1))))
  vdiffr::expect_doppelganger("prior-gamma-2", function()test_prior(prior("gamma", list(2, 2), list(1, 3))))
  vdiffr::expect_doppelganger("prior-gamma-3", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    test_prior(prior("gamma", list(shape = 2, "rate"  = 3),   list(0, 3)))
    test_prior(prior("gamma", list(shape = 2, "scale" = 1/3), list(0, 3)))
  })
})

test_that("Inverse-gamma prior distribution works", {

  vdiffr::expect_doppelganger("prior-invgamma-1", function()test_prior(prior("invgamma", list(1.5, 0.15)), skip_moments = TRUE))
  vdiffr::expect_doppelganger("prior-invgamma-2", function()test_prior(prior("invgamma", list(3, 2), list(1, 3))))

})

test_that("Exponential prior distribution works", {

  vdiffr::expect_doppelganger("prior-exp-1", function()test_prior(prior("exp", list(1.5))))
  vdiffr::expect_doppelganger("prior-exp-2", function()test_prior(prior("exp", list(2), list(1, 3))))
  vdiffr::expect_doppelganger("prior-exp-3", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    test_prior(prior("exp", list("rate"  = 3),   list(1, 3)))
    test_prior(prior("exp", list("scale" = 1/3), list(1, 3)))
  })

})

test_that("Beta prior distribution works", {

  vdiffr::expect_doppelganger("prior-beta-1", function()test_prior(prior("beta", list(1, 1))))
  vdiffr::expect_doppelganger("prior-beta-2", function()test_prior(prior("beta", list(25, 15), list(0.5, 1))))

})

test_that("Beta prior distribution works", {

  vdiffr::expect_doppelganger("prior-bernoulli-1", function()test_prior(prior("bernoulli", list(.50))))
  vdiffr::expect_doppelganger("prior-bernoulli-2", function()test_prior(prior("bernoulli", list(.66))))

})

test_that("Uniform prior distribution works", {

  vdiffr::expect_doppelganger("prior-uniform-1", function()test_prior(prior("uniform", list(0, 1))))
  vdiffr::expect_doppelganger("prior-uniform-2", function()test_prior(prior("uniform", list(1, 5))))

})

test_that("Spike prior distribution works", {

  vdiffr::expect_doppelganger("prior-point-1", function()test_prior(prior("point", list(1))))

})

test_that("spike and slab prior distribution works", {

  vdiffr::expect_doppelganger("prior-spike-and-slab-1",   function()test_prior(
    prior_spike_and_slab(
    prior("gamma", list(2, 2), list(0, Inf)),
    prior_inclusion = prior("beta", list(2, 1)))
  ))

})

test_that("PET & PEESE prior distribution works", {

  vdiffr::expect_doppelganger("prior-PET-1",   function()test_prior(prior_PET("normal", list(0, 1))))
  vdiffr::expect_doppelganger("prior-PEESE-1", function()test_prior(prior_PEESE("gamma", list(1, 1))))

})

test_that("One-sided weigthfunction prior distribution works", {

  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-1", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05), c(1, 1)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-2", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-3", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05), c(5, .5)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-4", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05, 0.10), c(5, 5, 1)))))

  # non-monotonic one-sided requires sampling
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-5", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05, 0.60), c(1, 1), c(1, 1))), skip_moments = TRUE))
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided-6", function()test_weightfunction(prior_weightfunction("one.sided", list(c(.05, 0.10, 0.60), c(1, 1, 1), c(1, 1))), skip_moments = TRUE))
})

test_that("Two-sided weigthfunction prior distribution works", {

  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided-1", function()test_weightfunction(prior_weightfunction("two.sided", list(c(.05), c(1, 1)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided-2", function()test_weightfunction(prior_weightfunction("two.sided", list(c(.05, 0.10), c(1, 1, 1)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided-3", function()test_weightfunction(prior_weightfunction("two.sided", list(c(.05), c(5, .5)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided-4", function()test_weightfunction(prior_weightfunction("two.sided", list(c(.05, 0.10), c(5, 5, 1)))))

})

test_that("One-sided.fixed weigthfunction prior distribution works", {

  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided.fixed-1", function()test_weightfunction(prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-one.sided.fixed-2", function()test_weightfunction(prior_weightfunction("one.sided.fixed", list(c(.05, 0.10), c(1, .2, .5)))))

})

test_that("Two-sided.fixed weigthfunction prior distribution works", {

  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided.fixed-1", function()test_weightfunction(prior_weightfunction("two.sided.fixed", list(c(.05), c(1, .5)))))
  vdiffr::expect_doppelganger("prior-weigthfunction-two.sided.fixed-2", function()test_weightfunction(prior_weightfunction("two.sided.fixed", list(c(.05, 0.10), c(1, .2, .5)))))

})

test_that("Vector prior distribution works", {
  # TODO
})

test_that("Orthonormal prior distribution works", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate Cauchy sampling does not exactly match across OSes

  p1   <- prior_factor("mnormal", list(mean = 0, sd = 0.5), contrast = "orthonormal")
  p2   <- prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "orthonormal")
  p3   <- prior_factor("point", list(0), contrast = "orthonormal")
  p1.5 <- p1.3 <- p1.2 <- p1
  p2.9 <- p2
  p3.3 <- p3.5 <- p3
  p1.2$parameters$K <- 2
  p1.3$parameters$K <- 3
  p1.5$parameters$K <- 5
  p2.9$parameters$K <- 9
  p3.3$parameters$K <- 3
  p3.5$parameters$K <- 5

  vdiffr::expect_doppelganger("prior-orthonormal-1-2", function()test_orthonormal(p1.2))
  vdiffr::expect_doppelganger("prior-orthonormal-1-3", function()test_orthonormal(p1.3))
  vdiffr::expect_doppelganger("prior-orthonormal-1-5", function()test_orthonormal(p1.5))
  vdiffr::expect_doppelganger("prior-orthonormal-2-9", function()test_orthonormal(p2.9, skip_moments = TRUE))
  vdiffr::expect_doppelganger("prior-orthonormal-3-3", function()test_orthonormal(p3.3))
  vdiffr::expect_doppelganger("prior-orthonormal-3-5", function()test_orthonormal(p3.5))

})

test_that("Treatment prior distribution works", {

  vdiffr::expect_doppelganger("prior-treatment-1", function()test_prior(prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")))
  vdiffr::expect_doppelganger("prior-treatment-2", function()test_prior(prior_factor("beta", list(alpha = 2, beta = 3), contrast = "treatment")))
  vdiffr::expect_doppelganger("prior-treatment-3", function()test_prior(prior_factor("spike", list(location = 1), contrast = "treatment")))

})

test_that("Independent prior distribution works", {

  vdiffr::expect_doppelganger("prior-independent-1", function()test_prior(prior_factor("gamma", list(shape = 2, rate = 3), contrast = "independent")))
  vdiffr::expect_doppelganger("prior-independent-2", function()test_prior(prior_factor("uniform", list(a = -0.5, b = 1), contrast = "independent")))
  vdiffr::expect_doppelganger("prior-independent-3", function()test_prior(prior_factor("spike", list(location = 1), contrast = "independent")))

})

test_that("Meandif prior distribution works", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes

  p1   <- prior_factor("mnormal", list(mean = 0, sd = 0.25), contrast = "meandif")
  p2   <- prior_factor("point", list(0), contrast = "orthonormal")
  p1.5 <- p1.3 <- p1.2 <- p1
  p2.5 <- p2.3 <- p2
  p1.2$parameters$K <- 2
  p1.3$parameters$K <- 3
  p1.5$parameters$K <- 5
  p2.3$parameters$K <- 3
  p2.5$parameters$K <- 5

  vdiffr::expect_doppelganger("prior-meandif-1-2", function()test_meandif(p1.2))
  vdiffr::expect_doppelganger("prior-meandif-1-3", function()test_meandif(p1.3))
  vdiffr::expect_doppelganger("prior-meandif-1-5", function()test_meandif(p1.5))
  vdiffr::expect_doppelganger("prior-meandif-2-3", function()test_meandif(p2.3))
  vdiffr::expect_doppelganger("prior-meandif-2-5", function()test_meandif(p2.5))

})
