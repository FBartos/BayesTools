skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Distributions of linear functionals of R-squared budget priors
# ============================================================================ #
#
# PURPOSE:
#   The density, distribution function, quantile function, random generation
#   and variance of a linear functional X of the coefficients of one term of an
#   R-squared budget prior (laws "exact", "normal" and "laplace"), and the
#   density and distribution function of the unit signal v = W * phi, from
#   R/priors-R2-density.R. The densities are deterministic nested
#   stats::integrate() quadratures; they are checked against closed forms,
#   independent numerical formulas, 18-digit mpmath references (literals below,
#   computed by nested tanh-sinh quadrature of the definitions), and
#   simulation written independently of the package generator.
#
# DEPENDENCIES:
#   - None (pure R, no JAGS)
#
# TAGS: @unit, @priors, @R2
# ============================================================================ #

r2d_spec     <- BayesTools:::.bt_R2_functional_spec
r2d_pdf      <- BayesTools:::.bt_R2_functional_pdf
r2d_cdf      <- BayesTools:::.bt_R2_functional_cdf
r2d_quantile <- BayesTools:::.bt_R2_functional_quantile
r2d_rng      <- BayesTools:::.bt_R2_functional_rng
r2d_variance <- BayesTools:::.bt_R2_functional_variance
r2d_signal_spec <- BayesTools:::.bt_R2_signal_spec
r2d_signal_pdf  <- BayesTools:::.bt_R2_signal_pdf
r2d_signal_cdf  <- BayesTools:::.bt_R2_signal_cdf

# One functional of every kind: the exact law with one and with three columns,
# one and several units, and both kernels. 'a != A' keeps the allocation share
# a separate mixing level (a = A absorbs it exactly).
r2d_cases <- list(
  exact_d1_n1 = list(law = "exact", a = 0.7, b = 2, n_units = 1, d = 1, scale = 1),
  exact_d3_n1 = list(law = "exact", a = 1.5, b = 2.5, n_units = 1, d = 3, scale = 1),
  exact_d1_n3 = list(law = "exact", a = 1.5, b = 2.5, alpha = 0.8, A = 2.4, n_units = 3, d = 1, scale = 1),
  exact_d3_n3 = list(law = "exact", a = 1.5, b = 2.5, alpha = 0.8, A = 2.4, n_units = 3, d = 3, scale = 1),
  normal_n1   = list(law = "normal", a = 1.5, b = 2.5, n_units = 1, d = 1, scale = 1),
  normal_n4   = list(law = "normal", a = 1.5, b = 2.5, alpha = 0.8, A = 2.4, n_units = 4, d = 1, scale = 0.8),
  laplace_n1  = list(law = "laplace", a = 0.6, b = 3, n_units = 1, d = 1, scale = 1),
  laplace_n4  = list(law = "laplace", a = 0.6, b = 3, alpha = 0.4, A = 3, n_units = 4, d = 1, scale = 1)
)
r2d_cheap  <- c("exact_d1_n1", "exact_d3_n1", "exact_d1_n3", "normal_n1", "laplace_n1")
r2d_nested <- c("exact_d3_n3", "normal_n4", "laplace_n4")
r2d_make   <- function(name) do.call(r2d_spec, r2d_cases[[name]])

# ---------------------------------------------------------------------------
# Independent references

# Closed form of Y = sqrt(W) * S with a fair sign S (the exact law, one column).
r2d_f_exact <- function(y, a, b) abs(y)^(2 * a - 1) / (beta(a, b) * (1 + y^2)^(a + b))
r2d_F_exact <- function(y, a, b) 1 / 2 + sign(y) * stats::pbeta(y^2 / (1 + y^2), a, b) / 2

# BetaPrime(a, b) density and distribution function.
r2d_betaprime_pdf <- function(q, a, b) q^(a - 1) * (1 + q)^(-a - b) / beta(a, b)

# Gauss hypergeometric function at z <= 0 by Euler's integral (C > B > 0).
r2d_hyp2f1 <- function(a1, b1, c1, z){
  stats::integrate(function(t) t^(b1 - 1) * (1 - t)^(c1 - b1 - 1) * (1 - z * t)^(-a1), 0, 1,
                   rel.tol = 1e-11, subdivisions = 500L)$value / beta(b1, c1 - b1)
}

# Density of v = W * phi, W ~ BetaPrime(a, b), phi ~ Beta(alpha, A - alpha).
r2d_signal_pdf_reference <- function(q, a, b, alpha, A){
  vapply(q, function(value){
    value^(-b - 1) * beta(alpha + b, A - alpha) * r2d_hyp2f1(a + b, alpha + b, A + b, -1 / value) /
      (beta(a, b) * beta(alpha, A - alpha))
  }, numeric(1))
}

# P(v <= q) = E[F_W(q / phi)] with F_W(w) = pbeta(w / (1 + w), a, b).
r2d_signal_cdf_reference <- function(q, a, b, alpha, A){
  vapply(q, function(value){
    stats::integrate(function(phi){
      stats::dbeta(phi, alpha, A - alpha) * stats::pbeta(value / (phi + value), a, b)
    }, 0, 1, rel.tol = 1e-11, subdivisions = 500L)$value
  }, numeric(1))
}

# Simulation from the definitions with other generators than the package's:
# gamma ratios for W, the first coordinate of a full Dirichlet vector for the
# share, a normalized normal vector for the uniform direction, a difference of
# two exponentials for the unit-variance Laplace law.
r2d_simulate <- function(case, n){

  weight <- stats::rgamma(n, case$a) / stats::rgamma(n, case$b)
  if(case$n_units > 1){
    rest  <- Reduce(`+`, lapply(seq_len(case$n_units - 1), function(i){
      stats::rgamma(n, (case$A - case$alpha) / (case$n_units - 1))
    }))
    first  <- stats::rgamma(n, case$alpha)
    weight <- weight * first / (first + rest)
  }
  shape <- switch(case$law,
    exact   = {
      z <- matrix(stats::rnorm(n * case$d), n, case$d)
      z[, 1] / sqrt(rowSums(z^2))
    },
    normal  = stats::rnorm(n),
    laplace = (stats::rexp(n) - stats::rexp(n)) / sqrt(2))
  case$scale * sqrt(weight) * shape
}

# Gauss-Legendre nodes and weights on (-1, 1) (Golub-Welsch).
r2d_gauss_legendre <- function(n){

  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1)] <- b
  jacobi[cbind(i + 1, i)] <- b
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(x = decomposition$values, w = 2 * decomposition$vectors[1, ]^2)
}

# The mass of a density between lo and hi by Gauss-Legendre on the log scale.
r2d_mass <- function(spec, lo, hi, nodes = 10){

  rule <- r2d_gauss_legendre(nodes)
  t    <- (log(lo) + log(hi)) / 2 + rule$x * (log(hi) - log(lo)) / 2
  sum(rule$w * r2d_pdf(spec, exp(t)) * exp(t)) * (log(hi) - log(lo)) / 2
}

# ---------------------------------------------------------------------------

test_that("the exact law of one unit and one column equals its closed form", {

  x      <- c(0.01, 0.3, 1, 2.5, 8)
  shapes <- list(c(0.3, 2), c(1 / 2, 1.5), c(1.7, 3.2), c(3, 0.8))
  problems <- expectation_problems({
    for(shape in shapes){
      a    <- shape[1]
      b    <- shape[2]
      info <- sprintf("a = %g, b = %g", a, b)
      spec <- r2d_spec("exact", a, b)
      expect_equal_each(r2d_pdf(spec, c(x, -x)), r2d_f_exact(c(x, -x), a, b), tolerance = 1e-7, info = info)
      expect_equal_each(r2d_cdf(spec, c(-x, x)), r2d_F_exact(c(-x, x), a, b), tolerance = 1e-7, info = info)
      expect_equal_each(r2d_cdf(spec, x, lower.tail = FALSE), 1 - r2d_F_exact(x, a, b), tolerance = 1e-7, info = info)
      # a deterministic scale divides the argument and the density
      scaled <- r2d_spec("exact", a, b, scale = 2.5)
      expect_equal_each(r2d_pdf(scaled, x), r2d_f_exact(x / 2.5, a, b) / 2.5, tolerance = 1e-7, info = info)
      expect_equal_each(r2d_cdf(scaled, -x), r2d_F_exact(-x / 2.5, a, b), tolerance = 1e-7, info = info)
      # the density integrates to 1 (the pieces avoid the singularity at 0 and the heavy tail)
      mass <- 2 * (stats::integrate(function(y) r2d_pdf(spec, y), 0, 1, rel.tol = 1e-10)$value +
                     stats::integrate(function(y) r2d_pdf(spec, y), 1, Inf, rel.tol = 1e-10)$value)
      expect_equal(mass, 1, tolerance = 1e-7, info = info)
    }
  })
  expect_identical(problems, character())
})

test_that("the unit signal with a = A is BetaPrime(alpha, b) and with a != A follows the hypergeometric form", {

  q <- c(0.01, 0.2, 1, 5, 40)

  # BetaPrime(A, b) * Beta(alpha, A - alpha) = BetaPrime(alpha, b)
  absorbed <- list(list(alpha = 0.4, A = 2, b = 1.5, n_units = 3), list(alpha = 1.2, A = 3, b = 0.7, n_units = 2),
                   list(alpha = 0.5, A = 1, b = 2, n_units = 4))
  problems <- expectation_problems({
    for(case in absorbed){
      spec <- r2d_signal_spec(a = case$A, b = case$b, alpha = case$alpha, A = case$A, n_units = case$n_units)
      info <- paste(names(case), unlist(case), collapse = ", ")
      expect_equal_each(r2d_signal_pdf(spec, q), r2d_betaprime_pdf(q, case$alpha, case$b),
                        tolerance = 1e-7, info = info)
      expect_equal_each(r2d_signal_cdf(spec, q), stats::pbeta(q / (1 + q), case$alpha, case$b),
                        tolerance = 1e-7, info = info)
      expect_equal_each(r2d_signal_cdf(spec, q, lower.tail = FALSE),
                        stats::pbeta(q / (1 + q), case$alpha, case$b, lower.tail = FALSE),
                        tolerance = 1e-7, info = info)
    }
  })
  expect_identical(problems, character())

  # general a: f(q) = q^(-b-1) B(alpha+b, A-alpha) 2F1(a+b, alpha+b; A+b; -1/q) / {B(a,b) B(alpha, A-alpha)}
  general <- list(list(a = 1.2, b = 2.5, alpha = 0.7, A = 2.1, n_units = 3),
                  list(a = 0.4, b = 1.5, alpha = 0.3, A = 4, n_units = 2),
                  list(a = 3, b = 0.8, alpha = 1, A = 1.5, n_units = 5))
  problems <- expectation_problems({
    for(case in general){
      spec <- r2d_signal_spec(case$a, case$b, case$alpha, case$A, case$n_units)
      info <- paste(names(case), unlist(case), collapse = ", ")
      expect_equal_each(r2d_signal_pdf(spec, q), r2d_signal_pdf_reference(q, case$a, case$b, case$alpha, case$A),
                        tolerance = 1e-6, info = info)
      reference <- r2d_signal_cdf_reference(q, case$a, case$b, case$alpha, case$A)
      expect_equal_each(r2d_signal_cdf(spec, q), reference, tolerance = 1e-7, info = info)
      expect_equal_each(r2d_signal_cdf(spec, q, lower.tail = FALSE), 1 - reference, tolerance = 1e-7, info = info)
    }
  })
  expect_identical(problems, character())

  # a unit signal is a probability distribution on [0, Inf)
  spec <- r2d_signal_spec(1.2, 2.5, 0.7, 2.1, 3)
  mass <- stats::integrate(function(v) r2d_signal_pdf(spec, v), 0, 1, rel.tol = 1e-9)$value +
    stats::integrate(function(v) r2d_signal_pdf(spec, v), 1, Inf, rel.tol = 1e-9)$value
  expect_equal(mass, 1, tolerance = 1e-7)
  expect_identical(r2d_signal_cdf(spec, c(-1, 0, Inf)), c(0, 0, 1))
  expect_identical(r2d_signal_cdf(spec, c(-1, 0, Inf), lower.tail = FALSE), c(1, 1, 0))
  expect_identical(r2d_signal_pdf(spec, c(-1, Inf)), c(0, 0))
})

test_that("the CDF matches the empirical CDF of the package generator, which matches an independent simulation", {

  n      <- 2e5
  probs  <- c(0.03, 0.2, 0.45, 0.7, 0.93)
  problems <- expectation_problems({
    for(i in seq_along(r2d_cases)){
      case <- r2d_cases[[i]]
      spec <- r2d_make(names(r2d_cases)[i])
      set.seed(100 + i)
      draws       <- r2d_rng(spec, n)
      independent <- r2d_simulate(case, n)
      # evaluation points from the independent simulation; conditional on them,
      # each empirical CDF value has the binomial standard error of its probability
      q    <- unname(stats::quantile(independent, probs))
      cdf  <- r2d_cdf(spec, q)
      se   <- sqrt(cdf * (1 - cdf) / n)
      info <- names(r2d_cases)[i]
      expect_true(all(abs(vapply(q, function(v) mean(draws <= v), numeric(1)) - cdf) < 4 * se + 1e-6),
                  info = paste(info, "package generator against the CDF"))
      # the two simulations differ in the draws and in the evaluation points
      expect_true(all(abs(vapply(q, function(v) mean(draws <= v), numeric(1)) - probs) < 4 * sqrt(2) * se + 1e-6),
                  info = paste(info, "package generator against the independent simulation"))
      expect_true(all(abs(cdf - probs) < 4 * se + 1e-6), info = paste(info, "CDF against the independent simulation"))
    }
  })
  expect_identical(problems, character())
})

test_that("pdf and CDF agree with mpmath references of the definitions", {

  # Values computed with mpmath (18 digits) by nested tanh-sinh quadrature of the
  # defining integrals (the base law at x / sqrt(s) mixed over every Beta share),
  # with the closed form of the exact base and Tricomi's U function for the
  # normal base: columns x, density, upper tail P(X > x).
  references <- list(
    normal_n1 = list(spec = r2d_spec("normal", 1.5, 2.5), x = c(0.05, 0.7, 3),
                     pdf = c(0.663829202195612, 0.266315582945131, 0.00880766153198311),
                     tail = c(0.466388633160438, 0.171159012198181, 0.00798521244400873)),
    laplace_n1 = list(spec = r2d_spec("laplace", 0.6, 3), x = c(0.1, 1, 5),
                      pdf = c(1.19946319228238, 0.0575613448097679, 0.000224823845834764),
                      tail = c(0.29483027104291, 0.0308811826921386, 0.00027022714676033)),
    normal_n4 = list(spec = r2d_spec("normal", 1.5, 2.5, 0.8, 2.4, 4), x = c(0.05, 0.7, 3),
                     pdf = c(1.47537325885351, 0.160161647449485, 0.0019205824276094),
                     tail = c(0.412345864081657, 0.070092421752388, 0.0015492096467166)),
    laplace_n4 = list(spec = r2d_spec("laplace", 0.6, 3, 0.4, 3, 4), x = c(0.1, 1, 5),
                      pdf = c(0.937972166113609, 0.00841529770979923, 9.4903290736001e-6),
                      tail = c(0.114679704687328, 0.00330404025408136, 1.01145743468301e-5)),
    exact_d3_n3 = list(spec = r2d_spec("exact", 1.5, 2.5, 0.8, 2.4, 3, d = 3), x = c(0.05, 0.7, 3),
                       pdf = c(1.86531466475651, 0.0925889940385075, 0.00013021994799116),
                       tail = c(0.389221290779206, 0.0240708564199508, 8.30126778463854e-5)),
    exact_d2_n3 = list(spec = r2d_spec("exact", 0.4, 3, 0.3, 0.9, 3, d = 2), x = c(0.05, 0.7, 3),
                       pdf = c(2.02012729221705, 0.0261681394901664, 1.33145095277142e-5),
                       tail = c(0.213381616977946, 0.0058014145473967, 7.10979886918043e-6))
  )
  problems <- expectation_problems({
    for(name in names(references)){
      reference <- references[[name]]
      expect_equal_each(r2d_pdf(reference$spec, reference$x), reference$pdf, tolerance = 1e-7, info = name)
      expect_equal_each(r2d_cdf(reference$spec, reference$x, lower.tail = FALSE), reference$tail,
                        tolerance = 1e-7, info = name)
    }
  })
  expect_identical(problems, character())
})

test_that("pdf is symmetric, equals the derivative of the CDF, integrates to the CDF, and the CDF is monotone", {

  h <- 1e-3
  problems <- expectation_problems({
    for(name in names(r2d_cases)){
      spec   <- r2d_make(name)
      nested <- name %in% r2d_nested
      # a nested (two-level) evaluation costs about 0.1 s: a single point, no grid
      x <- if(nested) 0.5 else c(0.2, 1.1)
      expect_identical(r2d_pdf(spec, -x), r2d_pdf(spec, x), info = name)
      expect_equal_each(r2d_cdf(spec, -x), 1 - r2d_cdf(spec, x), tolerance = 1e-12, info = name)
      expect_equal_each((r2d_cdf(spec, x + h) - r2d_cdf(spec, x - h)) / (2 * h), r2d_pdf(spec, x),
                        tolerance = 2e-5, info = name)
      if(!nested){
        expect_true(all(diff(r2d_cdf(spec, c(-5, -1, -0.1, 0.1, 1, 5))) > 0), info = name)
      }
      expect_identical(r2d_cdf(spec, c(-Inf, 0, Inf)), c(0, 1 / 2, 1), info = name)
      expect_identical(r2d_pdf(spec, c(-Inf, Inf)), c(0, 0), info = name)
    }
  })
  expect_identical(problems, character())

  # two-level cases: the density integrates to the difference of the upper tails
  # (a 10-node Gauss-Legendre rule on the log scale, accurate to about 1e-9 here,
  # at about 0.1 s per node and case)
  problems <- expectation_problems({
    for(name in r2d_nested){
      spec <- r2d_make(name)
      tails <- r2d_cdf(spec, c(0.3, 3), lower.tail = FALSE)
      expect_equal(r2d_mass(spec, 0.3, 3), tails[1] - tails[2], tolerance = 1e-7, info = name)
    }
  })
  expect_identical(problems, character())
})

test_that("the density integrates to 1 and the variance matches the second moment", {

  # one-level cases (and the absorbed share, a = A) are cheap enough to integrate
  cases <- list(
    exact_d3_n1   = r2d_make("exact_d3_n1"),
    normal_n1     = r2d_make("normal_n1"),
    laplace_n1    = r2d_make("laplace_n1"),
    exact_absorbed  = r2d_spec("exact", 2.4, 3.5, 0.8, 2.4, 3, d = 3, scale = 1.5),
    laplace_absorbed = r2d_spec("laplace", 2.4, 3.5, 0.8, 2.4, 3)
  )
  problems <- expectation_problems({
    for(name in names(cases)){
      spec <- cases[[name]]
      integral <- function(f){
        2 * (stats::integrate(function(x) f(x) * r2d_pdf(spec, x), 0, 1, rel.tol = 1e-9)$value +
               stats::integrate(function(x) f(x) * r2d_pdf(spec, x), 1, Inf, rel.tol = 1e-9)$value)
      }
      expect_equal(integral(function(x) 1), 1, tolerance = 1e-7, info = name)
      expect_equal(integral(function(x) x^2), r2d_variance(spec), tolerance = 1e-7, info = name)
    }
  })
  expect_identical(problems, character())

  # E[W] = a / (b - 1), E[phi] = alpha / A, E[U1^2] = 1 / d, scale^2
  expect_equal(r2d_variance(r2d_spec("exact", 1.5, 2.5, 0.8, 2.4, 3, d = 3, scale = 2)),
               4 * (1.5 / 1.5) * (0.8 / 2.4) / 3, tolerance = 1e-14)
  expect_equal(r2d_variance(r2d_spec("normal", 1.5, 2.5, 0.8, 2.4, 3, d = 3)), (1.5 / 1.5) * (0.8 / 2.4),
               tolerance = 1e-14)
  expect_equal(r2d_variance(r2d_spec("laplace", 3, 6)), 3 / 5, tolerance = 1e-14)
  expect_identical(r2d_variance(r2d_spec("exact", 1.5, 1)), Inf)
  expect_identical(r2d_variance(r2d_spec("normal", 1.5, 0.7, 0.8, 2.4, 3)), Inf)
  expect_identical(r2d_variance(r2d_spec("laplace", 0.2, 1)), Inf)
})

test_that("the quantile function inverts the CDF", {

  probs <- c(0.01, 0.2, 0.45, 0.6, 0.97, 0.9995)
  problems <- expectation_problems({
    for(name in r2d_cheap){
      spec <- r2d_make(name)
      quantiles <- r2d_quantile(spec, probs)
      expect_equal_each(r2d_cdf(spec, quantiles), probs, tolerance = 1e-8, info = name)
      expect_equal_each(r2d_quantile(spec, 1 - probs), -quantiles, tolerance = 1e-9, info = name)
      x <- c(-5, -1, -0.3, 0.3, 1, 5)
      expect_equal_each(r2d_quantile(spec, r2d_cdf(spec, x)), x, tolerance = 1e-6, info = name)
    }
    # a nested (two-level) case costs about 0.1 s per tail evaluation: one probability
    spec <- r2d_make("exact_d3_n3")
    expect_equal(r2d_cdf(spec, r2d_quantile(spec, 0.9)), 0.9, tolerance = 1e-8)
  })
  expect_identical(problems, character())

  spec <- r2d_make("exact_d1_n1")
  expect_identical(r2d_quantile(spec, c(0, 0.5, 1)), c(-Inf, 0, Inf))
  expect_identical(r2d_quantile(spec, NA_real_), NA_real_)
  # a deterministic scale multiplies the quantiles
  expect_equal(r2d_quantile(r2d_spec("exact", 0.7, 2, scale = 3), 0.8), 3 * r2d_quantile(spec, 0.8),
               tolerance = 1e-10)
  # the exact quantiles of the closed form
  expect_equal(r2d_quantile(spec, c(0.1, 0.8)),
               c(-1, 1) * sqrt(stats::qbeta(c(0.8, 0.6), 0.7, 2) / (1 - stats::qbeta(c(0.8, 0.6), 0.7, 2))),
               tolerance = 1e-8)
})

test_that("the density at 0 follows from the shapes", {

  eps <- 1e-9
  inv_half <- function(p, q) beta(p - 1 / 2, q) / beta(p, q)                   # E[S^(-1/2)], S ~ Beta(p, q)
  w_half   <- function(a, b) beta(a - 1 / 2, b + 1 / 2) / beta(a, b)           # E[W^(-1/2)], W ~ BetaPrime(a, b)

  # singular, critical and vanishing exact densities
  expect_identical(r2d_pdf(r2d_spec("exact", 0.3, 2), 0), Inf)
  expect_equal(r2d_pdf(r2d_spec("exact", 1 / 2, 2), 0), 1 / beta(1 / 2, 2), tolerance = 1e-12)
  expect_equal(r2d_pdf(r2d_spec("exact", 1 / 2, 2), eps), 1 / beta(1 / 2, 2), tolerance = 1e-7)
  expect_identical(r2d_pdf(r2d_spec("exact", 1.5, 2), 0), 0)
  expect_identical(r2d_pdf(r2d_spec("exact", 1.5, 2, scale = 2), c(0, 1e-3))[1], 0)
  # one column of a multi-column term: the projection U1 sits on the threshold
  # exactly, so the density is finite unless another factor ties with it
  expect_equal(r2d_pdf(r2d_spec("exact", 1.5, 2.5, d = 3), 0), 1 / beta(1 / 2, 1) * w_half(1.5, 2.5),
               tolerance = 1e-12)
  expect_equal(r2d_pdf(r2d_spec("exact", 1.5, 2.5, d = 3), eps), r2d_pdf(r2d_spec("exact", 1.5, 2.5, d = 3), 0),
               tolerance = 1e-6)
  expect_identical(r2d_pdf(r2d_spec("exact", 1 / 2, 2.5, d = 3), 0), Inf)
  expect_identical(r2d_pdf(r2d_spec("exact", 1.5, 2.5, 0.3, 2.4, 3, d = 1), 0), Inf)
  # one share at the threshold
  threshold <- r2d_spec("exact", 1.5, 2.5, 1 / 2, 2.4, 3, d = 1)
  expect_equal(r2d_pdf(threshold, 0), 1 / beta(1 / 2, 1.9) * w_half(1.5, 2.5), tolerance = 1e-12)
  expect_equal(r2d_pdf(threshold, eps), r2d_pdf(threshold, 0), tolerance = 1e-6)
  # kernels: k(0) * E[W^(-1/2)] * E[phi^(-1/2)], infinite when a share or W has shape <= 1/2
  expect_equal(r2d_pdf(r2d_spec("normal", 1.5, 2.5), 0), w_half(1.5, 2.5) / sqrt(2 * pi), tolerance = 1e-12)
  expect_equal(r2d_pdf(r2d_spec("normal", 1.5, 2.5), eps), w_half(1.5, 2.5) / sqrt(2 * pi), tolerance = 1e-7)
  expect_equal(r2d_pdf(r2d_spec("laplace", 1.5, 2.5, 1.2, 2.4, 4), 0),
               w_half(1.5, 2.5) * inv_half(1.2, 1.2) * sqrt(2) / 2, tolerance = 1e-12)
  expect_equal(r2d_pdf(r2d_spec("laplace", 1.5, 2.5, 1.2, 2.4, 4), eps),
               w_half(1.5, 2.5) * inv_half(1.2, 1.2) * sqrt(2) / 2, tolerance = 1e-6)
  expect_identical(r2d_pdf(r2d_spec("normal", 1 / 2, 2.5), 0), Inf)
  expect_identical(r2d_pdf(r2d_spec("laplace", 1.5, 2.5, 1 / 2, 2.4, 4), 0), Inf)

  # the unit signal: the density of BetaPrime(a, b) at 0 is b for a = 1
  expect_identical(r2d_signal_pdf(r2d_signal_spec(0.7, 2), 0), Inf)
  expect_identical(r2d_signal_pdf(r2d_signal_spec(1.4, 2), 0), 0)
  expect_equal(r2d_signal_pdf(r2d_signal_spec(1, 2), 0), 2, tolerance = 1e-12)
  expect_equal(r2d_signal_pdf(r2d_signal_spec(1, 2), 1e-10), 2, tolerance = 1e-8)
  # with a share: f(0) = b (A - 1) / (alpha - 1) for a = 1 and alpha > 1
  with_share <- r2d_signal_spec(1, 2, 3, 5, 3)
  expect_equal(r2d_signal_pdf(with_share, 0), 2 * 4 / 2, tolerance = 1e-12)
  expect_equal(r2d_signal_pdf(with_share, 1e-10), 4, tolerance = 1e-7)
  expect_identical(r2d_signal_pdf(r2d_signal_spec(1, 2, 1, 3, 3), 0), Inf)
  expect_identical(r2d_signal_pdf(r2d_signal_spec(1.5, 2, 3, 5, 3), 0), 0)
})

test_that("absorbing the share at a = A is exact", {

  # BetaPrime(A, b) * Beta(alpha, A - alpha) = BetaPrime(alpha, b): the pdf of the
  # absorbed (one-level) path is the limit of the two-level path as a -> A
  x <- c(0.05, 1.5)
  problems <- expectation_problems({
    for(law in c("exact", "normal", "laplace")){
      absorbed <- r2d_spec(law, 2.4, 2.5, 0.8, 2.4, 3, d = 3)
      nearby   <- r2d_spec(law, 2.4 * (1 + 1e-7), 2.5, 0.8, 2.4, 3, d = 3)
      expect_equal_each(r2d_pdf(absorbed, x), r2d_pdf(nearby, x), tolerance = 1e-5, info = law)
      expect_equal_each(r2d_cdf(absorbed, x, lower.tail = FALSE), r2d_cdf(nearby, x, lower.tail = FALSE),
                        tolerance = 1e-5, info = law)
    }
  })
  expect_identical(problems, character())
})

test_that("random generation follows the definition, including scale and the Beta projection", {

  set.seed(20261111)
  # deterministic structure: one unit, fair sign (exact, d = 1)
  draws <- r2d_rng(r2d_spec("exact", 2, 3, scale = 2), 1e4)
  expect_length(draws, 1e4)
  expect_true(all(is.finite(draws)))
  expect_true(abs(mean(draws > 0) - 1 / 2) < 4 * sqrt(1 / 4 / 1e4))
  # X^2 = scale^2 W phi Z^2: the mean of the draws of X^2 is the variance, with the
  # standard error sqrt((E[W^2] E[phi^2] E[Z^4] / (E[W] E[phi] E[Z^2])^2 - 1) / n)
  spec <- r2d_spec("normal", 2, 7, 0.5, 2, 3, scale = 1.5)
  n    <- 4e5
  draws <- r2d_rng(spec, n)
  relative_variance <- (2 * 3 / (6 * 5)) / (2 / 6)^2 * ((0.5 * 1.5 / (2 * 3)) / (0.5 / 2)^2) * 3 - 1
  expect_true(abs(mean(draws^2) / r2d_variance(spec) - 1) < 4 * sqrt(relative_variance / n))
  expect_length(r2d_rng(spec, 0), 0)
})

test_that("invalid specifications and arguments are refused with R2 conditions", {

  bad <- list(
    quote(r2d_spec("gauss", 1, 1)),
    quote(r2d_spec(c("exact", "normal"), 1, 1)),
    quote(r2d_spec(1, 1, 1)),
    quote(r2d_spec("exact", 0, 1)),
    quote(r2d_spec("exact", -1, 1)),
    quote(r2d_spec("exact", NA_real_, 1)),
    quote(r2d_spec("exact", Inf, 1)),
    quote(r2d_spec("exact", c(1, 2), 1)),
    quote(r2d_spec("exact", "1", 1)),
    quote(r2d_spec("exact", 1, 0)),
    quote(r2d_spec("exact", 1, -2)),
    quote(r2d_spec("exact", 1, NaN)),
    quote(r2d_spec("exact", 1, 1, scale = 0)),
    quote(r2d_spec("exact", 1, 1, scale = Inf)),
    quote(r2d_spec("exact", 1, 1, n_units = 0)),
    quote(r2d_spec("exact", 1, 1, n_units = 1.5)),
    quote(r2d_spec("exact", 1, 1, n_units = NA_integer_)),
    quote(r2d_spec("exact", 1, 1, d = 0)),
    quote(r2d_spec("exact", 1, 1, d = 2.5)),
    quote(r2d_spec("exact", 1, 1, alpha = 1, A = 2, n_units = 3, d = c(1, 2))),
    quote(r2d_spec("exact", 1, 1, alpha = 2, A = 2, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, alpha = 3, A = 2, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, alpha = 0, A = 2, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, alpha = -1, A = 2, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, alpha = 1, n_units = 3)),
    quote(r2d_spec("exact", 1, 1, alpha = 1, A = Inf, n_units = 3)),
    quote(r2d_quantile(r2d_spec("exact", 1, 1), 1.5)),
    quote(r2d_quantile(r2d_spec("exact", 1, 1), -0.1))
  )
  refused <- vapply(bad, function(call){
    inherits(tryCatch(eval(call), error = identity), "BayesTools_R2_invalid")
  }, logical(1))
  expect_identical(which(!refused), integer())

  # alpha and A are not needed (and ignored) for one unit
  expect_s3_class(r2d_spec("exact", 1, 1, alpha = 5, A = 2, n_units = 1), "BayesTools_R2_functional")
  expect_s3_class(r2d_signal_spec(1, 1), "BayesTools_R2_signal")
})

test_that("a quantile beyond the double range is a classed unavailability, not a wrong value", {

  # P(|X| > x) ~ x^(-2 b) for b = 0.5: the probability 1e-320 lies beyond x = exp(700)
  spec  <- r2d_spec("exact", 1, 0.5)
  error <- tryCatch(r2d_quantile(spec, 1e-320), error = identity)
  expect_s3_class(error, "BayesTools_R2_density_unavailable")
  expect_s3_class(error, "BayesTools_R2_unavailable")
  expect_s3_class(error, "BayesTools_R2_condition")
  expect_match(conditionMessage(error), "rejected by diagnostics", fixed = TRUE)
})

test_that("two-level laws through the log-product density match the nested integrals", {

  specs <- list(
    .bt_R2_functional_spec("exact", a = 1, b = 2, alpha = 1, A = 2, n_units = 2, d = 2, scale = 1),
    .bt_R2_functional_spec("exact", a = 1, b = 3, alpha = 0.95, A = 1.2, n_units = 2, d = 2, scale = 1),
    .bt_R2_functional_spec("normal", a = 0.3, b = 0.6, alpha = 0.2, A = 0.9, n_units = 3, scale = 1),
    .bt_R2_functional_spec("laplace", a = 1.5, b = 1, alpha = 0.5, A = 2, n_units = 4, scale = 1)
  )
  x <- c(1e-3, 0.3, 1, 4, 50)
  problems <- expectation_problems({
    for(spec in specs){
      for(what in c("pdf", "tail")){
        product <- vapply(x, function(z) .bt_R2_density_log(spec, z, what, "product"), numeric(1))
        nested <- vapply(x, function(z) .bt_R2_density_log(spec, z, what, "nested"), numeric(1))
        expect_equal_each(product, nested, tolerance = 1e-9)
      }
    }
  })
  expect_identical(problems, character())

  # the log-product density integrates to one, and its quadrature agrees with
  # stats::integrate() where it is used as the fallback
  first <- c(1, 0.5)
  second <- c(0.5, 0.5)
  mass <- stats::integrate(function(u) exp(.bt_R2_density_log_h(u, first, second, FALSE, 1e-10)), -Inf, 0,
    rel.tol = 1e-9)$value
  expect_equal(mass, 1, tolerance = 1e-7)
})
