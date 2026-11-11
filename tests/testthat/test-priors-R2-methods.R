skip_if_not_test_profile("unit")

# The empirical CDF of draws against a CDF at a few points, within 4 Monte
# Carlo standard errors.
.R2_expect_ecdf <- function(draws, cdf, q){
  p <- cdf(q)
  empirical <- vapply(q, function(v) mean(draws <= v), numeric(1))
  expect_true(all(abs(empirical - p) <= 4 * sqrt(p * (1 - p) / length(draws))))
}

.R2_methods_priors <- function(law = "exact", contrast = "orthonormal", n_per = c(8, 8, 8), alpha = c(x = 1, f = 1),
                               a = 1, b = 2, convention = "common"){
  set.seed(4)
  data <- data.frame(f = factor(rep(c("a", "b", "c"), times = n_per)), x = stats::rnorm(sum(n_per)))
  budget <- prior_R2(alpha = alpha, a = if(convention == "linked") NULL else a, b = b, law = law,
    convention = convention, concentration = "inherited", sigma = "sigma", contrast = contrast)
  info <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))
  info$prior_list
}

test_that("marginal methods of a scalar R2 term are the implied functional law", {

  priors <- .R2_methods_priors(b = 3)
  prior <- priors$mu_x
  spec <- .bt_R2_term_functional(prior, 1)
  # exact law, two units: share Beta(1, 1); the scale is the inverse SD of x
  expect_equal(spec$scale, abs(prior$map[1, 1]))
  expect_equal(c(spec$alpha, spec$A, spec$n_units), c(1, 2, 2))
  q <- c(-1, -0.2, 0.3, 2)
  expect_equal(mpdf(prior, q), .bt_R2_functional_pdf(spec, q))
  expect_equal(pdf(prior, q), mpdf(prior, q))
  expect_equal(mcdf(prior, q), .bt_R2_functional_cdf(spec, q))
  expect_equal(cdf(prior, q) + ccdf(prior, q), rep(1, 4), tolerance = 1e-12)
  expect_equal(mcdf(prior, mquant(prior, c(0.1, 0.9))), c(0.1, 0.9), tolerance = 1e-8)
  expect_equal(lpdf(prior, q), log(pdf(prior, q)))
  expect_identical(mean(prior), 0)
  # Var = scale^2 E[W] E[phi] = scale^2 (a / (b - 1)) (alpha / A)
  expect_equal(var(prior), spec$scale^2 * (1 / 2) * (1 / 2), tolerance = 1e-12)
  expect_true(is.nan(var(.R2_methods_priors(b = 1)$mu_x)))
  expect_true(is.nan(mean(.R2_methods_priors(b = 0.4)$mu_x)))
  expect_equal(range(prior), mquant(prior, c(0.01, 0.99)))

  draws <- rng(prior, 4000)
  expect_length(draws, 4000)
  .R2_expect_ecdf(draws, function(v) mcdf(prior, v), stats::quantile(draws, c(0.1, 0.3, 0.7, 0.9), names = FALSE))
})

test_that("factor level marginals use the declared design rows under the exact law", {

  prior <- .R2_methods_priors(contrast = "meandif")$mu_f
  rows <- attr(prior, "factor_design")
  spec <- .bt_R2_term_functional(prior, rows[1, ])
  expect_equal(spec$scale, sqrt(sum((t(prior$map) %*% rows[1, ])^2)))
  expect_equal(spec$d, 2L)
  expect_equal(mpdf(prior, c(0.1, 0.5)), .bt_R2_functional_pdf(spec, c(0.1, 0.5)))
  levels <- rng(prior, 3000)
  expect_identical(dim(levels), c(3000L, 3L))
  expect_equal(rowSums(levels), rep(0, 3000), tolerance = 1e-12)
  .R2_expect_ecdf(levels[, 2], function(v) mcdf(prior, v), c(-0.5, -0.1, 0.2))
  raw <- rng(prior, 10, transform_factor_samples = FALSE)
  expect_identical(dim(raw), c(10L, 2L))
  expect_error(quant(prior, 0.5), "Joint quantile functions are unavailable for factor priors", fixed = TRUE)
})

test_that("unavailable marginals refuse with classed conditions", {

  # unequal cell counts: the levels have different marginals
  unequal <- .R2_methods_priors(n_per = c(4, 8, 12))$mu_f
  expect_error(mpdf(unequal, 0.2), class = "BayesTools_R2_density_unavailable")
  expect_identical(dim(rng(unequal, 5)), c(5L, 3L))

  # coordinatewise kernels have no closed level marginal for three levels
  kernel <- .R2_methods_priors(law = "normal")$mu_f
  expect_error(mpdf(kernel, 0.2), class = "BayesTools_R2_unavailable")

  # treatment coefficients of a kernel are single coded coefficients
  treatment <- .R2_methods_priors(law = "laplace", contrast = "treatment")$mu_f
  expect_true(all(is.finite(mpdf(treatment, c(-0.3, 0.3)))))

  unbound <- prior_R2(alpha = c(x = 1), a = 1, b = 1, law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma")
  stub <- .bt_R2_term_stub(unbound, "x", "x", "continuous")
  expect_error(mpdf(stub, 0), class = "BayesTools_R2_unbound")
})

test_that("R2 term priors print their budget", {

  priors <- .R2_methods_priors(convention = "linked")
  expect_identical(print(priors$mu_x, silent = TRUE), "R2 exact: R2 ~ Beta(2, 2), allocation 1 of 2")
  expect_identical(print(priors$mu_x, silent = TRUE, short_name = TRUE), "R2 exact(2, 2; 1 of 2)")
  expect_output(print(priors$mu_f), "R2 exact: R2 ~ Beta(2, 2), allocation 1 of 2", fixed = TRUE)
})

test_that("prior densities of R2 terms come from the deterministic marginal", {

  priors <- .R2_methods_priors()
  scalar <- density(priors$mu_x, x_range = c(-2, 2), n_points = 21)
  expect_s3_class(scalar, "density.prior")
  expect_equal(scalar$y, mpdf(priors$mu_x, scalar$x))
  levels <- density(priors$mu_f, x_range = c(-2, 2), n_points = 21)
  expect_equal(levels$y, mpdf(priors$mu_f, levels$x))
})

test_that("joint prior draws keep the budget shared across terms", {

  priors <- .R2_methods_priors()
  draws <- .generate_prior_sample_matrix(priors, 500L, seed = 1)
  expect_true(all(c("mu_R2", "mu_R2_allocation[1]", "mu_x", "mu_f[1]", "mu_f[2]") %in% colnames(draws)))
  # the exact signal of each term is W times its share, from the same draw
  W <- draws[, "mu_R2"] / (1 - draws[, "mu_R2"])
  signal_x <- draws[, "mu_x"]^2 * priors$mu_x$metric[1, 1]
  expect_equal(signal_x, W * draws[, "mu_R2_allocation[1]"], tolerance = 1e-10)
  beta_f <- draws[, c("mu_f[1]", "mu_f[2]")]
  signal_f <- rowSums((beta_f %*% priors$mu_f$metric) * beta_f)
  expect_equal(signal_f, W * draws[, "mu_R2_allocation[2]"], tolerance = 1e-10)
})

test_that("a single-term budget has the marginal of its unit without an allocation", {

  set.seed(4)
  data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 8)), x = stats::rnorm(24))
  budget <- prior_R2(alpha = c(f = 1), a = 1, b = 3, law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma", contrast = "orthonormal")
  prior <- JAGS_formula(~ f, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))$prior_list$mu_f
  expect_identical(prior$n_units, 1L)
  expect_true(all(is.finite(mpdf(prior, c(-0.5, 0.5)))))
  expect_equal(var(prior), .bt_R2_functional_variance(.bt_R2_marginal_spec(prior)))
  expect_length(range(prior), 2L)
})

test_that("public tail and log-density methods keep representable tails", {

  set.seed(4)
  data <- data.frame(x = stats::rnorm(24))
  budget <- prior_R2(alpha = c(x = 1), a = 1, b = 2, law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma")
  prior <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))$prior_list$mu_x
  spec <- .bt_R2_marginal_spec(prior)
  far <- 1e5 * spec$scale
  expect_equal(ccdf(prior, far), .bt_R2_functional_cdf(spec, far, lower.tail = FALSE))
  expect_gt(ccdf(prior, far), 0)
  expect_equal(mccdf(prior, far), ccdf(prior, far))
  expect_true(is.finite(lpdf(prior, 1e100 * spec$scale)))
  expect_equal(mlpdf(prior, 1e100 * spec$scale), .bt_R2_functional_log_pdf(spec, 1e100 * spec$scale))
})
