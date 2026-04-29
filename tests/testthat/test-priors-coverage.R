# ============================================================================ #
# TEST FILE: Prior Distribution Coverage Tests
# ============================================================================ #
#
# PURPOSE:
#   Targeted tests to ensure code coverage for edge cases and error paths
#   in R/priors.R that are not covered by main test-priors.R
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#
# SKIP CONDITIONS:
#   - None (fast, pure R tests)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @coverage, @edge-cases, @priors
# ============================================================================ #


# ============================================================================ #
# SECTION: prior() constructor errors
# ============================================================================ #

test_that("prior() rejects unknown distribution names", {
  expect_error(prior("unknown_dist", list(0, 1)),
               "The specified distribution name")
})


# ============================================================================ #
# SECTION: prior_factor() contrast validation
# ============================================================================ #

test_that("prior_factor() requires multivariate prior for orthonormal/meandif contrasts", {
  # orthonormal contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "orthonormal"),
               "contrasts require multivariate prior")

  # meandif contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "meandif"),
               "contrasts require multivariate prior")

  # bernoulli is not a valid multivariate distribution
  expect_error(prior_factor("bernoulli", list(0.5), contrast = "orthonormal"),
               "contrasts require multivariate prior")
})


test_that("prior_factor() requires univariate prior for treatment contrast", {
  expect_error(prior_factor("mnormal", list(0, 1, 2), contrast = "treatment"),
               "contrasts require univariate prior")
})


# ============================================================================ #
# SECTION: spike_and_slab prior construction and helpers
# ============================================================================ #

test_that("prior_spike_and_slab() works with factor priors", {
  p_factor <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  p_factor$parameters[["K"]] <- 2

  p_ss <- prior_spike_and_slab(
    prior_parameter = p_factor,
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_true(is.prior.spike_and_slab(p_ss))
  expect_s3_class(p_ss, "prior.spike_and_slab")
})


test_that(".set_spike_and_slab_variable_attr() sets attributes correctly", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  p_ss2 <- BayesTools:::.set_spike_and_slab_variable_attr(p_ss, "test_attr", "test_value")
  expect_true(is.prior.spike_and_slab(p_ss2))

  # Error when not spike_and_slab
  expect_error(BayesTools:::.set_spike_and_slab_variable_attr(prior("normal", list(0, 1)), "attr", "val"),
               "only works with spike_and_slab priors")
})


test_that(".get_spike_and_slab_variable() requires spike_and_slab prior", {
  expect_error(BayesTools:::.get_spike_and_slab_variable(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


test_that(".get_spike_and_slab_inclusion() requires spike_and_slab prior", {
  expect_error(BayesTools:::.get_spike_and_slab_inclusion(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


# ============================================================================ #
# SECTION: prior_mixture() construction
# ============================================================================ #

test_that("prior_mixture() creates factor_mixture with spike component", {
  p1 <- prior("spike", list(0))
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture() handles prior_none in factor mixtures", {
  p1 <- prior_none()
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture() creates bias_mixture for PET/PEESE/weightfunction", {
  p_pet <- prior_PET("normal", list(0, 1))
  p_wf  <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  p_mix <- prior_mixture(list(p_pet, p_wf), components = c("a", "b"))

  expect_s3_class(p_mix, "prior.bias_mixture")
})


# ============================================================================ #
# SECTION: Distribution parameter validation
# ============================================================================ #

test_that("uniform prior requires a < b", {
  expect_error(prior("uniform", list(a = 5, b = 1)),
               "lower than")
})


# ============================================================================ #
# SECTION: rng() function with sample_components
# ============================================================================ #

test_that("rng() spike_and_slab returns component indicators", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  set.seed(1)
  components <- rng(p_ss, 10, sample_components = TRUE)

  expect_true(all(components %in% c(0, 1)))
  expect_length(components, 10)
})


test_that("rng() mixture returns component indicators", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )

  set.seed(1)
  components <- rng(p_mix, 10, sample_components = TRUE)

  expect_true(all(components %in% 1:2))
  expect_length(components, 10)
})


test_that("rng() factor_mixture with transform_factor_samples", {
  p_mix <- prior_mixture(
    list(
      prior("spike", list(0)),
      prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    ),
    components = c("null", "alt")
  )
  for (i in seq_along(p_mix)) {
    p_mix[[i]]$parameters[["K"]] <- 2
  }

  set.seed(1)
  samples <- rng(p_mix, 10, transform_factor_samples = FALSE)
  expect_true(is.matrix(samples))

  samples2 <- rng(p_mix, 10, transform_factor_samples = TRUE)
  expect_true(is.matrix(samples2))
  expect_equal(ncol(samples2), 3)  # K+1 columns
})


test_that("rng() orthonormal prior with transform_factor_samples", {
  p <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p$parameters[["K"]] <- 2

  set.seed(1)
  samples <- rng(p, 10, transform_factor_samples = TRUE)

  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), 3)  # K+1 columns
})


# ============================================================================ #
# SECTION: cdf() function edge cases
# ============================================================================ #

test_that("cdf() not implemented for spike_and_slab", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_error(cdf(p_ss, 0), "No cdfs are implemented for spike and slab")
})


test_that("cdf() handles truncated priors correctly", {
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))

  expect_true(cdf(p, 0) > 0)
  expect_true(cdf(p, -3) == 0)  # Below truncation
  expect_true(cdf(p, 3) >= 1 - 1e-6)  # Above truncation
})


# ============================================================================ #
# SECTION: ccdf() function edge cases
# ============================================================================ #

test_that("ccdf() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(ccdf(p_ss, 0), "No ccdf are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(ccdf(p_mix, 0), "No ccdf are implemented for prior mixtures")
})


test_that("ccdf() handles truncated priors correctly", {
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))
  expect_true(ccdf(p, 0) > 0)
  expect_true(ccdf(p, 3) == 0)  # Above truncation
})


# ============================================================================ #
# SECTION: lpdf() function edge cases
# ============================================================================ #

test_that("lpdf() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(lpdf(p_ss, 0), "No lpdf are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(lpdf(p_mix, 0), "No lpdf are implemented for prior mixtures")
})


# ============================================================================ #
# SECTION: quant() function edge cases
# ============================================================================ #

test_that("quant() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(quant(p_ss, 0.5), "No quant(ile)? functions? are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(quant(p_mix, 0.5), "No quant(ile)? functions? are implemented for prior mixtures")
})


test_that("quant() handles truncated priors with optimization", {
  p <- prior("normal", list(0, 1), truncation = list(0.5, 2))

  q <- quant(p, 0.5)
  expect_true(q > 0.5 && q < 2)

  q_low <- quant(p, 0.01)
  q_high <- quant(p, 0.99)
  expect_true(q_low >= 0.5)
  expect_true(q_high <= 2)
})


# ============================================================================ #
# SECTION: Multivariate distribution functions (mcdf, mccdf, mlpdf, mquant)
# ============================================================================ #

test_that("mcdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  cdf_val <- mcdf(p_orth, 0)
  expect_true(cdf_val >= 0 && cdf_val <= 1)

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  cdf_val2 <- mcdf(p_md, 0)
  expect_true(cdf_val2 >= 0 && cdf_val2 <= 1)
})


test_that("mccdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  ccdf_val <- mccdf(p_orth, 0)
  expect_true(ccdf_val >= 0 && ccdf_val <= 1)

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  ccdf_val2 <- mccdf(p_md, 0)
  expect_true(ccdf_val2 >= 0 && ccdf_val2 <= 1)
})


test_that("mlpdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  lpdf_val <- mlpdf(p_orth, 0)
  expect_true(is.finite(lpdf_val))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  lpdf_val2 <- mlpdf(p_md, 0)
  expect_true(is.finite(lpdf_val2))
})


test_that("mquant() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  q <- mquant(p_orth, 0.5)
  expect_true(is.numeric(q))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  q2 <- mquant(p_md, 0.5)
  expect_true(is.numeric(q2))
})


# ============================================================================ #
# SECTION: S3 dispatch and generic functions
# ============================================================================ #

test_that("pdf() S3 dispatch works for prior objects", {
  p <- prior("normal", list(0, 1))
  expect_true(is.numeric(pdf(p, 0)))
})


# ============================================================================ #
# SECTION: mean() function edge cases
# ============================================================================ #

test_that("mean() works for spike_and_slab priors", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(1, 1)),
    prior_inclusion = prior("point", list(0.5))
  )

  m <- mean(p_ss)
  expect_true(is.numeric(m))
  expect_equal(m, 0.5)  # mean(normal(1,1)) * 0.5
})


test_that("mean() handles truncated distributions and undefined moments", {
  # Truncated normal
  p <- prior("normal", list(0, 1), truncation = list(0, Inf))
  m <- mean(p)
  expect_true(m > 0)

  # Truncated t with df <= 1 returns NaN
  p_t <- prior("t", list(0, 1, 1), truncation = list(-1, 1))
  m_t <- mean(p_t)
  expect_true(is.nan(m_t))
})


test_that("mean() returns NaN for multivariate t with df <= 1", {
  p_mt <- prior_factor("mt", list(0, 1, 1), contrast = "orthonormal")
  p_mt$parameters[["K"]] <- 2

  m <- mean(p_mt)
  expect_true(is.nan(m))
})


# ============================================================================ #
# SECTION: var() function edge cases
# ============================================================================ #

test_that("var() S3 dispatch works for numeric vectors", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(var(x), stats::var(x))
})


test_that("var() works for spike_and_slab priors", {
  # spike_and_slab with beta inclusion
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  v <- var(p_ss)
  expect_true(is.numeric(v))
  expect_true(v > 0)

  # spike_and_slab with point inclusion
  p_ss2 <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0.5))
  )
  v2 <- var(p_ss2)
  expect_true(is.numeric(v2))
})


test_that("var() returns NaN for distributions with undefined variance", {
  # t with df <= 2 returns NaN for variance
  p_t <- prior("t", list(0, 1, 2), truncation = list(-1, 1))
  v <- var(p_t)
  expect_true(is.nan(v))

  # invgamma with shape <= 2 returns NaN
  p_ig <- prior("invgamma", list(2, 1), truncation = list(0.1, 10))
  v_ig <- var(p_ig)
  expect_true(is.nan(v_ig))
})


test_that("var() works for orthonormal and meandif priors", {
  # orthonormal with mpoint returns 0
  p_mp <- prior_factor("mpoint", list(0), contrast = "orthonormal")
  p_mp$parameters[["K"]] <- 2
  v <- var(p_mp)
  expect_equal(v, 0)

  # orthonormal with mt and df <= 2 returns NaN
  p_mt <- prior_factor("mt", list(0, 1, 2), contrast = "orthonormal")
  p_mt$parameters[["K"]] <- 2
  v_mt <- var(p_mt)
  expect_true(is.nan(v_mt))

  # meandif with mnormal returns positive variance
  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2
  v_md <- var(p_md)
  expect_true(v_md > 0)
})


test_that("var() not implemented for mixture priors", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(var(p_mix), "No var is implemented for prior mixtures")
})
