context("Prior distribution coverage tests")

# Targeted coverage tests for priors.R uncovered lines


test_that("Unknown distribution name error (line 115)", {
  expect_error(prior("unknown_dist", list(0, 1)),
               "The specified distribution name")
})


test_that("prior_factor orthonormal/meandif with non-vector prior error (lines 328, 334, 336)", {
  # orthonormal contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "orthonormal"),
               "contrasts require multivariate prior")

  # meandif contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "meandif"),
               "contrasts require multivariate prior")
})


test_that("prior_factor orthonormal/meandif with unsupported distribution (line 346)",
{
  # bernoulli is not a valid multivariate distribution
  expect_error(prior_factor("bernoulli", list(0.5), contrast = "orthonormal"),
               "contrasts require multivariate prior")
})


test_that("prior_factor treatment contrast requires univariate (line 358)", {
  # treatment contrast with multivariate dist should fail
  # mnormal requires 3 params, so it fails earlier - use a valid multivariate prior
  expect_error(prior_factor("mnormal", list(0, 1, 2), contrast = "treatment"),
               "contrasts require univariate prior")
})


test_that("prior_spike_and_slab with factor prior (lines 402-408)", {
  # spike_and_slab with factor prior as variable
  p_factor <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  p_factor$parameters[["K"]] <- 2

  p_ss <- prior_spike_and_slab(
    prior_parameter = p_factor,
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_true(is.prior.spike_and_slab(p_ss))
  expect_s3_class(p_ss, "prior.spike_and_slab")
})


test_that(".set_spike_and_slab_variable_attr (lines 486-497)", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  # Set an attribute on the variable component
  p_ss2 <- BayesTools:::.set_spike_and_slab_variable_attr(p_ss, "test_attr", "test_value")
  expect_true(is.prior.spike_and_slab(p_ss2))

  # Error when not spike_and_slab
  expect_error(BayesTools:::.set_spike_and_slab_variable_attr(prior("normal", list(0, 1)), "attr", "val"),
               "only works with spike_and_slab priors")
})


test_that(".get_spike_and_slab_variable error (line 459)", {
  expect_error(BayesTools:::.get_spike_and_slab_variable(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


test_that(".get_spike_and_slab_inclusion error (lines 471, 476)", {
  expect_error(BayesTools:::.get_spike_and_slab_inclusion(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


test_that("prior_mixture with factor prior containing spike (lines 522, 538)", {
  # Mixture of factor priors where one is a spike
  p1 <- prior("spike", list(0))
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture with prior_none factor (line 576)", {
  # Mixture with prior_none that should be converted to factor spike
  p1 <- prior_none()
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture bias mixture (lines 585, 600)", {
  # PET/PEESE/weightfunction mixture
  p_pet <- prior_PET("normal", list(0, 1))
  p_wf  <- prior_weightfunction("one.sided", list(steps = c(0.05), alpha = c(1, 1)))

  p_mix <- prior_mixture(list(p_pet, p_wf), components = c("a", "b"))

  expect_s3_class(p_mix, "prior.bias_mixture")
})


test_that("Uniform prior with a > b error (line 843)", {
  expect_error(prior("uniform", list(a = 5, b = 1)),
               "lower than")
})


test_that("rng spike_and_slab sample_components (line 1155)", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  set.seed(1)
  components <- rng(p_ss, 10, sample_components = TRUE)

  expect_true(all(components %in% c(0, 1)))
  expect_length(components, 10)
})


test_that("rng mixture sample_components (line 1190)", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )

  set.seed(1)
  components <- rng(p_mix, 10, sample_components = TRUE)

  expect_true(all(components %in% 1:2))
  expect_length(components, 10)
})


test_that("rng factor_mixture (lines 1221, 1236, 1240, 1244-1247)", {
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
  # Default: transform_factor_samples = FALSE
  samples <- rng(p_mix, 10, transform_factor_samples = FALSE)
  expect_true(is.matrix(samples))

  # With transform_factor_samples = TRUE
  samples2 <- rng(p_mix, 10, transform_factor_samples = TRUE)
  expect_true(is.matrix(samples2))
  expect_equal(ncol(samples2), 3)  # K+1 columns
})


test_that("rng orthonormal/meandif transform (line 1284)", {
  p <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p$parameters[["K"]] <- 2

  set.seed(1)
  samples <- rng(p, 10, transform_factor_samples = TRUE)

  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), 3)  # K+1 columns
})


test_that("cdf with truncation for spike_and_slab error (lines 1329, 1333)", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_error(cdf(p_ss, 0), "No cdfs are implemented for spike and slab")
})


test_that("cdf with truncation for simple prior (lines 1363, 1365, 1367, 1369)", {
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))

  # Test cdf at various points
  expect_true(cdf(p, 0) > 0)
  expect_true(cdf(p, -3) == 0)  # Below truncation
  expect_true(cdf(p, 3) >= 1 - 1e-6)  # Above truncation (approx 1)
})


test_that("ccdf errors and truncation (lines 1385, 1389, 1419-1425)", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_error(ccdf(p_ss, 0), "No ccdf are implemented for spike and slab")

  # Mixture error
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(ccdf(p_mix, 0), "No ccdf are implemented for prior mixtures")

  # ccdf with truncation
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))
  expect_true(ccdf(p, 0) > 0)
  expect_true(ccdf(p, 3) == 0)  # Above truncation
})


test_that("lpdf spike_and_slab and mixture errors (lines 1442, 1446, 1496, 1498)", {
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


test_that("quant spike_and_slab and mixture errors (lines 1528, 1532)", {
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


test_that("quant with non-default truncation optimization (lines 1561, 1582, 1584)", {
  # Truncated prior that requires optimization in quant
  p <- prior("normal", list(0, 1), truncation = list(0.5, 2))

  q <- quant(p, 0.5)
  expect_true(q > 0.5 && q < 2)

  # Also test edge quantiles
  q_low <- quant(p, 0.01)
  q_high <- quant(p, 0.99)
  expect_true(q_low >= 0.5)
  expect_true(q_high <= 2)
})


test_that("mcdf for orthonormal/meandif (lines 1681, 1685, 1689, 1711-1722)", {
  # orthonormal prior
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  cdf_val <- mcdf(p_orth, 0)
  expect_true(cdf_val >= 0 && cdf_val <= 1)

  # meandif prior
  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  cdf_val2 <- mcdf(p_md, 0)
  expect_true(cdf_val2 >= 0 && cdf_val2 <= 1)
})


test_that("mccdf for orthonormal/meandif (lines 1760-1801)", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  ccdf_val <- mccdf(p_orth, 0)
  expect_true(ccdf_val >= 0 && ccdf_val <= 1)

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  ccdf_val2 <- mccdf(p_md, 0)
  expect_true(ccdf_val2 >= 0 && ccdf_val2 <= 1)
})


test_that("mlpdf for orthonormal/meandif (lines 1840, 1844, 1870, 1874)", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  lpdf_val <- mlpdf(p_orth, 0)
  expect_true(is.finite(lpdf_val))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  lpdf_val2 <- mlpdf(p_md, 0)
  expect_true(is.finite(lpdf_val2))
})


test_that("mquant for orthonormal/meandif (lines 1933, 1937, 1963, 1967)", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  q <- mquant(p_orth, 0.5)
  expect_true(is.numeric(q))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  q2 <- mquant(p_md, 0.5)
  expect_true(is.numeric(q2))
})


test_that("pdf.default passes through to stats::dnorm for vectors", {
  # pdf.default calls stats::dnorm for numeric vectors, not an error
  # This tests the generic S3 dispatch working correctly
  p <- prior("normal", list(0, 1))
  expect_true(is.numeric(pdf(p, 0)))
})


test_that("mean.prior for spike_and_slab (line 2123)", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(1, 1)),
    prior_inclusion = prior("point", list(0.5))  # Fixed inclusion probability
  )

  m <- mean(p_ss)
  expect_true(is.numeric(m))
  expect_equal(m, 0.5)  # mean(normal(1,1)) * 0.5
})


test_that("mean.prior for truncated distributions (lines 2148, 2153)", {
  # Truncated normal
  p <- prior("normal", list(0, 1), truncation = list(0, Inf))
  m <- mean(p)
  expect_true(m > 0)

  # Truncated t with df <= 1 should return NaN
  p_t <- prior("t", list(0, 1, 1), truncation = list(-1, 1))
  m_t <- mean(p_t)
  expect_true(is.nan(m_t))
})


test_that("mean.prior for orthonormal/meandif with mt df<=1 (lines 2181, 2185, 2189)", {
  p_mt <- prior_factor("mt", list(0, 1, 1), contrast = "orthonormal")  # df = 1
  p_mt$parameters[["K"]] <- 2

  m <- mean(p_mt)
  expect_true(is.nan(m))
})


test_that("var dispatches to stats::var for vectors", {
  # var.default calls stats::var for numeric vectors
  x <- c(1, 2, 3, 4, 5)
  expect_equal(var(x), stats::var(x))
})


test_that("var.prior for spike_and_slab (lines 2276-2291)", {
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


test_that("var.prior for truncated distributions (lines 2316, 2321)", {
  # t with df <= 2 should return NaN for variance
  p_t <- prior("t", list(0, 1, 2), truncation = list(-1, 1))
  v <- var(p_t)
  expect_true(is.nan(v))

  # invgamma with shape <= 2 should return NaN
  p_ig <- prior("invgamma", list(2, 1), truncation = list(0.1, 10))
  v_ig <- var(p_ig)
  expect_true(is.nan(v_ig))
})


test_that("var.prior for orthonormal/meandif (lines 2350-2368)", {
  # orthonormal with mpoint
  p_mp <- prior_factor("mpoint", list(0), contrast = "orthonormal")
  p_mp$parameters[["K"]] <- 2
  v <- var(p_mp)
  expect_equal(v, 0)

  # orthonormal with mt and df <= 2
  p_mt <- prior_factor("mt", list(0, 1, 2), contrast = "orthonormal")
  p_mt$parameters[["K"]] <- 2
  v_mt <- var(p_mt)
  expect_true(is.nan(v_mt))

  # meandif with mnormal
  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2
  v_md <- var(p_md)
  expect_true(v_md > 0)
})


test_that("var.prior for mixture error", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(var(p_mix), "No var is implemented for prior mixtures")
})
