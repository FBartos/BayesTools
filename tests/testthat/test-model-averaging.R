context("Model-averaging functions")

test_that("Model-averaging functions work", {

  expect_equal(compute_inference(c(1,1), c(1, 1))$prior_probs, c(0.5, 0.5))
  expect_equal(compute_inference(c(1,1), c(1, 1))$post_probs,  c(0.5, 0.5))
  expect_equal(compute_inference(c(1,1), c(1, 1))$BF,          Inf)
  expect_equal(attr(compute_inference(c(1,1), c(1, 1)), "is_null"), c(FALSE, FALSE))

  expect_equal(compute_inference(c(1,4),   c(1, 1))$prior_probs,    c(0.2, 0.8))
  expect_equal(compute_inference(c(1,1,3), c(1, 1, 1))$prior_probs, c(0.2, 0.2, 0.6))
  expect_equal(compute_inference(c(1,1,4), c(1, 1, 1), c(F, T, F), conditional = TRUE)$prior_probs, c(0.2, 0, 0.8))

  expect_equal(compute_inference(c(1,4),   c(1, 1))$post_probs,    c(0.2, 0.8))
  expect_equal(compute_inference(c(1,1,3), c(1, 1, 1))$post_probs, c(0.2, 0.2, 0.6))
  expect_equal(compute_inference(c(1,1,4), c(1, 1, 1), c(F, T, F), conditional = TRUE)$post_probs, c(0.2, 0, 0.8))
  expect_equal(attr(compute_inference(c(1,1,4), c(1, 1, 1), c(2)), "is_null"), c(F, T, F))

  # automatically tests inclusion_bf as well
  expect_equal(compute_inference(c(1,1), c(1, 1), 1)$BF, 1)
  expect_equal(compute_inference(c(1,1), c(1, 2), c(F, T))$BF, exp(1-2))
  expect_equal(compute_inference(c(1,1,1), c(1, 1, 1), c(F, T, F))$BF, 1)
  expect_equal(compute_inference(c(1,1,1), c(1, 2, 1), c(F, T, F))$BF, exp(1-2))

  # and check BF formatting
  expect_equivalent(format_BF(c(0, 1, 2, Inf)), c(0, 1, 2, Inf))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), BF01 = TRUE), 1/c(0, 1, 2, Inf))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), logBF = TRUE), log(c(0, 1, 2, Inf)))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), BF01 = TRUE, logBF = TRUE), log(1/c(0, 1, 2, Inf)))
  expect_equal(attr(format_BF(1), "name"), "BF")
  expect_equal(attr(format_BF(1, logBF = TRUE), "name"), "log(BF)")
  expect_equal(attr(format_BF(1, BF01 = TRUE, logBF = TRUE), "name"), "log(1/BF)")

  # additional BF checks
  expect_equal(inclusion_BF(prior_probs = c(.5, .5), post_probs = c(.5, .5), is_null = c(T, F)), 1)
  expect_equal(inclusion_BF(prior_probs = c(.5, .5), post_probs = c(.75, .25), is_null = c(T, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(.25, .25, .25, .25), post_probs = c(.75, 0, .25, 0), is_null = c(T, T, F, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(.25, .25, .25, .25), post_probs = c(.65, .10, .20, 0.05), is_null = c(T, T, F, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(1, 0), post_probs = c(1, 0), is_null = c(T, F)), 0)
  expect_equal(inclusion_BF(prior_probs = c(1, 0), post_probs = c(1, 0), is_null = c(F, T)), Inf)

  # test the marglik versions of BF
  temp_prior_probs <- 1:6/sum(1:6)
  temp_margliks    <- -2:3
  temp_post_probs  <- bridgesampling::post_prob(temp_margliks, prior_prob = temp_prior_probs)
  expect_equal(
    inclusion_BF(prior_probs = temp_prior_probs, post_probs = temp_post_probs, is_null = rep(c(T, F), 3)),
    inclusion_BF(prior_probs = temp_prior_probs, margliks = temp_margliks, is_null = rep(c(T, F), 3))
  )

  # check for over/underflow
  temp_prior_probs <- 1:6/sum(1:6)
  temp_margliks    <- c(-2:2, 100)
  temp_post_probs  <- bridgesampling::post_prob(temp_margliks, prior_prob = temp_prior_probs)
  expect_true(is.infinite(inclusion_BF(prior_probs = temp_prior_probs, post_probs = temp_post_probs, is_null = rep(c(T, F), 3))))
  expect_false(is.infinite(inclusion_BF(prior_probs = temp_prior_probs, margliks = temp_margliks, is_null = rep(c(T, F), 3))))

  # additional omega mapping checks
  expect_equal(weightfunctions_mapping(prior_list = list(
    prior_none(),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),        prior_weights = 1/2),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),  prior_weights = 1/2)
  )), list(NULL, c(2, 1, 1), c(3, 2, 1)))
})
