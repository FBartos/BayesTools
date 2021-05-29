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
})
