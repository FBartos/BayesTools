skip_if_not_test_profile("unit")

test_that("model-averaging compatibility guards report metadata mismatches", {

  fake_fit <- structure(list(), class = "null_model")

  simple_a <- prior("normal", list(0, 1))
  simple_b <- prior("normal", list(0, 1))
  attr(simple_b, "interaction") <- TRUE
  attr(simple_b, "interaction_terms") <- c("x", "y")

  expect_error(
    .mix_posteriors.simple(
      fits = list(fake_fit, fake_fit),
      priors = list(simple_a, simple_b),
      parameter = "theta",
      post_probs = c(0.5, 0.5),
      seed = 1,
      n_samples = 10
    ),
    "non-matching prior factor type specifications",
    fixed = TRUE
  )

  factor_a <- prior_factor(
    "normal",
    list(0, 1),
    contrast = "treatment"
  )
  factor_b <- prior_factor(
    "normal",
    list(0, 1),
    contrast = "treatment"
  )
  attr(factor_a, "levels") <- 3
  attr(factor_b, "levels") <- 3
  attr(factor_a, "level_names") <- c("A", "B", "C")
  attr(factor_b, "level_names") <- c("X", "Y", "Z")

  expect_error(
    .mix_priors.factor(
      priors = list(factor_a, factor_b),
      parameter = "fac",
      seed = 1,
      n_samples = 10
    ),
    "non-matching prior factor type specifications",
    fixed = TRUE
  )

  expect_error(
    .mix_posteriors.factor(
      fits = list(fake_fit, fake_fit),
      priors = list(factor_a, factor_b),
      parameter = "fac",
      post_probs = c(0.5, 0.5),
      seed = 1,
      n_samples = 10
    ),
    "non-matching prior factor type specifications",
    fixed = TRUE
  )
})
