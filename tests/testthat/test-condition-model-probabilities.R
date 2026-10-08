skip_if_not_test_profile("unit")

test_that("bias posterior masks use their declared event family", {
  priors <- list(bias = prior_PET("normal", list(0, 1)))
  draws <- matrix(c(1, 2), 2L, dimnames = list(NULL, "bias_indicator"))
  expect_identical(.condition_event_label_posterior_mask(priors, draws, "bias"), c(TRUE, TRUE))
  expect_identical(.condition_event_label_posterior_mask(list(bias = prior_none()), draws, "bias"), c(FALSE, FALSE))
  p <- prior_mixture(list(prior_none(), priors$bias), is_null = c(TRUE, FALSE))
  for(rule in c("AND", "OR")){
    event <- .condition_event(list(bias = p), c("bias", "PET"), rule)
    expect_identical(.condition_event_posterior_mask(event, list(bias = p), draws), c(FALSE, TRUE))
  }
  expect_identical(.condition_event_label_posterior_mask(priors, draws, "PET"), c(TRUE, TRUE))
  expect_identical(.condition_event_label_posterior_mask(priors, draws, "PEESE"), c(FALSE, FALSE))
  expect_error(.condition_event_label_posterior_mask(list(bias = prior("normal", list(0, 1))), draws, "bias"),
    "The parameter 'bias' is not a conditional parameter.", fixed = TRUE)
})

test_that("declared positive option products survive event underflow", {
  priors <- list(a = prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(1e-200))),
    b = prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(1e-200))))
  event <- .condition_event(priors, c("a", "b"), "AND")
  options <- .condition_event_model_options(priors, event)
  expect_length(options$prior_lists, 1L)
  expect_identical(options$weights, 1)
  expect_identical(options$event_probability, 0)
  expect_equal(options$log_event_probability, 2 * log(1e-200), tolerance = 1e-12)
  expect_identical(options$log_weights, 0)
  context <- .prior_density_conditional_context(priors, c("a", "b"), c("a", "b"))
  expect_identical(context$model_weights, 1)
  expect_identical(context$model_log_weights, 0)
})

test_that("exact zero supplied event options stay impossible", {
  priors <- list(a = prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(0))),
    b = prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(1e-200))))
  options <- .condition_event_model_options(priors, .condition_event(priors, c("a", "b"), "AND"))
  expect_length(options$prior_lists, 0L)
  expect_identical(options$weights, numeric())
  expect_identical(options$event_probability, 0)
  expect_identical(options$log_event_probability, -Inf)
})
