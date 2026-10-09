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

test_that("within-fit branch event logs preserve extreme finite raw weights", {
  cases <- list(maximum = rep(.Machine$double.xmax, 2L), tiny = c(1e200, 1e-200), ordinary = c(1, 1))
  for(w in cases){
    p <- prior_mixture(list(prior("point", list(0), prior_weights = w[1L]),
      prior("normal", list(0, 1), prior_weights = w[2L])), is_null = c(TRUE, FALSE))
    original <- p
    event <- .condition_event(list(theta = p), "theta", "AND")
    result <- .condition_event_model_options(list(theta = p), event)
    reference <- log(w[2L]) - max(log(w)) - log(sum(exp(log(w) - max(log(w)))))
    expect_equal(result$log_event_probability, reference, tolerance = 1e-12)
    expect_equal(result$event_probability, exp(reference), tolerance = 1e-15)
    expect_length(result$prior_lists, 1L)
    expect_identical(unname(result$weights), 1)
    expect_identical(unname(result$log_weights), 0)
    .model_probability_validate(result$weights, result$log_weights, result$model_probability_declaration, normalized = TRUE)
    context <- .prior_density_conditional_context(list(theta = p), "theta", "theta")
    expect_identical(context$model_weights, 1)
    expect_identical(p, original)
    options <- .prior_density_condition_component(p)
    expect_true(is.finite(options[[2L]]$log_probability))
    bias <- prior_mixture(list(prior_PET("normal", list(0, 1), prior_weights = w[1L]),
      prior_PEESE("normal", list(0, 1), prior_weights = w[2L])))
    bias_result <- .condition_event_model_options(list(bias = bias), .condition_event(list(bias = bias), "PEESE", "AND"))
    expect_equal(bias_result$log_event_probability, reference, tolerance = 1e-12)
    expect_identical(bias_result$weights, 1)
  }
  rare <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(1e-200)))
  result <- .condition_event_model_options(list(a = rare, b = rare), .condition_event(list(a = rare, b = rare), c("a", "b"), "AND"))
  expect_equal(result$log_event_probability, 2 * log(1e-200), tolerance = 1e-12)
  expect_identical(result$weights, 1)
  expect_identical(result$event_probability, 0)
  zero <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(0)))
  result <- .condition_event_model_options(list(theta = zero), .condition_event(list(theta = zero), "theta", "AND"))
  expect_identical(result$log_event_probability, -Inf)
  expect_length(result$prior_lists, 0L)
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
