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

test_that("bias condition labels treat an expanded prior list as one branch each", {

  prior_list <- list(
    prior_none(),
    prior_PET("normal", list(0, 1)),
    prior_PEESE("normal", list(0, 1))
  )

  values <- BayesTools:::.condition_event_bias_label_values(
    prior_list,
    c("PET", "PEESE")
  )

  expect_equal(nrow(values), 3L)
  expect_equal(as.logical(values[, "PET"]), c(FALSE, TRUE, FALSE))
  expect_equal(as.logical(values[, "PEESE"]), c(FALSE, FALSE, TRUE))

  event <- list(
    conditional = "PET",
    families = list(bias = list(name = "bias", type = "bias", labels = "PET")),
    conditional_rule = "AND"
  )
  conditioned <- BayesTools:::.bias_prior_list_for_condition(prior_list, event)
  expect_equal(
    vapply(conditioned, BayesTools:::.prior_model_weight, numeric(1)),
    c(0, 1, 0)
  )
})

test_that("plot_models overall diamond uses the requested factor column", {

  skip_if_not_installed("ggplot2")

  make_summary <- function(parameter, mean, lo, hi) {
    out <- data.frame(
      Mean = mean,
      "0.025" = lo,
      "0.975" = hi,
      check.names = FALSE
    )
    class(out) <- c("BayesTools_table", "BayesTools_runjags_summary", "data.frame")
    attr(out, "parameters") <- parameter
    out
  }

  samples <- structure(
    cbind(
      `fac[A]` = c(0, 1, 2, 3),
      `fac[B]` = c(10, 20, 30, 40)
    ),
    class = c("mixed_posteriors.factor", "mixed_posteriors"),
    levels = 2,
    orthonormal = FALSE,
    meandif = FALSE,
    prior_list = list(
      prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  colnames(samples) <- c("fac[A]", "fac[B]")

  model <- list(
    fit = structure(list(), prior_list = list(fac = attr(samples, "prior_list")[[1]])),
    fit_summary = rbind(
      make_summary("fac[A]", 1.5, 0, 3),
      make_summary("fac[B]", 25, 10, 40)
    ),
    inference = list(prior_prob = 1, post_prob = 1, inclusion_BF = 1)
  )
  attr(model$fit_summary, "parameters") <- c("fac[A]", "fac[B]")

  plot_a <- BayesTools:::.plot_models.simple(
    models_summary = list(model$fit_summary),
    models_inference = list(model$inference),
    total_inference = structure(list(), is_null = FALSE),
    total_samples = samples,
    prior_list = attr(samples, "prior_list"),
    parameter = "fac[A]",
    par_name = "fac[A]",
    plot_type = "ggplot",
    prior = FALSE,
    conditional = FALSE,
    order = NULL,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )
  plot_b <- BayesTools:::.plot_models.simple(
    models_summary = list(model$fit_summary),
    models_inference = list(model$inference),
    total_inference = structure(list(), is_null = FALSE),
    total_samples = samples,
    prior_list = attr(samples, "prior_list"),
    parameter = "fac[B]",
    par_name = "fac[B]",
    plot_type = "ggplot",
    prior = FALSE,
    conditional = FALSE,
    order = NULL,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  built_a <- ggplot2::ggplot_build(plot_a)
  built_b <- ggplot2::ggplot_build(plot_b)
  overall_a <- built_a$data[[length(built_a$data) - 1L]]
  overall_b <- built_b$data[[length(built_b$data) - 1L]]

  expect_equal(
    sort(unique(overall_a$x)),
    sort(unique(c(
      unname(stats::quantile(samples[, "fac[A]"], .025)),
      mean(samples[, "fac[A]"]),
      unname(stats::quantile(samples[, "fac[A]"], .975))
    ))),
    tolerance = 1e-12
  )
  expect_equal(
    sort(unique(overall_b$x)),
    sort(unique(c(
      unname(stats::quantile(samples[, "fac[B]"], .025)),
      mean(samples[, "fac[B]"]),
      unname(stats::quantile(samples[, "fac[B]"], .975))
    ))),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(sort(unique(overall_a$x)), sort(unique(overall_b$x)))))
})
