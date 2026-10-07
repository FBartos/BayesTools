skip_if_not_test_profile("unit")

test_that("affine point coordinates avoid collapsed translations", {

  arguments <- list(posterior = c(1, 2, 3), prior = prior("normal", list(0, 1)),
                    parameter = "theta", density_method = "normal",
                    logBF = TRUE, columns = "all")
  baseline <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 0")))
  shifted <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta + 1e100 = 1e100")))
  expect_equal(as.numeric(shifted$BF), as.numeric(baseline$BF), tolerance = 1e-12)
  shifted_diagnostics <- attr(shifted, "numerical_diagnostics")[[1L]]
  expect_equal(shifted_diagnostics$evaluation_value, 0)
  expect_equal(shifted_diagnostics$requested_value, 1e100)
  scaled <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "2 * theta = 100")))
  direct <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 50")))
  expect_equal(as.numeric(scaled$BF), as.numeric(direct$BF), tolerance = 1e-12)
  scaled_diagnostics <- attr(scaled, "numerical_diagnostics")[[1L]]
  direct_diagnostics <- attr(direct, "numerical_diagnostics")[[1L]]
  expect_equal(scaled_diagnostics$log_prior_height,
    direct_diagnostics$log_prior_height - log(2), tolerance = 1e-12)
  expect_equal(scaled_diagnostics$log_posterior_height,
    direct_diagnostics$log_posterior_height - log(2), tolerance = 1e-12)
  divided <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta / 49 = 1")))
  diagnostics <- attr(divided, "numerical_diagnostics")[[1L]]
  expect_equal(diagnostics$log_prior_height, stats::dnorm(49, log = TRUE) + log(49), tolerance = 1e-12)
  expect_equal(diagnostics$log_posterior_height, stats::dnorm(49, 2, 1, log = TRUE) + log(49), tolerance = 1e-12)
  expect_equal(diagnostics$evaluation_value, 49)
  expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta - theta = 0"))),
               class = "BayesTools_point_mass_at_null")
  expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta * 1e-200 * 1e-200 = 0"))),
               class = "BayesTools_hypothesis_numerical_unavailable")
  for(expression in c("theta * (1e-200 * 1e-200) = 0", "theta * exp(-1000) = 0")){
    expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = expression))),
      class = "BayesTools_hypothesis_numerical_unavailable")
  }
  constant_function <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "sqrt(4) * theta = 100")))
  expect_equal(as.numeric(constant_function$BF), as.numeric(direct$BF), tolerance = 1e-12)
})

test_that("constant divisors preserve strict and inclusive atom boundaries", {

  draws <- data.frame(theta = c(48, 49, 50))
  forms <- c("theta / 49 >= 1", "theta / 49 > 1", "theta / -49 >= -1", "theta / -49 > -1")
  expected <- list(c(FALSE, TRUE, TRUE), c(FALSE, FALSE, TRUE),
                   c(TRUE, TRUE, FALSE), c(TRUE, FALSE, FALSE))
  observed <- lapply(forms, function(text){
    side <- hypothesis_parse(text)$statements[[1L]]$left
    .hypothesis_draw_region_indicator(side, draws)
  })
  expect_identical(observed, expected)
  expect_error(.hypothesis_linear_coefficients(parse(text = "theta * 1e-200 * 1e-200")[[1L]],
    "theta", draws), class = "BayesTools_hypothesis_numerical_unavailable")
})

test_that("formula measure classes cover declared laws and propagate missing metadata", {

  for(reason in c("state_dependent_map", "nonlinear_map", "incompatible_prior_recipes",
                 "missing_multiplier_law", "numerical_scale_unavailable")){
    expect_error(.bt_formula_density_stop("Unavailable.", reason = reason),
      class = "BayesTools_formula_measure_unavailable")
  }
  for(reason in c("unknown_target", "missing_source_coordinates", "missing_prior_list", "arbitrary_program_error")){
    condition <- tryCatch(.bt_formula_density_stop("Missing.", reason = reason), error = identity)
    expect_false(inherits(condition, "BayesTools_formula_measure_unavailable"))
  }
})

.inference_component_fixture <- function(multiplier_gate = FALSE, return_fit = FALSE){

  mixture <- prior_mixture(list(prior("uniform", list(if(multiplier_gate) 1 else 0, if(multiplier_gate) 2 else 1)),
    prior("uniform", list(10, 11))), is_null = c(FALSE, FALSE))
  slope <- if(multiplier_gate) prior("point", list(2)) else mixture
  attr(slope, "multiply_by") <- if(multiplier_gate) "sigma" else -3
  formula <- JAGS_formula(~ x, "mu", data.frame(x = c(10, 20, 30)),
    list(intercept = prior("point", list(5)), x = slope), formula_scale = TRUE)
  values <- cbind(mu_intercept = rep(5, 4), mu_x = if(multiplier_gate) rep(2, 4) else c(.2, .8, 10.2, 10.8))
  extra <- list()
  if(multiplier_gate){
    extra <- list(sigma = mixture)
    values <- cbind(values, sigma = c(1.2, 1.8, 10.2, 10.8), sigma_indicator = c(1, 1, 2, 2))
  }else values <- cbind(values, mu_x_indicator = c(1, 1, 2, 2))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(values)),
    c(formula$prior_list, extra), list(mu = formula$formula_design), list(mu = formula$formula_scale))
  if(return_fit) return(fit)
  marginal_posterior(as_mixed_posteriors(fit, c("mu_intercept", "mu_x")),
    "mu_intercept", formula = ~ x, at = list(x = 1), prior_samples = TRUE)$intercept
}

test_that("direct and compiled combinations retain canonical component supports", {

  for(gate in c(FALSE, TRUE)){
    producer <- .inference_component_fixture(gate)
    expected <- if(gate) list(c(7, 9), c(25, 27)) else list(c(2, 5), c(-28, -25))
    components <- .posterior_components_get(producer)
    expect_identical(lapply(components$supports, `[[`, "bounds"), expected)
    parameter <- "mu_intercept"
    quantity <- .as_hypothesis_quantities(producer, NULL,
      hypothesis_parse("mu_intercept + 0 = 2")$statements, parameter)[[1L]]
    direct <- .hypothesis_linear_point_marginal(quantity,
      hypothesis_parse("mu_intercept + 0 = 2")$statements[[1L]]$left)
    expect_identical(.posterior_components_get(direct$posterior)$index, components$index)
    expect_identical(lapply(.posterior_components_get(direct$posterior)$supports, `[[`, "bounds"), expected)
    fixed <- rep(0, length(producer))
    class(fixed) <- c("marginal_posterior.simple", "marginal_posterior", "numeric")
    fixed <- .bt_meta_assign(fixed, list(linear_weights = .bt_meta_get(producer, "linear_weights") * 0,
      atoms = .posterior_atoms_new(locations = matrix(0), mass = 1, column_names = "value"),
      prior_context = .bt_meta_get(producer, "prior_context"),
      condition = .bt_meta_get(producer, "condition")))
    levels <- list(A = fixed, B = producer)
    class(levels) <- c("marginal_posterior.factor", "marginal_posterior", "list")
    attr(levels, "parameter") <- parameter
    compiled <- hypothesis_linear_target(levels, "mu_intercept[B] - mu_intercept[A] = 2", parameter)
    expect_equal(as.numeric(compiled$posterior), as.numeric(producer), tolerance = 1e-14)
    expect_identical(lapply(.posterior_components_get(compiled$posterior)$supports, `[[`, "bounds"), expected)
    expect_identical(.bt_meta_get(compiled$posterior, "linear_weight_space"), "formula_contribution")
    boundary <- if(gate) 7 else 2
    baseline <- hypothesis_BF(producer, hypothesis = paste("mu_intercept =", boundary), parameter = parameter)
    direct_result <- hypothesis_BF(producer, hypothesis = paste("mu_intercept + 0 =", boundary), parameter = parameter)
    target <- hypothesis_linear_target(levels, paste("mu_intercept[B] - mu_intercept[A] =", boundary), parameter)
    compiled_result <- hypothesis_BF(target$posterior, hypothesis = target$hypothesis, parameter = target$parameter)
    expect_equal(as.numeric(direct_result$BF), as.numeric(baseline$BF), tolerance = 1e-13)
    expect_equal(as.numeric(compiled_result$BF), as.numeric(baseline$BF), tolerance = 1e-13)
  }
})

test_that("public missing formula targets and priors propagate outside the measure family", {

  fit <- .inference_component_fixture(return_fit = TRUE)
  missing_target <- tryCatch(JAGS_formula_prior_density(fit, "mu", "missing"), error = identity)
  expect_true(inherits(missing_target, "BayesTools_formula_prior_density_unavailable"))
  expect_false(inherits(missing_target, "BayesTools_formula_measure_unavailable"))
  expect_identical(missing_target$reason, "unknown_target")
  attr(fit, "prior_list") <- NULL
  missing_priors <- tryCatch(JAGS_formula_prior_density(fit, "mu", "mu_x"), error = identity)
  expect_true(inherits(missing_priors, "condition"))
  expect_false(inherits(missing_priors, "BayesTools_formula_measure_unavailable"))
})

test_that("repeated component gates and declared weight spaces require alignment", {

  producer <- .inference_component_fixture()
  context <- .bt_meta_get(producer, "prior_context")
  weights <- .bt_meta_get(producer, "linear_weights")[1L, ]
  components <- .posterior_components_get(producer)
  second <- .posterior_components_set(producer, .posterior_components_new(
    rev(components$index), components$supports, components$keys))
  expect_error(.hypothesis_linear_components(list(producer, second), context,
    weights, list(a = 0, b = 1), length(producer)), class = "BayesTools_linear_target_unavailable")
  second <- .bt_meta_set(producer, "linear_weight_space", "coefficient")
  levels <- list(A = producer, B = second)
  class(levels) <- c("marginal_posterior.factor", "marginal_posterior", "list")
  expect_error(hypothesis_linear_target(levels, "mu_intercept[A] + mu_intercept[B] = 2", "mu_intercept"),
    class = "BayesTools_linear_target_unavailable")
})
