skip_if_not_test_profile("unit")

test_that("named BF slices use the visible vector positions for every carrier", {

  values <- .format_BF_from_log(c(a = 710, b = 2, a = 750),
    bound_operator = c(">", NA, ">"), diagnostics = list(list(id = 1), list(id = 2), list(id = 3)))
  attr(values, "log_prior_height") <- c(10, 20, 30)
  attr(values, "log_posterior_height") <- c(720, 22, 780)
  for(index in list(c("b", "a", "a", "unknown"), c(3L, 1L, 1L, NA_integer_),
                   c(TRUE, FALSE, TRUE), -2L)){
    positions <- stats::setNames(seq_along(values), names(values))[index]
    selected <- values[index]
    expect_identical(as.numeric(selected), as.numeric(values)[positions])
    expect_identical(.BF_carrier_log(selected), .BF_carrier_log(values)[positions])
    expect_identical(attr(selected, "bound_operator"), attr(values, "bound_operator")[positions])
    expect_identical(attr(selected, "numerical_diagnostics"), attr(values, "numerical_diagnostics")[positions])
    expect_identical(attr(selected, "log_prior_height"), c(10, 20, 30)[positions])
    expect_identical(attr(selected, "log_posterior_height"), c(720, 22, 780)[positions])
  }
  selected <- values[c("b", "a", "a", "unknown")]
  table <- data.frame(parameter = seq_along(selected), BF = as.numeric(selected))
  table$BF <- selected
  class(table) <- c("BayesTools_table", "data.frame")
  attr(table, "type") <- c("string", "BF")
  expect_identical(as.numeric(update(table, logBF = TRUE)$BF), c(2, 710, 710, NA_real_))
  expect_identical(as.numeric(update(table, BF01 = TRUE)$BF), exp(-c(2, 710, 710, NA_real_)))
  expect_identical(.BF_carrier_log(values[]), .BF_carrier_log(values))
})

test_that("point inference retains finite log ordinates and orientation", {

  arguments <- list(posterior = c(1, 2, 3), prior = prior("normal", list(0, 1)),
                    parameter = "theta", density_method = "normal",
                    logBF = TRUE, columns = "all")
  implicit <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 40")))
  explicit <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 40 vs theta != 40")))
  reverse <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta != 40 vs theta = 40")))
  expect_equal(c(implicit$BF, explicit$BF, reverse$BF), c(-78, 78, -78), tolerance = 1e-12)
  expect_equal(attr(explicit, "raw_log_BF"), 78, tolerance = 1e-12)
  expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 1e100"))),
               class = "BayesTools_hypothesis_numerical_unavailable")
})

test_that("nonfinite Normal arithmetic refuses on raw and declared point readers", {

  expect_error(hypothesis_BF(c(1, 2, 3), prior("normal", list(0, 1e200)),
    "theta = 1e200", parameter = "theta", density_method = "normal"),
    class = "BayesTools_hypothesis_numerical_unavailable")
  expect_error(hypothesis_BF(c(-1e300, 0, 1e300), prior("normal", list(0, 1)),
    "theta = 0", parameter = "theta", density_method = "normal"),
    class = "BayesTools_hypothesis_numerical_unavailable")
  law <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))), c(theta = 1), n_grid = 128)
  make_posterior <- function(values){
    class(values) <- c("marginal_posterior.simple", "marginal_posterior", "numeric")
    .bt_meta_assign(values, list(prior_density = law, atoms = posterior_atom_attribute()))
  }
  levels <- list(unavailable = make_posterior(c(-1e300, 0, 1e300)), valid = make_posterior(c(-1, 0, 1)))
  class(levels) <- c("marginal_posterior.factor", "marginal_posterior", "list")
  result <- Savage_Dickey_BF(levels, normal_approximation = TRUE, silent = TRUE)
  expect_true(is.na(result$unavailable))
  expect_equal(as.numeric(result$valid), 1, tolerance = 1e-14)
  expect_s3_class(attr(result$unavailable, "numerical_diagnostics"), "BayesTools_hypothesis_numerical_unavailable")
})

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

test_that("affine Boolean regions use exact declared probabilities", {

  result <- hypothesis_BF(posterior = c(-2, -1, 1, 2),
    prior = prior("normal", list(0, 1)), parameter = "theta",
    hypothesis = "(theta + 1e100 > 1e100) & (theta / 49 <= 1)",
    seed = 337, columns = "all")
  expect_equal(result$prior, (stats::pnorm(49) - .5) / (1 - (stats::pnorm(49) - .5)), tolerance = 1e-14)
  expect_true(length(attr(result, "prior_numerical_diagnostics")) == 1L)
  expect_error(hypothesis_BF(c(-2, -1, 1, 2), prior("normal", list(0, 1)),
    "theta * 1e-200 * 1e-200 > 0", parameter = "theta"),
    class = "BayesTools_hypothesis_region_numerical_unavailable")
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

test_that("live BF carriers preserve direction and invalidate changed values", {

  values <- .format_BF_from_log(c(710, 750, 800), BF = rep(Inf, 3))
  inverse <- .format_BF_from_log(c(710, 750, 800), BF01 = TRUE, BF = rep(Inf, 3))
  expect_equal(as.numeric(inverse)[1L] / exp(-710), 1, tolerance = 1e-14)
  expect_equal(attr(values[c(3, 1, 1)], "canonical_log_BF"), c(800, 710, 710))
  expect_null(attr(values * 2, "canonical_log_BF"))
  expect_null(attr(abs(values), "canonical_log_BF"))
  changed <- values
  changed[] <- rep(Inf, 3)
  expect_null(attr(changed, "canonical_log_BF"))
})

test_that("declared atoms share the normalized region and Boolean event", {

  law <- .prior_linear_combination_density(list(theta = prior_mixture(list(
    prior("point", list(49)), prior("normal", list(49, 1))), is_null = c(TRUE, FALSE))),
    c(theta = 1), n_grid = 1024)
  forms <- c("theta / 49 >= 1", "theta / 49 > 1", "theta / -49 >= -1", "theta / -49 > -1",
             "!(theta / 49 < 1)", "(theta / 49 >= 1) | (theta / 49 > 1)")
  forms <- c(forms, "theta * (1 / 49) >= 1", "theta * (1 / 49) > 1")
  masses <- vapply(forms, function(text){
    .hypothesis_prior_density_prob(law, hypothesis_parse(text)$statements[[1L]]$left, "theta")
  }, numeric(1))
  expect_equal(unname(masses), c(.75, .25, .75, .25, .75, .75, .75, .25), tolerance = 1e-14)
})

test_that("supported deterministic masses and sampled fallback retain their variance", {

  law <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))), c(theta = 1), n_grid = 128)
  quantity <- .hypothesis_quantity_from_draws(data.frame(theta = c(-2, -1, 1, 2)),
    data.frame(theta = c(-3, -2, -1, 1)), "theta", parameter = "theta", prior_density = law)
  scalar <- hypothesis_parse("theta > 0")$statements[[1L]]$left
  boolean <- hypothesis_parse("theta > -1 & theta < 1")$statements[[1L]]$left
  sampled <- hypothesis_parse("abs(theta) < 2")$statements[[1L]]$left
  masses <- lapply(list(scalar, boolean, sampled), function(side) .hypothesis_region_mass(quantity, side, TRUE))
  expect_equal(as.numeric(masses[[1L]]), .5, tolerance = 1e-14)
  expect_equal(as.numeric(masses[[2L]]), stats::pnorm(1) - stats::pnorm(-1), tolerance = 1e-14)
  expect_identical(attr(masses[[3L]], "route"), "sampled")
  scalar_arguments <- list(quantity, scalar, TRUE)
  sampled_arguments <- list(quantity, sampled, TRUE)
  if("mass" %in% names(formals(.hypothesis_region_log_mass_mc_var))){
    scalar_arguments$mass <- masses[[1L]]
    sampled_arguments$mass <- masses[[3L]]
  }
  expect_identical(do.call(.hypothesis_region_log_mass_mc_var, scalar_arguments), 0)
  expect_equal(do.call(.hypothesis_region_log_mass_mc_var, sampled_arguments),
    stats::var(c(0, 0, 1, 1)) / (4 * .5^2), tolerance = 1e-14)
})

test_that("log-bound tables bind row identity and never revive replaced carriers", {

  result <- hypothesis_BF(c(1, 2, 3), prior("normal", list(0, 1)),
    c("theta = 40", "theta = 41"), parameter = "theta", density_method = "normal")
  selected <- result[c(2, 1, 1), , drop = FALSE]
  expect_equal(as.numeric(update(selected, logBF = TRUE)$BF), c(-80, -78, -78), tolerance = 1e-12)
  expect_equal(attr(selected, "raw_log_BF"), c(-80, -78, -78), tolerance = 1e-12)
  replaced <- result
  replaced$BF <- c(Inf, Inf)
  expect_true(all(is.infinite(as.numeric(update(replaced, logBF = TRUE)$BF))))
  replaced <- result
  replaced$BF[] <- c(0, 0)
  expect_null(attr(replaced$BF, "canonical_log_BF"))
  expect_true(all(is.infinite(as.numeric(update(replaced, logBF = TRUE)$BF))))
  for(exported in list(as.data.frame(selected), data.frame(selected))){
    expect_null(attr(exported$BF10, "canonical_log_BF"))
    expect_identical(class(exported), "data.frame")
  }
  expect_null(attr(result[, "Null", drop = FALSE], "raw_log_BF"))
})

test_that("huge finite prior logs retain eligibility with a genuine log ordinate", {

  gaussian <- prior("normal", list(0, 1))
  attr(gaussian, "multiply_by") <- "sigma"
  law <- .prior_linear_combination_density(list(theta = gaussian, sigma = prior("lognormal", list(0, 40))),
    c(theta = 1), n_grid = 128)
  log_height <- prior_density_ordinate(law, 0)$log_density
  expect_equal(log_height, 799.0810614667953, tolerance = 2e-13)
  posterior <- c(-1, 0, 1)
  class(posterior) <- c("marginal_posterior.simple", "marginal_posterior")
  posterior <- .bt_meta_assign(posterior, list(prior_density = law, atoms = posterior_atom_attribute()))
  expect_equal(attr(Savage_Dickey_BF(posterior, normal_approximation = TRUE, silent = TRUE), "log_BF"),
    800, tolerance = 2e-13)
  unit <- .posterior_log_ordinate_attribute(0, log_height, "unit point part", "precomputed")
  expect_equal(.posterior_ordinate_from_attribute(unit, 0)$log_y, log_height)
  expect_error(.hypothesis_log_ratio(-5e199, -5e199, normal = TRUE),
    class = "BayesTools_hypothesis_numerical_unavailable")
  expect_equal(.hypothesis_log_ratio(799.0811, 799.0811, normal = TRUE)$log_BF, 0)
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

test_that("nonlocal constructor diagnoses loss of derived tau", {

  for(distribution in c("moment", "invmoment")){
    parameters <- list(mode = 1e-200)
    if(distribution == "invmoment") parameters$df <- 1
    expect_error(prior(distribution, parameters),
      "The supplied 'mode' yields a derived 'tau' outside representable positive range.", fixed = TRUE)
  }
})

test_that("marginal estimate exports remove hidden logs and preserve numeric values", {

  posterior <- c(1, 2, 3)
  class(posterior) <- c("marginal_posterior.simple", "marginal_posterior", "numeric")
  law <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))), c(theta = 1), n_grid = 128)
  posterior <- .bt_meta_assign(posterior, list(prior_density = law, atoms = posterior_atom_attribute()))
  BF <- Savage_Dickey_BF(posterior, null_hypothesis = 40, normal_approximation = TRUE, silent = TRUE)
  result <- marginal_estimates_table(list(theta = posterior), list(theta = BF), "theta", logBF = TRUE)
  expect_equal(as.numeric(result$inclusion_BF), -78, tolerance = 1e-12)
  for(exported in list(as.data.frame(result), data.frame(result))){
    expect_null(attr(exported$inclusion_BF, "canonical_log_BF"))
    expect_null(attr(exported, "raw_log_BF"))
    expect_equal(as.numeric(exported$inclusion_BF), -78, tolerance = 1e-12)
  }
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
