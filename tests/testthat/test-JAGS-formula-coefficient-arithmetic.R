skip_if_not_test_profile("unit")

.formula_ratio_test_fit <- function(numerator = 1e-300, denominator = 2e-300,
                                    x_sd = 1e100, dynamic = FALSE,
                                    dynamic_denominator = TRUE,
                                    values = c(-1, 1, 2)){

  priors <- list(intercept = prior("point", list(0)), x = prior("point", list(0)),
    y = prior("point", list(0)), `x:y` = prior("normal", list(0, 1)))
  attr(priors$x, "multiply_by") <- if(dynamic && dynamic_denominator) "v" else denominator
  attr(priors$`x:y`, "multiply_by") <- if(dynamic) "u" else numerator
  compiled <- JAGS_formula(~ 1 + x * y, "mu",
    data.frame(x = c(-1, 0, 1) * x_sd, y = c(0, 1, 2)), priors, formula_scale = TRUE)
  columns <- names(compiled$prior_list)
  interaction <- columns[grepl("__xXx__", columns, fixed = TRUE)]
  draws <- matrix(0, length(values), length(columns), dimnames = list(NULL, columns))
  draws[, interaction] <- values
  if(dynamic){
    draws <- cbind(draws, u = rep(numerator, length(values)), v = rep(denominator, length(values)))
    compiled$prior_list <- c(compiled$prior_list,
      list(u = prior("normal", list(0, 1)), v = prior("normal", list(0, 1))))
  }
  .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
}

test_that("original static ratios preserve representable nonzero coefficients and laws", {

  fit <- .formula_ratio_test_fit()
  transform <- JAGS_formula_coefficient_transform(fit, "mu")
  interaction <- transform$source_names[grepl("__xXx__", transform$source_names, fixed = TRUE)]
  expected <- transform$basis_matrix["mu_x", interaction] * .5
  expect_identical(unname(expected), -5e-101)
  expect_identical(unname(transform$matrix["mu_x", interaction]), unname(expected))
  expect_identical(unname(transform$prior_recipes$mu_x$weights[interaction]), unname(expected))
  expect_true(interaction %in% transform$dependencies$source[transform$dependencies$target == "mu_x"])
  expect_identical(unname(transform_scale_samples(fit)[, "mu_x"]), c(-1, 1, 2) * unname(expected))
  law <- JAGS_formula_prior_density(fit, "mu", target = "mu_x")
  ordinate <- prior_density_ordinate(law, 0)
  expect_identical(ordinate$behavior, "regular")
  expect_identical(ordinate$point_mass, 0)
  expect_equal(ordinate$log_density, dnorm(0, 0, abs(expected), log = TRUE), tolerance = 1e-14)
  mixed <- as_mixed_posteriors(fit, "mu_x", transform_scaled = TRUE)$mu_x
  expect_length(.posterior_atoms_get(mixed)$mass, 0L)
  expect_null(posterior_metadata(mixed, "measure_unavailable"))
})

test_that("dynamic ratios use bounded alternate grouping and retain zero denominator checks", {

  for(x_sd in c(1e100, 1e-100)){
    fit <- .formula_ratio_test_fit(dynamic = TRUE, x_sd = x_sd,
      numerator = if(x_sd > 1) 1e-300 else 1e300,
      denominator = if(x_sd > 1) 2e-300 else 2e300)
    descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
    expect_true(all(is.na(descriptor$matrix["mu_x", ])))
    interaction <- descriptor$source_names[grepl("__xXx__", descriptor$source_names, fixed = TRUE)]
    expected <- c(-1, 1, 2) * .5 * unname(descriptor$basis_matrix["mu_x", interaction])
    expect_identical(unname(transform_scale_samples(fit)[, "mu_x"]), expected)
  }
  fit <- .formula_ratio_test_fit(dynamic = TRUE, values = c(0, 1))
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  raw <- as.matrix(fit)
  raw[, "v"] <- c(0, 2e-300)
  error <- tryCatch(.bt_apply_formula_coefficient_transform(raw, descriptor, "mu_x"), error = identity)
  expect_s3_class(error, "BayesTools_formula_transform_unavailable")
  expect_identical(error$reason, "zero_multiplier_denominator")
  expect_identical(unname(.bt_apply_formula_coefficient_transform(as.matrix(fit), descriptor, "mu_x")[, "mu_x"]),
    c(0, -5e-101))
})

test_that("static overflow grouping reverses and unresolved ratios refuse", {

  fit <- .formula_ratio_test_fit(numerator = 1e300, denominator = 2e300, x_sd = 1e-100)
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  interaction <- descriptor$source_names[grepl("__xXx__", descriptor$source_names, fixed = TRUE)]
  expect_identical(unname(descriptor$matrix["mu_x", interaction]),
    unname(descriptor$basis_matrix["mu_x", interaction] * .5))
  for(dynamic in c(FALSE, TRUE)){
    request <- function(){
      fit <- .formula_ratio_test_fit(numerator = 1e-300, denominator = 1, dynamic = dynamic)
      if(dynamic) transform_scale_samples(fit) else JAGS_formula_coefficient_transform(fit, "mu")
    }
    error <- tryCatch(request(), error = identity)
    expect_s3_class(error, "BayesTools_formula_transform_unavailable")
    expect_identical(error$reason, "nonfinite_transform")
    expect_identical(error$target, "mu_x")
    expect_true(nzchar(error$source))
    expect_true(length(error$observed) > 0L)
  }
  fit <- .formula_ratio_test_fit(numerator = 1e-300, denominator = 1e300,
    x_sd = 1e-100, dynamic = TRUE, values = 1e300)
  expect_error(transform_scale_samples(fit), class = "BayesTools_formula_transform_unavailable")
})

test_that("ordinary coefficient arithmetic retains multiplication-first bytes", {

  multipliers <- list(x = list(type = "constant", value = 7),
    y = list(type = "constant", value = 3))
  for(basis in c(.1, pi, 1e-200)){
    term <- .bt_formula_coefficient_term("x", "y", basis, multipliers, numeric(),
      c(x = "identity", y = "identity"))
    expect_identical(term$weight, basis * 7 / 3)
  }
  multipliers$x$value <- multipliers$y$value <- 0
  term <- .bt_formula_coefficient_term("x", "y", .1, multipliers, numeric(),
    c(x = "identity", y = "identity"))
  expect_identical(term$weight, .1)
  multipliers$x$value <- 1e-100
  multipliers$y$value <- 1e20
  term <- .bt_formula_coefficient_term("x", "y", 1e-200, multipliers, numeric(),
    c(x = "identity", y = "identity"))
  expect_identical(term$weight, 1e-200 * 1e-100 / 1e20)
  expect_true(term$weight != 0 && abs(term$weight) < .Machine$double.xmin)
  fit <- .formula_ratio_test_fit(numerator = 7, denominator = 3, x_sd = 10,
    dynamic = TRUE, values = c(.1, pi, 1e-100))
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  interaction <- descriptor$source_names[grepl("__xXx__", descriptor$source_names, fixed = TRUE)]
  expect_identical(unname(transform_scale_samples(fit)[, "mu_x"]),
    unname(descriptor$basis_matrix["mu_x", interaction]) * c(.1, pi, 1e-100) * 7 / 3)
})

test_that("unrepresentable dynamic recipe scales retain available per-draw values", {

  fit <- .formula_ratio_test_fit(numerator = 1e200, denominator = 1e300,
    dynamic = TRUE, dynamic_denominator = FALSE, values = c(-1, 1) * 1e200)
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  expect_true(all(is.na(descriptor$matrix["mu_x", ])))
  expect_identical(descriptor$targets$structural_status[descriptor$targets$target == "mu_x"], "dependent")
  expect_identical(descriptor$prior_recipes$mu_x$type, "unavailable")
  expect_identical(descriptor$prior_recipes$mu_x$reason, "numerical_scale_unavailable")
  expect_length(descriptor$prior_recipes$mu_x$weights, 0L)
  expect_identical(unname(transform_scale_samples(fit)[, "mu_x"]), c(1, -1))
  error <- tryCatch(JAGS_formula_prior_density(fit, "mu", target = "mu_x"), error = identity)
  expect_s3_class(error, "BayesTools_formula_measure_unavailable")
  expect_identical(error$reason, "numerical_scale_unavailable")
  leaf <- as_mixed_posteriors(fit, "mu_x", transform_scaled = TRUE)$mu_x
  expect_identical(as.numeric(leaf), c(1, -1))
  unavailable <- posterior_metadata(leaf, "measure_unavailable")
  expect_identical(sort(unavailable$measure), sort(c("prior_density", "atoms", "support")))
  expect_identical(unique(unavailable$reason), "numerical_scale_unavailable")
})

test_that("weighted logged dynamic outputs refuse while direct and identity outputs remain available", {

  make_fit <- function(mean = 20, dynamic = TRUE){
    slope <- prior("point", list(2))
    if(dynamic) attr(slope, "multiply_by") <- "sigma"
    formula <- ~ 1 + x
    attr(formula, "log(intercept)") <- TRUE
    compiled <- JAGS_formula(formula, "mu", data.frame(x = mean + c(-10, 0, 10)),
      list(intercept = prior("point", list(exp(5))), x = slope), formula_scale = TRUE)
    draws <- cbind(mu_intercept = rep(exp(5), 2), mu_x = rep(2, 2), sigma = c(-1, 1))
    priors <- c(compiled$prior_list, list(sigma = prior("normal", list(0, 1))))
    .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), priors,
      list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
  }
  fit <- make_fit()
  transform <- JAGS_formula_coefficient_transform(fit, "mu")
  expect_identical(transform$source_transforms[["mu_intercept"]], "log")
  expect_identical(transform$output_transforms[["mu_intercept"]], "exp")
  direct <- JAGS_formula_prior_density(fit, "mu", target = "mu_intercept")
  expect_identical(prior_density_ordinate(direct, -1)$behavior, "zero")
  expect_identical(prior_density_ordinate(direct, 0)$behavior, "zero")
  expect_identical(prior_density_ordinate(direct, 1)$behavior, "regular")
  error <- tryCatch(JAGS_formula_prior_density(fit, "mu", weights = c(mu_intercept = 1)), error = identity)
  expect_s3_class(error, "BayesTools_formula_prior_density_unavailable")
  expect_identical(error$reason, "nonlinear_map")
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  context <- posterior_metadata(mixed, "prior_context")
  list_error <- tryCatch(.prior_density_from_context(context, c(mu_intercept = 1)), error = identity)
  expect_s3_class(list_error, "BayesTools_formula_measure_unavailable")
  expect_identical(list_error$reason, "nonlinear_map")
  expect_s3_class(JAGS_formula_prior_density(fit, "mu", weights = c(mu_x = 1, mu_intercept = 0)), "prior_density")
  raw_context <- .prior_density_build_context(attr(fit, "prior_list", exact = TRUE), c("mu_intercept", "mu_x"))
  expect_identical(prior_density_ordinate(.prior_density_from_context(raw_context, c(mu_intercept = 1)), exp(5))$point_mass, 1)
  static <- make_fit(mean = 0, dynamic = FALSE)
  expect_identical(JAGS_formula_coefficient_transform(static, "mu")$targets$map_type[[1L]], "identity")
  expect_identical(prior_density_ordinate(JAGS_formula_prior_density(static, "mu", weights = c(mu_intercept = 1)), exp(5))$point_mass, 1)
  mixed <- as_mixed_posteriors(static, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  context <- posterior_metadata(mixed, "prior_context")
  expect_identical(prior_density_ordinate(.prior_density_from_context(context, c(mu_intercept = 1)), exp(5))$point_mass, 1)
})

test_that("mixed coefficient recipes refuse unrepresentable contribution conversion", {

  x_prior <- prior("normal", list(0, 1e100))
  y_prior <- prior("point", list(2))
  attr(x_prior, "multiply_by") <- 1e300
  attr(y_prior, "multiply_by") <- "sigma"
  compiled <- JAGS_formula(~ 1 + x + y, "mu",
    data.frame(x = c(-1, 0, 1) * 1e100, y = c(10, 20, 30)),
    list(intercept = prior("point", list(5)), x = x_prior, y = y_prior), formula_scale = TRUE)
  priors <- c(compiled$prior_list, list(sigma = prior("normal", list(0, 1))))
  draws <- cbind(mu_intercept = rep(5, 2), mu_x = c(-1, 1) * 1e100,
    mu_y = rep(2, 2), sigma = c(-1, 1))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), priors,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
  expect_identical(unname(transform_scale_samples(fit)[, "mu_x"]), c(-1, 1))
  expect_s3_class(JAGS_formula_prior_density(fit, "mu", target = "mu_x"), "prior_density")
  expect_s3_class(JAGS_formula_prior_density(fit, "mu", target = "mu_intercept"), "prior_density")
  error <- tryCatch(JAGS_formula_prior_density(fit, "mu", weights = c(mu_intercept = 1, mu_x = 1)), error = identity)
  expect_s3_class(error, "BayesTools_formula_measure_unavailable")
  expect_identical(error$reason, "numerical_scale_unavailable")
})

test_that("static subnormal ratio weights use value-aware coefficient arithmetic", {

  for(numerator in c(1e-200, -1e-200)){
    observed <- list()
    for(dynamic in c(FALSE, TRUE)){
      fit <- .formula_ratio_test_fit(numerator = numerator, denominator = 1e23,
        x_sd = 1e100, dynamic = dynamic, values = c(-1, 0, 1) * 1e200)
      transform <- JAGS_formula_coefficient_transform(fit, "mu")
      interaction <- transform$source_names[grepl("__xXx__", transform$source_names, fixed = TRUE)]
      expected <- transform$basis_matrix["mu_x", interaction] * c(-1, 0, 1) * 1e200 * numerator / 1e23
      value <- .bt_apply_formula_coefficient_transform(as.matrix(fit), transform, "mu_x")[, "mu_x"]
      expect_identical(unname(value[2L]), 0)
      expect_equal(unname(value[c(1L, 3L)] / expected[c(1L, 3L)]), c(1, 1), tolerance = 1e-14)
      observed[[length(observed) + 1L]] <- unname(value)
      if(!dynamic){
        expect_identical(transform$targets$map_type[transform$targets$target == "mu_x"], "affine")
        expect_error(JAGS_formula_prior_density(fit, "mu", target = "mu_x"), class = "BayesTools_formula_measure_unavailable")
      }
    }
    expect_identical(observed[[1L]], observed[[2L]])
  }
  fit <- .formula_ratio_test_fit(numerator = 1e-200, denominator = 1e23, values = 1)
  transform <- JAGS_formula_coefficient_transform(fit, "mu")
  value <- .bt_apply_formula_coefficient_transform(as.matrix(fit), transform, "mu_x")[, "mu_x"]
  expect_true(value != 0 && abs(value) < .Machine$double.xmin)
  fit <- .formula_ratio_test_fit(numerator = 1e-200, denominator = 1e23, values = 1e-200)
  expect_error(.bt_apply_formula_coefficient_transform(as.matrix(fit),
    JAGS_formula_coefficient_transform(fit, "mu"), "mu_x"), class = "BayesTools_formula_transform_unavailable")
})

test_that("compiled finite-density and fixed ratio owners retain truthful values and measures", {

  for(point in c(FALSE, TRUE)){
    priors <- list(intercept = prior("point", list(0)), x = prior("point", list(0)),
      y = prior("point", list(0)), `x:y` = if(point) prior("point", list(1e200)) else prior("lognormal", list(log(1e200), 1)))
    attr(priors$x, "multiply_by") <- 1e23
    attr(priors$`x:y`, "multiply_by") <- 1e-200
    compiled <- JAGS_formula(~ 1 + x * y, "mu", data.frame(x = c(-1, 0, 1) * 1e100, y = c(0, 1, 2)), priors, formula_scale = TRUE)
    columns <- names(compiled$prior_list)
    interaction <- columns[grepl("__xXx__", columns, fixed = TRUE)]
    draws <- matrix(0, 3L, length(columns), dimnames = list(NULL, columns))
    draws[, interaction] <- 1e200
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), compiled$prior_list,
      list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
    transform <- JAGS_formula_coefficient_transform(fit, "mu")
    expected <- unname(transform$basis_matrix["mu_x", interaction]) * 1e200 * 1e-200 / 1e23
    value <- .bt_apply_formula_coefficient_transform(as.matrix(fit), transform, "mu_x")[, "mu_x"]
    expect_equal(unname(value / expected), rep(1, 3L), tolerance = 1e-14)
    if(point){
      target <- transform$targets[transform$targets$target == "mu_x", ]
      expect_identical(target$structural_status, "structural")
      expect_equal(target$fixed_value / expected, 1, tolerance = 1e-14)
      expect_identical(prior_density_ordinate(JAGS_formula_prior_density(fit, "mu", target = "mu_x"), expected)$point_mass, 1)
      expect_error(JAGS_formula_prior_density(fit, "mu", target = "mu_x", context = list()),
        class = "BayesTools_formula_prior_density_unavailable")
    }else{
      expect_true(is.finite(lpdf(compiled$prior_list[[interaction]], 1e200)))
      condition <- tryCatch(JAGS_formula_prior_density(fit, "mu", target = "mu_x"), error = identity)
      expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
      expect_identical(condition$reason, "numerical_scale_unavailable")
      leaf <- as_mixed_posteriors(fit, "mu_x", transform_scaled = TRUE)$mu_x
      expect_equal(as.numeric(leaf) / expected, rep(1, 3L), tolerance = 1e-14)
      expect_identical(unique(posterior_metadata(leaf, "measure_unavailable")$reason), "numerical_scale_unavailable")
      expect_error(.bt_formula_measure_check(leaf, "prior_density"), class = "BayesTools_formula_measure_unavailable")
    }
    expect_s3_class(JAGS_formula_prior_density(fit, "mu", target = "mu_y"), "prior_density")
  }
})


test_that("weighted static recipes refuse lost nonzero products and retain cancellation", {

  info <- JAGS_formula(~ x, "mu", data.frame(x = c(-2, 0, 2)),
    list(intercept = prior("point", list(0)), x = prior("normal", list(0, 1))), formula_scale = TRUE)
  draws <- cbind(mu_intercept = 0, mu_x = c(-1, 0, 1))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), info$prior_list,
    list(mu = info$formula_design), list(mu = info$formula_scale))
  tiny <- .Machine$double.xmin * .Machine$double.eps
  condition <- expect_error(JAGS_formula_prior_density(fit, "mu", weights = c(mu_x = tiny)),
    class = "BayesTools_formula_measure_unavailable")
  expect_identical(condition$reason, "numerical_scale_unavailable")
  expect_null(conditionCall(condition))
  ordinary <- JAGS_formula_prior_density(fit, "mu", weights = c(mu_x = .5))
  expect_equal(prior_density_ordinate(ordinary, 0)$log_density,
    dnorm(0, sd = .25, log = TRUE), tolerance = 1e-14)
  zero <- JAGS_formula_prior_density(fit, "mu", weights = c(mu_intercept = 1))
  expect_identical(prior_density_ordinate(zero, 0)$point_mass, 1)
  shifted_info <- JAGS_formula(~ x, "mu", data.frame(x = c(0, 2, 4)),
    list(intercept = prior("point", list(0)), x = prior("normal", list(0, 1))), formula_scale = TRUE)
  shifted_fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), shifted_info$prior_list,
    list(mu = shifted_info$formula_design), list(mu = shifted_info$formula_scale))
  cancelled <- JAGS_formula_prior_density(shifted_fit, "mu", weights = c(mu_intercept = 1, mu_x = 2))
  expect_identical(prior_density_ordinate(cancelled, 0)$point_mass, 1)
})

test_that("ordinary fixed and draw products never certify lost nonzero values", {

  tiny <- .Machine$double.xmin * .Machine$double.eps
  make <- function(slope, values){
    info <- JAGS_formula(~ x, "mu", data.frame(x = c(-2, 0, 2)),
      list(intercept = prior("point", list(0)), x = slope), formula_scale = TRUE)
    draws <- cbind(mu_intercept = 0, mu_x = values)
    .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), info$prior_list,
      list(mu = info$formula_design), list(mu = info$formula_scale))
  }
  fixed <- make(prior("point", list(tiny)), rep(tiny, 3L))
  descriptor <- JAGS_formula_coefficient_transform(fixed, "mu")
  target <- descriptor$targets[descriptor$targets$target == "mu_x", , drop = FALSE]
  expect_identical(target$structural_status, "unavailable")
  expect_true(is.na(target$fixed_value))
  expect_identical(as.numeric(as_mixed_posteriors(fixed, "mu_x", transform_scaled = FALSE)$mu_x), rep(tiny, 3L))
  normal <- make(prior("normal", list(0, 1)), c(tiny, 2 * tiny, tiny))
  expect_error(as_mixed_posteriors(normal, "mu_x", transform_scaled = TRUE),
    class = "BayesTools_formula_transform_unavailable")
  expect_identical(as.numeric(as_mixed_posteriors(normal, "mu_x", transform_scaled = FALSE)$mu_x), c(tiny, 2 * tiny, tiny))
  ordinary <- make(prior("normal", list(0, 1)), c(-1, 0, 1))
  expect_identical(as.numeric(as_mixed_posteriors(ordinary, "mu_x", transform_scaled = TRUE)$mu_x), c(-.5, 0, .5))
  zero <- make(prior("point", list(0)), rep(0, 3L))
  expect_identical(as.numeric(as_mixed_posteriors(zero, "mu_x", transform_scaled = TRUE)$mu_x), rep(0, 3L))
})
