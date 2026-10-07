skip_if_not_test_profile("unit")

test_that("finite PET panels retain tiny physical partitions without extrapolation artifacts", {

  mu <- .petpeese_prior_components(prior("normal", list(.3, .7)))[[1L]]
  bias <- .petpeese_prior_components(prior("normal", list(0, .4)))[[1L]]
  q <- .30000000000001226
  value <- .petpeese_prior_probability_integral(mu, bias, .5, q)
  expect_equal(as.numeric(value), stats::pnorm(q, .3, sqrt(.7^2 + .2^2)), tolerance = 1e-7)
  diagnostics <- attr(value, "numerical_diagnostics")
  expect_lte(diagnostics$subdivisions, 200L)
  expect_identical(diagnostics$subdivisions, sum(vapply(diagnostics$panels, `[[`, integer(1), "subdivisions")))
  intervals <- do.call(rbind, lapply(diagnostics$panels, `[[`, "physical_interval"))
  expect_identical(intervals[, 1L], head(diagnostics$knots, -1L))
  expect_identical(intervals[, 2L], tail(diagnostics$knots, -1L))
  expect_true(any(intervals[, 2L] - intervals[, 1L] < 1e-14))
  expect_true(all(vapply(diagnostics$panels, function(panel){
    identical(panel$coordinate, "affine probability panel") &&
      identical(panel$jacobian, diff(panel$physical_interval)) &&
      panel$subdivisions <= panel$cap && panel$status %in% c("OK", "criterion_converged_at_cap")
  }, logical(1))))
  priors <- list(prior_PET("normal", list(0, .4), truncation = list(-Inf, Inf)))
  plotted <- .plot_data_prior_list.PETPEESE(priors, c(0, .5, 1), c(0, 1), NULL,
    3, 64, NULL, NULL, FALSE, list(prior("normal", list(.3, .7))))
  sd <- sqrt(.7^2 + (.4 * plotted$x)^2)
  expect_equal(plotted$y, rep(.3, 3), tolerance = 1e-6)
  expect_equal(plotted$y_lCI, .3 + stats::qnorm(.025) * sd, tolerance = 1e-6)
  expect_equal(plotted$y_uCI, .3 + stats::qnorm(.975) * sd, tolerance = 1e-6)
})

test_that("PET panels retain directly represented default primitive tails and strict truncation limits", {

  mu <- .petpeese_prior_components(prior("normal", list(0, 1)))[[1L]]
  bias <- .petpeese_prior_components(prior("normal", list(0, 1), list(0, Inf)))[[1L]]
  upper_p <- 2.9893409125836859e-17
  quantile <- .petpeese_prior_component_quantile(mu, 1, upper_p = upper_p)
  expect_equal(quantile, stats::qnorm(upper_p, lower.tail = FALSE), tolerance = 1e-12)
  q <- 7.6096083718688288
  value <- .petpeese_prior_probability_integral(mu, bias, 1, q)
  expect_equal(as.numeric(value), stats::pnorm(q / sqrt(2))^2, tolerance = 1e-7)
  diagnostics <- attr(value, "numerical_diagnostics")
  expect_lte(diagnostics$subdivisions, 200L)
  expect_equal(tail(diagnostics$panels, 1L)[[1L]]$physical_interval,
    c(.99999999999998623, 1), tolerance = 1e-15)
  for(scale in c(-1, 1)){
    targets <- c(-2, 0, 2)
    actual <- vapply(targets, function(q) as.numeric(.petpeese_prior_sum_cdf(mu, bias, scale, q)), numeric(1))
    reference <- if(scale > 0) stats::pnorm(targets / sqrt(2))^2 else{
      upper <- stats::pnorm(-targets / sqrt(2), lower.tail = FALSE)
      2 * upper - upper^2
    }
    expect_equal(actual, reference, tolerance = 1e-7)
  }
  truncated <- .petpeese_prior_components(prior("normal", list(0, 1), list(-.5, .5)))[[1L]]
  expect_error(.petpeese_prior_component_quantile(truncated, 1, upper_p),
    "truncated lower probability lost", fixed = TRUE)
  expect_error(.petpeese_prior_component_quantile(mu, 1, 0),
    "interior probability tail was lost", fixed = TRUE)
  expect_equal(.petpeese_prior_component_quantile(truncated, .75, .25),
    stats::qnorm(stats::pnorm(-.5) + .75 * (stats::pnorm(.5) - stats::pnorm(-.5))), tolerance = 1e-12)
})

test_that("posterior-only atomic ownership maps without inventing a prior promise", {

  pair <- .model_probability_pair(c(.25, .75), log(c(.25, .75)), "posterior", "ordinary")
  draws <- structure(c(3, 3, 5), class = c("mixed_posteriors", "mixed_posteriors.simple"), parameter = "theta")
  atoms <- .posterior_atoms_new(matrix(c(3, 5), ncol = 1), pair$probabilities,
    column_names = "theta", source = "model_probabilities",
    component_probabilities = pair$probabilities, component_log_probabilities = pair$logs,
    model_probability_declaration = pair$declaration)
  posterior_metadata(draws, "atoms") <- atoms
  mapped <- posterior_transform(draws, "lin", list(a = 1, b = 2))
  expect_identical(as.numeric(mapped), c(7, 7, 11))
  mapped_atoms <- posterior_metadata(mapped, "atoms")
  expect_identical(as.numeric(mapped_atoms$locations), c(7, 11))
  expect_identical(mapped_atoms$mass, pair$probabilities)
  expect_identical(mapped_atoms$component_log_probabilities, pair$logs)
  expect_null(posterior_metadata(mapped, "prior_density"))
  expect_null(posterior_metadata(mapped, "measure_unavailable"))
  bad <- atoms
  bad$model_probability_declaration$model_indices <- c(1L, 1L)
  expect_error(posterior_metadata(draws, "atoms") <- bad, "Model probability ownership", fixed = TRUE)
  promised <- draws
  promised <- .bt_meta_set(promised, "model_probabilities", list(
    prior = .model_probability_pair(c(.25, .75), log(c(.25, .75)), "prior", "ordinary"), posterior = pair))
  expect_error(posterior_transform(promised, "lin", list(a = 1, b = 2)),
    "source has no complete declared prior context", fixed = TRUE)
})

test_that("scalar omega declarations retain actual log owners and independent support", {

  priors <- list(prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5))),
    prior_weightfunction("one-sided", .1, wf_fixed(c(1, .25))))
  columns <- .weightfunction_mapping_info(priors)$names
  pair <- .model_probability_pair(c(1, 0), c(0, -1000), "posterior")
  laws <- .posterior_weightfunction_scalar_laws(priors, pair$probabilities, columns,
    posterior_pair = pair)
  expect_equal(laws$marginals[[1L]]$mass, 1)
  expect_identical(laws$marginals[[1L]]$component_log_probabilities, pair$logs)
  expect_identical(laws$marginals[[1L]]$model_probability_declaration, pair$declaration)
  expect_s3_class(laws$marginals[[2L]], "BayesTools_formula_atoms_unavailable")
  expect_identical(laws$marginals[[2L]]$diagnostics$log_posterior_probabilities, pair$logs)
  ordinary <- .model_probability_pair(c(.25, .75), log(c(.25, .75)), "posterior", "ordinary")
  laws <- .posterior_weightfunction_scalar_laws(priors, ordinary$probabilities, columns,
    posterior_pair = ordinary)
  expect_equal(laws$marginals[[1L]]$mass, 1)
  expect_equal(nrow(laws$marginals[[1L]]$locations), 1L)
  expect_equal(laws$marginals[[2L]]$mass, c(.25, .75))
  expect_equal(unname(laws$marginals[[2L]]$locations[, 1L]), c(.5, 1))
  fixed <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  continuous <- prior_weightfunction("one-sided", .05, wf_cumulative(c(1, 1)))
  probabilities <- c(.9, rep(.01, 8), .02 - 1e-16, 1e-16)
  partition_pair <- .model_probability_pair(probabilities, log(probabilities), "posterior", "ordinary")
  partition <- .posterior_weightfunction_scalar_laws(c(rep(list(fixed), 10), list(continuous)),
    probabilities, c("omega[0,0.05]", "omega[0.05,1]"), posterior_pair = partition_pair)$marginals[[2L]]
  if(inherits(partition, "BayesTools_formula_atoms_unavailable")){
    expect_identical(partition$reason, "numerical_model_probability_unavailable")
    expect_identical(partition$diagnostics$log_posterior_probabilities, partition_pair$logs)
  }else{
    expect_lt(sum(partition$mass), 1)
    expect_identical(partition$component_log_probabilities, partition_pair$logs)
  }
})

test_that("current declared coefficient and contribution recipes materialize their distinct laws", {

  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "scale"
  priors <- list(slope = slope, scale = prior("point", list(2)))
  context <- .prior_density_build_context(priors, c("slope", "scale"), n_grid = 512)
  raw <- structure(c(1, 2, 3), class = c("mixed_posteriors", "mixed_posteriors.simple"), parameter = "slope")
  posterior_metadata(raw, "prior_context") <- context
  contribution <- structure(c(2, 4, 6), class = class(raw), parameter = "slope")
  posterior_metadata(contribution, "prior_context") <- context
  posterior_metadata(contribution, "linear_weights") <- c(slope = 1, scale = 0)
  posterior_metadata(contribution, "linear_weight_space") <- "formula_contribution"
  set.seed(272)
  seed <- .Random.seed
  for(source in list(raw, contribution)){
    mapped <- posterior_transform(source, "lin", list(a = 10, b = 2))
    law <- posterior_metadata(mapped, "prior_density")
    expected_sd <- if(is.null(posterior_metadata(source, "linear_weights"))) 2 else 4
    expect_equal(prior_density_ordinate(law, 11)$log_density,
      stats::dnorm(11, 10, expected_sd, log = TRUE), tolerance = 1e-12)
    expect_identical(.Random.seed, seed)
  }
  expect_identical(attr(context$prior_list$slope, "multiply_by", exact = TRUE), "scale")
})

test_that("bias simplification preserves raw prior shares and owned posterior scalar masses", {

  # Synthetic declared-source carrier, without a fitted-model claim.
  priors <- list(prior_PET("point", list(2), prior_weights = 3),
    prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)), prior_weights = 7))
  columns <- c("PET", .weightfunction_prior_list_context(priors, merge = FALSE)$omega_names)
  pair <- .model_probability_pair(c(.25, .75), log(c(.25, .75)), "posterior", "ordinary")
  values <- rbind(c(2, 1, 1), c(2, 1, 1), c(0, 1, .5), c(0, 1, .5))
  colnames(values) <- columns
  class(values) <- c("mixed_posteriors", "mixed_posteriors.bias", "matrix", "array")
  attr(values, "prior_list") <- priors
  values <- .posterior_atoms_set(values, .posterior_atoms_new(
    rbind(c(2, 1, 1), c(0, 1, .5)), pair$probabilities, column_names = columns,
    component_probabilities = pair$probabilities, component_log_probabilities = pair$logs,
    model_probability_declaration = pair$declaration))
  values <- .posterior_weightfunction_declarations(values, priors, pair$probabilities, posterior_pair = pair)
  samples <- list(bias = values)
  for(parameter in c("omega", "PET")){
    selected <- .simplify_as_mixed_posterior_bias(samples, parameter)
    leaf <- selected[[parameter]]
    expect_identical(vapply(attr(leaf, "prior_list", exact = TRUE), .prior_model_weight, numeric(1)), c(3, 7))
    atoms <- posterior_metadata(leaf, "atoms")
    expect_identical(atoms$component_log_probabilities, pair$logs)
    expect_identical(atoms$model_probability_declaration, pair$declaration)
    if(parameter == "omega") expect_equal(atoms$marginals[[2L]]$mass, c(.25, .75))
  }
})

test_that("new transformed-prior refusals preserve optional existing model diagnostics", {

  law <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))), c(theta = 1))
  payload <- list(model_indices = 1:2, log_prior_probabilities = log(c(.25, .75)),
    log_posterior_probabilities = c(0, -1000), stage = "posterior")
  refusal <- errorCondition("The transformed prior law is unavailable.", call = NULL,
    class = c("BayesTools_formula_prior_density_unavailable", "BayesTools_formula_measure_unavailable"),
    reason = "numerical_model_probability_unavailable", diagnostics = payload)
  scalar <- structure(c(1, 2, 3), class = c("mixed_posteriors", "mixed_posteriors.simple"), parameter = "theta")
  scalar <- .bt_meta_set(scalar, "prior_density", law)
  scalar <- .bt_formula_measure_mark(scalar, "theta", "atoms", "Existing atom refusal.",
    cause = "numerical_model_probability_unavailable", diagnostics = payload)
  scalar <- .bt_formula_measure_mark(scalar, "theta", "support", "Existing support refusal.",
    cause = "unsupported_contribution_measure")
  matrix <- structure(matrix(c(1, 2, 3), ncol = 1, dimnames = list(NULL, "theta")), class = c("mixed_posteriors", "matrix", "array"))
  matrix <- .bt_meta_assign(matrix, .bt_meta_get_fields(scalar, .bt_meta_fields()))
  testthat::with_mocked_bindings({
    for(source in list(scalar, matrix, list(theta = scalar))){
      mapped <- posterior_transform(source, "lin", list(a = 1, b = 2))
      leaf <- if(is.list(mapped)) mapped$theta else mapped
      unavailable <- posterior_metadata(leaf, "measure_unavailable")
      expect_setequal(unavailable$measure, c("atoms", "support", "prior_density"))
      index <- match("prior_density", unavailable$measure)
      expect_identical(unavailable$cause[[index]], refusal$reason)
      expect_identical(unavailable$reason[[index]], conditionMessage(refusal))
      expect_identical(unavailable$diagnostics[[index]], payload)
      condition <- tryCatch(.bt_formula_measure_check(leaf, "prior_density", "theta"), error = identity)
      expect_identical(condition$reason, refusal$reason)
      expect_identical(condition$detail, conditionMessage(refusal))
      expect_identical(condition$diagnostics, payload)
    }
  }, .prior_density_output_transform = function(...) stop(refusal))
})

test_that("narrow huge finite pieces retain physical geometry and kernel scale", {

  bounds <- c(1e300, 1e300 * (1 + 1e-12))
  slope <- prior("normal", list(0, 1e-300))
  attr(slope, "multiply_by") <- "scale"
  law <- .prior_linear_combination_density(list(slope = slope,
    scale = prior("uniform", as.list(bounds))), c(slope = 1), n_grid = 8192)
  reference <- stats::integrate(function(probability){
    scale <- bounds[1L] + diff(bounds) * probability
    stats::dnorm(1, sd = 1e-300 * scale)
  }, 0, 1, rel.tol = 1e-12)$value
  result <- prior_density_ordinate(law, 1)
  expect_identical(result$behavior, "regular")
  expect_lt(abs(exp(result$log_density) / reference - 1), 1e-4)
  expect_true(all(result$provenance$integration$piece_evaluations <= 8192L))
})

test_that("transformed prior mixtures keep original knot and Jacobian correspondence", {

  priors <- list(prior("gamma", list(.5, 1)), prior("gamma", list(2, 1)))
  grid <- seq(0, 5, length.out = 101)
  data <- .plot_data_prior_list.simple(priors, grid, range(grid), NULL, 101, 1000,
    FALSE, FALSE, "exp_lin", list(a = 0, b = .5), FALSE)
  retained <- grid[-1L]
  expected <- (.5 * stats::dgamma(retained, .5) + .5 * stats::dgamma(retained, 2)) * 2 * sqrt(retained)
  expect_equal(data$density$x, sqrt(retained), tolerance = 1e-14)
  expect_equal(data$density$y, expected, tolerance = 1e-12)

  decreasing <- .plot_data_prior_list.simple(list(prior("uniform", list(0, 2))),
    seq(0, 1, length.out = 11), c(0, 1), NULL, 11, 100,
    FALSE, FALSE, "lin", list(a = 0, b = -1), FALSE)$density
  expect_equal(decreasing$x[1L], -1)
  expect_equal(decreasing$y[1L], .5)
  expect_error(.plot_data_prior_list.simple(priors, grid, range(grid), NULL, 101, 1000,
    FALSE, FALSE, "lin", list(a = 1e300, b = 1), FALSE),
    class = "BayesTools_nonmonotone_transformation")
})

test_that("finite inverse region loss refuses while genuine support events remain exact", {

  law <- .prior_linear_combination_density(list(theta = prior("gamma", list(.001, 1))),
    c(theta = 1e30), n_grid = 512)
  region <- .hypothesis_prior_region(quote(theta < 1e-300), "theta")
  expect_error(.prior_linear_density_region_probability(law, region),
    class = "BayesTools_numerical_condition")
  side <- hypothesis_parse("theta < 1e-300")$statements[[1L]]$left
  expect_error(.hypothesis_prior_density_prob(law, side, "theta"),
    class = "BayesTools_hypothesis_region")
  empty <- .hypothesis_prior_region(quote(theta < -1e-300), "theta")
  expect_equal(as.numeric(.prior_linear_density_region_probability(law, empty)), 0)
  ordinary <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))),
    c(theta = 2), n_grid = 512)
  expect_equal(as.numeric(.prior_linear_density_region_probability(ordinary,
    .hypothesis_prior_region(quote(theta < 1), "theta"))), stats::pnorm(.5))
})

test_that("PET-PEESE pair CDFs retain analytic Normal sums and actual subdivision counts", {

  mu <- .petpeese_prior_components(prior("normal", list(.3, .8)))[[1L]]
  bias <- .petpeese_prior_components(prior("normal", list(-.2, .4)))[[1L]]
  for(scale in c(-2, .5, 0)){
    value <- .petpeese_prior_sum_cdf(mu, bias, scale, .7)
    expected <- stats::pnorm(.7, mean = .3 - .2 * scale, sd = sqrt(.8^2 + (.4 * scale)^2))
    expect_lt(abs(as.numeric(value) - expected), 1e-7)
    if(scale != 0){
      diagnostic <- attr(value, "numerical_diagnostics")
      expect_lte(diagnostic$subdivisions, 200L)
      expect_equal(diagnostic$subdivisions,
        sum(vapply(diagnostic$panels, `[[`, integer(1), "subdivisions")))
      expect_true(all(vapply(diagnostic$panels, function(panel){
        panel$message == "OK" || panel$status == "criterion_converged_at_cap"
      }, logical(1))))
    }
  }
  expect_error(.petpeese_prior_scale("PEESE", 1e-200, 1),
    class = "BayesTools_numerical_condition")
  expect_error(.petpeese_prior_scale("PEESE", 1e200, 1),
    class = "BayesTools_numerical_condition")
  expect_identical(.petpeese_prior_scale("PEESE", 0, 1), 0)
})

test_that("public ensemble tables preserve duplicated and missing row denominator selection", {

  models <- lapply(0:3, function(value){
    list(marglik = .bt_marglik_manual_result(value), prior_weights = 1)
  })
  inference <- list(
    first = ensemble_inference(models[1:2], "first", list(first = c(TRUE, FALSE)))$first,
    second = ensemble_inference(models, "second", list(second = c(TRUE, TRUE, FALSE, FALSE)))$second
  )
  table <- ensemble_inference_table(inference, c("first", "second"))
  expect_identical(attr(table, "n_models"), c(2L, 4L))
  for(selection in list(c(2, 2), c("second", "second"), c(2, 1), c(FALSE, TRUE))){
    selected <- table[selection, , drop = FALSE]
    positions <- if(is.character(selection)) match(selection, rownames(table)) else seq_len(nrow(table))[selection]
    expect_identical(attr(selected, "n_models"), attr(table, "n_models")[positions])
    expect_equal(as.numeric(selected$inclusion_BF), as.numeric(table$inclusion_BF[positions]))
  }
  expect_identical(attr(table[c(NA_integer_, 2L), , drop = FALSE], "n_models"), c(NA_integer_, 4L))
  expect_identical(attr(table[integer(), , drop = FALSE], "n_models"), integer())
  expect_identical(attr(table["models"], "n_models"), attr(table, "n_models"))
  expect_true(any(grepl("2/4", capture.output(print(table[c(2, 2), , drop = FALSE])), fixed = TRUE)))
})

test_that("PET-PEESE required support transitions refuse while optional anchors record omissions", {

  mu <- .petpeese_prior_components(prior("normal", list(0, 1)))[[1L]]
  bias <- .petpeese_prior_components(prior("uniform", list(0, 1)))[[1L]]
  # Overflow in a declared finite support transition must propagate before
  # integration, even when quantile anchors are still representable.
  expect_error(.petpeese_prior_probability_integral(mu, bias, 1e308, -1e308),
    class = "BayesTools_numerical_condition")
  refused_cdf <- mu
  refused_cdf$cdf <- function(q){
    if(any(q == -1)) .petpeese_prior_numerical_stop("required CDF test", "The support-transition CDF is unavailable")
    stats::pnorm(q)
  }
  expect_error(.petpeese_prior_probability_integral(refused_cdf, bias, 1, 0),
    class = "BayesTools_numerical_condition")
  optional_bias <- bias
  optional_bias$quant <- function(p){
    if(any(p == .02)) .petpeese_prior_numerical_stop("optional quantile test", "The optional quantile is unavailable")
    stats::qunif(p)
  }
  optional_bias$partition_cache <- new.env(parent = emptyenv())
  result <- .petpeese_prior_probability_integral(mu, optional_bias, 1, 0)
  expect_lt(abs(as.numeric(result) - (stats::pnorm(-1) + stats::dnorm(0) - stats::dnorm(1))), 1e-7)
  omitted <- attr(result, "numerical_diagnostics")$optional_anchor_diagnostics
  expect_true(any(vapply(omitted, function(entry){
    entry$operation == "bias quantile 0.02" && inherits(entry$condition, "BayesTools_numerical_condition")
  }, logical(1))))
})
