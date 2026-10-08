skip_if_not_test_profile("unit")

test_that("owned structural zero splits agree without subtracting infinities", {
  owned <- function(value, log_value){
    pair <- .model_probability_pair(value, log_value, "prior", "raw")
    .set_prior_model_probability(prior("normal", list(0, 1)), pair$probabilities, pair$logs, pair$declaration)
  }
  result <- .model_probability_split_prior(owned(0, -Inf), owned(0, -Inf), .5)
  expect_identical(.prior_model_weight(result), 0)
  expect_identical(.prior_model_log_weight(result), -Inf)
  expect_identical(attr(result, "model_probability_declaration")$stage, "component")
  finite <- .model_probability_split_prior(owned(1, 0), owned(.5, log(.5)), .5)
  expect_identical(.prior_model_weight(finite), .5)
  expect_identical(.prior_model_log_weight(finite), log(.5))
  for(child in list(owned(1, 0), owned(.25, log(.25)))){
    expect_error(.model_probability_split_prior(owned(0, -Inf), child, .5), "contradicts the requested parent split", fixed = TRUE)
  }
  expect_identical(.prior_model_weight(.model_probability_split_prior(prior("normal", list(0, 1)), prior("normal", list(0, 1)), .5)), .5)
  expect_identical(.prior_model_weight(.model_probability_split_prior(owned(0, -Inf), prior("normal", list(0, 1)), .5)), 0)
  tiny <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(1e-200)), prior_weights = 1e-200)
  condition <- expect_error(.plot_prior_spike_and_slab_components(tiny), class = "BayesTools_numerical_range_limit")
  expect_s3_class(condition, "BayesTools_numerical_condition")
  expect_identical(condition$operation, "probability split")
  expect_identical(condition$requested_scale, "natural")
  expect_null(conditionCall(condition))
  pair <- .model_probability_pair(1e-200, log(1e-200), "prior", "raw")
  parent <- .set_prior_model_probability(prior("normal", list(0, 1)), pair$probabilities, pair$logs, pair$declaration)
  retained <- .model_probability_split_prior(parent, prior("normal", list(0, 1)), 1e-200)
  expect_identical(.prior_model_weight(retained), 0)
  expect_equal(.prior_model_log_weight(retained), 2 * log(1e-200), tolerance = 1e-12)
})

test_that("plain identical scalar point targets exclude only unrelated zero weights", {
  weights <- c(1e-300, 1e100)
  a <- lapply(weights, function(w) prior("point", list(2), prior_weights = w))
  b <- lapply(weights, function(w) prior("normal", list(0, 1), prior_weights = w))
  context <- .prior_density_model_mixture_context(list(a = a, b = b), c("a", "b"))
  density <- .prior_density_from_context(context, c(a = 1, b = 0))
  expect_identical(density$points, data.frame(x = 2, p = 1))
  expect_null(density$density)
  expect_error(.prior_density_from_context(context, c(a = 1, b = 1)), class = "BayesTools_formula_measure_unavailable")
  unequal <- context; unequal$prior_list$a[[2L]] <- prior("point", list(3), prior_weights = weights[2L])
  expect_error(.prior_density_from_context(unequal, c(a = 1, b = 0)), class = "BayesTools_formula_measure_unavailable")
  multiplier <- context
  for(i in seq_along(a)) attr(multiplier$prior_list$a[[i]], "multiply_by") <- "b"
  expect_error(.prior_density_from_context(multiplier, c(a = 1, b = 0)), class = "BayesTools_formula_measure_unavailable")
  expect_identical(.prior_density_from_context(context, c(a = 0, b = 0))$points, data.frame(x = 0, p = 1))
  compiled <- JAGS_formula(~x, "mu", data.frame(x = c(10, 20, 30)),
    list(intercept = prior("point", list(2)), x = prior("normal", list(0, 1))), formula_scale = TRUE)
  primitive <- .prior_density_context(compiled$prior_list, c("mu_intercept", "mu_x"), list(mu = compiled$formula_scale))
  transformed <- .prior_density_model_mixture_context(list(mu_intercept = a, mu_x = b), c("mu_intercept", "mu_x"))
  transformed$transforms <- primitive$transforms
  expect_error(.prior_density_from_context(transformed, c(mu_intercept = 1, mu_x = 0)), class = "BayesTools_formula_measure_unavailable")
})
source(testthat::test_path("common-functions.R"))

test_that("model probability refit requirements retain exact typed diagnostics", {
  priors <- list(prior("normal", list(0, 1)), prior("point", list(0)))
  context <- .prior_density_model_mixture_context(list(theta = priors), "theta")
  legacy <- context; legacy$schema_version <- 1L
  malformed <- context; malformed$n_grid <- 0
  pair <- .prior_model_probability_pair(priors)
  partially_owned <- priors
  partially_owned[[1L]] <- .set_prior_model_probability(priors[[1L]], pair$probabilities[[1L]],
    pair$logs[[1L]], within(pair$declaration, model_indices <- model_indices[1L]))
  cases <- list(
    list(read = function() .model_probability_validate(.5, log(.5), NULL),
      message = "Model probability ownership is missing or malformed. Recompute or refit with the current BayesTools version."),
    list(read = function() .model_probability_context_validate(malformed),
      message = "The model prior context is incomplete or malformed. Recompute or refit with the current BayesTools version."),
    list(read = function() .model_probability_context_validate(legacy),
      message = "Model prior contexts require current probability ownership. Recompute or refit with the current BayesTools version."),
    list(read = function() .prior_model_probability_pair(partially_owned),
      message = "The model prior distributions have incomplete probability ownership. Recompute or refit with the current BayesTools version.")
  )
  for(case in cases){
    condition <- expect_error(case$read(), case$message, fixed = TRUE, class = "BayesTools_refit_required")
    expect_identical(conditionMessage(condition), case$message)
    expect_null(conditionCall(condition))
  }
})

.model_probability_test_model <- function(prior_object, evidence, weight = 1) {
  values <- if(is.prior.point(prior_object)) rep(prior_object$parameters$location, 20) else seq(-1, 1, length.out = 20)
  values <- cbind(theta = values)
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(values)), sample = 20,
    summary.pars = list(mutate = NULL), monitor = "theta"), class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- list(theta = prior_object)
  list(fit = attach_test_parameter_map(fit), marglik = bridgesampling_object(evidence), prior_weights = weight)
}

.model_probability_weightfunction_model <- function(prior_object, evidence = 0, weight = 1) {
  omega <- if(prior_object$weights$type == "fixed") prior_object$weights$omega else c(1, .5)
  values <- matrix(rep(omega, each = 20), 20)
  colnames(values) <- paste0("omega[", seq_along(omega), "]")
  if(prior_object$weights$type != "fixed") values[, 2] <- seq(.1, .9, length.out = 20)
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(values)), sample = 20,
    summary.pars = list(mutate = NULL), monitor = colnames(values)), class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- list(bias = prior_object)
  list(fit = attach_test_parameter_map(fit), marglik = bridgesampling_object(evidence), prior_weights = weight)
}

test_that("joint numerical refusal retains independently declared omega scalars", {
  first <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  second <- prior_weightfunction("one-sided", .1, wf_fixed(c(1, .25)))
  mixed <- mix_posteriors(list(.model_probability_weightfunction_model(first),
    .model_probability_weightfunction_model(second, -1000)), "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  atoms <- posterior_metadata(mixed$bias, "atoms")
  expect_s3_class(atoms, "BayesTools_posterior_atoms")
  expect_error(posterior_atoms_free(mixed$bias), class = "BayesTools_formula_measure_unavailable")
  plot <- tryCatch(.plot_data_samples.weightparameter(list(omega = mixed$bias), 1L, 128L), error = identity)
  expect_false(inherits(plot, "condition"))
  if(!inherits(plot, "condition")){
    expect_equal(plot$points1$x, 1, tolerance = 0)
    expect_equal(plot$points1$y, 1, tolerance = 0)
    expect_null(plot$density)
  }
  unknown <- tryCatch(.plot_data_samples.weightparameter(list(omega = mixed$bias), 2L, 128L), error = identity)
  expect_s3_class(unknown, "BayesTools_formula_measure_unavailable")
  if(!is.null(atoms)){
    expect_false(atoms$joint_declared)
    expect_s3_class(atoms$joint_unavailable, "BayesTools_formula_measure_unavailable")
    expect_identical(atoms$joint_unavailable$diagnostics$log_posterior_probabilities, c(0, -1000))
    expect_null(conditionCall(atoms$joint_unavailable))
    expect_identical(dim(atoms$locations), c(0L, 3L))
    expect_equal(as.numeric(.posterior_atoms_for_column(atoms, 1L)$locations[, 1L]), 1, tolerance = 0)
    expect_equal(.posterior_atoms_for_column(atoms, 1L)$mass, 1, tolerance = 0)
    expect_null(.posterior_atoms_for_column(atoms, 2L))
    expect_error(.posterior_atoms_get(mixed$bias), class = "BayesTools_formula_measure_unavailable")
    expect_identical(.posterior_atoms_get(mixed$bias, allow_partial = TRUE), atoms)
    transformed <- posterior_transform(mixed$bias, "lin", list(a = 2, b = 3))
    expect_identical(dim(transformed), dim(mixed$bias))
    expect_equal(as.numeric(transformed), 2 + 3 * as.numeric(mixed$bias), tolerance = 0)
    transformed_atoms <- posterior_metadata(transformed, "atoms")
    expect_false(transformed_atoms$joint_declared)
    expect_equal(as.numeric(.posterior_atoms_for_column(transformed_atoms, 1L)$locations[, 1L]), 5, tolerance = 0)
    expect_error(.posterior_atoms_get(transformed), class = "BayesTools_formula_measure_unavailable")
    design <- matrix(c(2, 0, 0), 1L, dimnames = list("selected", colnames(mixed$bias)))
    scalar <- .posterior_atoms_linear_transform(atoms, design)
    expect_true(scalar$joint_declared)
    expect_equal(as.numeric(scalar$locations[, 1L]), 2, tolerance = 0)
    design[1L, 2L] <- 1
    expect_s3_class(.posterior_atoms_linear_transform(atoms, design), "BayesTools_formula_measure_unavailable")
    simplified <- .simplify_as_mixed_posterior_bias(mixed, "omega")
    expect_false(posterior_metadata(simplified$omega, "atoms")$joint_declared)
    numeric_copy <- as.numeric(mixed$bias[, 1L])
    posterior_metadata(numeric_copy, "atoms") <- atoms
    expect_error(.posterior_atoms_get(numeric_copy), class = "BayesTools_formula_measure_unavailable")
    malformed <- atoms; malformed$joint_declared <- c(FALSE, FALSE)
    expect_error(.posterior_atoms_from_attribute(malformed))
    malformed <- atoms; malformed$marginals <- rep(list(NULL), 3L); names(malformed$marginals) <- colnames(atoms$locations)
    expect_error(.posterior_atoms_from_attribute(malformed), "certified scalar marginals", fixed = TRUE)
  }
  equal <- mix_posteriors(list(.model_probability_weightfunction_model(first),
    .model_probability_weightfunction_model(first, -1000)), "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_true(.posterior_atoms_get(equal$bias)$joint_declared)
  expect_equal(.posterior_atoms_get(equal$bias)$mass, 1, tolerance = 0)
})

test_that("posterior support uses declared posterior logs and fixed inclusion boundaries", {
  narrow <- .model_probability_test_model(prior("uniform", list(1, 2)), 0)
  broad <- .model_probability_test_model(prior("uniform", list(-20, 20)), -Inf)
  narrow$fit$mcmc <- coda::mcmc.list(coda::mcmc(cbind(theta = seq(1.1, 1.9, length.out = 20))))
  # Reattach the unchanged declared geometry after replacing this synthetic fixture's draws.
  narrow$fit <- attach_test_parameter_map(narrow$fit)
  for(evidence in c(-Inf, -1000)){
    broad$marglik <- bridgesampling_object(evidence)
    mixed <- mix_posteriors(list(narrow, broad), "theta", list(c(FALSE, FALSE)), seed = 3, n_samples = 20)
    expected <- if(is.infinite(evidence)) c(1, 2) else c(-20, 20)
    expect_identical(posterior_metadata(mixed$theta, "support")$bounds, expected)
    expect_length(attr(mixed$theta, "prior_list"), 2L)
  }
  for(probability in c(0, 1, .5, 1 - 1e-9)){
    p <- prior_spike_and_slab(prior("uniform", list(1, 2)), prior("point", list(probability)))
    support <- .posterior_support_from_prior(p)
    expected <- if(probability == 0) c(0, 0) else if(probability == 1) c(1, 2) else c(0, 2)
    expect_identical(support$bounds, expected)
  }
  first <- narrow; first$prior_weights <- 1e-300
  second <- broad; second$prior_weights <- 1e100; second$marglik <- bridgesampling_object(0)
  mixed <- mix_posteriors(list(first, second), "theta", list(c(FALSE, FALSE)), seed = 3, n_samples = 20)
  context <- .prior_density_model_mixture_context(list(theta = attr(mixed$theta, "prior_list")), "theta")
  expect_identical(context$model_weights[1L], 0)
  expect_true(is.finite(context$model_log_weights[1L]))
  supports <- .posterior_components_supports(context, data.frame(.model = 1:2), c(theta = 1))
  expect_identical(supports[[1L]]$bounds, c(1, 2))
  expect_identical(supports[[2L]]$bounds, c(-20, 20))
})

test_that("weightfunction null locations are exactly one", {
  wf <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  expect_error(weightfunctions_mapping(list(wf, prior("point", list(1 - 1e-9)))))
  expect_silent(weightfunctions_mapping(list(wf, prior("point", list(1)))))
  expect_silent(weightfunctions_mapping(list(wf, prior_none())))
  expect_false(.is_prior_weightfunction_null(prior("point", list(1 - 1e-9))))
  expect_true(.is_prior_weightfunction_null(prior("point", list(1))))
})

test_that("compiled factor point mixtures retain the declared factor shape", {
  make <- function(location, contrast, levels = c("a", "b", "c")){
    data <- data.frame(f = factor(rep(levels, 2L), levels = levels))
    compiled <- JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
      f = prior_factor("point", list(location), contrast = contrast)))
    draws <- .generate_prior_sample_matrix(compiled$prior_list, 8L, seed = 17)
    fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), sample = 8L,
      summary.pars = list(mutate = NULL), monitor = colnames(draws)), class = c("runjags", "BayesTools_fit", "list"))
    attr(fit, "prior_list") <- compiled$prior_list
    attr(fit, "formula_design") <- list(mu = compiled$formula_design)
    list(fit = attach_test_parameter_map(fit), marglik = bridgesampling_object(0), prior_weights = 1)
  }
  for(contrast in c("treatment", "independent")){
    models <- list(make(0, contrast), make(2, contrast))
    mixed <- mix_posteriors(models, "mu_f", list(c(FALSE, FALSE)), seed = 1, n_samples = 20L)
    expect_true(all(mixed$mu_f %in% c(0, 2)))
    expect_identical(attr(mixed$mu_f, "factor_design"), attr(attr(models[[1L]]$fit, "prior_list")$mu_f, "factor_design"))
    atoms <- .posterior_atoms_get(mixed$mu_f)
    expect_true(atoms$joint_declared)
    expect_equal(sort(atoms$mass), c(.5, .5), tolerance = 0)
    levels <- marginal_posterior(mixed, "mu_f", use_formula = FALSE, prior_samples = FALSE)
    expect_identical(names(levels), c("a", "b", "c"))
    if(contrast == "treatment"){
      expect_identical(as.numeric(levels[[1L]]), rep(0, 20L))
      expect_equal(.posterior_atoms_get(levels[[1L]])$mass, 1, tolerance = 0)
    }
    incompatible <- list(models[[1L]], make(2, contrast, c("u", "v", "w")))
    expect_error(mix_posteriors(incompatible, "mu_f", list(c(FALSE, FALSE)), seed = 1, n_samples = 20L))
  }
})

test_that("empty allocated absent sources retain the complete ordered declaration", {
  ordered <- ordered_plot_test_fixture(prior("point", list(4)), allocation = c(.25, .75), levels = c("a", "b", "c"))$fit
  compiled <- JAGS_formula(~1, "mu", data.frame(f = factor(c("a", "b", "c"))), list(intercept = prior("point", list(0))))
  draws <- .generate_prior_sample_matrix(compiled$prior_list, 8L, seed = 17)
  absent <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), sample = 8L,
    summary.pars = list(mutate = NULL), monitor = colnames(draws)), class = c("runjags", "BayesTools_fit", "list"))
  attr(absent, "prior_list") <- compiled$prior_list
  attr(absent, "formula_design") <- list(mu = compiled$formula_design)
  absent <- attach_test_parameter_map(absent)
  for(evidence in c(-Inf, -1000)){
    models <- list(list(fit = ordered, marglik = bridgesampling_object(evidence), prior_weights = 1),
      list(fit = absent, marglik = bridgesampling_object(0), prior_weights = 1))
    mixed <- mix_posteriors(models, "mu_f", list(c(FALSE, TRUE)), seed = 1, n_samples = 20L)
    source <- posterior_metadata(mixed$mu_f, "ordered_source")
    expect_identical(dim(source$primitives), c(20L, 0L))
    expect_null(.bt_ordered_source_validate(source))
    expect_identical(source$model, rep(2L, 20L))
    expect_length(source$models, 2L)
    expect_true(is.prior.ordered(source$models[[1L]]$prior))
    expect_identical(source$model_log_probabilities[1L], evidence)
    expect_true(all(as.numeric(mixed$mu_f) == 0))
    if(is.finite(evidence)) expect_error(posterior_atoms_free(mixed$mu_f), class = "BayesTools_formula_measure_unavailable") else
      expect_equal(sum(.posterior_atoms_get(mixed$mu_f)$mass), 1, tolerance = 0)
    if(is.finite(evidence)){
      levels <- marginal_posterior(mixed, "mu_f", use_formula = FALSE, prior_samples = FALSE)
      extracted <- .plot_data_marginal_level_samples(list(reference = levels[[1L]]), "reference")[[1L]]
      expect_false(posterior_atoms_free(extracted))
      density <- .plot_data_marginal_samples.den(extracted, 64L, NULL, NULL, NULL)
      expect_equal(density$points1$x, 0, tolerance = 0)
      expect_equal(density$points1$y, 1, tolerance = 0)
      expect_error(posterior_atoms_free(mixed$mu_f), class = "BayesTools_formula_measure_unavailable")
    }
  }
  source <- posterior_metadata(as_mixed_posteriors(ordered, "mu_f")$mu_f, "ordered_source")
  invalid <- source; invalid$primitives <- invalid$primitives[, 0L, drop = FALSE]
  expect_match(.bt_ordered_source_validate(invalid), "allocated ordered primitive", fixed = TRUE)
  invalid <- source; invalid$primitives <- invalid$primitives[, setdiff(colnames(invalid$primitives), invalid$models[[1L]]$total_names), drop = FALSE]
  expect_false(is.null(.bt_ordered_source_validate(invalid)))
})

# This oracle centers evidence before combining it with the original raw logs.
.model_probability_reference <- function(weights, evidence) {
  raw <- log(weights)
  raw <- raw - max(raw)
  prior <- raw - log(sum(exp(raw)))
  evidence <- evidence - max(evidence[is.finite(evidence)])
  score <- prior + evidence
  score <- score - max(score)
  list(prior = prior, posterior = score - log(sum(exp(score))))
}

test_that("raw model odds survive four hundred decades", {
  weights <- c(1e-300, 1e100)
  evidence <- c(1000, 0)
  reference <- .model_probability_reference(weights, evidence)
  inference <- compute_inference(weights, evidence, is_null = c(FALSE, TRUE))
  expect_equal(attr(inference, "log_BF"), 1000, tolerance = 1e-12)
  expect_equal(inference$post_probs, exp(reference$posterior), tolerance = 1e-12)
  expect_true(all(inference$post_probs > 0))
  expect_equal(log(inference$post_probs[[1L]]) - log(inference$post_probs[[2L]]),
    1000 - 400 * log(10), tolerance = 1e-12)
  expect_equal(attr(inference, "log_prior_probs"), reference$prior, tolerance = 1e-12)
  expect_equal(attr(inference, "log_post_probs"), reference$posterior, tolerance = 1e-12)
})

test_that("common huge evidence preserves ordinary model odds", {
  inference <- compute_inference(c(9, 1), rep(1e300, 2), c(FALSE, TRUE))
  expected_prior <- c(1, 1 / 9) / sum(c(1, 1 / 9))
  expect_identical(inference$prior_probs, expected_prior)
  expect_equal(inference$post_probs, expected_prior, tolerance = 1e-14)
  expect_equal(attr(inference, "log_BF"), 0, tolerance = 1e-12)
})

test_that("subnormal derived odds use the original raw declaration", {
  weights <- c(1e-200, 1e120)
  reference <- .model_probability_reference(weights, c(740, 0))
  inference <- compute_inference(weights, c(740, 0), c(FALSE, TRUE))
  expect_equal(inference$post_probs, exp(reference$posterior), tolerance = 1e-12)
  expect_equal(attr(inference, "log_BF"), 740, tolerance = 1e-12)
  raw_subnormal <- compute_inference(c(1e-320, 1), c(740, 0), c(FALSE, TRUE))
  expect_equal(raw_subnormal$post_probs,
    exp(.model_probability_reference(c(1e-320, 1), c(740, 0))$posterior), tolerance = 1e-12)
})

test_that("tiny positive models retain failure policy and conditional meaning", {
  expect_error(compute_inference(c(1e-300, 1e100), c(NA, 0)), class = "BayesTools_marglik_failure")
  expect_error(compute_inference(c(1e-300, 1e100), c(Inf, 0)), "Infinite positive")
  expect_warning(dropped <- compute_inference(c(1e-300, 1e100), c(NA, 0), on_failure = "drop"), "Dropped model")
  expect_identical(dropped$prior_probs, c(0, 1))
  expect_equal(attr(dropped, "marglik_failure")$original_log_prior_prob, -400 * log(10), tolerance = 1e-12)
  expect_warning(zero <- compute_inference(c(1e-300, 1e100), c(NA, 0), on_failure = "zero"), "Assigned zero evidence")
  expect_true(is.finite(attr(zero, "log_prior_probs")[[1L]]))
  conditional <- compute_inference(c(1e-300, 1e100), c(0, 0), c(FALSE, TRUE), conditional = TRUE)
  expect_identical(conditional$prior_probs, c(1, 0))
  expect_identical(conditional$post_probs, c(1, 0))
  expect_identical(attr(conditional, "log_prior_probs"), c(0, -Inf))
})

test_that("finite evidence range loss is a typed refusal", {
  expect_error(compute_inference(c(1, 1), c(1e308, -1e308)), class = "BayesTools_numerical_condition")
})

test_that("public natural probability zero retains its original meaning", {
  expect_true(is.na(inclusion_BF(c(0, 1), margliks = c(1000, 0), is_null = c(FALSE, TRUE))))
  inference <- compute_inference(c(0, 1), c(NA, 0))
  expect_identical(inference$prior_probs, c(0, 1))
  expect_identical(inference$post_probs, c(0, 1))
  expect_null(attr(inference, "marglik_failure"))
})

test_that("tiny posterior model laws do not become structural zero", {
  models <- list(.model_probability_test_model(prior("point", list(0)), 0),
    .model_probability_test_model(prior("normal", list(0, 1)), -1000))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_true(all(is.finite(mixed$theta)))
  expect_error(posterior_atoms_free(mixed$theta), class = "BayesTools_formula_measure_unavailable")
  descriptive <- marginal_posterior(mixed, "theta", use_formula = FALSE, prior_samples = FALSE)
  expect_true(all(is.finite(descriptive)))
  expect_error(posterior_atoms_free(descriptive), class = "BayesTools_formula_measure_unavailable")
  original <- tryCatch(posterior_atoms_free(descriptive), error = identity)
  plotting <- tryCatch(plot_marginal(list(theta = descriptive), "theta",
    prior = FALSE, plot_type = "ggplot"), error = identity)
  expect_identical(class(plotting), class(original))
  expect_identical(conditionMessage(plotting), conditionMessage(original))
  expect_identical(plotting$reason, original$reason)
  expect_identical(plotting$detail, original$detail)
  expect_identical(plotting$diagnostics, original$diagnostics)
  expect_null(conditionCall(plotting))
})

test_that("unrepresentable prior model laws retain usable mixed draws", {
  models <- list(.model_probability_test_model(prior("normal", list(0, 1)), 1000, 1e-300),
    .model_probability_test_model(prior("normal", list(0, 2)), 0, 1e100))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_true(all(is.finite(mixed$theta)))
  expect_error(marginal_posterior(mixed, "theta", use_formula = FALSE, prior_samples = TRUE),
    class = "BayesTools_formula_measure_unavailable")
  expect_true(posterior_atoms_free(mixed$theta))
})

test_that("ordinary visible probabilities and seeded allocations are preserved", {
  weights <- c(2, 3, 5)
  evidence <- log(c(.25, 1.5, .5))
  prior <- (weights / max(weights)) / sum(weights / max(weights))
  scores <- log(prior) + evidence
  post <- exp(scores - max(scores)) / sum(exp(scores - max(scores)))
  inference <- compute_inference(weights, evidence)
  expect_identical(inference$prior_probs, prior)
  expect_identical(inference$post_probs, post)
  expect_identical(attr(inference, "model_probability_declaration")$prior$route, "ordinary")
  expect_identical(attr(inference, "model_probability_declaration")$posterior$route, "ordinary")
  models <- lapply(seq_along(weights), function(i){
    .model_probability_test_model(prior("normal", list(0, i)), evidence[[i]], weights[[i]])
  })
  set.seed(47)
  counts <- as.integer(stats::rmultinom(1L, 40L, post)[, 1L])
  indices <- unlist(lapply(counts[counts > 0], function(count) sample(20, count, replace = TRUE)), use.names = FALSE)
  components <- rep(seq_along(counts), counts)
  values <- seq(-1, 1, length.out = 20)[indices]
  caller <- .Random.seed
  mixed <- mix_posteriors(models, "theta", list(rep(FALSE, 3)), seed = 47, n_samples = 40)
  expect_identical(as.numeric(mixed$theta), values)
  expect_identical(.bt_meta_get(mixed$theta, "draw_index"), indices)
  expect_identical(.bt_meta_get(mixed$theta, "component"), components)
  expect_identical(.Random.seed, caller)
})

test_that("paired probability owners reject inconsistent zero and metadata", {
  pair <- .model_probability_pair(c(0, 1), c(-800, 0), "prior")
  expect_silent(.model_probability_validate(pair$probabilities, pair$logs, pair$declaration))
  expect_error(.model_probability_pair(c(0, 1), c(-10, 0), "prior"), "ownership")
  expect_error(.model_probability_pair(c(1e-20, 1), c(-Inf, 0), "prior"), "ownership")
  declaration <- pair$declaration
  declaration$eta <- Inf
  expect_error(.model_probability_validate(pair$probabilities, pair$logs, declaration), "ownership")
  declaration <- pair$declaration
  declaration$route <- "guessed"
  expect_error(.model_probability_validate(pair$probabilities, pair$logs, declaration), "ownership")
  inference <- compute_inference(c(1e-300, 1e100), c(1000, 0))
  attr(inference, "log_prior_probs") <- c(-10, 0)
  expect_error(.model_probability_inference_get(inference, "prior"), "ownership")
  inference <- compute_inference(c(1, 1), c(0, 0))
  owner <- attr(inference, "model_probability_declaration")
  owner$prior$stage <- "posterior"
  attr(inference, "model_probability_declaration") <- owner
  expect_error(.model_probability_inference_get(inference, "prior"), "stage")
})

test_that("final subnormal displays certify their designated exponential rounding", {
  minimum <- .Machine$double.xmin * .Machine$double.eps
  logs <- c(log(minimum), log(2 * minimum), log(1e-320), log(.Machine$double.xmin) - .1)
  visible <- exp(logs)
  expect_true(all(visible > 0 & visible < .Machine$double.xmin))
  pair <- .model_probability_pair(visible, logs, "prior", "stabilized")
  expect_silent(.model_probability_validate(visible, logs, pair$declaration))
  changed <- visible
  changed[[1L]] <- 2 * minimum
  expect_error(.model_probability_validate(changed, logs, pair$declaration), "ownership")
  changed_logs <- logs
  changed_logs[[1L]] <- logs[[1L]] + 1
  expect_error(.model_probability_validate(visible, changed_logs, pair$declaration), "ownership")
  changed <- visible
  changed[[1L]] <- 0
  expect_error(.model_probability_validate(changed, logs, pair$declaration), "ownership")
  changed_owner <- pair$declaration
  changed_owner$route <- "ordinary"
  expect_error(.model_probability_validate(visible, logs, changed_owner), "ownership")
  expect_silent(.model_probability_pair(0, log(minimum) - 1, "posterior", "stabilized"))
  raw <- .Machine$double.xmin * (1 - (1:32) * .Machine$double.eps)
  expect_true(any(raw != exp(log(raw))))
  expect_silent(.model_probability_pair(raw, log(raw), "prior", "raw"))
  expect_error(.model_probability_pair(raw, log(raw), "prior", "stabilized"), "ownership")
})

test_that("prior reset copy merge and split retain their correct owner", {
  inference <- compute_inference(c(1e-300, 1e100), c(1000, 0))
  pair <- .model_probability_inference_get(inference, "prior")
  priors <- lapply(seq_len(2), function(i){
    declaration <- pair$declaration
    declaration$model_indices <- as.integer(i)
    .set_prior_model_probability(prior("normal", list(0, 1)), pair$probabilities[[i]], pair$logs[[i]], declaration)
  })
  reset <- .set_prior_model_weight(priors[[1L]], .5)
  expect_null(attr(reset, "model_log_prior_weights", exact = TRUE))
  expect_null(attr(reset, "model_probability_declaration", exact = TRUE))
  expect_identical(.prior_model_log_weight(reset), log(.5))
  copied <- .prior_density_copy_parent_attributes(prior("point", list(0)), priors[[1L]])
  expect_identical(.prior_model_log_weight(copied), pair$logs[[1L]])
  independent <- .prior_density_copy_parent_attributes(priors[[2L]], priors[[1L]])
  expect_identical(.prior_model_log_weight(independent), pair$logs[[2L]])
  zero <- .marginal_posterior_zero_vector_prior(priors[[1L]], 2)
  expect_identical(.prior_model_log_weight(zero), pair$logs[[1L]])
  merged <- .simplify_prior_list(priors)
  expect_length(merged, 1L)
  expect_equal(.prior_model_log_weight(merged[[1L]]), 0, tolerance = 1e-12)
  split <- .model_probability_split_prior(priors[[1L]], prior("normal", list(0, 1)), .25)
  expect_identical(.prior_model_weight(split), 0)
  expect_equal(.prior_model_log_weight(split), pair$logs[[1L]] + log(.25), tolerance = 1e-12)
  context <- .prior_density_model_mixture_context(list(theta = priors), "theta")
  expect_identical(context$schema_version, 2L)
  expect_true(is.finite(context$model_log_weights[[1L]]))
  expect_error(.prior_density_from_context(context, c(theta = 1)), class = "BayesTools_formula_measure_unavailable")
  context$model_log_weights[[1L]] <- -10
  expect_error(.model_probability_context_validate(context), "ownership")
  context$schema_version <- 1L
  expect_error(.model_probability_context_validate(context), "Recompute or refit")
})

test_that("posterior atom extremes preserve exact continuous and point controls", {
  for(evidence in list(c(0, -1000), c(-1000, 0))){
    models <- list(.model_probability_test_model(prior("normal", list(0, 1)), evidence[[1L]]),
      .model_probability_test_model(prior("normal", list(0, 2)), evidence[[2L]]))
    mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
    expect_true(posterior_atoms_free(mixed$theta))
  }
  models <- list(.model_probability_test_model(prior("point", list(2)), 0),
    .model_probability_test_model(prior("point", list(2)), -1000))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  atoms <- posterior_metadata(mixed$theta, "atoms")
  expect_identical(atoms$mass, 1)
  expect_identical(as.numeric(atoms$locations), 2)
  models[[2L]] <- .model_probability_test_model(prior("point", list(3)), -1000)
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  condition <- tryCatch(posterior_atoms_free(mixed$theta), error = identity)
  expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
  expect_identical(condition$reason, "numerical_model_probability_unavailable")
  expect_identical(condition$diagnostics$model_indices, 1:2)
  expect_equal(condition$diagnostics$log_posterior_probabilities, c(0, -1000), tolerance = 1e-12)
  models[[2L]]$prior_weights <- 0
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_identical(posterior_metadata(mixed$theta, "atoms")$mass, 1)
})

test_that("structured refusal diagnostics survive row and output transforms", {
  models <- list(.model_probability_test_model(prior("point", list(0)), 0),
    .model_probability_test_model(prior("normal", list(0, 1)), -1000))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  payload <- posterior_metadata(mixed$theta, "measure_unavailable")$diagnostics[[1L]]
  subset <- .bt_draws_subset_rows(mixed$theta, c(1, 4, 7))
  expect_identical(posterior_metadata(subset, "measure_unavailable")$diagnostics[[1L]], payload)
  transformed <- posterior_transform(subset, transformation = "lin", transformation_arguments = list(a = 1, b = 2))
  expect_identical(posterior_metadata(transformed, "measure_unavailable")$diagnostics[[1L]], payload)
  malformed <- posterior_metadata(subset, "measure_unavailable")
  malformed$diagnostics[[1L]]$model_indices <- c(1L, 1L)
  expect_error(posterior_metadata(subset, "measure_unavailable") <- malformed, "invalid")
  malformed <- posterior_metadata(subset, "measure_unavailable")
  malformed$diagnostics[[1L]]$log_posterior_probabilities <- c(0, NA_real_)
  expect_error(posterior_metadata(subset, "measure_unavailable") <- malformed, "invalid")
})

test_that("actual model table consumers keep logs hidden and validate their owners", {
  models <- list(.model_probability_test_model(prior("normal", list(0, 1)), 1000, 1e-300),
    .model_probability_test_model(prior("normal", list(0, 2)), 0, 1e100))
  inference <- ensemble_inference(models, "theta", list(c(FALSE, TRUE)))
  table <- ensemble_inference_table(inference, "theta", logBF = TRUE)
  expect_equal(as.numeric(table$inclusion_BF), 1000, tolerance = 1e-12)
  expect_false(any(c("log_prior_probs", "log_post_probs", "model_probability_declaration") %in% names(as.data.frame(table))))
  model_inference <- models_inference(models)
  expect_equal(attr(model_inference[[1L]]$inference, "log_prior_prob"), -400 * log(10), tolerance = 1e-12)
  expect_equal(attr(model_inference[[1L]]$inference, "inclusion_log_BF"), 1000, tolerance = 1e-12)
  broken <- model_inference[[1L]]
  attr(broken$inference, "log_prior_prob") <- -10
  expect_error(model_summary_table(broken), "ownership")
  broken <- inference
  attr(broken$theta, "log_prior_probs") <- NULL
  expect_error(ensemble_inference_table(broken, "theta"), "ownership")
})

test_that("equal declared point laws remain available under extreme model odds", {
  models <- list(.model_probability_test_model(prior("point", list(2)), 0, 1e-300),
    .model_probability_test_model(prior("point", list(2)), 0, 1e100))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_silent(marginal <- marginal_posterior(mixed, "theta", use_formula = FALSE, prior_samples = TRUE))
  expect_identical(posterior_metadata(marginal, "atoms")$mass, 1)
  expect_identical(prior_density_ordinate(posterior_metadata(marginal, "prior_density"), 2)$behavior, "point_mass")
  priors <- attr(mixed$theta, "prior_list", exact = TRUE)
  expect_silent(.plot_data_prior_list.simple(priors, NULL, c(1, 3), NULL, 64, 100, FALSE, FALSE, NULL, NULL, FALSE))
  context <- .prior_density_model_mixture_context(list(theta = priors), "theta")
  context$n_grid <- 0
  condition <- tryCatch(.prior_density_from_context(context, c(theta = 1)), error = identity)
  expect_s3_class(condition, "error")
  expect_false(inherits(condition, "BayesTools_formula_measure_unavailable"))
})

test_that("ordinary context normalization retains its original direct route", {
  priors <- list(prior("normal", list(0, 1), prior_weights = 9),
    prior("normal", list(0, 2), prior_weights = 1))
  context <- .prior_density_model_mixture_context(list(theta = priors), "theta")
  expect_identical(context$model_weights, c(9, 1) / sum(c(9, 1)))
})

test_that("raw huge prior plot weights are normalized before omission", {
  priors <- list(prior("normal", list(0, 1), prior_weights = .Machine$double.xmax),
    prior("normal", list(1, 1), prior_weights = .Machine$double.xmax))
  data <- .plot_data_prior_list.simple(priors, seq(-5, 6, length.out = 128), c(-5, 6), NULL,
    128, 1000, FALSE, FALSE, NULL, NULL, FALSE)
  expect_true(length(data$density$x) > 0L)
  expect_equal(data$density$y, .5 * stats::dnorm(data$density$x) + .5 * stats::dnorm(data$density$x, 1), tolerance = 1e-12)
})

test_that("fixed weightfunction declarations provide complete joint points", {
  first <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  second <- prior_weightfunction("one-sided", .1, wf_fixed(c(1, .25)))
  models <- list(.model_probability_weightfunction_model(first), .model_probability_weightfunction_model(second))
  mixed <- mix_posteriors(models, "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  atoms <- posterior_metadata(mixed$bias, "atoms")
  expect_identical(unname(atoms$locations), rbind(c(1, .5, .5), c(1, 1, .25)))
  expect_identical(colnames(atoms$locations), colnames(mixed$bias))
  expect_equal(atoms$mass, c(.5, .5), tolerance = 1e-12)
  models[[2L]]$marglik <- bridgesampling_object(-1000)
  mixed <- mix_posteriors(models, "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_error(posterior_atoms_free(mixed$bias), class = "BayesTools_formula_measure_unavailable")
  models[[2L]] <- .model_probability_weightfunction_model(first, -1000)
  mixed <- mix_posteriors(models, "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  expect_identical(posterior_metadata(mixed$bias, "atoms")$mass, 1)
  for(prior_object in list(prior_weightfunction("one-sided", .05, wf_cumulative(c(1, 1))),
      prior_weightfunction("one-sided", .05, wf_independent(prior("gamma", list(2, 1)))))){
    models <- list(.model_probability_weightfunction_model(prior_object),
      .model_probability_weightfunction_model(prior_object, -1000))
    mixed <- mix_posteriors(models, "bias", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
    expect_length(posterior_metadata(mixed$bias, "atoms")$mass, 0L)
  }
})

test_that("fixed weightfunction invalid declarations and old model owners are strict", {
  prior_object <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  model <- .model_probability_weightfunction_model(prior_object)
  invalid <- model
  priors <- attr(invalid$fit, "prior_list", exact = TRUE)
  priors$bias$weights$omega[[1L]] <- .999
  attr(invalid$fit, "prior_list") <- priors
  condition <- tryCatch(mix_posteriors(list(invalid), "bias", list(FALSE), seed = 47, n_samples = 20), error = identity)
  expect_s3_class(condition, "error")
  expect_false(inherits(condition, "BayesTools_formula_measure_unavailable"))
  models <- list(.model_probability_test_model(prior("normal", list(0, 1)), 0),
    .model_probability_test_model(prior("normal", list(0, 2)), 0))
  mixed <- mix_posteriors(models, "theta", list(c(FALSE, FALSE)), seed = 47, n_samples = 20)
  old <- mixed$theta
  metadata <- attr(old, "bayestools_meta", exact = TRUE)
  metadata$model_probabilities <- NULL
  attr(old, "bayestools_meta") <- metadata
  expect_error(posterior_atoms_free(old), class = "BayesTools_refit_required")
})

test_that("weightfunction and PET plot consumers stabilize unsafe raw normalization", {
  priors <- list(prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)), prior_weights = .Machine$double.xmax),
    prior_weightfunction("one-sided", .05, wf_fixed(c(1, .25)), prior_weights = .Machine$double.xmax))
  context <- .weightfunction_prior_list_context(priors)
  expect_equal(context$model_weights, c(.5, .5), tolerance = 1e-12)
  bias <- list(prior_PET("point", list(0), prior_weights = .Machine$double.xmax),
    prior_PET("point", list(2), prior_weights = .Machine$double.xmax))
  mu <- list(prior("point", list(0)), prior("point", list(0)))
  deterministic <- .plot_data_prior_list.PETPEESE_deterministic(bias, c(0, 1), 2, NULL, NULL, mu)
  expect_true(all(is.finite(deterministic$y)))
  expect_equal(deterministic$y_lCI, c(0, 0), tolerance = 1e-12)
  expect_equal(deterministic$y_uCI, c(0, 2), tolerance = 1e-12)
  sampled <- .plot_data_prior_list.PETPEESE_sampled(bias, c(0, 1), 2, 1000, NULL, NULL, mu)
  expect_true(all(is.finite(sampled$y)))
})

test_that("PET plot guards certify the joint mapped point target", {
  weights <- c(1e-300, 1e100)
  bias <- list(prior_PET("point", list(2), prior_weights = weights[[1L]]),
    prior_PEESE("point", list(2), prior_weights = weights[[2L]]))
  mu <- list(prior("point", list(0)), prior("point", list(0)))
  expect_error(.plot_data_prior_list.PETPEESE_deterministic(bias, .5, 1, NULL, NULL, mu),
    class = "BayesTools_formula_measure_unavailable")
  expect_error(.plot_data_prior_list.PETPEESE_sampled(bias, .5, 1, 20, NULL, NULL, mu),
    class = "BayesTools_formula_measure_unavailable")
  bias <- list(prior_PET("point", list(0), prior_weights = weights[[1L]]),
    prior_PET("point", list(0), prior_weights = weights[[2L]]))
  mu <- list(prior("point", list(2)), prior("point", list(0)))
  expect_error(.plot_data_prior_list.PETPEESE_deterministic(bias, .5, 1, NULL, NULL, mu),
    class = "BayesTools_formula_measure_unavailable")
  expect_error(.plot_data_prior_list.PETPEESE_sampled(bias, .5, 1, 20, NULL, NULL, mu),
    class = "BayesTools_formula_measure_unavailable")
  same <- list(prior_PET("point", list(2), prior_weights = weights[[1L]]),
    prior_PET("point", list(2), prior_weights = weights[[2L]]))
  mu <- list(prior("point", list(0)), prior("point", list(0)))
  result <- .plot_data_prior_list.PETPEESE_deterministic(same, .5, 1, NULL, NULL, mu)
  expect_equal(result$y, 1, tolerance = 1e-12)
  zero <- list(prior_PET("point", list(0), prior_weights = weights[[1L]]),
    prior_PEESE("point", list(0), prior_weights = weights[[2L]]))
  mu <- list(prior("point", list(7)), prior("point", list(7)))
  result <- .plot_data_prior_list.PETPEESE_deterministic(zero, c(0, .5, 1), 3, NULL, NULL, mu)
  expect_equal(result$y, rep(7, 3), tolerance = 1e-12)
  none <- list(.set_prior_model_weight(prior_none(), weights[[1L]]), zero[[1L]], zero[[2L]])
  mu <- rep(list(prior("point", list(7))), 3)
  result <- .plot_data_prior_list.PETPEESE_deterministic(none, c(0, .5, 1), 3,
    NULL, NULL, mu, effect_direction = "negative")
  expect_equal(result$y, rep(7, 3), tolerance = 1e-12)
  pair <- .prior_model_probability_pair(same)
  owned_mu <- lapply(seq_len(2), function(i){
    declaration <- pair$declaration
    declaration$model_indices <- as.integer(i)
    .set_prior_model_probability(prior("point", list(0)), pair$probabilities[[i]], pair$logs[[i]], declaration)
  })
  expect_silent(.model_probability_petpeese_plot_check(same, owned_mu))
  missing <- owned_mu
  attr(missing[[1L]], "model_log_prior_weights") <- NULL
  expect_error(.model_probability_petpeese_plot_check(same, missing), "ownership")
  old <- owned_mu
  attr(old[[1L]], "model_probability_declaration")$schema_version <- 0L
  expect_error(.model_probability_petpeese_plot_check(same, old), "ownership")
  contradictory <- lapply(seq_len(2), function(i){
    declaration <- .model_probability_pair(.5, log(.5), "prior", "raw")$declaration
    declaration$model_indices <- as.integer(i)
    .set_prior_model_probability(prior("point", list(0)), .5, log(.5), declaration)
  })
  expect_error(.model_probability_petpeese_plot_check(same, contradictory), "not aligned")
  continuous <- lapply(seq_len(2), function(i) prior_PET("normal", list(0, 1), prior_weights = weights[[i]]))
  expect_error(.plot_data_prior_list.PETPEESE_deterministic(continuous, 0, 1, NULL, NULL,
    rep(list(prior("point", list(7))), 2)), class = "BayesTools_formula_measure_unavailable")
})


test_that("single fixed weightfunctions declare one complete joint point", {

  fixed <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  fit <- .model_probability_weightfunction_model(fixed)$fit
  samples <- as_mixed_posteriors(fit, "bias")$bias
  atoms <- posterior_metadata(samples, "atoms")
  expect_true(atoms$joint_declared)
  expect_identical(unname(atoms$locations), matrix(c(1, .5), 1L))
  expect_identical(atoms$mass, 1)
  expect_identical(vapply(atoms$marginals, function(x) x$locations[1L, 1L], numeric(1)),
    stats::setNames(c(1, .5), colnames(samples)))
  continuous <- prior_weightfunction("one-sided", .05, wf_independent(prior("beta", list(2, 2))))
  sibling <- as_mixed_posteriors(.model_probability_weightfunction_model(continuous)$fit, "bias")$bias
  expect_identical(nrow(posterior_metadata(sibling, "atoms")$locations), 0L)
  expect_identical(posterior_metadata(sibling, "atoms")$marginals[[1L]]$mass, 1)
})



test_that("bias-mixture scalar declarations use the retained compiler omega map", {

  for(type in c("two-sided", "one-sided")){
    prior <- prior_mixture(list(prior_none(),
      prior_weightfunction(type, .05, wf_fixed(c(1, .5)))))
    draws <- .generate_prior_sample_matrix(list(bias = prior), 20L, seed = 17)
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), list(bias = prior))
    # Materialize the registered deterministic monitors from the declared
    # fixed weight branches; the generic raw-prior collector omits them.
    draws <- cbind(draws, JAGS_evaluate_deterministic(fit, draws, nodes = "omega"))
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), list(bias = prior))
    samples <- as_mixed_posteriors(fit, "bias", conditional = "omega")$bias
    context <- attr(samples, "omega_context", exact = TRUE)
    atoms <- posterior_metadata(samples, "atoms")
    expect_identical(ncol(samples), length(context$names))
    expect_gt(ncol(samples), 0L)
    expect_identical(names(atoms$marginals), colnames(samples))
    expect_true(all(vapply(atoms$marginals, Negate(is.null), logical(1))))
    reference <- which(context$mapping[[1L]] == 1L)
    expect_true(length(reference) > 0L)
    expect_true(all(vapply(atoms$marginals[reference], function(x)
      identical(as.numeric(x$locations), 1) && identical(x$mass, 1), logical(1))))
    for(i in 1:2){
      plot <- plot_posterior(list(bias = samples), "omega", individual = TRUE,
        prior = FALSE, show_figures = reference[[1L]], plot_type = "ggplot")
      expect_type(plot, "list")
      expect_identical(names(plot), colnames(samples)[reference[[1L]]])
      expect_s3_class(plot[[1L]], "ggplot")
    }
  }
})

