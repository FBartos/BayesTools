skip_if_not_test_profile(c("unit", "fixture"))

# TEST FILE: Selection kernel backend priors
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("prior_phacking validates geometry, form, alpha, and stores constants", {

  ph <- prior_phacking(
    target      = .025,
    source      = .25,
    destination = .005,
    form        = "linear",
    alpha       = prior("beta", list(2, 3)),
    prior_weights = 2
  )

  expect_true(is_prior_phacking(ph))
  expect_equal(ph$distribution, "phacking")
  expect_equal(ph$side, "one-sided")
  expect_equal(ph$target, .025)
  expect_equal(ph$source, .25)
  expect_equal(ph$destination, .005)
  expect_equal(ph$form, "linear")
  expect_equal(ph$q, 1L)
  expect_equal(ph$prior_weights, 2)

  expect_error(prior_phacking(destination = .05), "destination < target < source")
  expect_error(prior_phacking(source = .01), "destination < target < source")
  expect_error(prior_phacking(form = "exponential"), "should be one of")
  expect_error(prior_phacking(side = "two-sided"), "one-sided")
  expect_error(prior_phacking(alpha = prior("normal", list(0, 1))), "support within")
  expect_error(prior_phacking(alpha = prior("point", list(1))), "\\[0, 1\\)")
})

test_that("prior_bias validates composition objects", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking  <- prior_phacking()
  combined  <- prior_bias(selection = selection, phacking = phacking, prior_weights = 3)

  expect_true(is_prior_bias(combined))
  expect_identical(combined$selection, selection)
  expect_identical(combined$phacking, phacking)
  expect_equal(combined$prior_weights, 3)

  expect_error(prior_bias(), "At least one")
  expect_error(prior_bias(selection = phacking), "weightfunction")
  expect_error(prior_bias(phacking = selection), "p-hacking")
  expect_error(
    prior_bias(
      selection = prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .5))),
      phacking = phacking
    ),
    "same side"
  )
})

test_that("p-hacking null calibration matches source interval moments", {

  target <- .025
  source <- .25
  destination <- .005
  alpha <- c(.1, .4, .8)

  z_a <- stats::qnorm(1 - source)
  z_b <- stats::qnorm(1 - target)
  source_linear <- stats::integrate(
    function(z) ((z - z_a) / (z_b - z_a)) * stats::dnorm(z),
    lower = z_a,
    upper = z_b
  )$value
  source_quadratic <- stats::integrate(
    function(z) ((z - z_a) / (z_b - z_a))^2 * stats::dnorm(z),
    lower = z_a,
    upper = z_b
  )$value

  expect_equal(
    phack_pi_null(alpha, "linear", source, destination, target = target),
    alpha * source_linear,
    tolerance = 1e-12
  )
  expect_equal(
    phack_pi_null(alpha, "quadratic", source, destination, target = target),
    alpha * source_quadratic,
    tolerance = 1e-12
  )
  expect_equal(
    phack_alpha_from_pi_null(
      phack_pi_null(alpha, "linear", source, destination, target = target),
      "linear", source, destination, target = target
    ),
    alpha,
    tolerance = 1e-12
  )
  expect_equal(
    phack_alpha_from_pi_null(
      phack_pi_null(alpha, "quadratic", source, destination, target = target),
      "quadratic", source, destination, target = target
    ),
    alpha,
    tolerance = 1e-12
  )
  expect_error(
    phack_pi_null(1, "linear", source, destination, target = target),
    "lower than 1"
  )
  expect_error(
    phack_alpha_from_pi_null(1, "linear", source, destination, target = target),
    "too large"
  )
  expect_error(
    phack_alpha_from_pi_null(
      phack_backend_constants("linear", source, destination, target = target)$pi_null_per_alpha,
      "linear", source, destination, target = target
    ),
    "too large"
  )
})

test_that("phack_backend_constants returns z-domain power constants", {

  constants <- phack_backend_constants("quadratic", source = .25, destination = .005, target = .025)

  expect_equal(constants$form, "quadratic")
  expect_equal(constants$q, 2L)
  expect_equal(constants$phack_kind, 2L)
  expect_equal(constants$z_source, stats::qnorm(1 - c(.25, .025)))
  expect_equal(constants$z_destination, stats::qnorm(1 - c(.025, .005)))
  expect_true(constants$source_null_mass > 0)
  expect_true(constants$destination_null_mass > 0)
  expect_true(constants$beta_null_per_alpha > 0)
})

test_that("p-hacking source depletion and destination inflation preserve null mass", {

  geometries <- list(
    c(source = .25, destination = .005, target = .025),
    c(source = .10, destination = .001, target = .01)
  )

  for(form in c("linear", "quadratic")){
    for(geometry in geometries){
      constants <- phack_backend_constants(
        form,
        source      = geometry[["source"]],
        destination = geometry[["destination"]],
        target      = geometry[["target"]]
      )
      alpha <- .73
      beta_null <- alpha * constants$beta_null_per_alpha

      source_integral <- stats::integrate(
        function(z){
          ((z - constants$z_source[1]) / diff(constants$z_source))^constants$q * stats::dnorm(z)
        },
        lower = constants$z_source[1],
        upper = constants$z_source[2]
      )$value
      destination_integral <- stats::integrate(
        function(z){
          ((constants$z_destination[2] - z) / diff(constants$z_destination))^constants$q * stats::dnorm(z)
        },
        lower = constants$z_destination[1],
        upper = constants$z_destination[2]
      )$value

      expect_equal(constants$source_null_mass, constants$pi_null_per_alpha, tolerance = 1e-12)
      expect_equal(source_integral, constants$source_null_mass, tolerance = 1e-12)
      expect_equal(destination_integral, constants$destination_null_mass, tolerance = 1e-12)
      expect_equal(
        constants$source_null_mass,
        constants$beta_null_per_alpha * constants$destination_null_mass,
        tolerance = 1e-12
      )
      expect_equal(
        1 - alpha * constants$source_null_mass + beta_null * constants$destination_null_mass,
        1,
        tolerance = 1e-12
      )
    }
  }
})

test_that("selection_backend_spec compiles none, step, phack, and combined priors", {

  selection <- prior_weightfunction("one-sided", c(.025, .05), wf_fixed(c(1, .5, .25)))
  phacking <- prior_phacking(form = "linear")
  combined <- prior_bias(selection, phacking)

  none_spec <- selection_backend_spec(prior_none())
  expect_equal(none_spec$mode, "none")
  expect_equal(none_spec$step$breaks, c(0, 1))
  expect_equal(none_spec$prior_code, "")

  step_spec <- selection_backend_spec(selection)
  expect_equal(step_spec$mode, "step")
  expect_equal(step_spec$step$breaks, c(0, .025, .05, 1))
  expect_equal(step_spec$step$coefficient_ids, paste0("omega[", 1:3, "]"))
  expect_match(step_spec$prior_code, "omega\\[2\\] <- 0.5")
  expect_false(grepl("alpha", step_spec$prior_code, fixed = TRUE))
  expect_false(grepl("phack_z_source\\[1\\]", step_spec$transform_code, fixed = TRUE))

  two_sided <- prior_weightfunction(
    "two-sided", c(.05, .10), wf_fixed(c(1, .5, .25))
  )
  two_sided_spec <- selection_backend_spec(two_sided)
  expect_equal(two_sided_spec$mode, "step")
  expect_equal(two_sided_spec$step$breaks, c(0, .025, .05, .95, .975, 1))
  expect_equal(two_sided_spec$step$coefficient_ids, paste0("omega[", 1:5, "]"))
  expect_match(two_sided_spec$prior_code, "omega_local\\[2\\] <- 0.5")
  expect_match(two_sided_spec$prior_code, "omega\\[4\\] <- omega_local\\[2\\]")
  expect_match(two_sided_spec$prior_code, "omega\\[5\\] <- omega_local\\[1\\]")

  phack_spec <- selection_backend_spec(phacking)
  expect_equal(phack_spec$mode, "phack_power")
  expect_equal(phack_spec$step$breaks, c(0, 1))
  expect_equal(phack_spec$phacking$form, "linear")
  expect_equal(phack_spec$phacking$q, 1L)
  expect_equal(
    phack_spec$phacking$branch_beta_null_per_alpha,
    phack_backend_constants("linear", .25, .005, target = .025)$beta_null_per_alpha
  )
  expect_match(phack_spec$prior_code, "alpha ~ dbeta\\(1,1\\)")
  expect_match(phack_spec$prior_code, "phack_kind <- 1")
  expect_match(phack_spec$prior_code, "beta_null <- alpha \\*")
  expect_match(phack_spec$prior_code, "omega\\[1\\] <- 1")
  expect_true(all(c("omega", "alpha", "phack_kind", "pi_null") %in% phack_spec$monitor))
  expect_false("phack_z_source" %in% names(phack_spec$data))
  expect_true("phack_component_beta_null_per_alpha" %in% names(phack_spec$data))
  expect_equal(phack_spec$transform_code, "")

  combined_spec <- selection_backend_spec(combined)
  expect_equal(combined_spec$mode, "step_phack_power")
  expect_equal(combined_spec$branch_type, "combined")
  expect_equal(combined_spec$step$breaks, c(0, .025, .05, 1))
  expect_match(combined_spec$prior_code, "omega\\[3\\] <- 0.25")
  expect_match(combined_spec$prior_code, "alpha ~ dbeta\\(1,1\\)")
  expect_match(combined_spec$prior_code, "phack_kind <- 1")
})

test_that("selection_backend_spec rejects malformed global breaks", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))

  expect_error(selection_backend_spec(selection, global_breaks = c(0)), "at least 0 and 1")
  expect_error(selection_backend_spec(selection, global_breaks = c(0, .025, .025, 1)), "duplicate")
  expect_error(selection_backend_spec(selection, global_breaks = c(0, .05, .025, 1)), "monotonically")
  expect_error(selection_backend_spec(selection, global_breaks = c(.001, .025, 1)), "start at 0 and end at 1")
  expect_error(selection_backend_spec(selection, global_breaks = c(0, .02, 1)), "contain all step-selection")

  spec <- selection_backend_spec(selection, global_breaks = c(0, .025, .50, 1))
  expect_equal(spec$step$breaks, c(0, .025, .50, 1))
  expect_equal(spec$step$coefficient_ids, paste0("omega[", 1:3, "]"))
})

test_that("selection initial omega expands two-sided weights onto one-sided global bins", {

  selection <- prior_weightfunction("two-sided", c(.05, .10), wf_fixed(c(1, .2, .5)))

  expect_equal(
    BayesTools:::.selection_initial_omega(selection, c(0, .025, .05, .95, .975, 1)),
    c(1, .2, .5, .2, 1)
  )
})

test_that("selection prior means handle point, cumulative, and log-scale components", {

  cumulative <- prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)))
  log_scale <- prior_weightfunction(
    "one-sided",
    c(.025),
    wf_independent(prior("normal", list(0, .2)), scale = "log_omega")
  )

  expect_equal(BayesTools:::.selection_weightfunction_mean(cumulative), c(1, 5 / 6, 1 / 2))
  expect_equal(BayesTools:::.selection_weightfunction_mean(log_scale), c(1, exp(.2^2 / 2)), tolerance = 1e-8)
  expect_equal(BayesTools:::.selection_prior_mean(prior("point", list(3))), 3)
  expect_equal(BayesTools:::.selection_prior_mean(prior("beta", list(2, 6))), 2 / 8)
})

test_that("selection_backend_spec compiles mixtures with active identity transforms", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking_linear <- prior_phacking(form = "linear")
  phacking_quadratic <- prior_phacking(form = "quadratic")
  combined <- prior_bias(selection, phacking_quadratic)

  bias <- prior_mixture(list(
    prior_none(),
    selection,
    phacking_linear,
    combined
  ))

  spec <- selection_backend_spec(bias)

  expect_equal(spec$mode, "step_phack_power")
  expect_equal(spec$branch_type, c("none", "weightfunction", "phack", "combined"))
  expect_equal(spec$step$breaks, c(0, .025, 1))
  expect_equal(spec$phacking$branch_phack_kind, c(1L, 2L))
  expect_match(spec$prior_code, "bias_indicator ~ dcat\\(c\\(1, 1, 1, 1\\)\\)")
  expect_match(spec$prior_code, "omega_component_1\\[2\\] <- 1")
  expect_match(spec$prior_code, "omega_component_3\\[2\\] <- 1")
  expect_match(spec$prior_code, "alpha_component_1 <- 0")
  expect_match(spec$prior_code, "alpha_component_2 <- 0")
  expect_match(spec$prior_code, "phack_kind_component_1 <- 0")
  expect_match(spec$prior_code, "beta_null_component_1 <- 0")
  expect_match(spec$prior_code, "phack_kind_component_2 <- 0")
  expect_match(spec$prior_code, "phack_kind_component_3 <- 1")
  expect_match(spec$prior_code, "beta_null_component_3 <- alpha_component_3 \\*")
  expect_match(spec$prior_code, "phack_kind_component_4 <- 2")
  expect_match(spec$transform_code, "omega\\[2\\] <- omega_component_1\\[2\\] \\* equals\\(bias_indicator, 1\\)")
  expect_match(spec$transform_code, "alpha <- alpha_component_1 \\* equals\\(bias_indicator, 1\\)")
  expect_match(spec$transform_code, "beta_null <- beta_null_component_1 \\* equals\\(bias_indicator, 1\\)")
  expect_match(spec$transform_code, "phack_kind <- phack_kind_component_1 \\* equals\\(bias_indicator, 1\\)")
})

test_that("JAGS_add_priors emits p-hacking active and inactive identity branches", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear")
  bias <- prior_mixture(list(prior_none(), selection, phacking, prior_bias(selection, phacking)))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "omega_component_1\\[1\\] <- 1")
  expect_match(syntax, "omega_component_3\\[2\\] <- 1")
  expect_match(syntax, "alpha_component_1 <- 0")
  expect_match(syntax, "alpha_component_2 <- 0")
  expect_match(syntax, "phack_kind_component_1 <- 0")
  expect_match(syntax, "phack_kind_component_2 <- 0")
  expect_match(syntax, "alpha_component_3 ~ dbeta\\(1,1\\)")
  expect_match(syntax, "phack_kind_component_3 <- 1")
  expect_match(syntax, "alpha <- alpha_component_1 \\* equals\\(bias_indicator, 1\\)")
  expect_match(syntax, "phack_kind <- phack_kind_component_1 \\* equals\\(bias_indicator, 1\\)")

  expect_true(all(c("bias_indicator", "omega", "alpha", "phack_kind", "pi_null") %in% JAGS_to_monitor(list(bias = bias))))
})

test_that("phacking-only mixtures still emit active omega identity", {

  bias <- prior_mixture(list(
    prior_none(),
    prior_phacking(form = "linear")
  ))

  spec <- selection_backend_spec(bias)
  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_equal(spec$mode, "phack_power")
  expect_equal(spec$step$n_bins, 1L)
  expect_match(spec$prior_code, "omega_component_1\\[1\\] <- 1")
  expect_match(spec$transform_code, "omega\\[1\\] <- omega_component_1\\[1\\] \\* equals\\(bias_indicator, 1\\)")
  expect_match(syntax, "omega\\[1\\] <- omega_component_1\\[1\\] \\* equals\\(bias_indicator, 1\\)")
  expect_true(all(c("bias_indicator", "omega", "alpha", "phack_kind", "pi_null") %in% JAGS_to_monitor(list(bias = bias))))
})

test_that("direct selection-family JAGS helpers use active public names", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  phacking  <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias      <- prior_bias(selection, phacking)

  syntax <- JAGS_add_priors("model{}", list(bias = bias))
  expect_match(syntax, "eta\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax, "alpha ~ dbeta\\(2,3\\)")
  expect_false(grepl("_component_", syntax, fixed = TRUE))

  inits <- JAGS_get_inits(list(bias = bias), chains = 1, seed = 1)[[1]]
  expect_true(all(c("eta", "alpha") %in% names(inits)))
  expect_false(any(grepl("_component_", names(inits), fixed = TRUE)))

  monitor <- JAGS_to_monitor(list(bias = bias))
  expect_true(all(c("omega", "eta", "alpha", "phack_kind", "pi_null") %in% monitor))
  expect_false(any(grepl("_component_", monitor, fixed = TRUE)))

  phacking_syntax <- JAGS_add_priors("model{}", list(phacking = phacking))
  expect_match(phacking_syntax, "alpha ~ dbeta\\(2,3\\)")
  expect_false(grepl("_component_", phacking_syntax, fixed = TRUE))

  phacking_inits <- JAGS_get_inits(list(phacking = phacking), chains = 1, seed = 1)[[1]]
  expect_true("alpha" %in% names(phacking_inits))
  expect_false(any(grepl("_component_", names(phacking_inits), fixed = TRUE)))
})

test_that("bias posterior extraction recognizes composed selection and phacking branches", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear")
  bias <- prior_mixture(
    list(prior_none(), selection, phacking, prior_bias(selection, phacking)),
    is_null = c(TRUE, FALSE, FALSE, FALSE)
  )

  model <- matrix(
    c(
      1, 1.0, 1.0, 0.0, 0.000, 0,
      2, 1.0, 0.5, 0.0, 0.000, 0,
      3, 1.0, 1.0, 0.2, 0.010, 1,
      4, 1.0, 0.5, 0.4, 0.020, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(model) <- c("bias_indicator", "omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = bias)

  mixed_all <- as_mixed_posteriors(model, parameters = "bias")
  expect_equal(colnames(mixed_all$bias), c("omega[0,0.025]", "omega[0.025,1]", "pi_null"))

  mixed_omega <- as_mixed_posteriors(model, parameters = "bias", conditional = "omega")
  expect_equal(colnames(mixed_omega$bias), c("omega[0,0.025]", "omega[0.025,1]"))
  expect_equal(attr(mixed_omega$bias, "models_ind"), c(2, 4))

  mixed_phacking <- as_mixed_posteriors(model, parameters = "bias", conditional = "phacking")
  expect_equal(colnames(mixed_phacking$bias), "pi_null")
  expect_equal(attr(mixed_phacking$bias, "models_ind"), c(3, 4))

  mixed_bias <- as_mixed_posteriors(model, parameters = "bias", conditional = "bias")
  conditioned_prior_weights <- sapply(attr(mixed_bias$bias, "prior_list"), function(prior) prior$prior_weights)
  expect_equal(conditioned_prior_weights, c(0, 1, 1, 1))
})

test_that("direct p-hacking and bias priors unpack mixed posteriors", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear")
  bias <- prior_bias(selection, phacking)

  ph_model <- matrix(
    c(
      1.0, .2, .010, 1,
      1.0, .4, .020, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(ph_model) <- c("omega[1]", "alpha", "pi_null", "phack_kind")
  class(ph_model) <- c("matrix", "BayesTools_fit")
  attr(ph_model, "prior_list") <- list(ph = phacking)

  mixed_phacking <- as_mixed_posteriors(ph_model, parameters = "ph")
  expect_equal(colnames(mixed_phacking$ph), "pi_null")

  bias_model <- matrix(
    c(
      1.0, .5, .2, .010, 1,
      1.0, .5, .4, .020, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(bias_model) <- c("omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")
  class(bias_model) <- c("matrix", "BayesTools_fit")
  attr(bias_model, "prior_list") <- list(pub_bias = bias)

  mixed_bias <- as_mixed_posteriors(bias_model, parameters = "pub_bias")
  expect_equal(colnames(mixed_bias$pub_bias), c("omega[0,0.025]", "omega[0.025,1]", "pi_null"))

  mixed_bias_omega <- as_mixed_posteriors(bias_model, parameters = "pub_bias", conditional = "omega")
  expect_equal(colnames(mixed_bias_omega$pub_bias), c("omega[0,0.025]", "omega[0.025,1]"))

  mixed_bias_phacking <- as_mixed_posteriors(bias_model, parameters = "pub_bias", conditional = "phacking")
  expect_equal(colnames(mixed_bias_phacking$pub_bias), "pi_null")
})


test_that("p-hacking report_scale controls public summary coordinates", {
  ph_alpha <- prior_phacking(form = "linear", report_scale = "alpha")
  ph_pi <- prior_phacking(form = "linear", report_scale = "pi_null")

  model <- matrix(
    c(
      .2, .010,
      .4, .020
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(model) <- c("alpha", "pi_null")
  class(model) <- c("matrix", "BayesTools_fit")

  attr(model, "prior_list") <- list(ph = ph_alpha)
  expect_equal(colnames(as_mixed_posteriors(model, parameters = "ph")$ph), "alpha")
  expect_equal(colnames(.as_mixed_priors(list(ph = ph_alpha), seed = 1, n_samples = 10)$ph), "alpha")
  expect_match(print(ph_alpha, silent = TRUE), "^alpha\\[phacking:")

  attr(model, "prior_list") <- list(ph = ph_pi)
  expect_equal(colnames(as_mixed_posteriors(model, parameters = "ph")$ph), "pi_null")
  expect_equal(colnames(.as_mixed_priors(list(ph = ph_pi), seed = 1, n_samples = 10)$ph), "pi_null")
  expect_match(print(ph_pi, silent = TRUE), "^pi_null\\[phacking:")
})

test_that("bias priors sample direct and mixture prior draws", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias <- prior_bias(selection, phacking)

  set.seed(19)
  direct <- rng(bias, 20)
  expect_equal(colnames(direct), c("omega[1]", "omega[2]", "alpha", "pi_null"))
  expect_equal(direct[, "omega[1]"], rep(1, 20))
  expect_equal(direct[, "omega[2]"], rep(.5, 20))
  expect_true(all(direct[, "alpha"] > 0 & direct[, "alpha"] < 1))
  expect_true(all(direct[, "pi_null"] >= 0))

  mixed_direct <- .as_mixed_priors(list(pub_bias = bias), seed = 20, n_samples = 20)
  expect_s3_class(mixed_direct$pub_bias, "mixed_posteriors.bias")
  expect_equal(colnames(mixed_direct$pub_bias), c("omega[0,0.025]", "omega[0.025,1]", "pi_null"))

  phacking_only <- prior_bias(phacking = phacking)
  expect_equal(colnames(rng(phacking_only, 10)), c("alpha", "pi_null"))
  mixed_phacking <- .as_mixed_priors(list(pub_bias = phacking_only), seed = 21, n_samples = 10)
  expect_equal(colnames(mixed_phacking$pub_bias), "pi_null")

  bias_mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_bias(selection = selection, prior_weights = 1),
    prior_bias(phacking = phacking, prior_weights = 1),
    prior_PET("normal", list(0, .2), truncation = list(-Inf, Inf), prior_weights = 1)
  ))
  mixture_rng <- rng(bias_mixture, 100)
  expect_equal(colnames(mixture_rng), c("omega[1]", "omega[2]", "alpha", "pi_null", "PET"))
  expect_setequal(attr(mixture_rng, "components"), 1:4)
  expect_true(all(mixture_rng[attr(mixture_rng, "components") == 1, c("omega[1]", "omega[2]")] == 1))
  expect_true(all(mixture_rng[attr(mixture_rng, "components") == 1, c("alpha", "pi_null", "PET")] == 0))

  mixed <- .as_mixed_priors(list(pub_bias = bias_mixture), seed = 22, n_samples = 40)
  expect_s3_class(mixed$pub_bias, "mixed_posteriors.bias")
  expect_s3_class(mixed$pub_bias, "mixed_posteriors.mixture")
  expect_equal(colnames(mixed$pub_bias), c("omega[0,0.025]", "omega[0.025,1]", "pi_null", "PET"))
  expect_setequal(attr(mixed$pub_bias, "models_ind"), 1:4)
  expect_true(all(mixed$pub_bias[attr(mixed$pub_bias, "models_ind") == 1, c("omega[0,0.025]", "omega[0.025,1]")] == 1))
  expect_true(all(mixed$pub_bias[attr(mixed$pub_bias, "models_ind") == 1, c("pi_null", "PET")] == 0))
  expect_equal(
    mixed$pub_bias[attr(mixed$pub_bias, "models_ind") == 2, "omega[0.025,1]"],
    rep(.5, sum(attr(mixed$pub_bias, "models_ind") == 2))
  )
  expect_true(all(mixed$pub_bias[attr(mixed$pub_bias, "models_ind") == 3, "pi_null"] > 0))
  expect_true(all(mixed$pub_bias[attr(mixed$pub_bias, "models_ind") != 4, "PET"] == 0))
})

test_that("selection-family prior generics fail explicitly when scalar semantics are ambiguous", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear")
  bias <- prior_bias(selection, phacking)

  expect_error(cdf(phacking, .5), "p-hacking priors")
  expect_error(ccdf(phacking, .5), "p-hacking priors")
  expect_error(lpdf(phacking, .5), "p-hacking priors")
  expect_error(pdf(phacking, .5), "p-hacking priors")
  expect_error(quant(phacking, .5), "p-hacking priors")
  expect_error(mean(phacking), "p-hacking priors")
  expect_error(var(phacking), "p-hacking priors")
  expect_error(sd(phacking), "p-hacking priors")
  expect_error(density(phacking), "p-hacking priors")
  expect_error(range(phacking), "p-hacking priors")

  expect_error(cdf(bias, .5), "composed bias priors")
  expect_error(ccdf(bias, .5), "composed bias priors")
  expect_error(lpdf(bias, .5), "composed bias priors")
  expect_error(pdf(bias, .5), "composed bias priors")
  expect_error(quant(bias, .5), "composed bias priors")
  expect_error(mean(bias), "composed bias priors")
  expect_error(var(bias), "composed bias priors")
  expect_error(sd(bias), "composed bias priors")
  expect_error(density(bias), "composed bias priors")
  expect_error(range(bias), "composed bias priors")
})

test_that("diagnostics density handles p-hacking and composed bias priors", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear")
  bias <- prior_mixture(
    list(prior_none(), phacking, prior_bias(selection, phacking)),
    is_null = c(TRUE, FALSE, FALSE)
  )

  alpha_data <- matrix(seq(.05, .95, length.out = 30), ncol = 1)
  colnames(alpha_data) <- "alpha"
  attr(alpha_data, "chain") <- rep(1, nrow(alpha_data))
  attr(alpha_data, "parameter") <- "alpha"
  attr(alpha_data, "prior") <- bias

  omega_data <- matrix(seq(.05, .95, length.out = 30), ncol = 1)
  colnames(omega_data) <- "omega[0.025,1]"
  attr(omega_data, "chain") <- rep(1, nrow(omega_data))
  attr(omega_data, "parameter") <- "omega"
  attr(omega_data, "prior") <- bias

  direct_data <- matrix(seq(.001, .01, length.out = 30), ncol = 1)
  colnames(direct_data) <- "pi_null"
  attr(direct_data, "chain") <- rep(1, nrow(direct_data))
  attr(direct_data, "parameter") <- "phacking"
  attr(direct_data, "prior") <- phacking

  expect_silent(.diagnostics_plot_data_density(alpha_data, n_points = 32, xlim = NULL))
  expect_silent(.diagnostics_plot_data_density(omega_data, n_points = 32, xlim = NULL))
  expect_silent(.diagnostics_plot_data_density(direct_data, n_points = 32, xlim = NULL))
})


test_that("diagnostic prior bounds use weightfunction relative-weight support", {
  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, 1.5)))
  bias <- prior_bias(selection = selection)

  expect_equal(.diagnostics_prior_bounds(selection, "omega"), list(lower = 0, upper = 1.5))
  expect_equal(.diagnostics_prior_bounds(bias, "omega"), list(lower = 0, upper = 1.5))
})

test_that("summary tables handle ordinary mixtures next to selection kernels", {

  skip_if_not_test_profile("fixture")
  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_selection_kernel_summary")

  fit <- readRDS(file.path(temp_fits_dir, "fit_selection_kernel_summary.RDS"))

  table <- suppressWarnings(runjags_estimates_table(fit, remove_diagnostics = TRUE))
  expect_true(all(c("mu (inclusion)", "bias (inclusion)", "pi_null") %in% rownames(table)))
})
