skip_if_not_test_profile(c("unit", "visual"))

# ============================================================================ #
# TEST FILE: Marginal Prior Samplers
# ============================================================================ #
#
# These tests cover the prior-sampling helpers in R/marginal-distributions.R.
# They deliberately avoid JAGS fixtures so covr can exercise the helper
# branches during ordinary package coverage runs.

expect_model_counts <- function(x, expected) {
  expect_equal(
    as.integer(table(attr(x, "models_ind"))),
    as.integer(expected)
  )
}

expect_mixed_prior <- function(x, parameter, extra_class) {
  expect_s3_class(x, "mixed_posteriors")
  expect_s3_class(x, extra_class)
  expect_equal(attr(x, "parameter"), parameter)
  expect_false(is.null(attr(x, "prior_list")))
}

make_treatment_prior <- function(mean = 0, sd = 0.2, prior_weights = 1, levels = 3) {
  x <- prior_factor("normal", list(mean, sd), contrast = "treatment", prior_weights = prior_weights)
  attr(x, "levels") <- levels
  x
}

make_independent_prior <- function(mean = 1, sd = 0.2, prior_weights = 1, levels = 3) {
  x <- prior_factor("normal", list(mean, sd), contrast = "independent", prior_weights = prior_weights)
  attr(x, "levels") <- levels
  x
}

make_meandif_prior <- function(mean = 0, sd = 0.2, prior_weights = 1, levels = 3) {
  x <- prior_factor("mnormal", list(mean = mean, sd = sd), contrast = "meandif", prior_weights = prior_weights)
  attr(x, "levels") <- levels
  x
}

test_that(".mix_priors helpers preserve mixture structure and prior attributes", {

  n_samples <- 400
  null_prior <- prior("spike", list(0), prior_weights = 1)
  normal_prior <- prior("normal", list(2, 0.25), prior_weights = 3)

  simple <- .mix_priors.simple(
    list(null_prior, normal_prior),
    parameter = "theta",
    seed = 11,
    n_samples = n_samples
  )

  expect_mixed_prior(simple, "theta", "mixed_posteriors.simple")
  expect_model_counts(simple, c(100, 300))
  expect_equal(simple[attr(simple, "models_ind") == 1], rep(0, 100))
  expect_equal(mean(simple[attr(simple, "models_ind") == 2]), 2, tolerance = 0.04)
  expect_equal(stats::sd(simple[attr(simple, "models_ind") == 2]), 0.25, tolerance = 0.06)

  set.seed(11)
  expected_simple <- c(
    rng(null_prior, 100, transform_factor_samples = FALSE),
    rng(normal_prior, 300, transform_factor_samples = FALSE)
  )
  expect_equal(as.numeric(simple), as.numeric(expected_simple))

  vector <- .mix_priors.vector(
    list(
      prior("point", list(0), prior_weights = 1),
      prior("mnormal", list(mean = 1, sd = 0.2, K = 2), prior_weights = 3)
    ),
    parameter = "beta",
    seed = 12,
    n_samples = n_samples
  )

  expect_mixed_prior(vector, "beta", "mixed_posteriors.vector")
  expect_equal(dim(vector), c(n_samples, 2))
  expect_equal(colnames(vector), c("beta[1]", "beta[2]"))
  expect_model_counts(vector, c(100, 300))
  expect_equal(unname(vector[attr(vector, "models_ind") == 1, ]), matrix(0, nrow = 100, ncol = 2))
  expect_equal(unname(colMeans(vector[attr(vector, "models_ind") == 2, ])), c(1, 1), tolerance = 0.04)

  treatment <- .mix_priors.factor(
    list(null_prior, make_treatment_prior(mean = 0, sd = 0.2, prior_weights = 3)),
    parameter = "fac",
    seed = 13,
    n_samples = n_samples
  )

  expect_mixed_prior(treatment, "fac", "mixed_posteriors.factor")
  expect_equal(dim(treatment), c(n_samples, 2))
  expect_equal(colnames(treatment), c("fac[2]", "fac[3]"))
  expect_true(isTRUE(attr(treatment, "treatment")))
  expect_equal(attr(treatment, "level_names"), 1:3)
  expect_model_counts(treatment, c(100, 300))
  expect_equal(unname(treatment[attr(treatment, "models_ind") == 1, ]), matrix(0, nrow = 100, ncol = 2))

  independent <- .mix_priors.factor(
    list(null_prior, make_independent_prior(mean = 1, sd = 0.2, prior_weights = 3)),
    parameter = "ind",
    seed = 14,
    n_samples = n_samples
  )

  expect_mixed_prior(independent, "ind", "mixed_posteriors.factor")
  expect_equal(dim(independent), c(n_samples, 3))
  expect_equal(colnames(independent), paste0("ind[", 1:3, "]"))
  expect_true(isTRUE(attr(independent, "independent")))
  expect_model_counts(independent, c(100, 300))
  expect_equal(unname(colMeans(independent[attr(independent, "models_ind") == 2, ])), c(1, 1, 1), tolerance = 0.05)

  meandif <- .mix_priors.factor(
    list(null_prior, make_meandif_prior(mean = 0, sd = 0.2, prior_weights = 3)),
    parameter = "md",
    seed = 15,
    n_samples = n_samples
  )

  expect_mixed_prior(meandif, "md", "mixed_posteriors.factor")
  expect_equal(dim(meandif), c(n_samples, 2))
  expect_true(isTRUE(attr(meandif, "meandif")))
  expect_model_counts(meandif, c(100, 300))
  expect_equal(unname(colMeans(meandif[attr(meandif, "models_ind") == 2, ])), c(0, 0), tolerance = 0.04)

  weightfunction <- .mix_priors.weightfunction(
    list(
      prior_none(prior_weights = 1),
      prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(2, 4)), prior_weights = 3)
    ),
    parameter = "omega",
    seed = 16,
    n_samples = n_samples
  )

  expect_mixed_prior(weightfunction, "omega", "mixed_posteriors.weightfunction")
  expect_equal(dim(weightfunction), c(n_samples, 2))
  expect_equal(colnames(weightfunction), c("omega[0,0.05]", "omega[0.05,1]"))
  expect_model_counts(weightfunction, c(100, 300))
  expect_equal(unname(weightfunction[attr(weightfunction, "models_ind") == 1, ]), matrix(1, nrow = 100, ncol = 2))
  expect_true(all(is.finite(weightfunction)))
  expect_true(all(weightfunction > 0))
})

test_that(".mix_priors dispatches supported prior families", {

  n_samples <- 240
  null_prior <- prior("spike", list(0), prior_weights = 1)

  out <- .mix_priors(
    list(
      theta = list(null_prior, prior("normal", list(2, 0.25), prior_weights = 3)),
      beta  = list(prior("point", list(0), prior_weights = 1), prior("mnormal", list(mean = 1, sd = 0.2, K = 2), prior_weights = 3)),
      fac   = list(null_prior, make_treatment_prior(mean = 0, sd = 0.2, prior_weights = 3)),
      omega = list(prior_none(prior_weights = 1), prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(2, 4)), prior_weights = 3))
    ),
    seed = 21,
    n_samples = n_samples
  )

  expect_equal(names(out), c("theta", "beta", "fac", "omega"))
  expect_mixed_prior(out$theta, "theta", "mixed_posteriors.simple")
  expect_mixed_prior(out$beta, "beta", "mixed_posteriors.vector")
  expect_mixed_prior(out$fac, "fac", "mixed_posteriors.factor")
  expect_mixed_prior(out$omega, "omega", "mixed_posteriors.weightfunction")
  expect_model_counts(out$theta, c(60, 180))
  expect_model_counts(out$beta, c(60, 180))
  expect_model_counts(out$fac, c(60, 180))
  expect_model_counts(out$omega, c(60, 180))
  expect_equal(attr(out$theta, "models_ind"), attr(out$beta, "models_ind"))
  expect_equal(attr(out$theta, "models_ind"), attr(out$fac, "models_ind"))
  expect_equal(attr(out$theta, "models_ind"), attr(out$omega, "models_ind"))

  formula_data <- data.frame(
    x = factor(c("A", "B", "A", "B"), levels = c("A", "B"))
  )
  formula_result <- JAGS_formula(
    formula = ~ x,
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior_factor("normal", list(0, 0.2), contrast = "treatment")
    )
  )
  formula_prior <- formula_result$prior_list$mu_x
  formula_prior[["prior_weights"]] <- 3

  formula_out <- .mix_priors(
    list(mu_x = list(null_prior, formula_prior)),
    seed = 22,
    n_samples = 120
  )

  expect_s3_class(formula_out$mu_x, "mixed_posteriors.formula")
  expect_equal(attr(formula_out$mu_x, "formula_parameter"), "mu")
  expect_equal(attr(formula_out$mu_x, "level_names"), c("A", "B"))

  expect_error(
    .mix_priors(list(bad = list(prior("normal", list(0, 1)), prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))))),
    "unsupported mixture of prior distributions"
  )
})

test_that(".as_mixed_priors helpers evaluate single-model prior families", {

  n_samples <- 300

  simple <- .as_mixed_priors.simple(prior("normal", list(1, 0.2)), "theta", seed = 31, n_samples = n_samples)
  expect_mixed_prior(simple, "theta", "mixed_posteriors.simple")
  expect_equal(length(simple), n_samples)
  expect_equal(mean(simple), 1, tolerance = 0.04)
  expect_equal(stats::sd(simple), 0.2, tolerance = 0.03)

  vector <- .as_mixed_priors.vector(prior("mnormal", list(mean = 1, sd = 0.2, K = 2)), "beta", seed = 32, n_samples = n_samples)
  expect_mixed_prior(vector, "beta", "mixed_posteriors.vector")
  expect_equal(dim(vector), c(n_samples, 2))
  expect_equal(unname(colMeans(vector)), c(1, 1), tolerance = 0.04)

  treatment <- .as_mixed_priors.factor(make_treatment_prior(mean = 0, sd = 0.2), "fac", seed = 33, n_samples = n_samples)
  expect_mixed_prior(treatment, "fac", "mixed_posteriors.factor")
  expect_equal(dim(treatment), c(n_samples, 2))
  expect_true(isTRUE(attr(treatment, "treatment")))

  independent <- .as_mixed_priors.factor(make_independent_prior(mean = 1, sd = 0.2), "ind", seed = 34, n_samples = n_samples)
  expect_mixed_prior(independent, "ind", "mixed_posteriors.factor")
  expect_equal(dim(independent), c(n_samples, 3))
  expect_true(isTRUE(attr(independent, "independent")))
  expect_equal(unname(colMeans(independent)), c(1, 1, 1), tolerance = 0.05)

  meandif <- .as_mixed_priors.factor(make_meandif_prior(mean = 0, sd = 0.2), "md", seed = 35, n_samples = n_samples)
  expect_mixed_prior(meandif, "md", "mixed_posteriors.factor")
  expect_equal(dim(meandif), c(n_samples, 2))
  expect_true(isTRUE(attr(meandif, "meandif")))

  weightfunction <- .as_mixed_priors.weightfunction(
    prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(2, 4))),
    "omega",
    seed = 36,
    n_samples = n_samples
  )
  expect_mixed_prior(weightfunction, "omega", "mixed_posteriors.weightfunction")
  expect_equal(dim(weightfunction), c(n_samples, 2))
  expect_true(all(is.finite(weightfunction)))
  expect_true(all(weightfunction > 0))

  spike_and_slab <- .as_mixed_priors.spike_and_slab(
    prior_spike_and_slab(
      prior("normal", list(2, 0.2)),
      prior_inclusion = prior("point", list(0.4))
    ),
    "slab",
    seed = 37,
    n_samples = n_samples
  )
  expect_mixed_prior(spike_and_slab, "slab", "mixed_posteriors.spike_and_slab")
  expect_lt(abs(mean(attr(spike_and_slab, "models_ind")) - 0.4), 0.10)
  expect_equal(spike_and_slab[attr(spike_and_slab, "models_ind") == 0], rep(0, sum(attr(spike_and_slab, "models_ind") == 0)))
  expect_equal(mean(spike_and_slab[attr(spike_and_slab, "models_ind") == 1]), 2, tolerance = 0.05)

  factor_spike_and_slab_prior <- prior_spike_and_slab(
    make_treatment_prior(mean = 0, sd = 0.2),
    prior_inclusion = prior("point", list(0.5))
  )
  factor_spike_and_slab <- .as_mixed_priors.spike_and_slab(
    factor_spike_and_slab_prior,
    "fac_slab",
    seed = 38,
    n_samples = n_samples
  )
  expect_mixed_prior(factor_spike_and_slab, "fac_slab", "mixed_posteriors.spike_and_slab")
  expect_equal(dim(factor_spike_and_slab), c(n_samples, 2))
  expect_true(all(factor_spike_and_slab[attr(factor_spike_and_slab, "models_ind") == 0, ] == 0))

  mixture <- .as_mixed_priors.mixture(
    prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(3, 0.2), prior_weights = 3)
      ),
      is_null = c(TRUE, FALSE)
    ),
    "mix",
    seed = 39,
    n_samples = n_samples
  )
  expect_mixed_prior(mixture, "mix", "mixed_posteriors.mixture")
  expect_model_counts(mixture, c(75, 225))
  expect_equal(mixture[attr(mixture, "models_ind") == 1], rep(0, 75))
  expect_equal(mean(mixture[attr(mixture, "models_ind") == 2]), 3, tolerance = 0.04)

  factor_mixture <- .as_mixed_priors.mixture(
    prior_mixture(
      list(
        {
          x <- prior_factor("normal", list(-1, 0.2), contrast = "treatment", prior_weights = 1)
          attr(x, "levels") <- 3
          x
        },
        make_treatment_prior(mean = 1, sd = 0.2, prior_weights = 3)
      ),
      is_null = c(FALSE, FALSE)
    ),
    "fac_mix",
    seed = 40,
    n_samples = n_samples
  )
  expect_mixed_prior(factor_mixture, "fac_mix", "mixed_posteriors.mixture")
  expect_equal(dim(factor_mixture), c(n_samples, 2))
  expect_model_counts(factor_mixture, c(75, 225))
  expect_equal(unname(colMeans(factor_mixture[attr(factor_mixture, "models_ind") == 1, ])), c(-1, -1), tolerance = 0.08)
  expect_equal(unname(colMeans(factor_mixture[attr(factor_mixture, "models_ind") == 2, ])), c(1, 1), tolerance = 0.08)
})

test_that(".as_mixed_priors applies conditional sampling rules", {

  spike_and_slab <- prior_spike_and_slab(
    prior("normal", list(2, 0.2)),
    prior_inclusion = prior("point", list(0.4))
  )
  mixture <- prior_mixture(
    list(
      prior("spike", list(0), prior_weights = 1),
      prior("normal", list(3, 0.2), prior_weights = 3)
    ),
    is_null = c(TRUE, FALSE)
  )

  conditional_and <- .as_mixed_priors(
    list(slab = spike_and_slab, mix = mixture),
    seed = 41,
    n_samples = 100,
    conditional = c("slab", "mix"),
    conditional_rule = "AND"
  )

  expect_equal(length(conditional_and$slab), 100)
  expect_true(all(attr(conditional_and$slab, "models_ind") == 1))
  expect_true(all(attr(conditional_and$mix, "models_ind") == 2))
  expect_false(isTRUE(attr(conditional_and, "all_alternative")))

  conditional_or <- .as_mixed_priors(
    list(slab = spike_and_slab, mix = mixture),
    seed = 42,
    n_samples = 100,
    conditional = c("slab", "mix"),
    conditional_rule = "OR"
  )

  expect_equal(length(conditional_or$slab), 100)
  expect_true(all(attr(conditional_or$slab, "models_ind") == 1 | attr(conditional_or$mix, "models_ind") == 2))

  expect_warning(
    nonconditional <- .as_mixed_priors(
      list(theta = prior("normal", list(0, 1))),
      seed = 43,
      n_samples = 20,
      conditional = "theta",
      conditional_rule = "AND"
    ),
    "not a conditional parameter"
  )
  expect_equal(length(nonconditional$theta), 20)
  expect_true(isTRUE(attr(nonconditional, "all_alternative")))

  expect_error(
    .as_mixed_priors(
      list(bad = prior_mixture(
        list(prior("normal", list(0, 1)), prior("normal", list(1, 1))),
        components = c("a", "b")
      )),
      seed = 44,
      n_samples = 20,
      conditional = "bad",
      conditional_rule = "AND"
    ),
    "conditional mixture posterior distributions"
  )
})

test_that("marginal prior sampler visuals show the expected mixture shapes", {
  skip_if_not_visual_tests()

  n_samples <- 2000
  null_prior <- prior("spike", list(0), prior_weights = 1)
  normal_prior <- prior("normal", list(2, 0.25), prior_weights = 3)

  simple <- .mix_priors.simple(list(null_prior, normal_prior), "theta", seed = 51, n_samples = n_samples)
  treatment <- .mix_priors.factor(
    list(null_prior, make_treatment_prior(mean = 0, sd = 0.2, prior_weights = 3)),
    "fac",
    seed = 52,
    n_samples = n_samples
  )
  weightfunction <- .mix_priors.weightfunction(
    list(
      prior_none(prior_weights = 1),
      prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(2, 4)), prior_weights = 3)
    ),
    "omega",
    seed = 53,
    n_samples = n_samples
  )

  vdiffr::expect_doppelganger("marginal-prior-samplers-mixtures", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mfrow = c(2, 3))
    graphics::hist(simple, breaks = 50, freq = FALSE, main = "simple mixture", xlab = "theta")
    graphics::curve(0.75 * stats::dnorm(x, 2, 0.25), add = TRUE, col = 2, lwd = 2)
    graphics::abline(v = 0, col = 4, lwd = 2)

    graphics::hist(treatment[, 1], breaks = 50, freq = FALSE, main = "treatment level 2", xlab = "fac[2]")
    graphics::curve(0.75 * stats::dnorm(x, 0, 0.2), add = TRUE, col = 2, lwd = 2)
    graphics::abline(v = 0, col = 4, lwd = 2)

    graphics::hist(treatment[, 2], breaks = 50, freq = FALSE, main = "treatment level 3", xlab = "fac[3]")
    graphics::curve(0.75 * stats::dnorm(x, 0, 0.2), add = TRUE, col = 2, lwd = 2)
    graphics::abline(v = 0, col = 4, lwd = 2)

    graphics::hist(weightfunction[, 1], breaks = 50, freq = FALSE, main = "omega[0,0.05]", xlab = "")
    graphics::abline(v = 1, col = 4, lwd = 2)

    graphics::hist(weightfunction[, 2], breaks = 50, freq = FALSE, main = "omega[0.05,1]", xlab = "")
    graphics::abline(v = 1, col = 4, lwd = 2)

    graphics::plot(
      treatment[, 1],
      treatment[, 2],
      pch = 16,
      cex = 0.35,
      main = "paired factor draws",
      xlab = "fac[2]",
      ylab = "fac[3]"
    )
    graphics::abline(0, 1, col = 2, lwd = 2)
  })
})

test_that("conditional marginal prior visuals remove the excluded components", {
  skip_if_not_visual_tests()

  spike_and_slab <- prior_spike_and_slab(
    prior("normal", list(2, 0.2)),
    prior_inclusion = prior("point", list(0.4))
  )
  mixture <- prior_mixture(
    list(
      prior("spike", list(0), prior_weights = 1),
      prior("normal", list(3, 0.2), prior_weights = 3)
    ),
    is_null = c(TRUE, FALSE)
  )

  conditional_and <- .as_mixed_priors(
    list(slab = spike_and_slab, mix = mixture),
    seed = 61,
    n_samples = 1000,
    conditional = c("slab", "mix"),
    conditional_rule = "AND"
  )
  conditional_or <- .as_mixed_priors(
    list(slab = spike_and_slab, mix = mixture),
    seed = 62,
    n_samples = 1000,
    conditional = c("slab", "mix"),
    conditional_rule = "OR"
  )

  or_probability <- 1 - (1 - 0.4) * (1 - 3 / 4)
  expect_lt(abs(mean(attr(conditional_or$slab, "models_ind") == 1) - 0.4 / or_probability), 0.05)
  expect_lt(abs(mean(attr(conditional_or$mix, "models_ind") == 2) - (3 / 4) / or_probability), 0.05)

  or_slab_density_weight <- mean(attr(conditional_or$slab, "models_ind") == 1)
  or_mix_density_weight  <- mean(attr(conditional_or$mix, "models_ind") == 2)

  vdiffr::expect_doppelganger("marginal-prior-samplers-conditional", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mfrow = c(2, 2))
    graphics::hist(conditional_and$slab, breaks = 50, freq = FALSE, main = "AND: spike-slab", xlab = "slab")
    graphics::curve(stats::dnorm(x, 2, 0.2), add = TRUE, col = 2, lwd = 2)
    graphics::hist(conditional_and$mix, breaks = 50, freq = FALSE, main = "AND: mixture", xlab = "mix")
    graphics::curve(stats::dnorm(x, 3, 0.2), add = TRUE, col = 2, lwd = 2)
    graphics::hist(conditional_or$slab, breaks = 50, freq = FALSE, main = "OR: spike-slab", xlab = "slab")
    graphics::abline(v = 0, col = 4, lwd = 2)
    graphics::curve(or_slab_density_weight * stats::dnorm(x, 2, 0.2), add = TRUE, col = 2, lwd = 2)
    graphics::hist(conditional_or$mix, breaks = 50, freq = FALSE, main = "OR: mixture", xlab = "mix")
    graphics::abline(v = 0, col = 4, lwd = 2)
    graphics::curve(or_mix_density_weight * stats::dnorm(x, 3, 0.2), add = TRUE, col = 2, lwd = 2)
  })
})
