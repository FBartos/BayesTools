skip_if_not_test_profile("unit")

test_that("homogeneous diagonal random effects have a shared SD label", {

  compiled <- JAGS_formula(
    ~ 1 + diag(1 + x | study, hom = TRUE), "mu",
    data.frame(study = factor(c("a", "a", "b", "b")), x = 1:4),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("gamma", list(2, 2)))
  )
  term <- compiled$formula_design$random_effects[[1L]]
  expect_identical(
    BayesTools:::.bt_random_effect_summary_sd_components(
      term, unique(term$sd_parameter_names)
    ),
    "shared"
  )
})

.random_effects_mean_variance_allocation_fit <- function(alpha = c(2, 3),
                                                          scale = "mean_variance",
                                                          sd = prior("gamma", list(2, 2))){

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    x = c(-1, 0, 1, 2)
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = scale,
        sd = sd,
        weights = prior("dirichlet", list(alpha = alpha))
      )
    )
  )
  allocation <- formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1]]
  samples <- matrix(
    c(
      2, 0.25, 0.75,
      2, 0.75, 0.25
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        allocation$source_node,
        paste0(allocation$weight_name, "[1]"),
        paste0(allocation$weight_name, "[2]")
      )
    )
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  attach_test_parameter_map(fit)
}

.random_effects_total_variance_allocation_fit <- function(alpha = c(2, 3)){

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = c("study", "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = alpha))
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1]]
  samples <- matrix(
    c(
      2, 0.25, 0.75,
      2, 0.75, 0.25
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        allocation$source_node,
        paste0(allocation$weight_name, "[1]"),
        paste0(allocation$weight_name, "[2]")
      )
    )
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  attach_test_parameter_map(fit)
}

.random_effects_gated_total_variance_allocation_fit <- function(){

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "allocation",
        terms = c(study = "study", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3))),
        inclusion = list(
          study = prior("spike", list(location = 0.5)),
          drug = prior("spike", list(location = 0.5))
        )
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1L]]
  indicators <- vapply(
    allocation$inclusion,
    `[[`,
    character(1),
    "indicator_name"
  )
  samples <- cbind(
    rep(2, 4),
    c(.25, .25, .25, .25),
    c(.75, .75, .75, .75),
    c(0, 1, 0, 1),
    c(0, 0, 1, 1)
  )
  colnames(samples) <- c(
    allocation$source_node,
    paste0(allocation$weight_name, "[1]"),
    paste0(allocation$weight_name, "[2]"),
    indicators
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  attach_test_parameter_map(fit)
}

test_that("random-effect summary posterior extracts mean-variance multipliers", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  multipliers <- random_effects_summary_posterior(fit, summary = "var_mult")
  multiplier_name <- "(mu) allocation: var_mult(x)"

  expect_s3_class(multipliers, "mixed_posteriors")
  expect_true(multiplier_name %in% names(multipliers))
  expect_s3_class(multipliers[[multiplier_name]], "marginal_posterior")
  expect_equal(unname(as.numeric(multipliers[[multiplier_name]])), c(1.5, 0.5), tolerance = 1e-12)
  multiplier_selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    multiplier_name,
    namespace = "mu"
  )
  expect_identical(
    parameter_transform(fit, multiplier_selection),
    list(type = "affine", offset = 0, scale = 2)
  )
  supplied_samples <- as.matrix(fit$mcmc[[1L]])
  supplied_samples[, grep("weight\\[1\\]$", colnames(supplied_samples))] <-
    c(.1, .2)
  supplied_samples[, grep("weight\\[2\\]$", colnames(supplied_samples))] <-
    c(.9, .8)
  supplied_draws <- parameter_draws(
    fit,
    multiplier_selection,
    model_samples = supplied_samples
  )
  expect_equal(
    unname(as.numeric(supplied_draws[[1L]][, 1L])),
    c(1.8, 1.6),
    tolerance = 1e-12
  )

  prior_density <- .bt_meta_get(multipliers[[multiplier_name]], "prior_density")
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(posterior_metadata(multipliers[[multiplier_name]], "support")$bounds, c(0, 2))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(0.5, 3, 2) / 2,
    tolerance = 1e-8
  )

  intercept_multiplier <- random_effects_summary_posterior(
    fit,
    summary = "var_mult",
    component = "intercept"
  )
  expect_equal(names(intercept_multiplier), "(mu) allocation: var_mult(intercept)")

  expect_error(
    random_effects_summary_posterior(fit, summary = "var_prop"),
    "Mean-variance SD-component allocations are returned by summary = \"var_mult\"",
    fixed = TRUE
  )
  expect_s3_class(
    plot_posterior(multipliers, multiplier_name, prior = TRUE, plot_type = "ggplot"),
    "ggplot"
  )
})

test_that("random-effect summary posterior extracts SD multipliers", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  multipliers <- random_effects_summary_posterior(fit, summary = "sd_mult")
  multiplier_name <- "(mu) allocation: sd_mult(x)"

  expect_true(multiplier_name %in% names(multipliers))
  expect_equal(
    unname(as.numeric(multipliers[[multiplier_name]])),
    sqrt(c(1.5, 0.5)),
    tolerance = 1e-12
  )

  prior_density <- .bt_meta_get(multipliers[[multiplier_name]], "prior_density")
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(posterior_metadata(multipliers[[multiplier_name]], "support")$bounds, c(0, sqrt(2)))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(0.5, 3, 2),
    tolerance = 1e-8
  )
})

test_that("full estimates summaries add SD multipliers to standard quantities", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  standard <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    return_samples = TRUE
  )
  full <- JAGS_estimates_table(
    fit,
    random_effects_summary = "full",
    return_samples = TRUE
  )

  expect_setequal(
    colnames(standard),
    c(
      "(mu) allocation: sd_common",
      "(mu) allocation: var_mult(intercept)",
      "(mu) allocation: var_mult(x)"
    )
  )
  expect_false(any(grepl(": sd_mult\\(", colnames(standard))))
  expect_false(any(grepl("(^|: )sd\\(", colnames(standard))))
  expect_false(any(grepl("var_common", colnames(standard), fixed = TRUE)))
  expect_true("(mu) allocation: var_mult(x)" %in% colnames(full))
  expect_true("(mu) allocation: sd_mult(x)" %in% colnames(full))
  expect_true("(mu) allocation: var_common" %in% colnames(full))
  expect_true("(mu) sd(x)" %in% colnames(full))
})

test_that("known group covariance retains its fitted kernel scale", {

  skip_if_not_installed("runjags")

  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")))
  kernel <- matrix(
    c(1, .25, .25, 1),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("s1", "s2"))
  )
  random_formula <- random_effects_formula(
    ~ 1 | study,
    group_covariance = random_group_covariance(kernel, scale = "none")
  )
  formula_result <- JAGS_formula(
    formula = random_formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  samples <- matrix(
    c(1, 2),
    ncol = 1L,
    dimnames = list(NULL, random_term$sd_parameter_names)
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

  standard <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    return_samples = TRUE
  )
  full <- JAGS_estimates_table(
    fit,
    random_effects_summary = "full",
    return_samples = TRUE
  )

  expect_identical(colnames(standard), "(mu) sd(intercept)")
  expect_true("(mu) sd(intercept)" %in% colnames(full))
  expect_true("(mu) var(intercept)" %in% colnames(full))

  simplified <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    simplify_names = TRUE,
    return_samples = TRUE
  )
  expect_identical(colnames(simplified), "(mu) sd")
  expect_equal(simplified[, 1L], samples[, 1L])
})

test_that("random-effect summary posterior extracts total-variance proportions", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_total_variance_allocation_fit()
  proportions <- random_effects_summary_posterior(fit, summary = "var_prop")
  proportion_name <- "(mu) allocation: var_prop(drug)"

  expect_true(proportion_name %in% names(proportions))
  expect_equal(
    unname(as.numeric(proportions[[proportion_name]])),
    c(0.75, 0.25),
    tolerance = 1e-12
  )

  prior_density <- .bt_meta_get(proportions[[proportion_name]], "prior_density")
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(posterior_metadata(proportions[[proportion_name]], "support")$bounds, c(0, 1))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 0.5),
    stats::dbeta(0.5, 3, 2),
    tolerance = 1e-7
  )

  expect_error(
    random_effects_summary_posterior(fit, summary = "var_mult"),
    "Variance-multiplier summaries are created only",
    fixed = TRUE
  )
})

test_that("gated total-variance summaries use realized totals and proportions", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_gated_total_variance_allocation_fit()
  catalog <- parameter_catalog(fit)
  total <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, "(mu) allocation: sd_total")
  )
  expect_equal(
    unname(as.numeric(total[[1L]][, 1L])),
    c(0, 1, sqrt(3), 2),
    tolerance = 1e-12
  )

  study <- parameter_draws(
    fit,
    parameter_catalog_resolve(
      catalog,
      "(mu) allocation: var_prop(study)"
    )
  )
  drug <- parameter_draws(
    fit,
    parameter_catalog_resolve(
      catalog,
      "(mu) allocation: var_prop(drug)"
    )
  )
  expect_equal(unname(as.numeric(study[[1L]][, 1L])), c(NA, 1, 0, .25))
  expect_equal(unname(as.numeric(drug[[1L]][, 1L])), c(NA, 0, 1, .75))

  proportions <- random_effects_summary_posterior(
    fit,
    summary = "var_prop"
  )
  expect_equal(
    unname(as.numeric(proportions[[
      "(mu) allocation: var_prop(drug)"
    ]])),
    c(0, 1, .75)
  )
  # the gated proportion's prior is the canonical mixed measure: with both
  # gates at 0.5 and weights Dirichlet(2, 3), drug is 0 (study only active),
  # 1 (drug only) or Beta(3, 2) (both), each with probability 1/3 given an
  # active component
  gated_prior <- .bt_meta_get(proportions[["(mu) allocation: var_prop(drug)"]], "prior_density")
  expect_equal(gated_prior$points$x, c(0, 1))
  expect_equal(gated_prior$points$p, c(1, 1) / 3, tolerance = 1e-12)
  expect_equal(
    prior_density_ordinate(gated_prior, .5)$log_density,
    log(stats::dbeta(.5, 3, 2) / 3),
    tolerance = 1e-12
  )

  estimates <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE
  )
  expect_equal(
    estimates["(mu) allocation: sd_total", "Mean"],
    mean(c(0, 1, sqrt(3), 2)),
    tolerance = 1e-12
  )
  expect_equal(
    estimates["(mu) allocation: var_prop(drug)", "Mean"],
    mean(c(0, 1, .75)),
    tolerance = 1e-12
  )
  component_names <- c("(mu) study: sd(intercept)", "(mu) drug: sd(intercept)")
  expect_false(any(component_names %in% rownames(estimates)))
  full <- JAGS_estimates_table(
    fit,
    random_effects_summary = "full",
    remove_diagnostics = TRUE
  )
  expect_equal(full[component_names, "Mean"], c(.5, sqrt(3) / 2),
               tolerance = 1e-12)
  conditional <- JAGS_estimates_table(
    fit,
    conditional = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE
  )
  expect_false(any(component_names %in% rownames(conditional)))
})


test_that("random-effect summary posteriors declare atoms from structure", {

  skip_if_not_installed("runjags")

  # Summaries whose canonical prior density has no point mass declare no atoms.
  multipliers <- random_effects_summary_posterior(
    .random_effects_mean_variance_allocation_fit(),
    summary = "var_mult"
  )
  for(multiplier in multipliers){
    atoms <- BayesTools:::.posterior_atoms_get(multiplier)
    expect_s3_class(atoms, "BayesTools_posterior_atoms")
    expect_identical(nrow(atoms$locations), 0L)
  }

  # Gate atoms take their masses from the inclusion-indicator draws: the
  # fixture's (study, drug) gates are (0, 0), (1, 0), (0, 1), (1, 1).
  fit <- .random_effects_gated_total_variance_allocation_fit()
  proportions <- random_effects_summary_posterior(fit, summary = "var_prop")
  # draw 1 (no active component) is undefined; drug is 0 in draw 2 (only study
  # active) and 1 in draw 3 (only drug active)
  drug_atoms <- BayesTools:::.posterior_atoms_get(
    proportions[["(mu) allocation: var_prop(drug)"]]
  )
  expect_equal(as.numeric(drug_atoms$locations[, 1L]), c(0, 1))
  expect_equal(drug_atoms$mass, c(1, 1) / 3)
  total <- random_effects_summary_posterior(fit, summary = "sd_total")
  total_atoms <- BayesTools:::.posterior_atoms_get(
    total[["(mu) allocation: sd_total"]]
  )
  expect_equal(as.numeric(total_atoms$locations[, 1L]), 0)
  expect_equal(total_atoms$mass, 1 / 4)
  expect_identical(
    unname(as.numeric(total[["(mu) allocation: sd_total"]]))[1L],
    0
  )

  # A scale prior with its own point mass is not a gate atom: the atom status
  # stays undeclared and plots stop instead of inferring point masses.
  spike_fit <- .random_effects_mean_variance_allocation_fit(
    sd = prior_mixture(
      list(prior("spike", list(0)), prior("gamma", list(2, 2))),
      is_null = c(TRUE, FALSE)
    )
  )
  common <- random_effects_summary_posterior(spike_fit, summary = "sd_common")
  expect_null(BayesTools:::.posterior_atoms_get(common[[1L]]))
  expect_error(
    plot_posterior(common, names(common)[1L], plot_type = "ggplot"),
    "Posterior atom status is unknown",
    fixed = TRUE
  )
})

test_that("parameter_mixed_posterior declares gate atoms and conditions on inclusion", {

  skip_if_not_installed("runjags")

  # (study, drug) gates are (0, 0), (1, 0), (0, 1), (1, 1) with the total SD 2
  # and weights (1/4, 3/4); both gates have prior probability 1/2
  fit <- .random_effects_gated_total_variance_allocation_fit()
  catalog <- parameter_catalog(fit)
  mixed <- function(name, conditional = FALSE){
    parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, name),
                              conditional = conditional)
  }
  atoms_of <- function(x){
    atoms <- posterior_metadata(x, "atoms")
    list(x = as.numeric(atoms$locations[, 1L]), mass = atoms$mass)
  }
  study_gate <- "mu__xRE_ALLOCx_allocation__include_study_indicator"
  drug_gate  <- "mu__xRE_ALLOCx_allocation__include_drug_indicator"

  # the component SD of study is 2 * gate * sqrt(1/4): its atom at 0 has the
  # posterior mass of the off gate (draws 1 and 3)
  study <- mixed("(mu) study: sd(intercept)")
  expect_s3_class(study, c("mixed_posteriors", "mixed_posteriors.simple",
                           "marginal_posterior.simple", "marginal_posterior"),
                  exact = TRUE)
  expect_identical(attr(study, "parameter"), "(mu) study: sd(intercept)")
  expect_equal(as.numeric(study), c(0, 1, 0, 1))
  expect_identical(atoms_of(study), list(x = 0, mass = .5))
  expect_equal(posterior_metadata(study, "prior_density")$points$p, .5)
  expect_identical(posterior_metadata(study, "support"), catalog$quantities$support[[
    match("(mu) study: sd(intercept)", catalog$quantities$canonical_name)
  ]])
  expect_true(posterior_metadata(study, "condition")$averaged)
  expect_false(posterior_atoms_free(study))
  study_included <- mixed("(mu) study: sd(intercept)", conditional = TRUE)
  expect_equal(as.numeric(study_included), c(1, 1))
  expect_true(posterior_atoms_free(study_included))
  expect_equal(nrow(posterior_metadata(study_included, "prior_density")$points), 0L)
  expect_equal(posterior_metadata(study_included, "prior_density")$density$mass, 1)
  expect_identical(
    posterior_metadata(study_included, "condition")[c("conditional", "conditional_rule", "averaged")],
    list(conditional = study_gate, conditional_rule = "AND", averaged = FALSE)
  )
  expect_identical(atoms_of(mixed("(mu) study: var(intercept)")), list(x = 0, mass = .5))

  # the total SD is 0 without an active component (draw 1); conditioning on
  # any active component keeps draws 2-4 and removes the atom
  total <- mixed("(mu) allocation: sd_total")
  expect_equal(as.numeric(total), c(0, 1, sqrt(3), 2))
  expect_identical(atoms_of(total), list(x = 0, mass = .25))
  total_included <- mixed("(mu) allocation: sd_total", conditional = TRUE)
  expect_equal(as.numeric(total_included), c(1, sqrt(3), 2))
  expect_true(posterior_atoms_free(total_included))
  expect_equal(posterior_metadata(total_included, "prior_density")$density$mass, 1)
  expect_identical(
    posterior_metadata(total_included, "condition")[c("conditional", "conditional_rule")],
    list(conditional = c(study_gate, drug_gate), conditional_rule = "OR")
  )

  # the proportion of drug is undefined in draw 1, 0 in draw 2 and 1 in draw
  # 3: atoms 1/3 each over the defined draws; conditional on its own gate
  # (draws 3 and 4) the atom at 1 keeps its renormalized mass, prior 1/2 at 1
  # and 1/2 Beta(3, 2)
  proportion <- mixed("(mu) allocation: var_prop(drug)")
  expect_equal(as.numeric(proportion), c(0, 1, .75))
  expect_identical(atoms_of(proportion), list(x = c(0, 1), mass = c(1, 1) / 3))
  expect_identical(
    posterior_metadata(proportion, "undefined_draws"),
    c("(mu) allocation: var_prop(drug)" = "allocation_active")
  )
  proportion_included <- mixed("(mu) allocation: var_prop(drug)", conditional = TRUE)
  expect_equal(as.numeric(proportion_included), c(1, .75))
  expect_identical(atoms_of(proportion_included), list(x = 1, mass = .5))
  prior_included <- posterior_metadata(proportion_included, "prior_density")
  expect_equal(prior_included$points$x, 1)
  expect_equal(prior_included$points$p, .5, tolerance = 1e-12)
  expect_equal(
    prior_density_ordinate(prior_included, .5)$log_density,
    log(.5 * stats::dbeta(.5, 3, 2)),
    tolerance = 1e-12
  )

  # quantities without an inclusion gate cannot be conditioned
  multiplier_fit <- .random_effects_mean_variance_allocation_fit()
  multiplier <- parameter_catalog_resolve(
    parameter_catalog(multiplier_fit), "(mu) allocation: var_mult(x)", namespace = "mu"
  )
  expect_true(posterior_atoms_free(parameter_mixed_posterior(multiplier_fit, multiplier)))
  expect_error(
    parameter_mixed_posterior(multiplier_fit, multiplier, conditional = TRUE),
    paste0(
      "The inclusion event of '(mu) allocation: var_mult(x)' is unavailable: ",
      "the quantity has no inclusion gate. Use 'conditional = FALSE'."
    ),
    fixed = TRUE
  )
  expect_error(parameter_mixed_posterior(multiplier_fit, multiplier, conditional = NA),
               "conditional", fixed = TRUE)
  expect_error(parameter_mixed_posterior(list(), multiplier),
               "'fit' must be a 'BayesTools_fit' object.", fixed = TRUE)
})

test_that("parameter_mixed_posterior reads source point masses from the mixture indicator", {

  skip_if_not_installed("runjags")

  # a mean-variance SD-component allocation whose scale prior is a spike at
  # 0 or gamma(2, 2) with prior probability 1/2 each; the monitored mixture
  # indicator selects the spike in the first draw
  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")), x = c(-1, 0, 1, 2))
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "mean_variance",
        sd = prior_mixture(
          list(prior("spike", list(0)), prior("gamma", list(2, 2))),
          is_null = c(TRUE, FALSE)
        ),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  allocation <- formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1]]
  samples <- matrix(
    c(0, 1, 0.25, 0.75,
      2, 2, 0.75, 0.25),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      allocation$source_node,
      paste0(allocation$source_node, "_indicator"),
      paste0(allocation$weight_name, "[1]"),
      paste0(allocation$weight_name, "[2]")
    ))
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)
  catalog <- parameter_catalog(fit)

  for(name in c("(mu) allocation: sd_common", "(mu) sd(x)", "(mu) var(x)")){
    x <- parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, name, namespace = "mu"))
    atoms <- posterior_metadata(x, "atoms")
    expect_equal(as.numeric(atoms$locations[, 1L]), 0, info = name)
    expect_equal(atoms$mass, .5, info = name)
    expect_equal(posterior_metadata(x, "prior_density")$points$p, .5, info = name)
  }
  common <- random_effects_summary_posterior(fit, summary = "sd_common")
  expect_equal(posterior_metadata(common[[1L]], "atoms")$mass, .5)
})

test_that("catalog quantities declare their exact support and definedness", {

  skip_if_not_installed("runjags")

  support_of <- function(fit, name){
    quantities <- parameter_catalog(fit)$quantities
    quantities$support[[match(name, quantities$canonical_name)]]
  }
  definedness_of <- function(fit, name){
    quantities <- parameter_catalog(fit)$quantities
    quantities$definedness[[match(name, quantities$canonical_name)]]
  }

  # mean-variance multipliers of K = 2 components: var_mult = K w in [0, K]
  # and sd_mult = sqrt(K w) in [0, sqrt(K)]; the gamma scale prior is
  # unbounded above, so composite SDs and variances span [0, Inf)
  fit <- .random_effects_mean_variance_allocation_fit()
  expect_equal(support_of(fit, "(mu) allocation: var_mult(x)")$bounds, c(0, 2))
  expect_equal(support_of(fit, "(mu) allocation: sd_mult(x)")$bounds, c(0, sqrt(2)))
  expect_equal(support_of(fit, "(mu) sd(x)")$bounds, c(0, Inf))
  expect_true(support_of(fit, "(mu) sd(x)")$exact)
  expect_identical(definedness_of(fit, "(mu) sd(x)"), "always")

  # a truncated scale prior bounds the one-to-one quantities exactly and
  # leaves the composite hull inexact (shares reach zero)
  truncated <- .random_effects_mean_variance_allocation_fit(
    sd = prior("normal", list(0, 1), list(.5, 3))
  )
  expect_equal(support_of(truncated, "(mu) allocation: sd_common")$bounds, c(.5, 3))
  expect_equal(support_of(truncated, "(mu) allocation: var_common")$bounds, c(.25, 9))
  expect_true(support_of(truncated, "(mu) allocation: var_common")$exact)
  expect_false(support_of(truncated, "(mu) sd(x)")$exact)

  # original-scale SDs of a scaled random slope: with SD priors on [0.5, Inf)
  # the intercept SD sqrt(sd_0^2 + (m / s)^2 sd_1^2) is at least
  # 0.5 sqrt(1 + (m / s)^2) and the slope SD sd_1 / s at least 0.5 / s, so the
  # hull [0, Inf) is not their exact support; half-normal SDs reach every
  # value in (0, Inf)
  scaled_slope_fit <- function(sd){
    data <- data.frame(x = c(2, 4, 6, 8, 3, 7), g = factor(rep(c("a", "b", "c"), 2)))
    formula_result <- JAGS_formula(
      ~ x + random(1 + x | g, covariance = "diag"), "mu", data,
      prior_list   = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = sd),
      formula_scale = list(x = TRUE)
    )
    prior_list <- formula_result$prior_list
    columns <- unlist(lapply(names(prior_list), function(parameter){
      BayesTools:::.prior_linear_prior_columns(parameter, prior_list[[parameter]])
    }))
    .parameter_catalog_test_fit(
      coda::mcmc.list(coda::mcmc(matrix(1, 4, length(columns), dimnames = list(NULL, columns)))),
      prior_list,
      formula_design = list(mu = formula_result$formula_design),
      formula_scale  = list(mu = formula_result$formula_scale)
    )
  }
  bounded_below <- scaled_slope_fit(prior("normal", list(0, 1), list(.5, Inf)))
  half_normal   <- scaled_slope_fit(prior("normal", list(0, 1), list(0, Inf)))
  for(name in c("(mu) sd(intercept)", "(mu) sd(x)", "(mu) var(intercept)")){
    expect_equal(support_of(bounded_below, name)$bounds, c(0, Inf))
    expect_false(support_of(bounded_below, name)$exact)
    expect_true(support_of(half_normal, name)$exact)
  }

  # gated total-variance proportions lie in [0, 1] and are undefined where
  # no component is active; the inclusion indicators take the values 0 and 1
  gated <- .random_effects_gated_total_variance_allocation_fit()
  proportion <- "(mu) allocation: var_prop(drug)"
  expect_equal(support_of(gated, proportion)$bounds, c(0, 1))
  expect_identical(definedness_of(gated, proportion), "allocation_active")
  expect_identical(support_of(gated, "(mu) allocation: inclusion(drug)")$points, c(0, 1))

  # parameter_draws() carries the catalog metadata of the selected quantity
  draws <- parameter_draws(gated, parameter_catalog_resolve(parameter_catalog(gated), proportion))
  expect_identical(
    posterior_metadata(draws, "undefined_draws"),
    stats::setNames("allocation_active", proportion)
  )
  expect_equal(posterior_metadata(draws, "support")[[proportion]]$bounds, c(0, 1))
  values <- as.numeric(as.matrix(draws))
  table <- ensemble_estimates_table(
    list(p = BayesTools:::.bt_meta_set(values, "undefined_draws", "allocation_active")),
    parameters = "p"
  )
  expect_equal(table["p", "Mean"], mean(values, na.rm = TRUE))

  # random-effect summary posteriors use the catalog support
  multipliers <- random_effects_summary_posterior(fit, summary = "sd_mult")
  expect_identical(
    posterior_metadata(multipliers[["(mu) allocation: sd_mult(x)"]], "support"),
    support_of(fit, "(mu) allocation: sd_mult(x)")
  )
})

test_that("inherited gates zero descendant realized allocations", {

  allocation <- list(
    scale          = "total_variance",
    n_targets      = 2L,
    inclusion      = list(),
    parent_factors = list(list(inclusion_name = "parent_gate"))
  )
  weights <- matrix(
    c(0.25, 0.75, 0.4, 0.6),
    nrow = 2L,
    byrow = TRUE
  )
  samples <- cbind(parent_gate = c(0, 1))
  realized <- .bt_random_effect_summary_realized_allocation(
    allocation    = allocation,
    weights       = weights,
    model_samples = samples
  )

  expect_equal(realized[["total_fraction"]], c(0, 1))
  expect_true(all(is.na(realized[["proportions"]][1L, ])))
  expect_equal(realized[["proportions"]][2L, ], weights[2L, ])
})

test_that("random-effect summary posterior handles singular Dirichlet boundaries", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit(alpha = c(0.5, 2))
  multipliers <- random_effects_summary_posterior(
    fit,
    summary = "var_mult",
    component = "intercept",
    n_prior_points = 64
  )
  prior_density <- .bt_meta_get(multipliers[[1]], "prior_density")

  expect_equal(posterior_metadata(multipliers[[1]], "support")$bounds, c(0, 2))
  expect_true(all(is.finite(prior_density$density$x)))
  expect_true(all(is.finite(prior_density$density$y)))
  expect_identical(
    BayesTools:::.prior_linear_density_height(prior_density, 0),
    Inf
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(.5, .5, 2) / 2,
    tolerance = 1e-12
  )
  # the lower bound is a singular boundary (w^(-1/2) at 0)
  expect_identical(prior_density_ordinate(prior_density, 0)$behavior, "infinite")

  prior_plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density,
    n_points = 32
  )
  expect_true(all(is.finite(prior_plot_data$density$y)))
})

test_that("unit-scale Dirichlet summaries keep singular upper bounds off the grid", {

  skip_if_not_installed("runjags")

  # With scale 1, the unit reference value is the upper support bound. It is
  # singular when the complementary Dirichlet mass sum(alpha) - alpha_i < 1.
  for(alpha in list(c(0.5, 0.5), c(1, 0.6))){
    fit <- .random_effects_total_variance_allocation_fit(alpha = alpha)
    proportions <- random_effects_summary_posterior(
      fit,
      summary = "var_prop",
      n_prior_points = 64
    )
    for(index in 1:2){
      component <- c("study", "drug")[index]
      prior_density <- .bt_meta_get(proportions[[paste0("(mu) allocation: var_prop(", component, ")")]], "prior_density")
      alpha_i <- alpha[index]
      beta_i  <- sum(alpha) - alpha_i
      info <- paste0("alpha = c(", toString(alpha), "), ", component)

      expect_s3_class(prior_density, "prior_linear_density")
      expect_equal(
        posterior_metadata(proportions[[paste0("(mu) allocation: var_prop(", component, ")")]], "support")$bounds,
        c(0, 1),
        info = info
      )
      expect_true(all(is.finite(prior_density$density$y)), info = info)
      # a singular upper bound stays off the plotted values
      plot_data <- BayesTools:::.prior_linear_density_to_plot_data(prior_density, n_points = 32)
      expect_true(all(is.finite(plot_data$density$y)), info = info)
      expect_identical(
        prior_density_ordinate(prior_density, 1)$behavior,
        if(beta_i < 1) "infinite" else "regular",
        info = info
      )
      # Reference: the analytic Beta(alpha_i, sum(alpha) - alpha_i) margin.
      expect_equal(
        BayesTools:::.prior_linear_density_height(prior_density, 0.5),
        stats::dbeta(0.5, alpha_i, beta_i),
        tolerance = 1e-12,
        info = info
      )
      expect_identical(
        BayesTools:::.prior_linear_density_height(prior_density, 1),
        if(beta_i < 1) Inf else alpha_i
      )
    }
  }

  # A total-variance SD-component allocation exposes sqrt(w) with scale 1.
  fit <- .random_effects_mean_variance_allocation_fit(
    alpha = c(0.5, 0.5),
    scale = "total_variance"
  )
  multipliers <- random_effects_summary_posterior(
    fit,
    summary = "sd_mult",
    component = "x",
    n_prior_points = 64
  )
  prior_density <- .bt_meta_get(multipliers[[1L]], "prior_density")
  expect_equal(posterior_metadata(multipliers[[1L]], "support")$bounds, c(0, 1))
  expect_true(all(is.finite(prior_density$density$y)))
  expect_true(all(is.finite(
    BayesTools:::.prior_linear_density_to_plot_data(prior_density, n_points = 32)$density$y
  )))
  expect_identical(prior_density_ordinate(prior_density, 1)$behavior, "infinite")
  # Reference: density of sqrt(w), w ~ Beta(0.5, 0.5), is dbeta(x^2) * 2x.
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 0.5),
    stats::dbeta(0.25, 0.5, 0.5) * 2 * 0.5,
    tolerance = 1e-12
  )
  proportions <- random_effects_summary_posterior(
    fit,
    summary = "var_prop",
    n_prior_points = 64
  )
  expect_length(proportions, 2L)
  expect_true(all(vapply(proportions, function(x){
    all(is.finite(.bt_meta_get(x, "prior_density")$density$y))
  }, logical(1))))
})

test_that("scaled-Beta prior densities preserve square-root endpoint limits", {

  # sqrt(4 w), w ~ Beta(0.5, 1), the density route of SD multipliers
  prior_density <- BayesTools:::.bt_parameter_prior_density_transformed(
    prior("beta", list(alpha = .5, beta = 1)),
    list(type = "sqrt_scale", scale = 4),
    n_grid = 64,
    tail_prob = BayesTools:::.prior_linear_density_tail_prob()
  )

  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 0),
    1 / sqrt(4),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 2),
    2 * .5 / 2,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, -1),
    0
  )
})

test_that("raw random-slope SD rows keep their own level labels", {

  # Index-like levels: the level-renamed SD column `f[2]` (level 2) has the
  # same name as the backend coordinate `f[2]` (level 3).
  data <- data.frame(
    f = factor(rep(1:3, 4)),
    g = factor(rep(c("A", "B", "C", "D"), each = 3))
  )
  formula_result <- JAGS_formula(
    ~ 1 + diag(0 + f | g), "mu", data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  expect_identical(
    unname(random_term$sd_leaves$leaf_terms),
    c("f[2]", "f[3]")
  )
  # Distinct constant draws identify each SD: level 2 = 0.5, level 3 = 1.5.
  posterior <- matrix(
    rep(c(0, 0.5, 1.5), 4),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(posterior)),
      sample = nrow(posterior),
      summary.pars = list(mutate = NULL),
      monitor = colnames(posterior)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

  standard <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE
  )
  raw <- JAGS_estimates_table(
    fit,
    random_effects_summary = "raw",
    remove_diagnostics = TRUE
  )
  # Reference: the standard table of the same fit.
  expect_equal(standard["(mu) sd(f[2])", "Mean"], 0.5)
  expect_equal(standard["(mu) sd(f[3])", "Mean"], 1.5)
  expect_equal(raw["(mu) g: sd(f[2])", "Mean"], standard["(mu) sd(f[2])", "Mean"])
  expect_equal(raw["(mu) g: sd(f[3])", "Mean"], standard["(mu) sd(f[3])", "Mean"])
  expect_identical(sum(grepl("sd(", rownames(raw), fixed = TRUE)), 2L)
})

test_that("raw group rows require their group level labels", {

  raw_names <- c("mu__xREx__g_xRE_Zx[2,1]", "mu__xREx__g_xRE_COEFx[1,1]")
  label <- function(group_levels){
    BayesTools:::.bt_random_effect_summary_raw_effect_matrix_names(
      names = raw_names,
      raw_names = raw_names,
      stem = "mu__xREx__g",
      matrix = "_xRE_Zx",
      label = "z",
      components = "intercept",
      group = "g",
      group_levels = group_levels,
      prefix = "(mu) "
    )
  }

  # Numeric group labels: row 2 is the group labelled "20", never "2".
  expect_identical(
    label(c("10", "20")),
    c("(mu) z(g[20], intercept)", "mu__xREx__g_xRE_COEFx[1,1]")
  )
  expect_error(
    label(NULL),
    paste0(
      "Random-effect summary metadata do not identify group 2 of grouping ",
      "factor 'g'. Refit the model with the current BayesTools version."
    ),
    fixed = TRUE
  )
  expect_error(
    label("10"),
    "do not identify group 2 of grouping factor 'g'",
    fixed = TRUE
  )
})

# A fit of an allocation over random intercepts of 'study', 'drug' (and
# 'site') with scale prior 'sd', Dirichlet 'alpha' and optional inclusion
# gates; the draws only name the monitored coordinates.
.allocation_prior_fit <- function(sd, alpha, inclusion = NULL){

  terms <- c("study", "drug", "site")[seq_along(alpha)]
  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug  = factor(c("a", "b", "a", "b")),
    site  = factor(c("x", "y", "y", "x"))
  )
  formula <- stats::as.formula(paste(
    "~ 1 +",
    paste0("random(1 | ", terms, ", name = '", terms, "', covariance = 'diag')",
           collapse = " + ")
  ))
  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "allocation",
        terms = stats::setNames(terms, terms),
        sd = sd,
        weights = prior("dirichlet", list(alpha = alpha)),
        inclusion = inclusion
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1L]]
  gates <- vapply(allocation$inclusion, `[[`, character(1), "indicator_name")
  columns <- c(
    allocation$source_node,
    if(is.prior.mixture(sd)) paste0(allocation$source_node, "_indicator"),
    paste0(allocation$weight_name, "[", seq_along(alpha), "]"),
    gates
  )
  samples <- matrix(1 / length(alpha), nrow = 4L, ncol = length(columns),
                    dimnames = list(NULL, columns))
  samples[, allocation$source_node] <- 1
  samples[, gates] <- 1
  if(is.prior.mixture(sd)){
    samples[, paste0(allocation$source_node, "_indicator")] <- 1
  }
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attach_test_parameter_map(fit)
}

# The fitted model's generative process of the allocated SDs (component i:
# T g_i sqrt(w_i); total: T sqrt(sum_j g_j w_j)) with w ~ Dirichlet(alpha)
# over all components.
.allocation_prior_draws <- function(n, draw_sd, alpha, gate_probability){

  eta <- matrix(stats::rgamma(n * length(alpha), shape = rep(alpha, each = n)),
                nrow = n)
  w <- eta / rowSums(eta)
  g <- matrix(stats::rbinom(n * length(alpha), 1, rep(gate_probability, each = n)),
              nrow = n)
  scale <- draw_sd(n)
  list(
    components = scale * g * sqrt(w),
    total      = scale * sqrt(rowSums(g * w)),
    gates      = g
  )
}

# Bin probabilities of a prior density (its region route) against draws of the
# generative process: each within 4 binomial standard errors; the atom at 0
# as well.
.expect_allocation_bins <- function(density, draws, edges, info){

  n <- length(draws)
  model_zero <- sum(density$points$p[density$points$x == 0])
  zero_se <- sqrt(max(model_zero * (1 - model_zero), 1 / n) / n)
  expect_lte(abs(mean(draws == 0) - model_zero), 4 * zero_se)
  for(k in seq_len(length(edges) - 1L)){
    lower <- edges[k]
    upper <- edges[k + 1L]
    region <- list(
      intervals = matrix(c(lower, upper), 1L),
      indicator = function(x) x > lower & x <= upper
    )
    probability <- as.numeric(BayesTools:::.prior_linear_density_region_probability(density, region))
    se <- sqrt(probability * (1 - probability) / n)
    expect_lte(abs(mean(draws > lower & draws <= upper) - probability), 4 * se)
  }
}

test_that("allocation-derived SDs, variances and totals have exact prior densities", {

  skip_if_not_installed("runjags")

  # the structure of RoBMA's BMA.mv_random_components: a two-component
  # allocation with Dirichlet(1, 1) weights, both components gated with
  # probability 1/2, and the scale prior 0.6 HN(0.25) + 0.4 HN(0.75)
  scale_prior <- prior_mixture(list(
    prior("normal", list(0, .25), list(0, Inf), prior_weights = 3),
    prior("normal", list(0, .75), list(0, Inf), prior_weights = 2)
  ))
  fit <- .allocation_prior_fit(
    sd = scale_prior,
    alpha = c(1, 1),
    inclusion = list(study = prior("spike", list(.5)), drug = prior("spike", list(.5)))
  )
  catalog <- parameter_catalog(fit)
  density_of <- function(name, conditional = FALSE){
    BayesTools:::.bt_parameter_prior_density_quantity(
      fit, parameter_catalog_resolve(catalog, name), n_grid = 1024L,
      tail_prob = BayesTools:::.prior_linear_density_tail_prob(),
      conditional = conditional
    )
  }
  # model formulas with an independent quadrature: the component SD is
  # (1/2) h(y) and the total (1/4) f_T(y) + (1/2) h(y), where
  # h(y) = int_0^1 f_T(y / sqrt(s)) / sqrt(s) ds (the Beta(1, 1) share)
  f_T <- function(y) .6 * 2 * stats::dnorm(y, 0, .25) + .4 * 2 * stats::dnorm(y, 0, .75)
  h <- function(y){
    stats::integrate(function(s) f_T(y / sqrt(s)) / sqrt(s), 0, 1, rel.tol = 1e-13)$value
  }
  for(name in c("(mu) study: sd(intercept)", "(mu) drug: sd(intercept)")){
    density <- density_of(name)
    expect_identical(attr(density, "adaptive_evaluation")$kind, "allocation_product")
    expect_equal(density$points$p[density$points$x == 0], .5, tolerance = 1e-12)
    for(y in c(1e-3, .06, .3, 1.5)){
      ordinate <- prior_density_ordinate(density, y)
      expect_true(ordinate$exact)
      expect_equal(exp(ordinate$log_density), h(y) / 2, tolerance = 1e-8)
    }
  }
  # RoBMA's fit gives 1.7780 and 2.3492 at 0.06 (MC 1.7746, se 0.0105)
  expect_equal(exp(prior_density_ordinate(density_of("(mu) study: sd(intercept)"), .06)$log_density),
               1.7780422878, tolerance = 1e-9)
  total <- density_of("(mu) allocation: sd_total")
  expect_equal(total$points$p[total$points$x == 0], .25, tolerance = 1e-12)
  expect_equal(exp(prior_density_ordinate(total, .06)$log_density),
               f_T(.06) / 4 + h(.06) / 2, tolerance = 1e-8)
  expect_equal(exp(prior_density_ordinate(total, .06)$log_density), 2.3492269474, tolerance = 1e-9)
  # variances by the Jacobian: f_V(v) = f(sqrt(v)) / (2 sqrt(v))
  variance <- density_of("(mu) study: var(intercept)")
  expect_equal(exp(prior_density_ordinate(variance, .0036)$log_density),
               h(.06) / 2 / (2 * .06), tolerance = 1e-8)
  expect_true(prior_density_ordinate(variance, .0036)$exact)
  # its plotted curve starts at the infinite density next to zero
  plotted <- BayesTools:::.prior_linear_density_to_plot_data(variance)$density
  expect_lt(min(plotted$x), 1e-3)
  expect_true(all(is.finite(plotted$y)))
  var_total <- density_of("(mu) allocation: var_total")
  expect_equal(exp(prior_density_ordinate(var_total, .0036)$log_density),
               (f_T(.06) / 4 + h(.06) / 2) / (2 * .06), tolerance = 1e-8)
  # the atoms keep their class at 0 with the continuous part's behavior
  at_zero <- prior_density_ordinate(variance, 0)
  expect_identical(at_zero$behavior, "point_mass")
  expect_identical(at_zero$provenance$continuous_behavior, "infinite")
  expect_identical(prior_density_ordinate(total, 0)$provenance$continuous_behavior, "regular")
  # conditional on the inclusion event the gate atoms drop out
  included <- density_of("(mu) study: sd(intercept)", conditional = TRUE)
  expect_identical(nrow(included$points), 0L)
  expect_equal(exp(prior_density_ordinate(included, .06)$log_density), h(.06), tolerance = 1e-8)
  total_included <- density_of("(mu) allocation: sd_total", conditional = TRUE)
  expect_equal(exp(prior_density_ordinate(total_included, .06)$log_density),
               (f_T(.06) / 4 + h(.06) / 2) / .75, tolerance = 1e-8)

  # bin probabilities against 1e6 draws of the generative process
  set.seed(20260925)
  draw_T <- function(n){
    abs(stats::rnorm(n, 0, ifelse(stats::runif(n) < .6, .25, .75)))
  }
  draws <- .allocation_prior_draws(1e6, draw_T, c(1, 1), c(.5, .5))
  edges <- c(0, .01, .03, .06, .12, .25, .5, 1, Inf)
  .expect_allocation_bins(density_of("(mu) study: sd(intercept)"), draws$components[, 1L], edges)
  .expect_allocation_bins(density_of("(mu) drug: sd(intercept)"), draws$components[, 2L], edges)
  .expect_allocation_bins(total, draws$total, edges)
  .expect_allocation_bins(variance, draws$components[, 1L]^2, edges^2)
  .expect_allocation_bins(var_total, draws$total^2, edges^2)
  .expect_allocation_bins(included, draws$components[draws$gates[, 1L] == 1, 1L], edges)
  .expect_allocation_bins(total_included, draws$total[rowSums(draws$gates) > 0], edges)
})

test_that("allocation prior densities cover partial sets, ungated components and mean-variance shares", {

  skip_if_not_installed("runjags")

  # three components with Dirichlet(1, 2, 3); study and drug gated (0.3,
  # 0.6), site always active: the total mixes the sets {site}, {study, site},
  # {drug, site} and the full set; T ~ gamma(2, 2)
  fit <- .allocation_prior_fit(
    sd = prior("gamma", list(2, 2)),
    alpha = c(1, 2, 3),
    inclusion = list(study = prior("spike", list(.3)), drug = prior("spike", list(.6)))
  )
  catalog <- parameter_catalog(fit)
  density_of <- function(name, conditional = FALSE){
    BayesTools:::.bt_parameter_prior_density_quantity(
      fit, parameter_catalog_resolve(catalog, name), n_grid = 1024L,
      tail_prob = BayesTools:::.prior_linear_density_tail_prob(),
      conditional = conditional
    )
  }
  names <- c("(mu) study: sd(intercept)", "(mu) drug: sd(intercept)",
             "(mu) site: sd(intercept)", "(mu) allocation: sd_total")
  densities <- lapply(names, density_of)
  expect_equal(densities[[1L]]$points$p, .7, tolerance = 1e-12)
  expect_equal(densities[[2L]]$points$p, .4, tolerance = 1e-12)
  expect_identical(nrow(densities[[3L]]$points), 0L)
  expect_identical(nrow(densities[[4L]]$points), 0L)
  for(density in densities){
    expect_true(prior_density_ordinate(density, .4)$exact)
  }
  # the total of the site-only set is T sqrt(W), W ~ Beta(3, 3) (probability
  # 0.7 * 0.4); with the full set's T and the two other partial sets
  f_T <- function(y) stats::dgamma(y, 2, 2)
  h <- function(y, a, b){
    stats::integrate(function(s) f_T(y / sqrt(s)) * stats::dbeta(s, a, b) / sqrt(s),
                     0, 1, rel.tol = 1e-13)$value
  }
  expect_equal(
    exp(prior_density_ordinate(densities[[4L]], .4)$log_density),
    .3 * .6 * f_T(.4) + .7 * .4 * h(.4, 3, 3) + .3 * .4 * h(.4, 4, 2) +
      .7 * .6 * h(.4, 5, 1),
    tolerance = 1e-8
  )

  set.seed(20260926)
  draws <- .allocation_prior_draws(1e6, function(n) stats::rgamma(n, 2, 2),
                                   c(1, 2, 3), c(.3, .6, 1))
  edges <- c(0, .05, .15, .3, .5, .8, 1.2, 2, Inf)
  for(i in 1:3){
    .expect_allocation_bins(densities[[i]], draws$components[, i], edges)
  }
  .expect_allocation_bins(densities[[4L]], draws$total, edges)
  .expect_allocation_bins(density_of("(mu) allocation: var_total"), draws$total^2, edges^2)
  .expect_allocation_bins(density_of("(mu) study: sd(intercept)", conditional = TRUE),
                          draws$components[draws$gates[, 1L] == 1, 1L], edges)

  # a mean-variance SD-component allocation: sd(x) = T sqrt(2 w_x), w_x ~
  # Beta(3, 2) with k = K = 2
  mean_variance <- .random_effects_mean_variance_allocation_fit()
  selection <- parameter_catalog_resolve(parameter_catalog(mean_variance), "(mu) sd(x)", namespace = "mu")
  density <- parameter_prior_density(mean_variance, selection)
  expect_identical(attr(density, "adaptive_evaluation")$kind, "allocation_product")
  expect_equal(
    exp(prior_density_ordinate(density, .4)$log_density),
    stats::integrate(function(s) f_T(.4 / sqrt(2 * s)) * stats::dbeta(s, 3, 2) / sqrt(2 * s),
                     0, 1, rel.tol = 1e-13)$value,
    tolerance = 1e-8
  )
  eta <- matrix(stats::rgamma(2e6, shape = rep(c(2, 3), each = 1e6)), ncol = 2L)
  sd_x <- stats::rgamma(1e6, 2, 2) * sqrt(2 * eta[, 2L] / rowSums(eta))
  .expect_allocation_bins(density, sd_x, edges)
})

test_that("SD components of nested allocations keep a plotting density refused for point tests", {

  skip_if_not_installed("runjags")

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug  = factor(c("a", "b", "a", "b")),
    x     = c(-1, 0, 1, 2)
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      root = random_variance_allocation(
        name = "root", terms = c(study = "study", drug = "drug"),
        sd = prior("gamma", list(2, 2))
      ),
      split = random_variance_allocation(
        name = "split", terms = "study", target = "sd_component",
        parent = allocation_ref("root", "study")
      )
    )
  )
  prior_list <- formula_result$prior_list
  columns <- unique(unlist(lapply(names(prior_list), function(parameter){
    BayesTools:::.prior_linear_prior_columns(parameter, prior_list[[parameter]])
  })))
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(matrix(.5, 4L, length(columns),
                                                  dimnames = list(NULL, columns)))),
         sample = 4L),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)
  catalog <- parameter_catalog(fit)
  nested <- catalog$quantities$canonical_name[grepl("sd(x)", catalog$quantities$canonical_name, fixed = TRUE)]
  expect_length(nested, 1L)
  density <- parameter_prior_density(fit, parameter_catalog_resolve(catalog, nested))
  expect_null(attr(density, "adaptive_evaluation"))
  ordinate <- prior_density_ordinate(density, .4)
  expect_identical(ordinate$behavior, "unknown")
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "nested variance allocation", fixed = TRUE)
  status <- prior_ordinate_status(density, .4)
  expect_identical(status$condition, "BayesTools_inexact_ordinate")
  expect_match(status$reason, "SD component of a nested variance allocation", fixed = TRUE)
  # the drug block, one Dirichlet share of the root, stays exact
  drug <- parameter_prior_density(
    fit, parameter_catalog_resolve(catalog, "(mu) drug: sd(intercept)")
  )
  expect_true(prior_density_ordinate(drug, .4)$exact)
})

test_that("gate-only allocations and inclusion indicators have exact priors and declared atoms", {

  skip_if_not_installed("runjags")

  # a gate-only allocation: the study SD is T g, T ~ gamma(2, 2), g ~
  # Bernoulli(0.4); the gate draws are (0, 1, 0, 1)
  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(allocation = random_variance_allocation(
      name = "allocation", terms = "study", sd = prior("gamma", list(2, 2)),
      inclusion = list(study = prior("spike", list(location = .4)))
    ))
  )
  allocation <- formula_result$formula_design$random_allocations[[1L]]
  gate <- vapply(allocation$inclusion, `[[`, character(1), "indicator_name")
  samples <- cbind(c(2, 1.5, 3, 1), c(0, 1, 0, 1))
  colnames(samples) <- c(allocation$source_node, gate)
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)
  catalog <- parameter_catalog(fit)
  resolve <- function(name) parameter_catalog_resolve(catalog, name)

  sd <- parameter_prior_density(fit, resolve("(mu) sd(intercept)"))
  expect_identical(attr(sd, "adaptive_evaluation")$kind, "allocation_product")
  expect_equal(sd$points$p[sd$points$x == 0], .6, tolerance = 1e-12)
  ordinate <- prior_density_ordinate(sd, .4)
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density), .4 * stats::dgamma(.4, 2, 2), tolerance = 1e-12)
  variance <- parameter_prior_density(fit, resolve("(mu) var(intercept)"))
  expect_equal(exp(prior_density_ordinate(variance, .16)$log_density),
               .4 * stats::dgamma(.4, 2, 2) / (2 * .4), tolerance = 1e-12)
  region <- list(intervals = matrix(c(0, 1), 1L), indicator = function(x) x > 0 & x <= 1)
  expect_equal(as.numeric(.prior_linear_density_region_probability(sd, region)),
               .4 * stats::pgamma(1, 2, 2), tolerance = 1e-10)

  mixed <- parameter_mixed_posterior(fit, resolve("(mu) sd(intercept)"))
  expect_equal(as.numeric(mixed), c(0, 1.5, 0, 1))
  atoms <- posterior_metadata(mixed, "atoms")
  expect_equal(as.numeric(atoms$locations[, 1L]), 0)
  expect_equal(atoms$mass, .5)
  included <- parameter_mixed_posterior(fit, resolve("(mu) sd(intercept)"), conditional = TRUE)
  expect_equal(as.numeric(included), c(1.5, 1))
  expect_true(posterior_atoms_free(included))
  expect_equal(
    exp(prior_density_ordinate(posterior_metadata(included, "prior_density"), .4)$log_density),
    stats::dgamma(.4, 2, 2), tolerance = 1e-12
  )

  # the inclusion indicator: Bernoulli(0.4) prior, atoms from the indicator
  inclusion <- parameter_prior_density(fit, resolve("(mu) allocation: inclusion(study)"))
  expect_equal(inclusion$points$x, c(0, 1))
  expect_equal(inclusion$points$p, c(.6, .4), tolerance = 1e-12)
  expect_identical(prior_density_ordinate(inclusion, 1)$behavior, "point_mass")
  indicator <- parameter_mixed_posterior(fit, resolve("(mu) allocation: inclusion(study)"))
  indicator_atoms <- posterior_metadata(indicator, "atoms")
  expect_equal(as.numeric(indicator_atoms$locations[, 1L]), c(0, 1))
  expect_equal(indicator_atoms$mass, c(.5, .5))
})

test_that("totals with an ungated component cannot be conditioned on inclusion", {

  skip_if_not_installed("runjags")

  # study and drug gated, site ungated: some component is always active, so
  # the inclusion event of the total is certain; a declared condition on the
  # study and drug gates would not describe the kept draws or the prior
  fit <- .allocation_prior_fit(
    sd = prior("gamma", list(2, 2)),
    alpha = c(1, 2, 3),
    inclusion = list(study = prior("spike", list(.5)), drug = prior("spike", list(.4)))
  )
  catalog <- parameter_catalog(fit)
  for(name in c("(mu) allocation: sd_total", "(mu) allocation: var_total")){
    expect_error(
      parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, name), conditional = TRUE),
      paste0(
        "The inclusion event of '", name, "' is unavailable: the quantity has ",
        "no inclusion gate. Use 'conditional = FALSE'."
      ),
      fixed = TRUE
    )
    total <- parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, name))
    expect_true(posterior_metadata(total, "condition")$averaged)
    expect_true(posterior_atoms_free(total))
  }
  # a gated component of the same allocation still conditions on its gate
  study <- parameter_mixed_posterior(
    fit, parameter_catalog_resolve(catalog, "(mu) study: sd(intercept)"), conditional = TRUE
  )
  expect_identical(
    posterior_metadata(study, "condition")$conditional,
    "mu__xRE_ALLOCx_allocation__include_study_indicator"
  )

  # with every component gated the event is the OR of the gates
  gated <- .allocation_prior_fit(
    sd = prior("gamma", list(2, 2)),
    alpha = c(1, 2),
    inclusion = list(study = prior("spike", list(.5)), drug = prior("spike", list(.4)))
  )
  total <- parameter_mixed_posterior(
    gated, parameter_catalog_resolve(parameter_catalog(gated), "(mu) allocation: sd_total"),
    conditional = TRUE
  )
  expect_identical(
    posterior_metadata(total, "condition")[c("conditional", "conditional_rule", "averaged")],
    list(
      conditional = c("mu__xRE_ALLOCx_allocation__include_study_indicator",
                      "mu__xRE_ALLOCx_allocation__include_drug_indicator"),
      conditional_rule = "OR",
      averaged = FALSE
    )
  )
})
