skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Random-effect evaluation
# ============================================================================ #
#
# PURPOSE:
#   JAGS_evaluate_formula() reconstruction of monitored and latent random effects.
# ============================================================================ #

test_that("JAGS_evaluate_formula evaluates monitored observed random effects", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(
      10, 1,  0.5, 0.1, -0.5, -0.2,
      20, 2,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expected <- cbind(
    c(9.4, 10.5, 10.3, 11.1),
    c(18.4, 20.6, 21.1, 22.8)
  )
  expect_equal(unname(prediction), expected)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = data.frame(x = 0, id = factor("c", levels = "c")),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "New random-effect level",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x +
        diag(1 + x | id) +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = transform(df, drug = factor(c("d1", "d1", "d2", "d2"))),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "do not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + random(1 + x | drug, name = "id", covariance = "diag"),
      parameter = "mu",
      data = transform(df, drug = id),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + us(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + diag(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )

  slope_only_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  slope_only_fit <- coda::mcmc(matrix(
    c(0, 1, 2),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[2,1]"
    ))
  ))
  attr(slope_only_fit, "formula_design") <- list(mu = slope_only_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = slope_only_fit,
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = transform(df, x = c(-1, NA, 1, 2)),
      prior_list = slope_only_result$prior_list,
      formula_target = "conditional"
    ),
    "missing predictor values",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula reconstructs observed random effects from latent samples", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  sd_prior <- prior("gamma", list(2, 2))
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  expect_true("mu__xREx__id_xRE_Zx" %in% formula_result$add_parameters)
  expect_false("mu__xREx__id_xRE_COEFx" %in% formula_result$add_parameters)

  posterior <- matrix(
    c(
      10, 1, 2, 3,  0.5, 0.1, -0.5, -0.2,
      20, 2, 4, 5,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[1,2]",
      "mu__xREx__id_xRE_Zx[2,1]",
      "mu__xREx__id_xRE_Zx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expected <- cbind(
    c(9.7, 11.0, 9.4, 9.8),
    c(19.4, 22.4, 18.1, 18.6)
  )
  expect_equal(unname(prediction), expected)

  incomplete_fit <- coda::mcmc(posterior[, colnames(posterior) != "mu__xREx__id_xRE_Zx[2,2]", drop = FALSE])
  attr(incomplete_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = incomplete_fit,
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "cannot be reconstructed from the posterior samples",
    fixed = TRUE,
    class = "BayesTools_refit_monitoring"
  )
})

test_that("JAGS_evaluate_formula rejects new levels for known group covariance", {

  df <- data.frame(
    id = factor(c("a", "b", "a", "c"), levels = c("a", "b", "c"))
  )
  K <- matrix(
    c(2, .4, .2,
      .4, 3, .5,
      .2, .5, 4),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  formula_result <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 1)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  posterior <- matrix(
    c(10, 2, .1, -.2, .3),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names, z_names))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = random_effects$formula,
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(prediction)),
    10 + 2 * c(.1, -.2, .1, .3),
    tolerance = 1e-12
  )
  for(new_levels in c("zero", "sample")){
    conditional_prediction <- JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional",
      new_levels = new_levels
    )
    expect_equal(conditional_prediction, prediction, tolerance = 1e-12)
  }

  new_data <- data.frame(
    id = factor(c("a", "d"), levels = c("a", "b", "c", "d"))
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = new_data,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "known group covariance",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = new_data,
      prior_list = formula_result$prior_list,
      formula_target = "conditional",
      new_levels = "sample"
    ),
    "known group covariance",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula reconstructs allocated random effects from simplex coordinates", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    drug = factor(c("a", "b", "a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  posterior <- matrix(
    c(
      10, 2, 1, 3,  0.1, 0.2, 0.3,  1, 2,
      20, 4, 2, 2,  0.5, 0.6, 0.7, -1, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[3,1]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  first_draw <- 10 +
    c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
    sqrt(3) * c(1, 2, 1, 2, 1, 2)
  second_sd <- 4 * sqrt(0.5)
  second_draw <- 20 +
    second_sd * c(0.5, 0.5, 0.6, 0.6, 0.7, 0.7) +
    second_sd * c(-1, 1, -1, 1, -1, 1)
  expect_equal(unname(prediction), unname(cbind(first_draw, second_draw)), tolerance = 1e-12)

  public_simplex <- matrix(
    c(
      5, 2, 0.25, 0.75,  0.1, 0.2, 0.3,  1, 2
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[3,1]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  public_fit <- coda::mcmc(public_simplex)
  attr(public_fit, "formula_design") <- list(mu = formula_result$formula_design)
  public_prediction <- JAGS_evaluate_formula(
    fit = public_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(public_prediction)),
    5 +
      c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
      sqrt(3) * c(1, 2, 1, 2, 1, 2),
    tolerance = 1e-12
  )

  boundary_simplex <- cbind(
    public_simplex,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 0,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 1
  )
  boundary_simplex[, "mu__xRE_ALLOCx_allocation__weight[1]"] <- 0
  boundary_simplex[, "mu__xRE_ALLOCx_allocation__weight[2]"] <- 1
  boundary_fit <- coda::mcmc(boundary_simplex)
  attr(boundary_fit, "formula_design") <- list(mu = formula_result$formula_design)
  boundary_prediction <- JAGS_evaluate_formula(
    fit = boundary_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(boundary_prediction)),
    5 + 2 * c(1, 2, 1, 2, 1, 2),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula reconstructs row-indexed external SD random effects", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  posterior <- matrix(
    c(
      10,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    fitted_rows = seq_len(nrow(df)),
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(prediction)),
    10 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
    ),
    tolerance = 1e-12
  )

  multi_posterior <- matrix(
    c(
      10,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4,

      -5,
      1, 10, 100, 1000,
      4, 1,
      -2, 5,
      7, -11
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, colnames(posterior))
  )
  multi_fit <- coda::mcmc(multi_posterior)
  attr(multi_fit, "formula_design") <- list(mu = formula_result$formula_design)
  multi_prediction <- JAGS_evaluate_formula(
    fit = multi_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    fitted_rows = seq_len(nrow(df)),
    formula_target = "conditional"
  )
  expect_equal(dim(multi_prediction), c(4L, 2L))
  expect_equal(
    unname(multi_prediction),
    unname(cbind(
      10 + c(
        2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
        4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
        6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
        8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
      ),
      -5 + c(
        1    * (sqrt(4 / 5) * -2 + sqrt(1 / 5) *   7),
        10   * (sqrt(4 / 5) * -2 + sqrt(1 / 5) * -11),
        100  * (sqrt(4 / 5) *  5 + sqrt(1 / 5) *   7),
        1000 * (sqrt(4 / 5) *  5 + sqrt(1 / 5) * -11)
      )
    )),
    tolerance = 1e-12
  )

  newdata <- data.frame(
    study = factor(c("s2", "s1"), levels = c("s1", "s2")),
    drug = factor(c("b", "a"), levels = c("a", "b"))
  )
  new_posterior <- posterior[, c(
    "mu_intercept",
    "tau[1]",
    "tau[4]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
    "mu__xREx__study_xRE_Zx[1,1]",
    "mu__xREx__study_xRE_Zx[2,1]",
    "mu__xREx__drug_xRE_Zx[1,1]",
    "mu__xREx__drug_xRE_Zx[2,1]"
  ), drop = FALSE]
  new_posterior[, "tau[1]"] <- 20
  new_posterior[, "tau[4]"] <- 10
  new_fit <- coda::mcmc(new_posterior)
  attr(new_fit, "formula_design") <- list(mu = formula_result$formula_design)
  new_prediction <- JAGS_evaluate_formula(
    fit = new_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list,
    fitted_rows = c(4L, 1L),
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(new_prediction)),
    10 + c(
      10 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4),
      20 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3)
    ),
    tolerance = 1e-12
  )

  sd_component_df <- data.frame(
    x = c(-1, 0, 1, 2),
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  sd_component_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = sd_component_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  sd_component_posterior <- matrix(
    c(
      5,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,2]"
    ))
  )
  sd_component_fit <- coda::mcmc(sd_component_posterior)
  attr(sd_component_fit, "formula_design") <- list(mu = sd_component_result$formula_design)
  sd_component_prediction <- JAGS_evaluate_formula(
    fit = sd_component_fit,
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = sd_component_df,
    prior_list = sd_component_result$prior_list,
    fitted_rows = seq_len(nrow(sd_component_df)),
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(sd_component_prediction)),
    5 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
    ),
    tolerance = 1e-12
  )

  missing_source_fit <- coda::mcmc(
    posterior[, colnames(posterior) != "tau[2]", drop = FALSE]
  )
  attr(missing_source_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = missing_source_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df[1:2, , drop = FALSE],
      prior_list = formula_result$prior_list,
      fitted_rows = 1:2,
      formula_target = "conditional"
    ),
    "missing values for prediction row(s): 2",
    fixed = TRUE
  )

  scalar_source_posterior <- posterior[, c(
    "mu_intercept",
    "tau[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
    "mu__xREx__study_xRE_Zx[1,1]",
    "mu__xREx__study_xRE_Zx[2,1]",
    "mu__xREx__drug_xRE_Zx[1,1]",
    "mu__xREx__drug_xRE_Zx[2,1]"
  ), drop = FALSE]
  colnames(scalar_source_posterior)[colnames(scalar_source_posterior) == "tau[1]"] <- "tau"
  scalar_source_fit <- coda::mcmc(scalar_source_posterior)
  attr(scalar_source_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = scalar_source_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df[1, , drop = FALSE],
      prior_list = formula_result$prior_list,
      fitted_rows = 1L,
      formula_target = "conditional"
    ),
    "missing values for prediction row(s): 1",
    fixed = TRUE
  )

  coefficient_names <- as.vector(BayesTools:::.bt_random_effect_coefficient_names(
    random_term = formula_result$formula_design$random_effects[[1]],
    n_groups = formula_result$formula_design$random_effects[[1]]$n_groups,
    n_columns = formula_result$formula_design$random_effects[[1]]$n_columns
  ))
  coefficient_only_posterior <- cbind(
    posterior[, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
    ), drop = FALSE],
    matrix(
      seq_along(coefficient_names),
      nrow = 1,
      dimnames = list(NULL, coefficient_names)
    )
  )
  coefficient_only_fit <- coda::mcmc(coefficient_only_posterior)
  attr(coefficient_only_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = coefficient_only_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      fitted_rows = seq_len(nrow(df)),
      formula_target = "conditional"
    ),
    "cannot be reconstructed from the posterior samples",
    fixed = TRUE,
    class = "BayesTools_refit_monitoring"
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = data.frame(
        study = factor("s3", levels = "s3"),
        drug = factor("a", levels = c("a", "b"))
      ),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "New random-effect level",
    fixed = TRUE
  )

  nested_df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | estimate, name = "estimate", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = nested_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  nested_posterior <- matrix(
    c(
      10,
      0.25, 0.75,
      2, 3, 4, 5,
      0.1, 0.2, -0.1, 0.4,
      1, -0.5
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "mu__xREx__estimate_xRE_Zx[1,1]",
      "mu__xREx__estimate_xRE_Zx[2,1]",
      "mu__xREx__estimate_xRE_Zx[3,1]",
      "mu__xREx__estimate_xRE_Zx[4,1]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]"
    ))
  )
  nested_fit <- coda::mcmc(nested_posterior)
  attr(nested_fit, "formula_design") <- list(mu = nested_result$formula_design)
  nested_prediction <- JAGS_evaluate_formula(
    fit = nested_fit,
    formula = ~ 1 +
      random(1 | estimate, name = "estimate", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = nested_df,
    prior_list = nested_result$prior_list,
    fitted_rows = seq_len(nrow(nested_df)),
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(nested_prediction)),
    10 + c(
      2 * (sqrt(0.25) *  0.1 + sqrt(0.75) *  1),
      3 * (sqrt(0.25) *  0.2 + sqrt(0.75) *  1),
      4 * (sqrt(0.25) * -0.1 + sqrt(0.75) * -0.5),
      5 * (sqrt(0.25) *  0.4 + sqrt(0.75) * -0.5)
    ),
    tolerance = 1e-12
  )

  slope_df <- data.frame(
    x = c(1, 2, 3, 4),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  slope_result <- JAGS_formula(
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  slope_newdata <- data.frame(
    x = c(2, 4),
    study = factor(c("s1", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b"), levels = c("a", "b"))
  )
  slope_posterior <- matrix(
    c(
      0, 0,
      2, 3,
      0.25, 0.75,
      1, 2,
      -1, 0.5,
      0, 0
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "tau[2]",
      "tau[4]",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  slope_fit <- coda::mcmc(slope_posterior)
  attr(slope_fit, "formula_design") <- list(mu = slope_result$formula_design)
  attr(slope_fit, "formula_scale") <- list(mu = slope_result$formula_scale)
  slope_prediction <- JAGS_evaluate_formula(
    fit = slope_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_result$prior_list,
    fitted_rows = c(2L, 4L),
    formula_target = "conditional"
  )
  scaled_x <- (slope_newdata$x - slope_result$formula_scale$mu_x$mean) /
    slope_result$formula_scale$mu_x$sd
  expect_equal(
    unname(drop(slope_prediction)),
    c(
      2 * sqrt(0.25) * (1 + scaled_x[1] * 2),
      3 * sqrt(0.25) * (-1 + scaled_x[2] * 0.5)
    ),
    tolerance = 1e-12
  )
  slope_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$base_tau * data$x[seq_len(n_rows)]
    },
    inputs = "base_tau"
  )
  slope_values_result <- JAGS_formula(
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(slope_value_source),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  slope_values_posterior <- cbind(
    slope_posterior[, !colnames(slope_posterior) %in% c("tau[1]", "tau[2]"), drop = FALSE],
    "base_tau" = 1
  )
  slope_values_fit <- coda::mcmc(slope_values_posterior)
  attr(slope_values_fit, "formula_design") <- list(mu = slope_values_result$formula_design)
  attr(slope_values_fit, "formula_scale") <- list(mu = slope_values_result$formula_scale)
  slope_values_prediction <- JAGS_evaluate_formula(
    fit = slope_values_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_values_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(slope_values_prediction)),
    c(
      2 * sqrt(0.25) * (1 + scaled_x[1] * 2),
      4 * sqrt(0.25) * (-1 + scaled_x[2] * 0.5)
    ),
    tolerance = 1e-12
  )
  slope_values_stale_fit <- coda::mcmc(cbind(
    slope_values_posterior,
    "tau[1]" = 999,
    "tau[2]" = 999
  ))
  attr(slope_values_stale_fit, "formula_design") <- list(mu = slope_values_result$formula_design)
  attr(slope_values_stale_fit, "formula_scale") <- list(mu = slope_values_result$formula_scale)
  slope_values_stale_prediction <- JAGS_evaluate_formula(
    fit = slope_values_stale_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_values_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(slope_values_stale_prediction)),
    unname(drop(slope_values_prediction)),
    tolerance = 1e-12
  )

  cs_df <- data.frame(
    f = factor(c("a", "b", "a", "b")),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  cs_result <- JAGS_formula(
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = cs_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      study = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  cs_posterior <- matrix(
    c(
      0,
      2, 3,
      0.25, 0.75,
      0,
      1, 2,
      -1, 0.5,
      0, 0
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_rho_z",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  cs_fit <- coda::mcmc(cs_posterior)
  attr(cs_fit, "formula_design") <- list(mu = cs_result$formula_design)
  cs_prediction <- JAGS_evaluate_formula(
    fit = cs_fit,
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = cs_df[1:2, , drop = FALSE],
    prior_list = cs_result$prior_list,
    fitted_rows = 1:2,
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(cs_prediction)),
    c(
      2 * sqrt(0.25) * 1,
      3 * sqrt(0.25) * 2
    ),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula reconstructs correlated latent random effects", {

  df <- data.frame(
    x = c(0, 1, 2, -1),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 2),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  expect_true("mu__xREx__study_xRE_Zx" %in% formula_result$add_parameters)
  expect_true("mu__xREx__study_xRE_CORx_L" %in% formula_result$add_parameters)
  expect_false("mu__xREx__study_xRE_COEFx" %in% formula_result$add_parameters)

  L22 <- sqrt(1 - 0.5^2)
  posterior <- matrix(
    c(
      0, 0,
      2, 3,
      1, 0, 0.5, L22,
      1, 2, -1, 0.25
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_intercept",
      "mu__xREx__study_x",
      "mu__xREx__study_xRE_CORx_L[1,1]",
      "mu__xREx__study_xRE_CORx_L[1,2]",
      "mu__xREx__study_xRE_CORx_L[2,1]",
      "mu__xREx__study_xRE_CORx_L[2,2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  coef_a_x <- 3 * (0.5 * 1 + L22 * 2)
  coef_b_x <- 3 * (0.5 * -1 + L22 * 0.25)
  expect_equal(
    unname(drop(prediction)),
    c(2, 2 + coef_a_x, -2 + 2 * coef_b_x, -2 - coef_b_x),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula reconstructs point-SD latent random effects", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(1, 0.25, -0.5),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expect_equal(unname(drop(prediction)), c(1.5, 0))
})

test_that("JAGS_evaluate_formula scales random-slope prediction data", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  scale_info <- formula_result$formula_scale$mu_x
  posterior <- matrix(
    c(0, 0, 2, 1, -1),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  newdata <- data.frame(
    x = c(2, 4),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expected_x <- (newdata$x - scale_info$mean) / scale_info$sd
  expect_equal(unname(drop(prediction)), expected_x * c(2, -2))
})

test_that("structured random-effect prediction indexes raw data when fixed predictors are scaled", {

  df <- data.frame(
    time = rep(c(1, 2, 3), 2),
    id = factor(rep(c("a", "b"), each = 3), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))

  for (structure in c("cs", "hcs", "ar1", "car", "har")) {
    structured_formula <- stats::as.formula(
      if (identical(structure, "car")) {
        "~ 1 + time + car(0 + time | id)"
      } else {
        paste0("~ 1 + time + ", structure, "(time | id)")
      }
    )
    formula_result <- JAGS_formula(
      formula = structured_formula,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        time = prior("normal", list(0, 1))
      ),
      formula_scale = list(time = TRUE),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          cor = rho_prior,
          monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
        )
      )
    )
    random_term <- formula_result$formula_design$random_effects[[1]]
    coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
      random_term = random_term,
      n_groups = length(random_term$group_levels),
      n_columns = random_term$n_columns
    )
    target_column_name <- if (identical(structure, "car")) "time_2" else "time2"
    target_column <- match(target_column_name, random_term$column_names)
    expect_false(is.na(target_column), info = structure)

    posterior_names <- c("mu_intercept", "mu_time", as.vector(coefficient_names))
    posterior <- matrix(
      0,
      nrow = 1,
      ncol = length(posterior_names),
      dimnames = list(NULL, posterior_names)
    )
    posterior[, coefficient_names[match("b", random_term$group_levels), target_column]] <- 42
    fit <- coda::mcmc(posterior)
    attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
    attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

    prediction <- JAGS_evaluate_formula(
      fit = fit,
      formula = structured_formula,
      parameter = "mu",
      data = data.frame(time = 2, id = factor("b", levels = c("a", "b"))),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    )
    expect_equal(unname(drop(prediction)), 42, info = structure)
  }
})

test_that("transform_scale_samples leaves structured random-effect SDs on the fitted scale", {

  df <- data.frame(
    time = rep(c(1, 2, 3), 2),
    id = factor(rep(c("a", "b"), each = 3), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))

  for (structure in c("hcs", "har")) {
    formula_result <- JAGS_formula(
      formula = stats::as.formula(paste0("~ 1 + time + ", structure, "(time | id)")),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        time = prior("normal", list(0, 1))
      ),
      formula_scale = list(time = TRUE),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          cor = rho_prior
        )
      )
    )
    random_term <- formula_result$formula_design$random_effects[[1]]
    sd_cols <- unique(random_term$sd_parameter_names)
    posterior <- matrix(
      seq_along(sd_cols),
      nrow = 1,
      dimnames = list(NULL, sd_cols)
    )
    transformed <- BayesTools:::.bt_transform_scale_posterior(
      posterior,
      list(mu = formula_result$formula_scale)
    )
    expect_equal(
      unname(transformed[1, sd_cols]),
      seq_along(sd_cols),
      info = structure
    )
  }
})

test_that("JAGS_evaluate_formula rejects duplicate requested random-effect block names", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id) + random(1 | id, name = "id", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "Random-effect block names in the supplied formula must be unique.",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula rejects malformed fitted random-effect metadata without structure", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  malformed_design <- formula_result$formula_design
  malformed_design$random_effects[[1]]$structure <- NULL
  malformed_design$random_effects[[1]]$covariance <- "diag"

  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "formula_design") <- list(mu = malformed_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "missing canonical 'random_term\\$structure'"
  )

  malformed_design <- formula_result$formula_design
  malformed_design$random_effects[[1]]$homogeneous_sd <- NULL
  attr(fit, "formula_design") <- list(mu = malformed_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )
})

test_that("JAGS_evaluate_formula sums multiple monitored random-effect blocks", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    drug = factor(c("d1", "d2", "d1", "d2"), levels = c("d1", "d2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2))),
      drug_slope = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  posterior <- matrix(
    c(
      1, 2,
      0.5, -0.5,
      0.1, 0.2
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_xRE_COEFx[1,1]",
      "mu__xREx__study_xRE_COEFx[2,1]",
      "mu__xREx__drug_slope_xRE_COEFx[1,1]",
      "mu__xREx__drug_slope_xRE_COEFx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expect_equal(unname(drop(prediction)), c(1.5, 3.7, 4.7, 7.1))

  study_prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expect_equal(unname(drop(study_prediction)), c(1.5, 3.5, 4.5, 6.5))
})
