skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Random-effect scaling and back-transformation
# ============================================================================ #
#
# PURPOSE:
#   Predictor scaling in random effects and the back-transformation of random-effect
#   SDs and correlations.
# ============================================================================ #

test_that("JAGS_formula scales predictors used only in random effects", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  scale_info <- formula_result$formula_scale$mu_x
  expect_equal(scale_info$mean, mean(df$x))
  expect_equal(scale_info$sd, stats::sd(df$x))
  expect_equal(
    unname(formula_result$data$mu__xREx__id_xRE_DATAx[, "x"]),
    unname((df$x - scale_info$mean) / scale_info$sd)
  )
  scale_leaves <- attr(formula_result$formula_scale, "random_effect_sd_leaves")
  expect_s3_class(scale_leaves[["__xREx__id"]], "BayesTools_random_effect_sd_leaves")
  expect_equal(scale_leaves[["__xREx__id"]]$leaf_terms, c(mu__xREx__id_x = "x"))

  posterior <- matrix(
    c(0, 2, 1, -1),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
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
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expected_x <- (newdata$x - scale_info$mean) / scale_info$sd
  expect_equal(unname(drop(prediction)), expected_x * c(2, -2))
})

test_that("random-effect grouping variables use raw data when predictors are scaled", {

  df <- data.frame(
    x = c(-1, 0.1, 0.2, 100),
    xpos = c(FALSE, TRUE, TRUE, TRUE),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 | xpos, name = "xpos", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      xpos = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )

  raw_group_map <- as.numeric(factor(df$xpos, levels = levels(as.factor(df$xpos))))
  expect_equal(formula_result$formula_design$random_effects[[1]]$group_map, raw_group_map)

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = formula_result$formula_design$random_effects[[1]],
    n_groups = 2,
    n_columns = 1
  )
  posterior <- matrix(
    c(0, 0, 10, 20),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + random(1 | xpos, name = "xpos", covariance = "diag"),
    parameter = "mu",
    data = data.frame(x = 0.1, xpos = TRUE),
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(unname(drop(prediction)), 20)
})

test_that("random-effect grouping factors preserve unused training levels", {

  # Declared levels without fitting rows stay fitted groups, but prediction
  # routes them through the new-level policy (maintainer decision M10 (a)).
  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b", "unused"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
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
  random_term <- formula_result$formula_design$random_effects[[1]]
  expect_equal(random_term$group_levels, c("a", "b", "unused"))
  expect_equal(random_term$group_observed_levels, c("a", "b"))
  expect_equal(random_term$n_groups, 3L)

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(0, 10, 20, 30),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  unused_data <- data.frame(
    id = factor(c("b", "unused"), levels = c("a", "b", "unused"))
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = unused_data,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    paste0(
      "New random-effect level(s) for block 'id' are not supported by ",
      "JAGS_evaluate_formula(): unused. Declared grouping level(s) unused ",
      "have no rows in the fitting data and are predicted as new levels. ",
      "Use 'new_levels = \"zero\"' or 'new_levels = \"sample\"' to predict ",
      "them."
    ),
    fixed = TRUE
  )
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = unused_data,
    prior_list = formula_result$prior_list,
    formula_target = "conditional",
    new_levels = "zero"
  )
  expect_equal(unname(drop(prediction)), c(20, 0))
})

test_that("random effects evaluate on non-mu formula parameters", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "sigma",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  expect_equal(
    names(formula_result$prior_list),
    c("sigma_intercept", "sigma__xREx__id_intercept")
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(1, 10, 20),
    nrow = 1,
    dimnames = list(NULL, c("sigma_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(sigma = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "sigma",
    data = data.frame(id = factor(c("a", "b"))),
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(unname(drop(prediction)), c(11, 21))
})

test_that("character random-factor predictors keep fitted levels for prediction", {

  df <- data.frame(
    f = c("b", "a", "c", "a", "b", "c"),
    id = factor(c("g1", "g1", "g1", "g2", "g2", "g2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
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
  random_term <- formula_result$formula_design$random_effects[[1]]
  expect_equal(random_term$xlevels$f, c("a", "b", "c"))

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(0, seq_along(coefficient_names)),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(f = c("c", "a"), id = factor(c("g1", "g2"))),
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expected <- c(
    posterior[1, coefficient_names[1, 2]],
    0
  )
  expect_equal(unname(drop(prediction)), unname(expected))

  prediction_extra_level <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(
      f = factor(c("c", "a"), levels = c("a", "b", "c", "unused")),
      id = factor(c("g1", "g2"))
    ),
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )
  expect_equal(unname(drop(prediction_extra_level)), unname(expected))
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data.frame(
        f = factor(c("unused", "a"), levels = c("a", "b", "c", "unused")),
        id = factor(c("g1", "g2"))
      ),
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "do not match the levels",
    fixed = TRUE
  )
})

.canonical_us_sd_leaves <- function(){

  leaf_names <- c(
    "mu__xREx__id_intercept",
    "mu__xREx__id_x"
  )
  out <- list(
    random_structure = "us",
    leaf_names_by_column = leaf_names,
    leaf_terms_by_column = c("intercept", "x"),
    leaf_names = leaf_names,
    leaf_terms = stats::setNames(c("intercept", "x"), leaf_names),
    column_names = c("(Intercept)", "x")
  )
  class(out) <- c("BayesTools_random_effect_sd_leaves", "list")
  list(`__xREx__id` = out)
}

test_that("transform_scale_samples unscales correlated random-effect SDs with covariance", {

  posterior <- matrix(
    c(1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  # the standardization of x with the fitted design of ~ x (required by
  # original-scale transforms); the random-effect SD structure under test is
  # set by hand
  formula_scale <- list(
    mu = formula_scale_for_test(~ x, list(x = list(mean = 5, sd = 1)))
  )
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

  transformed <- BayesTools:::.bt_transform_scale_posterior(posterior, formula_scale)

  expect_equal(
    unname(transformed[1, "mu__xREx__id_intercept"]),
    sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2),
    tolerance = 1e-12
  )
  expect_equal(unname(transformed[1, "mu__xREx__id_x"]), 2, tolerance = 1e-12)
  expected_cor <- (0.8 * 1 * 2 - 5 * 2^2) /
    (sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2) * 2)
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[2,1]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expected_L <- t(chol(matrix(c(1, expected_cor, expected_cor, 1), 2, 2)))
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_L[2,1]"]),
    expected_L[2, 1],
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_L[2,2]"]),
    expected_L[2, 2],
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples requires correlations for scaled correlated random-effect SDs", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_cols <- random_term$sd_parameter_names
  posterior <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, sd_cols)
  )

  expect_error(
    BayesTools:::.bt_transform_scale_posterior(posterior, list(mu = formula_result$formula_scale)),
    "requires random-effect correlation samples",
    fixed = TRUE
  )

  R_names <- outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column){
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  partial_R_posterior <- cbind(posterior, 1)
  colnames(partial_R_posterior)[ncol(partial_R_posterior)] <- R_names[1, 1]
  expect_error(
    BayesTools:::.bt_transform_scale_posterior(partial_R_posterior, list(mu = formula_result$formula_scale)),
    "correlation samples are incomplete",
    fixed = TRUE
  )

  # fitted-scale summaries read the catalog's correlation coordinates; the original-scale
  # table first transforms the block, which needs the monitored correlation samples
  expect_error(
    .parameter_catalog_random_summary_samples(
      model_samples = posterior,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design),
      mode = "standard",
      formula_scale = list(mu = formula_result$formula_scale),
      transform_scaled = FALSE
    ),
    "missing or invalid canonical correlation coordinates",
    fixed = TRUE
  )
  expect_error(
    .parameter_catalog_random_summary_samples(
      model_samples = posterior,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design),
      mode = "standard",
      formula_scale = list(mu = formula_result$formula_scale)
    ),
    "requires random-effect correlation samples",
    fixed = TRUE
  )
})

test_that("JAGS_estimates_table backtransforms random-effect correlations", {

  testthat::skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(4, 5, 6),
    id = factor(c("a", "a", "b"))
  )
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_names <- random_term$sd_parameter_names
  R_names <- outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column) {
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L)
  source_sd <- c(1, 2)
  source_cor <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  source_L <- t(chol(source_cor))
  posterior <- matrix(
    rep(c(0, 0, source_sd, as.vector(source_cor), as.vector(source_L)), 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu_intercept",
        "mu_x",
        sd_names,
        as.vector(R_names),
        as.vector(L_names)
      )
    )
  )

  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  semantic_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = FALSE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  original_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  original_table <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  source_cov <- diag(source_sd, nrow = 2) %*% source_cor %*%
    diag(source_sd, nrow = 2)
  expected_cov <- M %*% source_cov %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov[1, 2] / prod(expected_sd)

  # transform_scaled = FALSE keeps the fitted-scale SDs and correlation
  expect_equal(
    unname(semantic_samples[, "(mu) sd(intercept)"]),
    rep(source_sd[1], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(semantic_samples[, "(mu) sd(x)"]),
    rep(source_sd[2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(semantic_samples[, "(mu) cor(intercept,x)"]),
    rep(source_cor[1, 2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(intercept)"]),
    rep(expected_sd[1], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(x)"]),
    rep(expected_sd[2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) cor(intercept,x)"]),
    rep(expected_cor, 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_table["(mu) cor(intercept,x)", "Mean"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_lt(expected_cor, 0)
})

test_that("backtransformed random-effect correlations use the transformed LKJ primitives", {

  testthat::skip_if_not_installed("runjags")

  # Real fits always monitor the LKJ primitives, which the Cholesky evaluator
  # prefers over the monitored L; they must follow the original-scale L.
  df <- data.frame(
    x = c(4, 5, 6, 8),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_names <- random_term$sd_parameter_names
  R_names <- outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column) {
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L)
  u_names <- random_term$correlation$primitive_names
  expect_identical(u_names, "mu__xREx__id_xRE_CORx_lkj_u[1]")

  source_sd <- rbind(c(1, 2), c(0.5, 1.5))
  source_rho <- c(0.8, -0.3)
  rows <- lapply(seq_len(2), function(draw) {
    source_cor <- matrix(c(1, source_rho[draw], source_rho[draw], 1), 2, 2)
    source_L <- t(chol(source_cor))
    c(0, 0, source_sd[draw, ], as.vector(source_cor), as.vector(source_L),
      (source_rho[draw] + 1) / 2)
  })
  posterior <- do.call(rbind, rows)
  colnames(posterior) <- c(
    "mu_intercept", "mu_x", sd_names, as.vector(R_names),
    as.vector(L_names), u_names
  )

  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  scale_info <- formula_result$formula_scale$mu_x
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  expected_cor <- vapply(seq_len(2), function(draw) {
    source_cor <- matrix(c(1, source_rho[draw], source_rho[draw], 1), 2, 2)
    source_cov <- diag(source_sd[draw, ]) %*% source_cor %*% diag(source_sd[draw, ])
    expected_cov <- M %*% source_cov %*% t(M)
    expected_cov[1, 2] / sqrt(prod(diag(expected_cov)))
  }, numeric(1))

  original_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  expect_equal(
    unname(original_samples[, "(mu) cor(intercept,x)"]),
    expected_cor,
    tolerance = 1e-12
  )
  draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(parameter_catalog(fit), "(mu) cor(intercept,x)")
  )
  expect_equal(
    as.numeric(as.matrix(draws[[1]])),
    expected_cor,
    tolerance = 1e-12
  )
  expect_true(all(abs(expected_cor - source_rho) > 0.1))

  transformed <- transform_scale_samples(fit)
  expect_equal(
    unname(transformed[, u_names]),
    (expected_cor + 1) / 2,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[, "mu__xREx__id_xRE_CORx_L[2,1]"]),
    expected_cor,
    tolerance = 1e-12
  )
})

test_that("backtransformed LKJ coordinates are consistent and missing for undefined correlations", {

  testthat::skip_if_not_installed("runjags")

  # Monitored canonical partial correlations (cpc = 2u - 1) follow the
  # rewritten primitives. A zero slope SD (e.g., an excluded spike-and-slab
  # SD) leaves the original-scale correlation undefined: that draw is missing
  # instead of failing every correlation summary.
  df <- data.frame(
    x = c(4, 5, 6, 8),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1, include_primitives = TRUE),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_names <- random_term$sd_parameter_names
  R_names <- as.vector(outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column) {
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  ))
  L_names <- as.vector(BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L))
  u_names <- random_term$correlation$primitive_names
  cpc_names <- random_term$correlation$cpc_names
  expect_identical(cpc_names, "mu__xREx__id_xRE_CORx_lkj_cpc[1]")

  source_sd <- rbind(c(1, 2), c(0.5, 1.5), c(0.7, 0))
  source_rho <- c(0.8, -0.3, 0.4)
  rows <- lapply(seq_len(3), function(draw) {
    source_cor <- matrix(c(1, source_rho[draw], source_rho[draw], 1), 2, 2)
    c(0, 0, source_sd[draw, ], as.vector(source_cor),
      as.vector(t(chol(source_cor))), (source_rho[draw] + 1) / 2,
      source_rho[draw])
  })
  posterior <- do.call(rbind, rows)
  colnames(posterior) <- c(
    "mu_intercept", "mu_x", sd_names, R_names, L_names, u_names, cpc_names
  )
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  scale_info <- formula_result$formula_scale$mu_x
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  expected_cor <- vapply(seq_len(2), function(draw) {
    source_cor <- matrix(c(1, source_rho[draw], source_rho[draw], 1), 2, 2)
    source_cov <- diag(source_sd[draw, ]) %*% source_cor %*% diag(source_sd[draw, ])
    expected_cov <- M %*% source_cov %*% t(M)
    expected_cov[1, 2] / sqrt(prod(diag(expected_cov)))
  }, numeric(1))
  expected_cor <- c(expected_cor, NA_real_)

  transformed <- transform_scale_samples(fit)
  expect_equal(unname(transformed[, L_names[2]]), expected_cor, tolerance = 1e-12)
  expect_equal(unname(transformed[, u_names]), (expected_cor + 1) / 2,
               tolerance = 1e-12)
  expect_equal(unname(transformed[, cpc_names]), expected_cor, tolerance = 1e-12)

  original_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  expect_equal(
    unname(original_samples[, "(mu) cor(intercept,x)"]),
    expected_cor,
    tolerance = 1e-12
  )
  draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(parameter_catalog(fit), "(mu) cor(intercept,x)")
  )
  expect_equal(as.numeric(as.matrix(draws[[1]])), expected_cor, tolerance = 1e-12)
})

test_that("transform_prior_samples unscales correlated random slopes from LKJ primitives", {

  testthat::skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(4, 5, 6, 8),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 2),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_names <- random_term$sd_parameter_names
  R_names <- as.vector(outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column) {
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  ))
  L_names <- as.vector(BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L))
  u_names <- random_term$correlation$primitive_names
  posterior <- matrix(
    0.5,
    nrow = 2,
    ncol = 2 + length(sd_names) + 4 + 4 + 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x", sd_names, R_names, L_names, u_names))
  )
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  raw <- transform_prior_samples(fit, n_samples = 200, seed = 7, formula_scale = list())
  fit <- attach_test_parameter_map(fit)
  transformed <- transform_prior_samples(fit, n_samples = 200, seed = 7)
  expect_true(all(c(sd_names, R_names, L_names, u_names) %in% colnames(raw)))
  expect_true(all(c(sd_names, R_names, L_names, u_names) %in% colnames(transformed)))

  # Raw prior draws: the monitored L and R follow the sampled primitives.
  raw_L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(raw[, u_names, drop = FALSE], K = 2L)
  expect_equal(unname(raw[, L_names[2]]), raw_L[, 2, 1], tolerance = 1e-12)
  expect_equal(unname(raw[, R_names[2]]), raw_L[, 2, 1], tolerance = 1e-12)

  # Original-scale draws: covariance M Sigma M' from the same raw draws.
  scale_info <- formula_result$formula_scale$mu_x
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  expected <- t(vapply(seq_len(nrow(raw)), function(draw) {
    rho <- raw[draw, R_names[2]]
    sds <- raw[draw, sd_names]
    source_cov <- diag(sds) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sds)
    expected_cov <- M %*% source_cov %*% t(M)
    expected_sd <- sqrt(diag(expected_cov))
    c(expected_sd, expected_cov[1, 2] / prod(expected_sd))
  }, numeric(3)))
  expect_equal(unname(transformed[, sd_names]), expected[, 1:2], tolerance = 1e-10)
  expect_equal(unname(transformed[, R_names[2]]), expected[, 3], tolerance = 1e-10)
  expect_equal(unname(transformed[, L_names[2]]), expected[, 3], tolerance = 1e-10)
  expect_equal(unname(transformed[, u_names]), (expected[, 3] + 1) / 2, tolerance = 1e-10)
})

test_that("transform_scale_samples updates valid random-effect correlations draw-wise", {

  posterior <- matrix(
    c(
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6,
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  # the standardization of x with the fitted design of ~ x (required by
  # original-scale transforms); the random-effect SD structure under test is
  # set by hand
  formula_scale <- list(
    mu = formula_scale_for_test(~ x, list(x = list(mean = 5, sd = 1)))
  )
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

  transformed <- BayesTools:::.bt_transform_scale_posterior(posterior, formula_scale)

  expected_cor <- (0.8 * 1 * 2 - 5 * 2^2) /
    (sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2) * 2)
  expected_L <- t(chol(matrix(c(1, expected_cor, expected_cor, 1), 2, 2)))
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[2, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[2, "mu__xREx__id_xRE_CORx_L[2,1]"]),
    expected_L[2, 1],
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples clears invalid transformed random-effect correlations", {

  posterior <- matrix(
    c(
      0, 0, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6,
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  # the standardization of x with the fitted design of ~ x (required by
  # original-scale transforms); the random-effect SD structure under test is
  # set by hand
  formula_scale <- list(
    mu = formula_scale_for_test(~ x, list(x = list(mean = 5, sd = 1)))
  )
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

  transformed <- BayesTools:::.bt_transform_scale_posterior(posterior, formula_scale)

  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]))
  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_L[2,1]"]))
  expect_false(is.na(transformed[2, "mu__xREx__id_xRE_CORx_R[1,2]"]))

  all_invalid <- BayesTools:::.bt_transform_scale_posterior(
    posterior[1, , drop = FALSE],
    formula_scale
  )
  all_invalid_cor_cols <- grep("_xRE_CORx_[RL]", colnames(all_invalid), value = TRUE)
  expect_true(all(is.na(all_invalid[, all_invalid_cor_cols, drop = FALSE])))

})

test_that("factor point priors do not add unindexed columns when unscaling", {

  testthat::skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 4, 6, 2, 8, 5, 3, 7, 9),
    f = factor(rep(c("a", "b", "c"), 3L), levels = c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    formula = ~ x * f,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 5)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("spike", list(0.3), contrast = "treatment"),
      "x:f" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    formula_scale = list(x = TRUE)
  )
  expect_identical(
    attr(formula_result$formula_scale, "point_terms"),
    c(`mu_f[1]` = 0.3, `mu_f[2]` = 0.3)
  )

  # The spiked factor coefficients are monitored as indexed JAGS nodes.
  posterior <- cbind(
    mu_intercept = c(1, 2),
    mu_x = c(0.5, -0.25),
    `mu_f[1]` = 0.3,
    `mu_f[2]` = 0.3,
    `mu_x__xXx__f[1]` = c(0.2, 0.4),
    `mu_x__xXx__f[2]` = c(-0.6, 0.1)
  )
  scale_info <- formula_result$formula_scale$mu_x
  expected_f <- 0.3 - scale_info$mean / scale_info$sd *
    posterior[, c("mu_x__xXx__f[1]", "mu_x__xXx__f[2]")]

  transformed <- BayesTools:::.bt_transform_scale_posterior(
    posterior,
    formula_scale = list(mu = formula_result$formula_scale)
  )
  expect_identical(colnames(transformed), colnames(posterior))
  expect_equal(
    unname(transformed[, c("mu_f[1]", "mu_f[2]")]),
    unname(expected_f),
    tolerance = 1e-12
  )

  # Fits created before point terms were stored per coefficient carry one
  # unindexed factor point term.
  legacy_scale <- formula_result$formula_scale
  attr(legacy_scale, "point_terms") <- c(mu_f = 0.3)
  expect_identical(
    colnames(BayesTools:::.bt_transform_scale_posterior(posterior, formula_scale = list(mu = legacy_scale))),
    colnames(posterior)
  )

  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  table <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    remove_spike_0 = FALSE,
    remove_diagnostics = TRUE
  )
  expect_false("(mu) f" %in% rownames(table))
  expect_equal(
    unname(table[c("(mu) f[b]", "(mu) f[c]"), "Mean"]),
    unname(colMeans(expected_f)),
    tolerance = 1e-12
  )

  prior_samples <- transform_prior_samples(fit, n_samples = 16, seed = 3)
  expect_false("mu_f" %in% colnames(prior_samples))
  expect_true(all(c("mu_f[1]", "mu_f[2]") %in% colnames(prior_samples)))
})

test_that("transform_scale_samples keeps indexed random-effect correlations in one block", {

  df <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    f = factor(c("a", "b", "a", "b", "a", "b")),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(0 + f + x:f | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_cols <- random_term$sd_parameter_names
  expect_equal(
    unname(random_term$sd_leaves$leaf_terms[sd_cols]),
    c("f[b]", "f[a]__xXx__x", "f[b]__xXx__x")
  )

  source_sd <- c(1.1, 0.8, 1.2)
  source_cor <- matrix(0.25, nrow = 3, ncol = 3)
  diag(source_cor) <- 1
  cor_names <- outer(
    seq_along(sd_cols),
    seq_along(sd_cols),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  posterior <- matrix(
    c(source_sd, as.vector(source_cor)),
    nrow = 1,
    dimnames = list(NULL, c(sd_cols, as.vector(cor_names)))
  )

  transformed <- BayesTools:::.bt_transform_scale_posterior(
    posterior,
    list(mu = formula_result$formula_scale)
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- diag(3)
  M[1, 3] <- -scale_info$mean / scale_info$sd
  M[2, 2] <- 1 / scale_info$sd
  M[3, 3] <- 1 / scale_info$sd
  expected_cov <- M %*% diag(source_sd, nrow = 3) %*% source_cor %*%
    diag(source_sd, nrow = 3) %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov / tcrossprod(expected_sd)

  expect_equal(unname(transformed[1, sd_cols]), expected_sd, tolerance = 1e-12)
  expect_equal(
    unname(transformed[1, as.vector(cor_names)]),
    as.vector(expected_cor),
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples unscales random-factor correlations in column space", {

  df <- data.frame(
    x = 1:9,
    f = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 3))
  )
  formula_result <- NULL
  expect_no_warning(
    formula_result <- JAGS_formula(
      formula = ~ 1 + us(1 + x + f | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      formula_scale = list(x = TRUE),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          terms = list(
            f = prior("gamma", list(2, 2))
          ),
          contrasts = c(f = "orthonormal"),
          cor = prior_lkj(eta = 1)
        )
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  # the orthonormal factor term is one exchangeable unit: one SD, and an
  # identity correlation block of its two columns
  expect_equal(
    random_term$sd_parameter_names[3:4],
    rep("mu__xREx__id_f", 2L)
  )
  sd_cols <- unique(random_term$sd_parameter_names)
  R_names <- outer(
    seq_len(random_term$n_columns),
    seq_len(random_term$n_columns),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- outer(
    seq_len(random_term$n_columns),
    seq_len(random_term$n_columns),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_L[", row, ",", column, "]")
    })
  )
  source_sd <- c(1.5, 2, 3, 3)
  source_cor <- matrix(
    c(
      1, 0.2, 0.3, 0.4,
      0.2, 1, 0.5, 0.6,
      0.3, 0.5, 1, 0,
      0.4, 0.6, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  source_L <- t(chol(source_cor))
  posterior <- matrix(
    c(source_sd[1:3], as.vector(source_cor), as.vector(source_L)),
    nrow = 1,
    dimnames = list(NULL, c(sd_cols, as.vector(R_names), as.vector(L_names)))
  )

  transformed <- BayesTools:::.bt_transform_scale_posterior(
    posterior,
    list(mu = formula_result$formula_scale)
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- diag(random_term$n_columns)
  M[1, 2] <- -scale_info$mean / scale_info$sd
  M[2, 2] <- 1 / scale_info$sd
  expected_cov <- M %*% diag(source_sd, nrow = 4) %*% source_cor %*%
    diag(source_sd, nrow = 4) %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov / tcrossprod(expected_sd)
  expected_L <- t(chol(expected_cor))

  # the unit's columns are not rescaled: its SD and identity block are kept
  expect_identical(expected_sd[3], expected_sd[4])
  expect_identical(expected_cor[3, 4], 0)
  expect_equal(
    unname(transformed[1, sd_cols]),
    expected_sd[1:3],
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, as.vector(R_names)]),
    as.vector(expected_cor),
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, as.vector(L_names)]),
    as.vector(expected_L),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(source_cor[1, 4], expected_cor[1, 4])))
})

test_that("JAGS_formula refuses homogeneous random-effect SD blocks with a scaled slope", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  compile <- function(formula, random = prior_random(id = random_block(sd = prior("gamma", list(2, 2)))),
                      formula_scale = list(x = TRUE),
                      priors = list(intercept = prior("normal", list(0, 1)))){
    JAGS_formula(
      formula = formula, parameter = "mu", data = df,
      prior_list = priors,
      formula_scale = formula_scale, prior_random = random
    )
  }
  fixed_slope_priors <- list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))
  refusal <- "Cannot unscale homogeneous random-effect SD 'mu__xREx__id_sd'"

  # one shared SD for several columns, one of them a scaled slope
  for(formula in list(
    ~ 1 + id(1 + x | id),
    ~ 1 + random(1 + x | id, covariance = "id"),
    ~ 1 + diag(1 + x | id, hom = TRUE)
  )){
    expect_error(compile(formula), refusal, fixed = TRUE)
    # F314 (2), F318 (1): a classed refusal whose remedy fits JAGS_formula()
    condition <- tryCatch(compile(formula), error = identity)
    expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
    expect_identical(condition$reason, "homogeneous_scaled_block")
    expect_identical(condition$parameter, "mu")
    expect_identical(condition$coordinates, "mu__xREx__id_sd")
    expect_match(
      conditionMessage(condition),
      paste0("Use a heterogeneous random-effect SD structure (for example diag() or us()), ",
             "or do not standardize the slope with 'formula_scale'."),
      fixed = TRUE
    )
    expect_no_match(conditionMessage(condition), "fitted scale", fixed = TRUE)
  }
  # nothing is refused when x is not scaled, for a single column, for
  # heterogeneous blocks, and for an SD that comes from an external source
  expect_s3_class(compile(~ 1 + id(1 + x | id), formula_scale = NULL)$formula_design, "BayesTools_formula_design")
  expect_s3_class(compile(~ 1 + id(0 + x | id))$formula_design, "BayesTools_formula_design")
  for(formula in list(~ 1 + x + id(1 | id), ~ 1 + x + diag(1 + x | id))){
    expect_s3_class(compile(formula, priors = fixed_slope_priors)$formula_design, "BayesTools_formula_design")
  }
  expect_s3_class(compile(~ 1 + id(1 + x | id),
    prior_random(id = random_block(sd_source = random_sd_source("tau"))))$formula_design, "BayesTools_formula_design")
})

test_that("transform_scale_samples guards homogeneous random-effect SD scaling", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  slope_result <- JAGS_formula(
    formula = ~ 1 + id(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  slope_posterior <- matrix(
    2,
    nrow = 1,
    dimnames = list(NULL, "mu__xREx__id_sd")
  )
  slope_transformed <- BayesTools:::.bt_transform_scale_posterior(
    slope_posterior,
    list(mu = slope_result$formula_scale)
  )
  expect_equal(
    unname(slope_transformed[1, "mu__xREx__id_sd"]),
    2 / slope_result$formula_scale$mu_x$sd,
    tolerance = 1e-12
  )

  # JAGS_formula() no longer builds an `id(1 + x | id)` block with a scaled x,
  # so the transform-time guard is checked on hand-built leaves: the leaves of
  # the unscaled block beside the scale of the scaled single-column block
  intercept_slope_result <- JAGS_formula(
    formula = ~ 1 + id(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  hand_built_scale <- slope_result$formula_scale
  attr(hand_built_scale, "random_effect_sd_leaves") <- list(
    "__xREx__id" = intercept_slope_result$formula_design$random_effects[[1L]]$sd_leaves
  )
  expect_error(
    BayesTools:::.bt_transform_scale_posterior(
      slope_posterior,
      list(mu = hand_built_scale)
    ),
    "Cannot unscale homogeneous random-effect SD",
    fixed = TRUE
  )
  condition <- tryCatch(
    BayesTools:::.bt_transform_scale_posterior(
      slope_posterior,
      list(mu = hand_built_scale)
    ),
    error = identity
  )
  expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
  expect_identical(condition$reason, "homogeneous_scaled_block")
})
