skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Random-effect bridge reconstruction
# ============================================================================ #
#
# PURPOSE:
#   Reconstruction of formula and random-effect parameters for bridge sampling.
# ============================================================================ #

test_that("random-effect formulas are guarded in fixed-only downstream evaluators", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"))
  )
  fit <- coda::mcmc(matrix(
    c(0, 0),
    ncol = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  fixed_prior <- prior("normal", list(0, 1))
  attr(fixed_prior, "parameter") <- "mu"

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + (1 || id),
      parameter = "mu",
      data = df,
      prior_list = list(mu_intercept = fixed_prior)
    ),
    "needs the fitted formula design of parameter 'mu'",
    fixed = TRUE
  )

  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  formula_result$formula_design$formula_scale <- .bt_formula_scale_finalize(
    formula_result$formula_scale, formula_result$formula_design, formula_result$prior_list)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional"
    ),
    "cannot be reconstructed",
    fixed = TRUE
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1,
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "includes random effects, so JAGS_evaluate_formula() needs an explicit 'formula_target'",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(),
      formula_list = list(mu = ~ 1 + diag(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      ))
    ),
    "supplied formula-related inputs do not fully match",
    fixed = TRUE
  )

  bridge_incomplete_fit <- coda::mcmc(matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu__xREx__id_sd"))
  ))
  expect_error(
    JAGS_bridgesampling(
      fit = bridge_incomplete_fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(),
      formula_list = list(mu = ~ 1 + id(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      ))
    ),
    "supplied formula inputs cannot replace fitted replay metadata",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = NULL
    ),
    "requires posterior samples of standardized latent random effects",
    fixed = TRUE,
    class = "BayesTools_refit_monitoring"
  )
})

test_that("JAGS bridgesampling can reconstruct formula parameters from fitted design metadata", {

  df <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  samples <- c(mu_intercept = 10, mu_x = 2)

  reconstructed <- JAGS_marglik_parameters_formula(
    samples = samples,
    formula_list = NULL,
    formula_data_list = NULL,
    formula_prior_list = list(mu = formula_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = formula_result$formula_design)
  )
  expected_x <- as.vector(formula_result$formula_design$model_matrix[, "x"])
  expect_equal(reconstructed$mu, 10 + 2 * expected_x, tolerance = 1e-12)

  log_formula <- ~ 1 + x
  attr(log_formula, "log(intercept)") <- TRUE
  log_result <- JAGS_formula(
    formula = log_formula,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("gamma", list(2, 2)),
      x = prior("normal", list(0, 1))
    )
  )
  log_reconstructed <- JAGS_marglik_parameters_formula(
    samples = c(mu_intercept = exp(1), mu_x = 2),
    formula_list = NULL,
    formula_data_list = NULL,
    formula_prior_list = list(mu = log_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = log_result$formula_design)
  )
  expect_equal(log_reconstructed$mu, 1 + 2 * df$x, tolerance = 1e-12)
})

test_that("JAGS bridgesampling uses fitted formula metadata and errors on supplied mismatches", {

  df <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "prior_list") <- formula_result$prior_list
  formula_result$formula_design$formula_scale <- .bt_formula_scale_finalize(
    formula_result$formula_scale, formula_result$formula_design, formula_result$prior_list)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list()
    ),
    "posterior' does not contain",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      formula_list = list(mu = ~ 1 + x),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = prior_list),
      formula_scale_list = list(wrong_name = list(x = TRUE))
    ),
    "supplied formula-related inputs do not fully match",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      formula_list = list(mu = ~ 1 + x),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = prior_list),
      formula_scale_list = TRUE
    ),
    "supplied formula-related inputs do not fully match",
    fixed = TRUE
  )

  # the fitted formula priors and re-specified equal priors (without the
  # bookkeeping attributes JAGS_formula() adds) are dropped silently
  respecified_prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x         = prior("normal", list(0, 1))
  )
  for(supplied in list(formula_result$prior_list, respecified_prior_list)){
    expect_no_warning(
      expect_error(
        JAGS_bridgesampling(
          fit = fit,
          log_posterior = function(parameters, data) 0,
          data = list(),
          prior_list = supplied
        ),
        "posterior' does not contain",
        fixed = TRUE
      )
    )
  }

  # a different prior for a formula parameter is an error, not a silent override
  different_prior_list <- formula_result$prior_list
  different_prior_list[["mu_intercept"]] <- prior("normal", list(0, 5))
  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = different_prior_list
    ),
    paste0(
      "JAGS_bridgesampling() received formula priors in 'prior_list' that ",
      "differ from the fitted formula priors: mu_intercept. Formula priors are ",
      "taken from the fitted formula design; remove them from 'prior_list' or ",
      "refit the model with the intended priors."
    ),
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
})

test_that("formula prior duplicates are compared by their semantic fields", {

  df <- data.frame(
    x = c(1, 3, 7, 2, 5, 4),
    f = factor(c("a", "b", "c", "a", "b", "c"))
  )
  slope_prior <- prior("normal", list(0, 1))
  attr(slope_prior, "multiply_by") <- "sigma"
  user_priors <- list(
    intercept = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    ),
    x = slope_prior,
    f = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  formula_result <- JAGS_formula(~ x + f, "mu", df, user_priors)
  design_list <- list(mu = formula_result$formula_design)
  compare <- function(prior_list){
    BayesTools:::.bt_JAGS_bridge_non_formula_prior_list(
      prior_list          = prior_list,
      formula_design_list = design_list,
      check_duplicates    = TRUE
    )
  }
  respecified <- stats::setNames(user_priors, paste0("mu_", names(user_priors)))

  # bookkeeping attributes of the fitted priors (formula parameter, levels,
  # factor design, ...) and the model prior weight are ignored
  expect_length(compare(formula_result$prior_list), 0L)
  expect_length(compare(respecified), 0L)
  reweighted <- respecified
  reweighted$mu_x$prior_weights <- 3
  expect_length(compare(c(reweighted, list(sigma = prior("gamma", list(2, 2))))), 1L)

  # the distribution, parameters, truncation, contrast, multiply_by, and
  # mixture components are compared
  changed <- list(
    mu_x = prior("normal", list(0, 1)),
    mu_x = {p <- prior("normal", list(0, 1), truncation = list(0, Inf)); attr(p, "multiply_by") <- "sigma"; p},
    mu_x = {p <- prior("normal", list(0, 2)); attr(p, "multiply_by") <- "sigma"; p},
    mu_x = {p <- prior("normal", list(0, 1)); attr(p, "multiply_by") <- "tau"; p},
    mu_f = prior_factor("normal", list(0, 1), contrast = "independent"),
    mu_f = prior_factor("cauchy", list(0, 1), contrast = "treatment"),
    mu_intercept = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(2, 1))
    ),
    mu_intercept = prior("normal", list(0, 1))
  )
  for(i in seq_along(changed)){
    expect_error(
      compare(changed[i]),
      paste0("differ from the fitted formula priors: ", names(changed)[i], "."),
      fixed = TRUE
    )
  }

  # mixture priors (components, mixing weights, null/alternative labels) and
  # ordered priors (total, allocation, contrast, id)
  mixture <- function(weight = 2, sd = 3, is_null = c(FALSE, FALSE, TRUE)){
    prior_mixture(
      list(
        prior("normal", list(0, 1), prior_weights = weight),
        prior("normal", list(0, sd)),
        prior("point", list(0))
      ),
      is_null = is_null
    )
  }
  ordered <- function(sd = 1, alpha = c(2, 2), contrast = "cumulative", id = NULL){
    prior_ordered(
      prior("normal", list(0, sd), truncation = list(0, Inf)),
      allocation = prior("dirichlet", list(alpha = alpha)),
      contrast   = contrast,
      id         = id
    )
  }
  df$o <- factor(c("l", "m", "h", "l", "m", "h"), levels = c("l", "m", "h"))
  formula_result <- JAGS_formula(~ x + o, "mu", df, list(intercept = mixture(), x = slope_prior, o = ordered()))
  design_list <- list(mu = formula_result$formula_design)
  expect_length(compare(formula_result$prior_list), 0L)
  expect_length(compare(list(mu_intercept = mixture(), mu_o = ordered())), 0L)
  changed <- list(
    mu_intercept = mixture(weight = 1),
    mu_intercept = mixture(sd = 2),
    mu_intercept = mixture(is_null = c(FALSE, TRUE, TRUE)),
    mu_o         = ordered(sd = 2),
    mu_o         = ordered(alpha = c(1, 1)),
    mu_o         = ordered(contrast = "cumulative_levels"),
    mu_o         = ordered(id = "shared")
  )
  for(i in seq_along(changed)){
    expect_error(
      compare(changed[i]),
      paste0("differ from the fitted formula priors: ", names(changed)[i], "."),
      fixed = TRUE
    )
  }
})

test_that("formula random effects require prior_random with no legacy fallback", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"))
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "intercept|id" = prior("gamma", list(2, 2))
      )
    ),
    "Formula random effects require 'prior_random'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = ~ 1 + diag(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1))))
    ),
    "JAGS_fit() requires 'formula_random_prior_list'",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "intercept|id" = prior("gamma", list(2, 2))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      )
    ),
    "no longer supported",
    fixed = TRUE
  )
})

test_that("formula random effects expose bridge-ready stochastic coordinates", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
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
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = formula_result$formula_design)
  )
  expected_z <- c(
    "mu__xREx__id_xRE_Zx[1,1]",
    "mu__xREx__id_xRE_Zx[2,1]",
    "mu__xREx__id_xRE_Zx[1,2]",
    "mu__xREx__id_xRE_Zx[2,2]"
  )
  expect_equal(bridge$parameters, expected_z)
  expect_equal(bridge$bounds$lb, stats::setNames(rep(-Inf, 4), expected_z))
  expect_equal(bridge$bounds$ub, stats::setNames(rep( Inf, 4), expected_z))
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = NULL,
      add_bounds = list(lb = -Inf, ub = Inf),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(lb = -Inf, ub = Inf),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "names must be unique and match",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(
        lb = stats::setNames(0, expected_z[1]),
        ub = stats::setNames(Inf, expected_z[1])
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "conflict",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(
        lb = stats::setNames(-Inf, expected_z[1]),
        ub = stats::setNames(1 - 1e-12, expected_z[1])
      ),
      bridge_parameters = expected_z[1],
      bridge_bounds = list(
        lb = stats::setNames(-Inf, expected_z[1]),
        ub = stats::setNames(1, expected_z[1])
      )
    ),
    "conflict",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(
        lb = stats::setNames(-Inf, "wrong_name"),
        ub = stats::setNames( Inf, expected_z[1])
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "names must match",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = c(expected_z[1], expected_z[2]),
      add_bounds = list(
        lb = stats::setNames(c(-Inf, -Inf), c(expected_z[1], expected_z[1])),
        ub = stats::setNames(c( Inf,  Inf), c(expected_z[1], expected_z[2]))
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "unique",
    fixed = TRUE
  )
  merged_bridge <- BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
    add_parameters = expected_z[1],
    add_bounds = list(
      lb = stats::setNames(-Inf, expected_z[1]),
      ub = stats::setNames( Inf, expected_z[1])
    ),
    bridge_parameters = bridge$parameters,
    bridge_bounds = bridge$bounds
  )
  expect_equal(merged_bridge$add_parameters, expected_z)

  samples <- c(
    mu_intercept = 10,
    mu_x = 1,
    mu__xREx__id_intercept = 0.5,
    mu__xREx__id_x = 0.25,
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4
  )
  parameters <- JAGS_marglik_parameters_formula(
    samples = samples,
    formula_list = list(mu = formula_result$formula),
    formula_data_list = list(mu = formula_result$data),
    formula_prior_list = list(mu = formula_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = formula_result$formula_design)
  )

  expect_equal(parameters$mu, c(10.5, 12, 15.5, 17.5), tolerance = 1e-12)
  expect_equal(
    bayestools_reference_formula_random_log_prior(
      samples,
      list(mu = formula_result$formula_design)
    ),
    sum(stats::dnorm(samples[expected_z], log = TRUE)),
    tolerance = 1e-12
  )
})

test_that("formula random-effect bridge helpers handle LKJ and scalar rho blocks", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )

  lkj_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1.5)
      )
    )
  )
  lkj_u <- "mu__xREx__id_xRE_CORx_lkj_u[1]"
  expect_true(lkj_u %in% lkj_result$add_parameters)
  lkj_bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = lkj_result$formula_design)
  )
  expect_true(lkj_u %in% lkj_bridge$parameters)
  expect_equal(lkj_bridge$bounds$lb[[lkj_u]], 0)
  expect_equal(lkj_bridge$bounds$ub[[lkj_u]], 1)

  bridge_fitted <- list(mu = lkj_result$formula_design)
  bridge_changed <- bridge_fitted
  bridge_changed$mu$random_effects[[1]]$homogeneous_sd <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = bridge_fitted[["mu"]],
      rebuilt = bridge_changed[["mu"]]
    ),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )

  bridge_changed <- bridge_fitted
  bridge_changed$mu$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = bridge_fitted[["mu"]],
      rebuilt = bridge_changed[["mu"]]
    ),
    "missing canonical 'random_term\\$correlation'"
  )

  malformed_lkj <- lkj_result$formula_design
  malformed_lkj$random_effects[[1]]$structure <- NULL
  malformed_lkj$random_effects[[1]]$covariance <- "us"
  expect_error(
    BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
      list(mu = malformed_lkj)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  lkj_samples <- c(
    mu_intercept = 0,
    mu_x = 0,
    mu__xREx__id_intercept = 0.5,
    mu__xREx__id_x = 0.25,
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4,
    "mu__xREx__id_xRE_CORx_lkj_u[1]" = 0.75
  )
  lkj_parameters <- JAGS_marglik_parameters_formula(
    samples = lkj_samples,
    formula_list = list(mu = lkj_result$formula),
    formula_data_list = list(mu = lkj_result$data),
    formula_prior_list = list(mu = lkj_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = lkj_result$formula_design)
  )
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(0.75, K = 2)
  z <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
  coef <- sweep(z %*% t(L), 2, c(0.5, 0.25), "*")
  expected_lkj <- rowSums(coef[lkj_result$formula_design$random_effects[[1]]$group_map, ] *
    lkj_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(lkj_parameters$mu, unname(expected_lkj), tolerance = 1e-12)
  expect_equal(
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = lkj_result$formula_design)
    ),
    sum(stats::dnorm(lkj_samples[grep("xRE_Zx", names(lkj_samples), value = TRUE)], log = TRUE)) +
      BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(0.75, K = 2, eta = 1.5),
    tolerance = 1e-12
  )
  expect_error(
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = malformed_lkj)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  missing_correlation_lkj <- lkj_result$formula_design
  missing_correlation_lkj$random_effects[[1]]$correlation <- NULL
  expect_error(
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = missing_correlation_lkj)
    ),
    "missing canonical 'random_term\\$correlation'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples,
      missing_correlation_lkj$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation'"
  )
  expect_error(
    .parameter_catalog_random_summary_samples(
      model_samples = matrix(
        unname(lkj_samples),
        nrow = 1,
        dimnames = list(NULL, names(lkj_samples))
      ),
      prior_list = lkj_result$prior_list,
      formula_design = list(mu = missing_correlation_lkj),
      mode = "standard"
    ),
    "missing canonical 'random_term\\$correlation'"
  )

  missing_eta_lkj <- lkj_result$formula_design
  missing_eta_lkj$random_effects[[1]]$correlation$eta <- NULL
  expect_error(
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = missing_eta_lkj)
    ),
    "missing canonical 'random_term\\$correlation\\$eta'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples,
      malformed_lkj$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$structure'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples[names(lkj_samples) != lkj_u],
      lkj_result$formula_design$random_effects[[1]]
    ),
    "missing random-effect correlation coordinates",
    fixed = TRUE
  )

  ar_result <- JAGS_formula(
    formula = ~ 1 + x + ar1(idx | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_equal(ar_result$formula_design$random_effects[[1]]$correlation$sample_name, "mu__xREx__id_rho_z")

  ar_rho_bridge <- BayesTools:::.bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    list(mu = ar_result$formula_design)
  )
  expect_equal(ar_rho_bridge$parameters, "mu__xREx__id_rho_z")
  expect_equal(ar_rho_bridge$bounds$lb[["mu__xREx__id_rho_z"]], -Inf)
  expect_equal(ar_rho_bridge$bounds$ub[["mu__xREx__id_rho_z"]],  Inf)

  ar_samples <- c(
    mu_intercept = 0,
    mu_x = 0,
    mu__xREx__id_sd = 0.5,
    mu__xREx__id_rho_z = atanh(0.25),
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4
  )
  ar_parameters <- JAGS_marglik_parameters_formula(
    samples = ar_samples,
    formula_list = list(mu = ar_result$formula),
    formula_data_list = list(mu = ar_result$data),
    formula_prior_list = list(mu = ar_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = ar_result$formula_design)
  )
  ar_L <- t(chol(matrix(c(1, 0.25, 0.25, 1), nrow = 2)))
  ar_coef <- 0.5 * (z %*% t(ar_L))
  expected_ar <- rowSums(ar_coef[ar_result$formula_design$random_effects[[1]]$group_map, ] *
    ar_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(ar_parameters$mu, unname(expected_ar), tolerance = 1e-12)
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      c(ar_samples, mu__xREx__id_rho = 0.99),
      ar_result$formula_design$random_effects[[1]]
    ),
    0.25,
    tolerance = 1e-12
  )

  ar_saturated_samples <- ar_samples
  ar_saturated_samples[["mu__xREx__id_rho_z"]] <- 20
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      ar_saturated_samples,
      ar_result$formula_design$random_effects[[1]]
    ),
    -Inf
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_saturated_samples,
      ar_result$formula_design$random_effects[[1]]
    ),
    "missing or invalid scalar correlation coordinates",
    fixed = TRUE
  )

  missing_rho_scale <- ar_result$formula_design
  missing_rho_scale$random_effects[[1]]$correlation$rho_scale <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_samples,
      missing_rho_scale$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation\\$rho_scale'"
  )

  missing_rho_bounds <- ar_result$formula_design
  missing_rho_bounds$random_effects[[1]]$correlation$bounds <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_samples,
      missing_rho_bounds$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation\\$bounds'"
  )

  car_df <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  car_result <- JAGS_formula(
    formula = ~ 1 + car(time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_equal(car_result$formula_design$random_effects[[1]]$correlation$sample_name, "mu__xREx__id_rho_z")
  expect_equal(car_result$formula_design$random_effects[[1]]$correlation$time_values, c(0, 0.5, 2))

  car_samples <- c(
    mu_intercept = 0,
    mu__xREx__id_sd = 0.5,
    mu__xREx__id_rho_z = atanh(0.25),
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[1,3]" = 3,
    "mu__xREx__id_xRE_Zx[2,1]" = 4,
    "mu__xREx__id_xRE_Zx[2,2]" = 5,
    "mu__xREx__id_xRE_Zx[2,3]" = 6
  )
  car_parameters <- JAGS_marglik_parameters_formula(
    samples = car_samples,
    formula_list = list(mu = car_result$formula),
    formula_data_list = list(mu = car_result$data),
    formula_prior_list = list(mu = car_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = car_result$formula_design)
  )
  car_R <- 0.25^abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
  car_L <- t(chol(car_R))
  car_z <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)
  car_coef <- 0.5 * (car_z %*% t(car_L))
  expected_car <- rowSums(car_coef[car_result$formula_design$random_effects[[1]]$group_map, ] *
    car_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(car_parameters$mu, unname(expected_car), tolerance = 1e-12)

  car_logit <- car_result
  car_logit$formula_design$random_effects[[1]]$correlation$rho_scale <- "logit"
  car_logit$formula_design$random_effects[[1]]$correlation$sample_name <- "mu__xREx__id_rho_logit"
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      c(car_samples[names(car_samples) != "mu__xREx__id_rho_z"], mu__xREx__id_rho_logit = stats::qlogis(0.25)),
      car_logit$formula_design$random_effects[[1]]
    ),
    0.25,
    tolerance = 1e-12
  )
  car_logit_saturated_samples <- c(
    car_samples[names(car_samples) != "mu__xREx__id_rho_z"],
    mu__xREx__id_rho_logit = 1000
  )
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      car_logit_saturated_samples,
      car_logit$formula_design$random_effects[[1]]
    ),
    -Inf
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      car_logit_saturated_samples,
      car_logit$formula_design$random_effects[[1]]
    ),
    "missing or invalid scalar correlation coordinates",
    fixed = TRUE
  )

  hcs_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2"), each = 3), levels = c("g1", "g2"))
  )
  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = hcs_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  hcs_random_term <- hcs_result$formula_design$random_effects[[1]]
  hcs_rho_bridge <- BayesTools:::.bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    list(mu = hcs_result$formula_design)
  )
  expect_equal(hcs_rho_bridge$parameters, "mu__xREx__id_rho_z")
  expect_equal(
    hcs_rho_bridge$bounds$lb[["mu__xREx__id_rho_z"]],
    atanh(-1 / 2),
    tolerance = 1e-12
  )
  expect_equal(hcs_rho_bridge$bounds$ub[["mu__xREx__id_rho_z"]], Inf)

  hcs_z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = hcs_random_term,
    n_groups = hcs_random_term$n_groups,
    n_columns = hcs_random_term$n_columns
  ))
  hcs_boundary_samples <- c(
    stats::setNames(rep(0, length(hcs_z_names)), hcs_z_names),
    mu__xREx__id_rho_z = atanh(-1 / 2)
  )
  expect_equal(
    bayestools_reference_formula_random_log_prior(
      hcs_boundary_samples,
      list(mu = hcs_result$formula_design)
    ),
    -Inf
  )
  hcs_interior_samples <- hcs_boundary_samples
  hcs_interior_samples[["mu__xREx__id_rho_z"]] <- 0
  expect_true(is.finite(
    bayestools_reference_formula_random_log_prior(
      hcs_interior_samples,
      list(mu = hcs_result$formula_design)
    )
  ))

  hcs_fixed_rho <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = hcs_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("point", list(location = 0.2))
      )
    )
  )
  hcs_fixed_term <- hcs_fixed_rho$formula_design$random_effects[[1]]
  expect_null(BayesTools:::.bt_JAGS_bridge_scalar_rho_parameter(hcs_fixed_term))
  expect_equal(hcs_fixed_term$correlation$sample_fixed, 0.2)
  hcs_fixed_z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = hcs_fixed_term,
    n_groups = hcs_fixed_term$n_groups,
    n_columns = hcs_fixed_term$n_columns
  ))
  hcs_fixed_samples <- stats::setNames(rep(0, length(hcs_fixed_z_names)), hcs_fixed_z_names)
  expect_true(is.finite(
    bayestools_reference_formula_random_log_prior(
      hcs_fixed_samples,
      list(mu = hcs_fixed_rho$formula_design)
    )
  ))
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      hcs_fixed_samples,
      hcs_fixed_term
    ),
    tanh(0.2),
    tolerance = 1e-12
  )
})
