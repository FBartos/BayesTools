skip_if_not_test_profile("unit")

.formula_prediction_data <- function(){
  data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
}

.formula_prediction_sd_prior <- function(){
  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.formula_prediction_result <- function(random_effects_compile = NULL,
                                       new_levels = NULL){
  JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = .formula_prediction_data(),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = .formula_prediction_sd_prior(),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = TRUE,
          correlation = FALSE
        ),
        new_levels = new_levels
      )
    ),
    random_effects_compile = random_effects_compile
  )
}

.formula_prediction_fit <- function(result){
  posterior <- matrix(
    c(
      10, 1, 2, 3,  0.5, 0.1, -0.5, -0.2,
      20, 2, 2, 3,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit
}

test_that("conditional target equals fixed target for fixed-only formulas", {

  df <- .formula_prediction_data()
  result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  posterior <- matrix(
    c(10, 1, 20, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "fixed"
  )
  conditional <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "conditional"
  )

  expect_equal(conditional, fixed)
})

test_that("formula prediction validates parameter names and seeds", {
  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)

  expect_error(
    JAGS_evaluate_formula(fit, parameter = NA_character_),
    "'parameter' argument cannot contain NA/NaN values",
    fixed = TRUE
  )
  expect_error(
    JAGS_predict_formula(fit, parameter = "mu", seed = .Machine$integer.max + 1),
    "'seed' must be equal or lower than",
    fixed = TRUE
  )
})

test_that("formula design metadata preserves scaling during prediction", {

  df <- data.frame(
    x = c(10, 20, 30),
    id = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  fixed_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE)
  )
  fixed_posterior <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  )
  fixed_fit <- coda::mcmc(fixed_posterior)
  attr(fixed_fit, "formula_design") <- list(mu = fixed_result$formula_design)

  scaled_x <- as.numeric(scale(df$x))
  fixed_prediction <- JAGS_evaluate_formula(
    fit = fixed_fit,
    parameter = "mu",
    formula_target = "fixed"
  )
  expect_equal(unname(drop(fixed_prediction)), 1 + 2 * scaled_x)

  random_result <- JAGS_formula(
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
        sd = .formula_prediction_sd_prior(),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- random_result$formula_design$random_effects[[1L]]
  coefficient_names <- as.vector(BayesTools:::.bt_random_effect_coefficient_names(
    random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  random_posterior <- matrix(
    c(1, 2, 0.5, -0.25),
    nrow = 1,
    dimnames = list(
      NULL,
      c("mu_intercept", "mu_x", coefficient_names)
    )
  )
  random_fit <- coda::mcmc(random_posterior)
  attr(random_fit, "formula_design") <- list(mu = random_result$formula_design)

  conditional_prediction <- JAGS_evaluate_formula(
    fit = random_fit,
    parameter = "mu",
    formula_target = "conditional"
  )
  expect_equal(
    unname(drop(conditional_prediction)),
    1 + 2 * scaled_x + c(0.5, -0.25, 0.5) * scaled_x
  )
})

test_that("formula_target fixed and conditional preserve explicit semantics", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  df <- .formula_prediction_data()

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = result$prior_list
    ),
    "silently dropping group-level contributions",
    fixed = TRUE
  )

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "fixed"
  )
  expect_equal(
    unname(fixed),
    cbind(10 + df$x, 20 + 2 * df$x)
  )

  conditional <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "conditional"
  )
  expect_equal(
    unname(conditional),
    cbind(
      c(9.4, 10.5, 10.3, 11.1),
      c(18.4, 20.6, 21.1, 22.8)
    )
  )

  default_prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu"
  )
  expect_equal(default_prediction$value, conditional)
  expect_equal(default_prediction$mean, fixed)
  expect_equal(default_prediction$random, conditional - fixed)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      new_levels = "zero"
    ),
    "'new_levels' can be used only",
    fixed = TRUE
  )
})

test_that("conditional target handles new levels by explicit policy", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "conditional"
    ),
    "New random-effect level",
    fixed = TRUE
  )

  zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional",
    new_levels = "zero"
  )
  expect_equal(
    unname(zero),
    cbind(c(10.5, 11, 12), c(20.6, 22, 24))
  )

  repeated_new <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("c", "c", "c"), levels = c("a", "b", "c"))
  )
  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = repeated_new,
    formula_target = "fixed"
  )
  set.seed(123)
  sampled <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = repeated_new,
    formula_target = "conditional",
    new_levels = "sample"
  )
  random <- sampled - fixed
  expect_equal(
    unname(random[3, ] - random[2, ]),
    unname(random[2, ] - random[1, ]),
    tolerance = 1e-12
  )
})

test_that("structured new-level sampling uses only requested column subsets", {

  factor_levels <- sprintf("level_%03d", seq_len(113L))
  group_levels  <- sprintf("group_%02d", seq_len(17L))
  df <- data.frame(
    f = factor(factor_levels, levels = factor_levels),
    id = factor(rep(group_levels, length.out = length(factor_levels)),
                levels = group_levels)
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  result <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  sd_name     <- unique(random_term$sd_parameter_names)
  rho_name    <- random_term$correlation$rho_name
  posterior <- cbind(
    mu_intercept = c(10, 20),
    sd = c(2, 3),
    rho = c(-0.005, 0.4)
  )
  colnames(posterior)[2:3] <- c(sd_name, rho_name)
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  new_data <- data.frame(
    f = factor(
      factor_levels[c(2L, 113L, 2L, 50L, 113L)],
      levels = factor_levels
    ),
    id = factor(
      c("new_1", "new_1", "new_1", "new_2", "new_2"),
      levels = c(group_levels, "new_1", "new_2")
    )
  )

  subset_cholesky <- BayesTools:::.bt_random_effect_structured_subset_cholesky
  prediction_subset_transform <-
    BayesTools:::.bt_random_effect_prediction_structured_subset_transform
  requested_columns <- list()
  testthat::local_mocked_bindings(
    .bt_random_effect_marginal_covariance_correlation_draws = function(...){
      stop("global correlation materialized", call. = FALSE)
    },
    .bt_random_effect_prediction_structured_subset_transform = function(
        structure, columns, latent, rho, coordinates){

      requested_columns[[length(requested_columns) + 1L]] <<- columns
      prediction_subset_transform(
        structure = structure,
        columns = columns,
        latent = latent,
        rho = rho,
        coordinates = coordinates
      )
    },
    .package = "BayesTools"
  )

  set.seed(812)
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    prior_list = result$prior_list,
    formula_target = "conditional",
    new_levels = "sample"
  )

  set.seed(812)
  expected_random <- matrix(0, nrow = nrow(new_data), ncol = nrow(posterior))
  for(draw in seq_len(nrow(posterior))){
    L_1 <- subset_cholesky(
      structure = "cs",
      columns = c(2L, 113L),
      rho = posterior[draw, rho_name],
      global_n_columns = length(factor_levels),
      column_coordinates = seq_along(factor_levels)
    )
    effect_1 <- drop(L_1 %*% stats::rnorm(2L)) * posterior[draw, sd_name]
    expected_random[1:3, draw] <- effect_1[c(1L, 2L, 1L)]

    L_2 <- subset_cholesky(
      structure = "cs",
      columns = c(50L, 113L),
      rho = posterior[draw, rho_name],
      global_n_columns = length(factor_levels),
      column_coordinates = seq_along(factor_levels)
    )
    effect_2 <- drop(L_2 %*% stats::rnorm(2L)) * posterior[draw, sd_name]
    expected_random[4:5, draw] <- effect_2
  }
  expected <- expected_random + matrix(
    posterior[, "mu_intercept"],
    nrow = nrow(new_data),
    ncol = nrow(posterior),
    byrow = TRUE
  )

  expect_equal(unname(prediction), unname(expected), tolerance = 1e-12)
  expect_equal(
    requested_columns,
    rep(list(c(2L, 113L), c(50L, 113L)), nrow(posterior))
  )
  expect_lt(max(lengths(requested_columns)), length(factor_levels))
})

test_that("structured new-level subset sampling covers HCS, AR1, HAR, and CAR", {

  factor_levels <- letters[1:5]
  time_values   <- c(0, 0.5, 2, 5, 9)
  factor_data <- data.frame(
    f = factor(rep(factor_levels, 2L), levels = factor_levels),
    id = factor(rep(c("old_1", "old_2"), each = length(factor_levels)))
  )
  car_data <- data.frame(
    time = rep(time_values, 2L),
    id = factor(rep(c("old_1", "old_2"), each = length(time_values)))
  )
  specifications <- list(
    hcs = list(formula = ~ 1 + hcs(f | id), data = factor_data),
    ar1 = list(formula = ~ 1 + ar1(f | id), data = factor_data),
    har = list(formula = ~ 1 + har(f | id), data = factor_data),
    car = list(formula = ~ 1 + car(time | id), data = car_data)
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  subset_calls <- character()
  subset_cholesky <- BayesTools:::.bt_random_effect_structured_subset_cholesky
  prediction_subset_transform <-
    BayesTools:::.bt_random_effect_prediction_structured_subset_transform
  testthat::local_mocked_bindings(
    .bt_random_effect_marginal_covariance_correlation_draws = function(...){
      stop("global correlation materialized", call. = FALSE)
    },
    .bt_random_effect_prediction_structured_subset_transform = function(
        structure, columns, latent, rho, coordinates){

      subset_calls <<- c(subset_calls, structure)
      prediction_subset_transform(
        structure = structure,
        columns = columns,
        latent = latent,
        rho = rho,
        coordinates = coordinates
      )
    },
    .package = "BayesTools"
  )

  for(structure in names(specifications)){
    specification <- specifications[[structure]]
    result <- JAGS_formula(
      formula = specification$formula,
      parameter = "mu",
      data = specification$data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          rho = prior("normal", list(0, 0.5))
        )
      )
    )
    random_term <- result$formula_design$random_effects[[1L]]
    sd_names    <- unique(random_term$sd_parameter_names)
    rho_name    <- random_term$correlation$rho_name
    posterior_values <- c(
      mu_intercept = 0,
      stats::setNames(seq_along(sd_names) + 1, sd_names),
      stats::setNames(0.35, rho_name)
    )
    posterior <- matrix(
      posterior_values,
      nrow = 1L,
      dimnames = list(NULL, names(posterior_values))
    )
    fit <- coda::mcmc(posterior)
    attr(fit, "formula_design") <- list(mu = result$formula_design)
    new_data <- if(identical(structure, "car")){
      data.frame(
        time = time_values[c(1L, 3L, 5L)],
        id = factor(rep("new", 3L), levels = c("old_1", "old_2", "new"))
      )
    }else{
      data.frame(
        f = factor(factor_levels[c(1L, 3L, 5L)], levels = factor_levels),
        id = factor(rep("new", 3L), levels = c("old_1", "old_2", "new"))
      )
    }

    set.seed(914)
    prediction <- JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      prior_list = result$prior_list,
      formula_target = "conditional",
      new_levels = "sample"
    )

    columns <- c(1L, 3L, 5L)
    coordinates <- if(identical(structure, "car")){
      time_values
    }else{
      seq_along(factor_levels)
    }
    sd_draws <- BayesTools:::.bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = length(factor_levels),
      posterior = posterior,
      prior_list = result$prior_list
    )
    L <- subset_cholesky(
      structure = structure,
      columns = columns,
      rho = posterior[1L, rho_name],
      global_n_columns = length(factor_levels),
      column_coordinates = coordinates
    )
    set.seed(914)
    expected <- drop(L %*% stats::rnorm(length(columns))) *
      sd_draws[1L, columns]

    expect_equal(
      unname(drop(prediction)),
      unname(expected),
      tolerance = 1e-12,
      info = structure
    )
  }

  expect_equal(subset_calls, names(specifications))
})

test_that("stored new-level policies are used when no override is supplied", {

  result <- .formula_prediction_result(
    new_levels = random_new_levels(method = "zero")
  )
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  implicit_zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional"
  )
  explicit_zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional",
    new_levels = "zero"
  )
  expect_equal(implicit_zero, explicit_zero)

  marginal_zero <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "marginal",
    marginal_method = "covariance"
  )
  expect_true(all(marginal_zero$vcov$samples[, 2:3, ] == 0))
  expect_true(all(marginal_zero$vcov$samples[, , 2:3] == 0))
})

test_that("JAGS_predict_formula composes fixed means with marginal covariance", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      blocks = "id"
    ),
    "'blocks' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      new_levels = "zero"
    ),
    "'new_levels' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "conditional",
      marginal_method = "sample"
    ),
    "'marginal_method' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "marginal"
    ),
    "explicit new-level policy",
    fixed = TRUE
  )

  prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "marginal",
    marginal_method = "covariance",
    new_levels = "sample"
  )
  direct_vcov <- random_effects_marginal_vcov(
    fit = fit,
    parameter = "mu",
    data = new_data,
    posterior_samples = as.matrix(fit),
    prior_list = result$prior_list,
    new_levels = "sample"
  )

  expect_s3_class(prediction, "BayesTools_formula_prediction")
  expect_equal(
    unname(prediction$value),
    unname(JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "fixed"
    ))
  )
  expect_equal(prediction$vcov$samples, direct_vcov$samples)
})

test_that("marginal target includes marginalized random-effect blocks", {

  result <- .formula_prediction_result(
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  random_term <- result$formula_design$random_effects[[1]]
  posterior <- matrix(
    c(10, 1, 2, 3),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      random_term$sd_parameter_names
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "conditional"
    ),
    "compiled as marginalized",
    fixed = TRUE
  )

  prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "marginal",
    marginal_method = "covariance"
  )
  expect_equal(prediction$vcov$metadata$included_blocks, "id")
  expect_equal(prediction$vcov$metadata$blocks$id$compile_mode, "marginalized")
})
