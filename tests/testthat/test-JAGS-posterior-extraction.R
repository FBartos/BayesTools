# ============================================================================ #
# TEST FILE: Posterior Density Extraction Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for posterior extraction helper functions including
#   .extract_posterior_samples and .remove_auxiliary_parameters.
#
# DEPENDENCIES:
#   - rjags, runjags, coda: For JAGS model handling
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"), skip_if_not_installed("runjags")
#   - Note: Creates mock objects, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates mock runjags objects for testing
#
# TAGS: @evaluation, @JAGS, @posterior-extraction
# ============================================================================ #

# Tests for posterior extraction helper functions
test_that(".extract_posterior_samples extracts samples correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  
  # Load runjags to ensure S3 methods are registered
  library(runjags)

  # Create a proper runjags object structure for testing
  # The runjags package has an as.mcmc method that handles mcmc.list objects
  set.seed(123)
  mcmc1 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc2 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc_list <- coda::mcmc.list(mcmc1, mcmc2)
  
  # Create a minimal runjags object
  fit <- structure(
    list(mcmc = mcmc_list),
    class = c("runjags", "list")
  )

  # Test matrix extraction (as_list = FALSE)
  # This calls coda::as.mcmc on the runjags object which returns an mcmc object
  samples_matrix <- BayesTools:::.extract_posterior_samples(fit, as_list = FALSE)
  # mcmc objects inherit from matrix
  expect_true(inherits(samples_matrix, "mcmc"))
  expect_equal(ncol(samples_matrix), 1)
  expect_true("mu" %in% colnames(samples_matrix))
  expect_equal(nrow(samples_matrix), 200)  # 100 samples x 2 chains merged

  # Test list extraction (as_list = TRUE)
  samples_list <- BayesTools:::.extract_posterior_samples(fit, as_list = TRUE)
  expect_true(inherits(samples_list, "mcmc.list"))
  expect_equal(length(samples_list), 2) # 2 chains
})


test_that(".remove_auxiliary_parameters removes invgamma support", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with invgamma support parameter
  model_samples <- matrix(rnorm(100), ncol = 2)
  colnames(model_samples) <- c("sigma", "inv_sigma")

  prior_list <- list(
    sigma = prior("invgamma", list(1, 1))
  )

  result <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)

  expect_false("inv_sigma" %in% colnames(result$model_samples))
  expect_true("sigma" %in% colnames(result$model_samples))
})


test_that(".remove_auxiliary_parameters removes indexed factor invgamma support", {
  model_samples <- matrix(rnorm(400), ncol = 4)
  colnames(model_samples) <- c("theta[1]", "theta[2]", "inv_theta[1]", "inv_theta[2]")

  theta_prior <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  theta_prior$parameters$K <- 2
  prior_list <- list(theta = theta_prior)

  result <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)

  expect_equal(colnames(result$model_samples), c("theta[1]", "theta[2]"))
})


test_that(".process_spike_and_slab handles conditional samples", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with spike and slab
  model_samples <- matrix(c(
    rnorm(50, 0, 1),  # mu values
    rep(1, 50)        # indicator (all in slab)
  ), ncol = 2)
  colnames(model_samples) <- c("mu", "mu_indicator")

  prior_list <- list(
    mu = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("spike", list(0.5))
    )
  )

  result <- BayesTools:::.process_spike_and_slab(
    model_samples, prior_list, "mu",
    conditional = TRUE, remove_inclusion = FALSE, warnings = NULL
  )

  expect_true("mu (inclusion)" %in% colnames(result$model_samples))
  expect_false("mu_indicator" %in% colnames(result$model_samples))
  expect_true(is.prior.simple(result$prior_list$mu))
})


test_that(".apply_parameter_transformations applies transformations", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples
  model_samples <- matrix(rnorm(100, 0, 1), ncol = 1)
  colnames(model_samples) <- "mu"

  prior_list <- list(
    mu = prior("normal", list(0, 1))
  )

  # Apply exp transformation
  transformations <- list(
    mu = list(fun = exp, arg = list())
  )

  result <- BayesTools:::.apply_parameter_transformations(
    model_samples, transformations, prior_list
  )

  expect_true(all(result[, "mu"] > 0))  # exp makes all values positive
  expect_equal(ncol(result), 1)
})


test_that(".rename_factor_levels renames treatment factors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with factor
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior_obj, "levels") <- 4  # 4 levels total (treatment has K-1 parameters for K levels)
  attr(prior_obj, "level_names") <- c("A", "B", "C", "D")  # Should be a vector, not a list
  
  prior_list <- list(group = prior_obj)

  result <- BayesTools:::.rename_factor_levels(model_samples, prior_list)

  expect_true("group[B]" %in% colnames(result))
  expect_true("group[C]" %in% colnames(result))
  expect_true("group[D]" %in% colnames(result))
  expect_false("group[1]" %in% colnames(result))
})


test_that(".rename_factor_levels keeps interaction level on the factor term", {

  model_samples <- matrix(rnorm(400), ncol = 4)
  colnames(model_samples) <- c(
    "mu_alloc[1]",
    "mu_alloc[2]",
    "mu_alloc__xXx__year[1]",
    "mu_alloc__xXx__year[2]"
  )

  alloc_prior <- prior_factor("normal", list(0, 1), contrast = "independent")
  attr(alloc_prior, "levels") <- 2
  attr(alloc_prior, "level_names") <- c("random", "systematic")

  interaction_prior <- prior_factor("normal", list(0, 1), contrast = "independent")
  attr(interaction_prior, "levels") <- 2
  attr(interaction_prior, "level_names") <- list(alloc = c("random", "systematic"))
  attr(interaction_prior, "interaction") <- TRUE

  prior_list <- list(
    mu_alloc = alloc_prior,
    mu_alloc__xXx__year = interaction_prior
  )

  renamed <- BayesTools:::.rename_factor_levels(model_samples, prior_list)

  expect_equal(
    colnames(renamed),
    c(
      "mu_alloc[random]",
      "mu_alloc[systematic]",
      "mu_alloc[random]__xXx__year",
      "mu_alloc[systematic]__xXx__year"
    )
  )
  expect_equal(
    format_parameter_names(colnames(renamed), formula_parameters = "mu"),
    c(
      "(mu) alloc[random]",
      "(mu) alloc[systematic]",
      "(mu) alloc[random]:year",
      "(mu) alloc[systematic]:year"
    )
  )
})


test_that(".rename_factor_levels handles multi-factor independent interactions", {

  model_samples <- matrix(rnorm(600), ncol = 6)
  colnames(model_samples) <- paste0("mu_a__xXx__year__xXx__b[", 1:6, "]")

  interaction_prior <- prior_factor("normal", list(0, 1), contrast = "independent")
  attr(interaction_prior, "levels") <- 6
  attr(interaction_prior, "level_names") <- list(
    a = c("a1", "a2"),
    b = c("b1", "b2", "b3")
  )
  attr(interaction_prior, "interaction") <- TRUE

  renamed <- BayesTools:::.rename_factor_levels(
    model_samples,
    list(mu_a__xXx__year__xXx__b = interaction_prior)
  )

  expect_equal(
    colnames(renamed),
    c(
      "mu_a[a1]__xXx__year__xXx__b[b1]",
      "mu_a[a2]__xXx__year__xXx__b[b1]",
      "mu_a[a1]__xXx__year__xXx__b[b2]",
      "mu_a[a2]__xXx__year__xXx__b[b2]",
      "mu_a[a1]__xXx__year__xXx__b[b3]",
      "mu_a[a2]__xXx__year__xXx__b[b3]"
    )
  )
  expect_equal(
    format_parameter_names(colnames(renamed), formula_parameters = "mu"),
    c(
      "(mu) a[a1]:year:b[b1]",
      "(mu) a[a2]:year:b[b1]",
      "(mu) a[a1]:year:b[b2]",
      "(mu) a[a2]:year:b[b2]",
      "(mu) a[a1]:year:b[b3]",
      "(mu) a[a2]:year:b[b3]"
    )
  )
})


test_that(".rename_factor_levels handles multi-factor treatment interactions", {

  model_samples <- matrix(rnorm(200), ncol = 2)
  colnames(model_samples) <- paste0("mu_a__xXx__year__xXx__b[", 1:2, "]")

  interaction_prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(interaction_prior, "levels") <- 3
  attr(interaction_prior, "level_names") <- list(
    a = c("a1", "a2"),
    b = c("b1", "b2", "b3")
  )
  attr(interaction_prior, "interaction") <- TRUE

  renamed <- BayesTools:::.rename_factor_levels(
    model_samples,
    list(mu_a__xXx__year__xXx__b = interaction_prior)
  )

  expect_equal(
    colnames(renamed),
    c(
      "mu_a[a2]__xXx__year__xXx__b[b2]",
      "mu_a[a2]__xXx__year__xXx__b[b3]"
    )
  )
  expect_equal(
    format_parameter_names(colnames(renamed), formula_parameters = "mu"),
    c(
      "(mu) a[a2]:year:b[b2]",
      "(mu) a[a2]:year:b[b3]"
    )
  )
})


test_that(".transform_factor_contrasts transforms orthonormal to differences", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with orthonormal contrasts
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- 4  # 4 levels total (orthonormal has K-1 parameters for K levels)
  attr(prior_obj, "level_names") <- c("A", "B", "C", "D")  # Should be a vector, not a list
  
  prior_list <- list(group = prior_obj)

  expect_silent(
    result <- BayesTools:::.transform_factor_contrasts(
      model_samples, prior_list, transform_factors = TRUE
    )
  )

  expect_message(
    BayesTools:::.transform_factor_contrasts(
      model_samples,
      prior_list,
      transform_factors = TRUE,
      transformations = list(group = list(fun = exp, arg = list()))
    ),
    "transformation was applied"
  )

  # Should have 4 columns after transformation (one per level)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), paste0("group[dif: ", c("A", "B", "C", "D"), "]"))
})


test_that(".transform_factor_contrasts handles multi-factor transformed interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )

  for(interaction_contrast in c("orthonormal", "meandif")){
    formula_result <- JAGS_formula(
      formula = ~ a * b,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
        b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
        "a:b"     = prior_factor("mnormal", list(0, 1), contrast = interaction_contrast)
      )
    )
    interaction_prior <- formula_result$prior_list$mu_a__xXx__b

    model_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
    colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

    transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
      model_samples,
      list(mu_a__xXx__b = interaction_prior),
      transform_factors = TRUE
    ))
    expected <- model_samples %*% t(attr(interaction_prior, "factor_design"))

    expect_equal(unname(transformed), unname(expected))
    expected_names <- paste0(
      "mu_a[dif: ",
      rep(c("a1", "a2"), times = 3),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2", "b3"), each = 2),
      "]"
    )
    expect_equal(
      colnames(transformed),
      expected_names
    )
    expect_equal(anyDuplicated(colnames(transformed)), 0L)
  }
})

test_that(".transform_factor_contrasts handles one-coefficient interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2"), levels = c("b1", "b2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  expect_equal(BayesTools:::.JAGS_prior_factor_names("mu_a__xXx__b", interaction_prior), "mu_a__xXx__b")

  model_samples <- matrix(seq_len(10), nrow = 10, ncol = 1)
  colnames(model_samples) <- "mu_a__xXx__b"

  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = interaction_prior),
    transform_factors = TRUE
  ))
  expected <- model_samples %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(unname(transformed), unname(expected))
  expect_equal(ncol(transformed), 4L)
  expect_equal(
    colnames(transformed),
    paste0(
      "mu_a[dif: ",
      rep(c("a1", "a2"), times = 2),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2"), each = 2),
      "]"
    )
  )
})

test_that(".transform_factor_contrasts ignores unrelated transformations", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2"), levels = c("b1", "b2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  model_samples <- matrix(seq_len(10), nrow = 5, ncol = 2)
  colnames(model_samples) <- c("mu_intercept", "mu_a__xXx__b")

  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(
      mu_intercept = formula_result$prior_list$mu_intercept,
      mu_a__xXx__b = interaction_prior
    ),
    transform_factors = TRUE,
    transformations = list(mu_intercept = list(fun = exp))
  ))

  expect_equal(unname(transformed[, -1, drop = FALSE]), unname(model_samples[, 2, drop = FALSE] %*% t(attr(interaction_prior, "factor_design"))))
})

test_that(".transform_factor_contrasts reconstructs multi-factor designs from metadata", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b
  expected_design <- attr(interaction_prior, "factor_design")

  model_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

  metadata_only_prior <- interaction_prior
  attr(metadata_only_prior, "factor_design") <- NULL
  attr(metadata_only_prior, "factor_cell_names") <- NULL
  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = metadata_only_prior),
    transform_factors = TRUE
  ))
  expect_equal(unname(transformed), unname(model_samples %*% t(expected_design)))

  unnamed_contrast_prior <- metadata_only_prior
  attr(unnamed_contrast_prior, "factor_contrasts") <- unname(attr(unnamed_contrast_prior, "factor_contrasts"))
  transformed_unnamed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = unnamed_contrast_prior),
    transform_factors = TRUE
  ))
  expect_equal(unname(transformed_unnamed), unname(model_samples %*% t(expected_design)))

  inferred_contrast_prior <- metadata_only_prior
  attr(inferred_contrast_prior, "factor_contrasts") <- NULL
  model_samples_full <- cbind(
    mu_a = seq_len(nrow(model_samples)),
    `mu_b[1]` = seq_len(nrow(model_samples)) + 10,
    `mu_b[2]` = seq_len(nrow(model_samples)) + 20,
    model_samples
  )
  transformed_inferred <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples_full,
    list(
      mu_a = formula_result$prior_list$mu_a,
      mu_b = formula_result$prior_list$mu_b,
      mu_a__xXx__b = inferred_contrast_prior
    ),
    transform_factors = TRUE
  ))
  expected_names <- paste0(
    "mu_a[dif: ",
    rep(c("a1", "a2"), times = 3),
    "]__xXx__b[dif: ",
    rep(c("b1", "b2", "b3"), each = 2),
    "]"
  )
  expect_equal(unname(transformed_inferred[, expected_names]), unname(model_samples %*% t(expected_design)))
})

test_that(".transform_factor_contrasts validates multi-factor metadata", {

  prior_obj <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- 3
  attr(prior_obj, "level_names") <- list(a = c("a1", "a2"), b = c("b1", "b2"))
  attr(prior_obj, "interaction") <- TRUE
  attr(prior_obj, "factor_terms") <- c("a", "b")
  attr(prior_obj, "factor_contrasts") <- c(a = "contr.orthonormal")

  model_samples <- matrix(seq_len(10), nrow = 5, ncol = 2)
  colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

  expect_error(
    BayesTools:::.factor_term_design_from_metadata(prior_obj),
    "incomplete"
  )

  missing_contrasts_prior <- prior_obj
  attr(missing_contrasts_prior, "factor_contrasts") <- NULL
  expect_error(
    BayesTools:::.factor_term_design_from_metadata(missing_contrasts_prior),
    "missing"
  )

  missing_levels_prior <- prior_obj
  attr(missing_levels_prior, "level_names") <- NULL
  attr(missing_levels_prior, "levels") <- NULL
  expect_error(
    BayesTools:::.factor_term_design_from_metadata(missing_levels_prior),
    "level names"
  )

  mismatch_prior <- prior_obj
  attr(mismatch_prior, "factor_contrasts") <- c(a = "contr.orthonormal", b = "contr.orthonormal")
  attr(mismatch_prior, "factor_design") <- diag(3)
  expect_error(
    BayesTools:::.transform_factor_contrasts(
      model_samples,
      list(mu_a__xXx__b = mismatch_prior),
      transform_factors = TRUE
    ),
    "has 3 coefficient columns"
  )
})

test_that("as_mixed_posteriors propagates multi-factor contrast metadata", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_a__xXx__b[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list

  mixed <- as_mixed_posteriors(fit, parameters = "mu_a__xXx__b")

  expect_equal(attr(mixed$mu_a__xXx__b, "factor_terms"), attr(interaction_prior, "factor_terms"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_contrasts"), attr(interaction_prior, "factor_contrasts"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_design"), attr(interaction_prior, "factor_design"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_cell_names"), attr(interaction_prior, "factor_cell_names"))

  transformed <- transform_factor_samples(mixed)$mu_a__xXx__b
  expect_equal(
    as.vector(transformed),
    as.vector(posterior %*% t(attr(interaction_prior, "factor_design")))
  )
})


test_that(".filter_parameters removes spike at 0 priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("point", list(0)),  # spike at 0
    tau = prior("normal", list(1, 1))
  )

  # With remove_spike_0 = TRUE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = TRUE)
  expect_true("delta" %in% result)
  expect_false("mu" %in% result)
  expect_false("tau" %in% result)

  # With remove_spike_0 = FALSE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = FALSE)
  expect_equal(length(result), 0)
})


test_that(".filter_parameters removes character specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = c("mu", "tau"), remove_spike_0 = FALSE)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
  expect_false("delta" %in% result)
})


test_that(".filter_parameters removes non-formula parameters when TRUE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with formula attributes
  prior_formula <- prior("normal", list(0, 1))
  attr(prior_formula, "parameter") <- "y"

  prior_list <- list(
    intercept = prior_formula,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = TRUE, remove_spike_0 = FALSE)
  expect_true("sigma" %in% result)
  expect_false("intercept" %in% result)
})


test_that(".filter_parameters removes formula-specific parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_formulas = "y", remove_spike_0 = FALSE)
  expect_true("intercept_y" %in% result)
  expect_true("slope_y" %in% result)
  expect_false("intercept_x" %in% result)
  expect_false("sigma" %in% result)
})


test_that(".filter_parameters keeps only specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("delta" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters keeps only specified formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_true("intercept_x" %in% result)
  expect_true("sigma" %in% result)
})


test_that(".filter_parameters combines keep_parameters and keep_formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  # Keep formula "y" and parameter "sigma"
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "sigma", keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_false("sigma" %in% result)
  expect_true("intercept_x" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is removed", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component (simulating bias)
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  # When bias is removed, PET should also be removed
  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains PEESE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PEESE component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PEESE("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PEESE" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains weightfunction", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with weightfunction component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("omega" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters keeps bias-related parameters when bias is kept", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only bias is kept, mu and tau should be removed, but PET should be kept
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "bias", remove_spike_0 = FALSE)
  expect_false("bias" %in% result)
  expect_false("PET" %in% result)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters handles non-mixture bias priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a single PET prior named bias (not a mixture)
  bias_prior <- prior_PET("normal", list(0, 1))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is not in keep list", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only mu is kept, bias should be removed along with PET
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_true("tau" %in% result)
})


test_that("helper functions work with runjags estimates extraction", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Test the helper functions with mock data (not full integration)
  # This tests that our refactored code correctly uses the helpers
  
  # Create mock posterior samples
  set.seed(123)
  model_samples <- matrix(rnorm(200), ncol = 2, dimnames = list(NULL, c("mu", "inv_sigma")))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("invgamma", list(1, 1))
  )
  
  # Test that remove_auxiliary_parameters helper works
  cleaned <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)
  
  # Should remove inv_sigma
  expect_false("inv_sigma" %in% colnames(cleaned$model_samples))
  expect_true("mu" %in% colnames(cleaned$model_samples))
  expect_equal(ncol(cleaned$model_samples), 1)
})
