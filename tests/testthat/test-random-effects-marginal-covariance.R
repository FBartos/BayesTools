skip_if_not_test_profile("unit")

.re_cov_sd_prior <- function(){
  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.re_cov_fixed_priors <- function(){
  list(intercept = prior("normal", list(0, 1)))
}

.re_cov_formula <- function(formula, data, prior_random,
                            random_effects_compile = NULL){
  JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = .re_cov_fixed_priors(),
    prior_random = prior_random,
    random_effects_compile = random_effects_compile
  )
}

.re_cov_term <- function(result, block_name){
  random_effects <- result$formula_design$random_effects
  block_names <- vapply(random_effects, `[[`, character(1), "block_name")
  random_effects[[match(block_name, block_names)]]
}

.re_cov_posterior <- function(values){
  matrix(
    as.numeric(values),
    nrow = 1,
    dimnames = list(NULL, names(values))
  )
}

.re_cov_posterior_draws <- function(...){
  draws <- list(...)
  posterior <- do.call(rbind, draws)
  colnames(posterior) <- names(draws[[1L]])
  posterior
}

.re_cov_output <- function(result, posterior, data = NULL, blocks = NULL,
                           new_levels = NULL, diagonal_only = FALSE,
                           fitted_rows = NULL){
  random_effects_marginal_vcov(
    result$formula_design,
    data = data,
    posterior_samples = posterior,
    prior_list = result$prior_list,
    blocks = blocks,
    new_levels = new_levels,
    fitted_rows = fitted_rows,
    diagonal_only = diagonal_only
  )
}

.re_cov_first <- function(out){
  unname(out$samples[1, , ])
}

.re_cov_dense_diagonal <- function(out){
  n_draws <- dim(out$samples)[1L]
  n_rows  <- dim(out$samples)[2L]
  unname(t(vapply(
    seq_len(n_draws),
    function(draw) diag(out$samples[draw, , ]),
    numeric(n_rows)
  )))
}

test_that("single-coefficient correlation Cholesky draws retain dimensions", {

  posterior <- matrix(0, nrow = 2L, ncol = 1L)
  cholesky <- array(1, dim = c(2L, 1L, 1L))
  random_term <- list(block_name = "study")

  expect_no_error(
    .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
      cholesky   = cholesky,
      random_term = random_term,
      n_columns  = 1L,
      posterior  = posterior
    )
  )
})

.re_cov_expand <- function(model_matrix, group_map, G){
  out <- matrix(0, nrow(model_matrix), nrow(model_matrix))
  for(rows in split(seq_len(nrow(model_matrix)), group_map)){
    Z <- model_matrix[rows, , drop = FALSE]
    out[rows, rows] <- Z %*% G %*% t(Z)
  }
  unname(out)
}

.re_cov_sd_values <- function(random_term, sd){
  sd_names <- random_term$sd_parameter_names
  if(length(unique(sd_names)) == 1L){
    return(stats::setNames(sd[1L], unique(sd_names)))
  }

  stats::setNames(sd, sd_names)
}

.re_cov_rho_sample <- function(random_term, rho){
  stats::setNames(atanh(rho), random_term$correlation$sample_name)
}

.re_cov_cholesky_values <- function(random_term, L){
  out <- numeric()
  for(row in seq_len(nrow(L))){
    for(column in seq_len(ncol(L))){
      name <- paste0(
        random_term$parameter_stem,
        "_xRE_CORx_L[", row, ",", column, "]"
      )
      out[name] <- L[row, column]
    }
  }
  out
}

.re_cov_compound_correlation <- function(K, rho){
  R <- matrix(rho, K, K)
  diag(R) <- 1
  R
}

.re_cov_ar1_correlation <- function(K, rho){
  rho^abs(outer(seq_len(K), seq_len(K), "-"))
}

.re_cov_car_correlation <- function(time, rho){
  rho^abs(outer(time, time, "-"))
}

.re_cov_structured_data <- function(){
  data.frame(
    id = factor(rep(c("a", "b"), each = 3), levels = c("a", "b")),
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    time = rep(c(0, 2, 5), 2)
  )
}

.re_cov_expect_structured <- function(formula, structure, sd, rho, R){
  result <- .re_cov_formula(
    formula = formula,
    data = .re_cov_structured_data(),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  sd <- if(length(sd) == 1L){
    rep(sd, random_term$n_columns)
  }else{
    sd
  }

  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, sd),
    .re_cov_rho_sample(random_term, rho)
  ))
  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    R * tcrossprod(sd)
  )

  expect_equal(random_term$structure, structure)
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out),
               tolerance = 1e-12)
  expect_equal(out$metadata$structures, stats::setNames(structure, "id"))
}

test_that("diag random intercept and slope covariance uses ZGZ' by group", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, c(2, 3)))
  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    diag(c(2, 3)^2)
  )

  expect_s3_class(out, "BayesTools_random_effects_marginal_vcov")
  expect_equal(dim(out$samples), c(1L, 4L, 4L))
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out),
               tolerance = 1e-12)
  expect_equal(dim(diagonal$samples), c(1L, 4L))
  expect_identical(diagonal$metadata$representation, "diagonal")
  expect_true(diagonal$metadata$diagonal_only)
  expect_false(diagonal$metadata$dense)
  expect_true(is.na(diagonal$metadata$dense_entries))
  expect_equal(diagonal$metadata$sample_entries, 4)
  expect_equal(out$metadata$blocks$id$n_columns, 2L)
  expect_false(out$metadata$blocks$id$row_varying_sd)

  new_data <- data.frame(
    id = factor(c("a", "b", "a"), levels = c("a", "b")),
    x = c(4, 5, 6)
  )
  out_new <- .re_cov_output(result, posterior, data = new_data)
  Z_new <- cbind("(Intercept)" = 1, x = new_data$x)
  expected_new <- .re_cov_expand(Z_new, c(1L, 2L, 1L), diag(c(2, 3)^2))
  expect_equal(.re_cov_first(out_new), expected_new, tolerance = 1e-12)
  expect_equal(out_new$metadata$data_source, "data")
  expect_equal(out_new$metadata$n_rows, 3L)

  new_level_data <- data.frame(
    id = factor(c("a", "c", "c", "b"), levels = c("a", "b", "c")),
    x = c(4, 5, 6, 7)
  )
  expect_error(
    .re_cov_output(result, posterior, data = new_level_data),
    "explicit new-level policy",
    fixed = TRUE
  )
  out_new_level <- .re_cov_output(
    result,
    posterior,
    data = new_level_data,
    new_levels = "sample"
  )
  Z_new_level <- cbind("(Intercept)" = 1, x = new_level_data$x)
  expected_new_level <- .re_cov_expand(
    Z_new_level,
    c(1L, 3L, 3L, 2L),
    diag(c(2, 3)^2)
  )
  expect_equal(.re_cov_first(out_new_level), expected_new_level, tolerance = 1e-12)
  expect_equal(out_new_level$metadata$blocks$id$n_groups, 3L)
  expect_equal(out_new_level$metadata$blocks$id$fitted_n_groups, 2L)
  expect_equal(out_new_level$metadata$blocks$id$group_levels, c("a", "b", "c"))
  expect_equal(out_new_level$metadata$blocks$id$new_group_levels, "c")
  expect_equal(out_new_level$metadata$blocks$id$new_level_rows, c(2L, 3L))
  expect_equal(out_new_level$metadata$blocks$id$new_levels$method, "sample")

  out_zero_level <- .re_cov_output(
    result,
    posterior,
    data = new_level_data,
    new_levels = "zero"
  )
  expected_zero_level <- expected_new_level
  expected_zero_level[2:3, ] <- 0
  expected_zero_level[, 2:3] <- 0
  expect_equal(.re_cov_first(out_zero_level), expected_zero_level, tolerance = 1e-12)
  expect_equal(out_zero_level$metadata$blocks$id$new_group_levels, "c")
  expect_equal(out_zero_level$metadata$blocks$id$new_level_rows, c(2L, 3L))
  expect_equal(out_zero_level$metadata$blocks$id$new_levels$method, "zero")
  diagonal_zero_level <- .re_cov_output(
    result,
    posterior,
    data = new_level_data,
    new_levels = "zero",
    diagonal_only = TRUE
  )
  expect_equal(
    unname(diagonal_zero_level$samples),
    matrix(diag(expected_zero_level), nrow = 1L),
    tolerance = 1e-12
  )

  result_stored_zero <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        new_levels = random_new_levels(method = "zero")
      )
    )
  )
  random_term_stored_zero <- .re_cov_term(result_stored_zero, "id")
  posterior_stored_zero <- .re_cov_posterior(
    .re_cov_sd_values(random_term_stored_zero, c(2, 3))
  )
  out_stored_zero <- .re_cov_output(
    result_stored_zero,
    posterior_stored_zero,
    data = new_level_data
  )
  expect_equal(.re_cov_first(out_stored_zero), expected_zero_level, tolerance = 1e-12)
})

test_that("known group covariance uses tau squared times ZKZ prime", {

  df <- data.frame(
    id = factor(c("b", "a", "c", "b"), levels = c("b", "a", "c"))
  )
  K <- matrix(
    c(4, 1, .5,
      1, 9, 2,
      .5, 2, 16),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  result <- .re_cov_formula(
    formula = random_effects,
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, .5))
  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)
  expected <- unname(.5^2 * K[as.character(df$id), as.character(df$id)])

  expect_equal(random_term$group_covariance$scale, "none")
  expect_equal(
    BayesTools:::.bt_random_effect_sd_summary_label("intercept", "id", random_term),
    "sd_multiplier(intercept | id)"
  )
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out),
               tolerance = 1e-12)
  expect_equal(diag(.re_cov_first(out)), .5^2 * c(9, 4, 16, 9))
  expect_equal(
    out$metadata$blocks$id$group_covariance$kernel,
    K[levels(df$id), levels(df$id)]
  )

  new_data <- data.frame(
    id = factor(c("c", "a"), levels = c("b", "a", "c"))
  )
  out_new <- .re_cov_output(result, posterior, data = new_data)
  expect_equal(
    .re_cov_first(out_new),
    unname(.5^2 * K[as.character(new_data$id), as.character(new_data$id)]),
    tolerance = 1e-12
  )

  new_level_data <- data.frame(
    id = factor(c("a", "d"), levels = c("a", "b", "c", "d"))
  )
  expect_error(
    .re_cov_output(result, posterior, data = new_level_data),
    "known group covariance",
    fixed = TRUE
  )
  expect_error(
    .re_cov_output(result, posterior, data = new_level_data, new_levels = "sample"),
    "known group covariance",
    fixed = TRUE
  )

  marginalized_result <- .re_cov_formula(
    formula = random_effects,
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  marginalized_term <- .re_cov_term(marginalized_result, "id")
  marginalized_posterior <- .re_cov_posterior(
    .re_cov_sd_values(marginalized_term, .5)
  )
  marginalized_out <- .re_cov_output(marginalized_result, marginalized_posterior)

  expect_equal(marginalized_term$compile_mode, "marginalized")
  expect_equal(.re_cov_first(marginalized_out), expected, tolerance = 1e-12)
  expect_equal(
    marginalized_out$metadata$blocks$id$group_covariance$kernel,
    K[levels(df$id), levels(df$id)]
  )
  expect_error(
    .re_cov_output(marginalized_result, marginalized_posterior, data = new_level_data),
    "known group covariance",
    fixed = TRUE
  )
})

test_that("marginal variance factors support diagonal known group covariance scaling", {

  df <- data.frame(
    id = factor(c("b", "a", "c"), levels = c("b", "a", "c"))
  )
  K <- diag(c(4, 9, 16))
  dimnames(K) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expected <- list(
    none = c(9, 4, 16),
    cor = c(1, 1, 1),
    cor0 = c(1, 1, 1),
    cov0 = c(9, 4, 16)
  )

  for(scale in names(expected)){
    random_effects <- random_effects_formula(
      ~ 1 | id,
      group_covariance = random_group_covariance(K, scale = scale)
    )
    result <- .re_cov_formula(
      formula = random_effects,
      data = df,
      prior_random = prior_random(
        id = random_block(sd = .re_cov_sd_prior())
      ),
      random_effects_compile = random_effects_compile(marginalized = "id")
    )
    factors <- random_effects_marginal_variance_factors(
      result$formula_design,
      require_one_to_one = TRUE
    )
    block <- factors$blocks$id

    expect_s3_class(factors, "BayesTools_random_effects_marginal_variance_factors")
    expect_equal(factors$included_blocks, "id")
    expect_equal(block$compile_mode, "marginalized")
    expect_equal(block$group_covariance$scale, scale)
    expect_equal(unname(block$row_multiplier), expected[[scale]], tolerance = 1e-12)
    expect_true(block$row_covariance_diagonal)
    expect_true(block$one_to_one)
  }
})

test_that("marginal variance factors default to marginalized blocks and accept named blocks", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"),
                      levels = c("e1", "e2", "e3", "e4"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | estimate, name = "estimate", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      study = random_block(sd = .re_cov_sd_prior()),
      estimate = random_block(sd = .re_cov_sd_prior())
    ),
    random_effects_compile = random_effects_compile(marginalized = "estimate")
  )

  default <- random_effects_marginal_variance_factors(result$formula_design)
  expect_equal(default$included_blocks, "estimate")
  expect_equal(default$skipped_blocks$block_name, "study")
  expect_equal(default$skipped_blocks$reason, "not marginalized")
  expect_equal(unname(default$blocks$estimate$row_multiplier), rep(1, 4))
  expect_true(default$blocks$estimate$row_covariance_diagonal)

  named <- random_effects_marginal_variance_factors(
    result$formula_design,
    blocks = "study",
    require_diagonal = FALSE
  )
  expect_equal(named$included_blocks, "study")
  expect_equal(named$blocks$study$compile_mode, "sampled")
  expect_equal(named$skipped_blocks$block_name, "estimate")
  expect_equal(named$skipped_blocks$reason, "not requested")
  expect_error(
    random_effects_marginal_variance_factors(
      result$formula_design,
      blocks = c("study", "study")
    ),
    "'blocks' must be unique",
    fixed = TRUE
  )
  expect_error(
    random_effects_marginal_variance_factors(
      result$formula_design,
      blocks = "missing"
    ),
    "Unknown random-effect block",
    fixed = TRUE
  )
})

test_that("marginal variance factors protect diagonal-only known covariance use", {

  df <- data.frame(
    id = factor(c("b", "a", "c"), levels = c("b", "a", "c"))
  )
  K <- matrix(
    c(4, 1, .5,
      1, 9, 2,
      .5, 2, 16),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  result <- .re_cov_formula(
    formula = random_effects,
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )

  unavailable <- tryCatch(
    random_effects_marginal_variance_factors(result$formula_design),
    error = identity
  )
  expect_s3_class(
    unavailable,
    "BayesTools_random_effects_marginal_variance_unavailable"
  )
  expect_match(
    conditionMessage(unavailable),
    "off-diagonal row covariance",
    fixed = TRUE
  )
  expect_identical(unavailable$block_name, "id")
  expect_identical(unavailable$reason, "non_diagonal_row_covariance")
  factors <- random_effects_marginal_variance_factors(
    result$formula_design,
    require_diagonal = FALSE
  )
  expected_row_covariance <- K[as.character(df$id), as.character(df$id)]

  expect_false(factors$blocks$id$row_covariance_diagonal)
  expect_equal(
    unname(factors$blocks$id$row_multiplier),
    unname(diag(expected_row_covariance)),
    tolerance = 1e-12
  )
  expect_gt(factors$blocks$id$max_off_diagonal, 0)
})

test_that("marginal variance factor diagonal checks are scale invariant", {

  df <- data.frame(
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  correlation <- matrix(
    c(1, .5,
      .5, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  scales <- c(1e-20, 1, 1e20)

  for(scale in scales){
    K <- scale * correlation
    random_effects <- random_effects_formula(
      ~ 1 | id,
      group_covariance = random_group_covariance(K, scale = "none")
    )
    result <- .re_cov_formula(
      formula = random_effects,
      data = df,
      prior_random = prior_random(
        id = random_block(sd = .re_cov_sd_prior())
      ),
      random_effects_compile = random_effects_compile(marginalized = "id")
    )

    unavailable <- tryCatch(
      random_effects_marginal_variance_factors(result$formula_design),
      error = identity
    )
    expect_s3_class(
      unavailable,
      "BayesTools_random_effects_marginal_variance_unavailable"
    )
    expect_identical(
      unavailable$reason,
      "non_diagonal_row_covariance",
      info = paste("scale =", scale)
    )

    factors <- random_effects_marginal_variance_factors(
      result$formula_design,
      require_diagonal = FALSE
    )
    expect_false(
      factors$blocks$id$row_covariance_diagonal,
      info = paste("scale =", scale)
    )
  }

  diagonal_status <- lapply(
    scales,
    function(scale){
      BayesTools:::.bt_random_effect_row_covariance_diagonal_status(
        scale * diag(c(1, 4))
      )
    }
  )
  expect_true(all(vapply(diagonal_status, `[[`, logical(1), "is_diagonal")))

  heterogeneous <- diag(c(1e20, 1e-8, 1e-8))
  heterogeneous[2L, 3L] <- heterogeneous[3L, 2L] <- 5e-9
  heterogeneous_status <-
    BayesTools:::.bt_random_effect_row_covariance_diagonal_status(heterogeneous)
  expect_false(heterogeneous_status$is_diagonal)
})

test_that("marginal variance factors reject repeated groups when required", {

  df <- data.frame(
    id = factor(c("a", "a", "b"), levels = c("a", "b"))
  )
  K <- diag(c(4, 9))
  dimnames(K) <- list(c("a", "b"), c("a", "b"))
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  result <- .re_cov_formula(
    formula = random_effects,
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )

  unavailable <- tryCatch(
    random_effects_marginal_variance_factors(
      result$formula_design,
      require_diagonal = FALSE,
      require_one_to_one = TRUE
    ),
    error = identity
  )
  expect_s3_class(
    unavailable,
    "BayesTools_random_effects_marginal_variance_unavailable"
  )
  expect_match(conditionMessage(unavailable), "one-to-one", fixed = TRUE)
  expect_identical(unavailable$block_name, "id")
  expect_identical(unavailable$reason, "repeated_groups")
})

test_that("marginal variance factor limitations signal stable conditions", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  row_indexed <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd_source = random_sd_source("tau", shape = "row"))
    )
  )
  multiple_columns <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )

  cases <- list(
    row_indexed_sd = row_indexed$formula_design,
    multiple_columns = multiple_columns$formula_design
  )
  messages <- c(
    row_indexed_sd = "row-indexed external SD sources",
    multiple_columns = "one random-effect column"
  )
  for(reason in names(cases)){
    unavailable <- tryCatch(
      random_effects_marginal_variance_factors(
        cases[[reason]],
        blocks = "id",
        require_diagonal = FALSE
      ),
      error = identity
    )
    expect_s3_class(
      unavailable,
      "BayesTools_random_effects_marginal_variance_unavailable"
    )
    expect_match(
      conditionMessage(unavailable),
      messages[[reason]],
      fixed = TRUE
    )
    expect_identical(unavailable$block_name, "id")
    expect_identical(unavailable$reason, reason)
  }
})

test_that("id covariance shares one SD across independent columns", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "id"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, 2))
  out <- .re_cov_output(result, posterior)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    diag(rep(2^2, random_term$n_columns))
  )

  expect_equal(random_term$structure, "id")
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
})

test_that("diag covariance uses draw-specific SD samples", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior_draws(
    .re_cov_sd_values(random_term, c(2, 3)),
    .re_cov_sd_values(random_term, c(4, 5))
  )
  out <- .re_cov_output(result, posterior)
  expected_1 <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    diag(c(2, 3)^2)
  )
  expected_2 <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    diag(c(4, 5)^2)
  )

  expect_equal(out$metadata$n_draws, 2L)
  expect_equal(unname(out$samples[1, , ]), expected_1, tolerance = 1e-12)
  expect_equal(unname(out$samples[2, , ]), expected_2, tolerance = 1e-12)
})

test_that("us covariance uses monitored LKJ Cholesky samples", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id"),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        cor = prior_lkj(eta = 2)
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  rho <- 0.5
  L <- matrix(c(1, 0, rho, sqrt(1 - rho^2)), 2, 2, byrow = TRUE)
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, c(2, 3)),
    .re_cov_cholesky_values(random_term, L)
  ))
  out <- .re_cov_output(result, posterior)
  diagonal <- testthat::with_mocked_bindings(
    .re_cov_output(result, posterior, diagonal_only = TRUE),
    .bt_random_effect_marginal_covariance_correlation_draws = function(...) {
      stop("dense correlation reconstruction must not be used", call. = FALSE)
    },
    .package = "BayesTools"
  )
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    R * tcrossprod(c(2, 3))
  )

  expect_equal(random_term$structure, "us")
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out),
               tolerance = 1e-12)
  expect_equal(out$metadata$blocks$id$correlation_type, "lkj")

  bad_posterior <- posterior
  bad_posterior[1, paste0(random_term$parameter_stem, "_xRE_CORx_L[2,2]")] <- NA_real_
  expect_error(
    .re_cov_output(result, bad_posterior),
    "must be finite",
    fixed = TRUE
  )
})

test_that("us covariance can reconstruct LKJ primitive samples", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3),
    z = c(3, 2, 1, 0)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x + z | id, name = "id"),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        cor = prior_lkj(eta = 2, include_primitives = TRUE)
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  u <- c(0.2, 0.6, 0.8)
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(matrix(u, nrow = 1L), K = 3L)
  if(length(dim(L)) == 3L){
    L <- L[1L, , ]
  }
  sd <- c(2, 3, 4)
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, sd),
    stats::setNames(u, random_term$correlation$primitive_names)
  ))
  out <- .re_cov_output(result, posterior)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    tcrossprod(L) * tcrossprod(sd)
  )

  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)

  duplicated_result <- result
  duplicated_names <- random_term$correlation$primitive_names
  duplicated_names[2L] <- duplicated_names[1L]
  duplicated_correlation <-
    duplicated_result$formula_design$random_effects[[1L]]$correlation
  duplicated_correlation$primitive_names <- duplicated_names
  names(duplicated_correlation$primitive_bounds$lb) <- duplicated_names
  names(duplicated_correlation$primitive_bounds$ub) <- duplicated_names
  duplicated_result$formula_design$random_effects[[1L]]$correlation <-
    duplicated_correlation

  expect_error(
    .re_cov_output(duplicated_result, posterior),
    "must define unique 'random_term$correlation$primitive_names'",
    fixed = TRUE
  )
})

test_that("cs covariance uses scalar-rho compound symmetry", {

  .re_cov_expect_structured(
    formula = ~ 1 + cs(f | id),
    structure = "cs",
    sd = 2,
    rho = 0.25,
    R = .re_cov_compound_correlation(3L, 0.25)
  )
})

test_that("hcs covariance uses heterogeneous SDs with compound symmetry", {

  .re_cov_expect_structured(
    formula = ~ 1 + hcs(f | id),
    structure = "hcs",
    sd = c(2, 3, 4),
    rho = 0.25,
    R = .re_cov_compound_correlation(3L, 0.25)
  )
})

test_that("ar1 covariance uses discrete autoregressive distances", {

  .re_cov_expect_structured(
    formula = ~ 1 + ar1(f | id),
    structure = "ar1",
    sd = 2,
    rho = 0.5,
    R = .re_cov_ar1_correlation(3L, 0.5)
  )
})

test_that("har covariance uses heterogeneous SDs with AR1 correlations", {

  .re_cov_expect_structured(
    formula = ~ 1 + har(f | id),
    structure = "har",
    sd = c(2, 3, 4),
    rho = 0.5,
    R = .re_cov_ar1_correlation(3L, 0.5)
  )
})

test_that("car covariance uses continuous-time distances", {

  .re_cov_expect_structured(
    formula = ~ 1 + car(time | id),
    structure = "car",
    sd = 2,
    rho = 0.5,
    R = .re_cov_car_correlation(c(0, 2, 5), 0.5)
  )
})

test_that("car one-hot covariance routes through the stable factor", {

  result <- .re_cov_formula(
    formula = ~ 1 + car(time | id),
    data = .re_cov_structured_data(),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, 2),
    .re_cov_rho_sample(random_term, 0.5)
  ))
  outputs <- testthat::with_mocked_bindings(
    list(
      full = .re_cov_output(result, posterior),
      diagonal = .re_cov_output(
        result,
        posterior,
        diagonal_only = TRUE
      )
    ),
    .bt_random_effect_marginal_variance_one_sparse = function(...) {
      stop("one-sparse shortcut must not be used", call. = FALSE)
    },
    .bt_random_effect_marginal_covariance_structured_one_hot = function(...) {
      stop("one-hot shortcut must not be used", call. = FALSE)
    },
    .bt_random_effect_marginal_covariance_correlation_draws = function(...) {
      stop("dense correlation reconstruction must not be used", call. = FALSE)
    },
    .package = "BayesTools"
  )
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    .re_cov_car_correlation(c(0, 2, 5), 0.5) * 2^2
  )

  expect_equal(.re_cov_first(outputs$full), expected, tolerance = 1e-12)
  expect_equal(
    unname(outputs$diagonal$samples),
    .re_cov_dense_diagonal(outputs$full),
    tolerance = 1e-12
  )
})

test_that("car covariance preserves upper-rho half-gap contrast variance", {

  result <- .re_cov_formula(
    formula = ~ 1 + car(time | id),
    data = data.frame(
      id = factor(c("a", "a")),
      time = c(1, 1.5)
    ),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  random_term$model_matrix <- matrix(
    c(1, -1),
    nrow = 1L,
    dimnames = list("contrast", random_term$column_names)
  )
  random_term$group_map <- 1L
  result$formula_design$random_effects[[1L]] <- random_term

  rho <- 1 - .Machine$double.eps / 2
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, 1),
    .re_cov_rho_sample(random_term, rho)
  ))
  outputs <- testthat::with_mocked_bindings(
    list(
      full = .re_cov_output(result, posterior),
      diagonal = .re_cov_output(
        result,
        posterior,
        diagonal_only = TRUE
      )
    ),
    .bt_random_effect_marginal_covariance_correlation_draws = function(...) {
      stop("dense correlation reconstruction must not be used", call. = FALSE)
    },
    .bt_random_effect_marginal_variance_one_sparse = function(...) {
      stop("one-sparse shortcut must not be used", call. = FALSE)
    },
    .bt_random_effect_marginal_covariance_structured_one_hot = function(...) {
      stop("one-hot shortcut must not be used", call. = FALSE)
    },
    .package = "BayesTools"
  )
  log_phi <- 0.5 * log(rho)
  expected <- (1 - exp(log_phi))^2 - expm1(2 * log_phi)

  expect_identical(dim(outputs$full$samples), c(1L, 1L, 1L))
  expect_identical(dim(outputs$diagonal$samples), c(1L, 1L))
  expect_gt(outputs$full$samples[1L, 1L, 1L], 0)
  expect_gt(outputs$diagonal$samples[1L, 1L], 0)
  expect_equal(
    unname(outputs$full$samples[1L, 1L, 1L]) / expected,
    1,
    tolerance = 1e-12
  )
  expect_equal(
    unname(outputs$diagonal$samples[1L, 1L]) / expected,
    1,
    tolerance = 1e-12
  )
})

test_that("car covariance validates every compact coordinate copy", {

  result <- .re_cov_formula(
    formula = ~ 1 + car(time | id),
    data = .re_cov_structured_data(),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, 1),
    .re_cov_rho_sample(random_term, 0.5)
  ))
  with_coordinates <- function(correlation_time, car_time){
    modified <- result
    modified$formula_design$random_effects[[1L]]$correlation$time_values <-
      correlation_time
    modified$formula_design$random_effects[[1L]]$car$time_values <- car_time
    modified
  }

  duplicate <- with_coordinates(c(0, 0, 5), c(0, 0, 5))
  expect_error(
    .re_cov_output(duplicate, posterior, diagonal_only = TRUE),
    "finite, unique, and strictly increasing",
    fixed = TRUE
  )

  unordered <- with_coordinates(c(0, 5, 2), c(0, 5, 2))
  expect_error(
    .re_cov_output(unordered, posterior),
    "finite, unique, and strictly increasing",
    fixed = TRUE
  )

  conflicting <- with_coordinates(c(0, 2, 5), c(0, 3, 5))
  expect_error(
    .re_cov_output(conflicting, posterior),
    "conflicting canonical CAR time coordinates",
    fixed = TRUE
  )

  nonfinite_group <- result
  nonfinite_group$formula_design$random_effects[[1L]]$group_map[1L] <- Inf
  expect_error(
    .re_cov_output(nonfinite_group, posterior),
    "valid 'random_term$group_map'",
    fixed = TRUE
  )
})

test_that("structured covariance uses fitted column order for supplied data", {

  result <- .re_cov_formula(
    formula = ~ 1 + cs(f | id),
    data = .re_cov_structured_data(),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, 2),
    .re_cov_rho_sample(random_term, 0.25)
  ))
  new_data <- data.frame(
    id = factor(c("c", "c", "a", "b"), levels = c("c", "b", "a")),
    f = factor(c("c", "a", "b", "c"), levels = c("c", "b", "a")),
    time = c(5, 0, 2, 5)
  )
  out <- .re_cov_output(
    result,
    posterior,
    data = new_data,
    new_levels = "sample"
  )
  Z_new <- matrix(0, nrow = nrow(new_data), ncol = random_term$n_columns)
  colnames(Z_new) <- random_term$column_names
  Z_new[cbind(seq_len(nrow(new_data)), c(3L, 1L, 2L, 3L))] <- 1
  expected <- .re_cov_expand(
    Z_new,
    c(3L, 3L, 1L, 2L),
    .re_cov_compound_correlation(3L, 0.25) * 2^2
  )

  expect_equal(out$metadata$blocks$id$column_names, random_term$column_names)
  expect_equal(out$metadata$blocks$id$group_levels, c("a", "b", "c"))
  expect_equal(out$metadata$blocks$id$n_groups, 3L)
  expect_equal(out$metadata$blocks$id$fitted_n_groups, 2L)
  expect_equal(
    matrix(as.numeric(out$metadata$blocks$id$model_matrix), nrow = nrow(Z_new)),
    matrix(as.numeric(Z_new), nrow = nrow(Z_new)),
    tolerance = 1e-12
  )
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
})

test_that("point SD and scalar-rho priors are materialized on the right scale", {

  df <- .re_cov_structured_data()
  posterior <- .re_cov_posterior(c(unrelated = 1))

  fisher_result <- .re_cov_formula(
    formula = ~ 1 + cs(f | id),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 2)),
        rho = prior("point", list(location = 0.2))
      )
    )
  )
  fisher_term <- .re_cov_term(fisher_result, "id")
  fisher_out <- .re_cov_output(fisher_result, posterior)
  fisher_expected <- .re_cov_expand(
    fisher_term$model_matrix,
    fisher_term$group_map,
    .re_cov_compound_correlation(3L, tanh(0.2)) * 2^2
  )
  expect_equal(.re_cov_first(fisher_out), fisher_expected, tolerance = 1e-12)

  raw_result <- .re_cov_formula(
    formula = ~ 1 + cs(f | id),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 2)),
        covariance = random_covariance(
          rho = prior("point", list(location = 0.25)),
          rho_scale = "rho"
        )
      )
    )
  )
  raw_term <- .re_cov_term(raw_result, "id")
  raw_out <- .re_cov_output(raw_result, posterior)
  raw_expected <- .re_cov_expand(
    raw_term$model_matrix,
    raw_term$group_map,
    .re_cov_compound_correlation(3L, 0.25) * 2^2
  )
  expect_equal(.re_cov_first(raw_out), raw_expected, tolerance = 1e-12)
})

test_that("sparse structured covariance expands row pairs without global matrices", {

  K <- 678L
  row_column <- c(1L, 50L, 678L, 4L, 300L, 50L, 2L, 678L)
  df <- data.frame(
    id = factor(rep(c("a", "b"), each = 4L)),
    index = factor(
      paste0("level_", row_column),
      levels = paste0("level_", seq_len(K))
    )
  )
  result <- .re_cov_formula(
    formula = ~ 1 + cs(index | id),
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 2)),
        rho = prior("point", list(location = 0.2))
      )
    )
  )
  term <- .re_cov_term(result, "id")
  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )

  out <- .re_cov_output(
    result,
    .re_cov_posterior(c(unrelated = 1))
  )
  fitted_column <- max.col(term$model_matrix, ties.method = "first")
  expected <- outer(term$group_map, term$group_map, "==") *
    tanh(0.2)^outer(fitted_column, fitted_column, function(x, y) x != y) * 4
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
})

test_that("scalar rho posterior columns use canonical sample precedence", {

  result <- .re_cov_formula(
    formula = ~ 1 + cs(f | id),
    data = .re_cov_structured_data(),
    prior_random = prior_random(
      id = random_block(
        sd = .re_cov_sd_prior(),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(random_term, 2),
    stats::setNames(0.9, random_term$correlation$rho_name),
    .re_cov_rho_sample(random_term, 0.25)
  ))
  out <- .re_cov_output(result, posterior)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    .re_cov_compound_correlation(3L, 0.25) * 2^2
  )

  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
})

test_that("crossed independent random intercept blocks are summed", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("d1", "d2", "d1", "d2"), levels = c("d1", "d2"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      study = random_block(sd = .re_cov_sd_prior()),
      drug = random_block(sd = .re_cov_sd_prior())
    )
  )
  study <- .re_cov_term(result, "study")
  drug <- .re_cov_term(result, "drug")
  posterior <- .re_cov_posterior(c(
    .re_cov_sd_values(study, 2),
    .re_cov_sd_values(drug, 3)
  ))
  study_expected <- .re_cov_expand(
    study$model_matrix,
    study$group_map,
    matrix(2^2, 1L, 1L)
  )
  drug_expected <- .re_cov_expand(
    drug$model_matrix,
    drug$group_map,
    matrix(3^2, 1L, 1L)
  )

  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)
  expect_equal(.re_cov_first(out), study_expected + drug_expected)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out))
  expect_equal(diagonal$metadata$included_blocks, c("study", "drug"))
  expect_equal(out$metadata$included_blocks, c("study", "drug"))
  expect_equal(out$metadata$skipped_blocks$block_name, character())
  expect_true(out$metadata$dense)
  expect_false(out$metadata$potentially_expensive)

  study_only <- .re_cov_output(result, posterior, blocks = "study")
  expect_equal(.re_cov_first(study_only), study_expected)
  expect_equal(study_only$metadata$included_blocks, "study")
  expect_equal(study_only$metadata$skipped_blocks$block_name, "drug")
  expect_equal(study_only$metadata$skipped_blocks$reason, "not requested")
  expect_error(
    .re_cov_output(result, posterior, blocks = c("study", "study")),
    "'blocks' must be unique",
    fixed = TRUE
  )
  expect_error(
    random_effects_marginal_vcov(
      result$formula_design,
      posterior_samples = posterior,
      prior_list = result$prior_list,
      unused = TRUE
    ),
    "Unused argument(s): unused.",
    fixed = TRUE
  )

  broken_result <- result
  broken_result$formula_design$random_effects[[2L]]$model_matrix <-
    broken_result$formula_design$random_effects[[2L]]$model_matrix[1:3, , drop = FALSE]
  broken_result$formula_design$random_effects[[2L]]$group_map <-
    broken_result$formula_design$random_effects[[2L]]$group_map[1:3]
  expect_error(
    .re_cov_output(broken_result, posterior),
    "produced dimensions",
    fixed = TRUE
  )

  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "prior_list") <- result$prior_list
  extracted <- random_effects_marginal_vcov(fit, parameter = "mu")
  expect_equal(extracted$samples, out$samples)
})

test_that("nested random-effect blocks are expanded and summed", {

  df <- data.frame(
    district = factor(c("d1", "d1", "d2", "d2"), levels = c("d1", "d2")),
    school = factor(c("s1", "s2", "s3", "s4"),
                    levels = c("s1", "s2", "s3", "s4"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | district / school, covariance = "diag"),
    data = df,
    prior_random = prior_random(sd = .re_cov_sd_prior())
  )
  random_effects <- result$formula_design$random_effects
  expect_length(random_effects, 2L)

  posterior_values <- numeric()
  expected <- matrix(0, nrow(df), nrow(df))
  for(block_i in seq_along(random_effects)){
    random_term <- random_effects[[block_i]]
    sd <- block_i + 1
    posterior_values <- c(
      posterior_values,
      .re_cov_sd_values(random_term, sd)
    )
    expected <- expected + .re_cov_expand(
      random_term$model_matrix,
      random_term$group_map,
      matrix(sd^2, 1L, 1L)
    )
  }

  out <- .re_cov_output(result, .re_cov_posterior(posterior_values))
  expect_equal(.re_cov_first(out), expected)

  posterior <- .re_cov_posterior_draws(
    posterior_values,
    posterior_values / 2
  )
  dense <- .re_cov_output(result, posterior)
  factors <- random_effects_marginal_factor_states(
    result$formula_design,
    posterior_samples = posterior,
    prior_list        = result$prior_list,
    blocks            = vapply(random_effects, `[[`, character(1), "block_name"),
    row_blocks        = list(seq_len(nrow(df)))
  )

  expect_s3_class(
    factors,
    "BayesTools_random_effects_marginal_factor_states"
  )
  expect_identical(factors$row_blocks, list(seq_len(nrow(df))))
  for(draw in seq_len(nrow(posterior))){
    reconstructed <- matrix(0, nrow(df), nrow(df))
    for(block in seq_along(factors$factor_plans)){
      plan  <- factors$factor_plans[[block]]
      state <- factors$factor_states[[draw]][[block]]
      G     <- tcrossprod(state$coefficient_factor)
      reconstructed <- reconstructed +
        outer(plan$group_map, plan$group_map, "==") *
        tcrossprod(plan$model_matrix %*% G, plan$model_matrix)
    }
    expect_equal(
      unname(reconstructed),
      unname(dense$samples[draw, , ]),
      tolerance = 1e-12
    )
  }
})

test_that("row-varying direct SD sources weight rows inside groups", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd_source = random_sd_source("tau", shape = "row"))
    )
  )
  random_term <- .re_cov_term(result, "id")
  tau_1 <- c(1, 2, 3, 4)
  tau_2 <- c(2, 1, 4, 3)
  posterior <- .re_cov_posterior_draws(
    stats::setNames(tau_1, paste0("tau[", 1:4, "]")),
    stats::setNames(tau_2, paste0("tau[", 1:4, "]"))
  )
  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)
  expected_1 <- .re_cov_expand(
    matrix(tau_1, ncol = 1L),
    random_term$group_map,
    matrix(1, 1L, 1L)
  )
  expected_2 <- .re_cov_expand(
    matrix(tau_2, ncol = 1L),
    random_term$group_map,
    matrix(1, 1L, 1L)
  )

  expect_true(out$metadata$blocks$id$row_varying_sd)
  expect_equal(out$metadata$n_draws, 2L)
  expect_equal(unname(out$samples[1, , ]), expected_1)
  expect_equal(unname(out$samples[2, , ]), expected_2)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out))
  expect_error(
    .re_cov_output(result, posterior, data = df),
    "requires an explicit 'fitted_rows' mapping",
    fixed = TRUE
  )
  reordered <- .re_cov_output(
    result,
    posterior,
    data = df[c(4, 2, 2), , drop = FALSE],
    fitted_rows = c(4, 2, 2)
  )
  expected_reordered_1 <- .re_cov_expand(
    matrix(tau_1[c(4, 2, 2)], ncol = 1L),
    c(2L, 1L, 1L),
    matrix(1, 1L, 1L)
  )
  expected_reordered_2 <- .re_cov_expand(
    matrix(tau_2[c(4, 2, 2)], ncol = 1L),
    c(2L, 1L, 1L),
    matrix(1, 1L, 1L)
  )
  expect_equal(unname(reordered$samples[1, , ]), expected_reordered_1)
  expect_equal(unname(reordered$samples[2, , ]), expected_reordered_2)
  expect_error(
    .re_cov_output(
      result,
      posterior,
      data = df[1, , drop = FALSE],
      fitted_rows = 5L
    ),
    "equal or lower than 4",
    fixed = TRUE
  )
})

test_that("row-varying SD sources with values functions support supplied data", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    tau_scale = c(1, 2, 3, 4)
  )
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$tau_total * data$tau_scale[seq_len(n_rows)]
    }
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd_source = random_sd_source(tau_source))
    )
  )
  posterior <- .re_cov_posterior(c(tau_total = 2))
  new_data <- data.frame(
    id = factor(c("a", "a", "b"), levels = c("a", "b")),
    tau_scale = c(5, 7, 11)
  )
  out <- .re_cov_output(result, posterior, data = new_data)
  tau <- 2 * new_data$tau_scale
  expected <- .re_cov_expand(
    matrix(tau, ncol = 1L),
    c(1L, 1L, 2L),
    matrix(1, 1L, 1L)
  )

  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
})

test_that("row-varying SD-component allocation weights columns and rows", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    x = c(0, 1, 2, 3)
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      random_variance_allocation(
        terms = "id",
        target = "sd_component",
        scale = "total_variance",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  random_term <- .re_cov_term(result, "id")
  tau <- c(2, 3, 5, 7)
  weights <- c(0.25, 0.75)
  weight_name <- random_term$sd_binding$allocations[[1L]]$weight_name
  posterior <- .re_cov_posterior(c(
    stats::setNames(tau, paste0("tau[", seq_along(tau), "]")),
    stats::setNames(weights, paste0(weight_name, "[", seq_along(weights), "]"))
  ))
  out <- .re_cov_output(result, posterior)
  diagonal <- .re_cov_output(result, posterior, diagonal_only = TRUE)

  weighted_design <- random_term$model_matrix *
    matrix(tau, nrow = nrow(random_term$model_matrix),
           ncol = ncol(random_term$model_matrix)) *
    matrix(sqrt(weights), nrow = nrow(random_term$model_matrix),
           ncol = ncol(random_term$model_matrix), byrow = TRUE)
  expected <- .re_cov_expand(
    weighted_design,
    random_term$group_map,
    diag(ncol(random_term$model_matrix))
  )

  expect_true(out$metadata$blocks$id$row_varying_sd)
  expect_equal(.re_cov_first(out), expected, tolerance = 1e-12)
  expect_equal(unname(diagonal$samples), .re_cov_dense_diagonal(out),
               tolerance = 1e-12)

  eta_name <- BayesTools:::.JAGS_prior_dirichlet_eta_name(weight_name)
  eta_posterior <- .re_cov_posterior(c(
    stats::setNames(tau, paste0("tau[", seq_along(tau), "]")),
    stats::setNames(c(1, 3), paste0(eta_name, "[", seq_along(weights), "]"))
  ))
  eta_out <- .re_cov_output(result, eta_posterior)
  expect_equal(.re_cov_first(eta_out), expected, tolerance = 1e-12)

  eta_draws <- .re_cov_posterior_draws(
    c(
      stats::setNames(tau, paste0("tau[", seq_along(tau), "]")),
      stats::setNames(c(1, 3), paste0(eta_name, "[", seq_along(weights), "]"))
    ),
    c(
      stats::setNames(tau, paste0("tau[", seq_along(tau), "]")),
      stats::setNames(c(3, 1), paste0(eta_name, "[", seq_along(weights), "]"))
    )
  )
  dense <- .re_cov_output(result, eta_draws)
  factors <- random_effects_marginal_factor_states(
    result$formula_design,
    posterior_samples = eta_draws,
    prior_list        = result$prior_list,
    blocks            = "id",
    row_blocks        = split(seq_len(nrow(df)), random_term$group_map)
  )
  plan <- factors$factor_plans[[1L]]
  for(draw in seq_len(nrow(eta_draws))){
    state <- factors$factor_states[[draw]][[1L]]
    Z     <- plan$model_matrix * state$row_scale
    G     <- tcrossprod(state$coefficient_factor)
    reconstructed <- outer(plan$group_map, plan$group_map, "==") *
      tcrossprod(Z %*% G, Z)
    expect_equal(
      unname(reconstructed),
      unname(dense$samples[draw, , ]),
      tolerance = 1e-12
    )
  }

  invalid_eta <- eta_draws[1L, , drop = FALSE]
  invalid_eta[, paste0(eta_name, "[", seq_along(weights), "]")] <- 0
  expect_error(
    random_effects_marginal_factor_states(
      result$formula_design,
      posterior_samples = invalid_eta,
      prior_list        = result$prior_list,
      blocks            = "id",
      row_blocks        = split(seq_len(nrow(df)), random_term$group_map)
    ),
    "finite, non-negative, and have a positive sum",
    fixed = TRUE
  )
})

test_that("marginalized blocks keep usable covariance metadata", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      study = random_block(sd = .re_cov_sd_prior())
    ),
    random_effects_compile = random_effects_compile(marginalized = "study")
  )
  random_term <- .re_cov_term(result, "study")
  posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, 2))
  out <- .re_cov_output(result, posterior)
  expected <- .re_cov_expand(
    random_term$model_matrix,
    random_term$group_map,
    matrix(2^2, 1L, 1L)
  )

  expect_equal(random_term$compile_mode, "marginalized")
  expect_equal(out$metadata$blocks$study$compile_mode, "marginalized")
  expect_equal(.re_cov_first(out), expected)
})

test_that("diagonal-only covariance scales with draws times rows", {

  n_draws <- 40000L
  n_rows  <- 82L
  df <- data.frame(
    id = factor(rep("a", n_rows))
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- matrix(
    2,
    nrow = n_draws,
    ncol = 1L,
    dimnames = list(NULL, random_term$sd_parameter_names)
  )

  out <- .re_cov_output(result, posterior, diagonal_only = TRUE)

  expect_equal(dim(out$samples), c(n_draws, n_rows))
  expect_true(all(out$samples == 4))
  expect_equal(out$metadata$sample_entries, as.numeric(n_draws) * n_rows)
  expect_lt(out$metadata$estimated_size_bytes, 40 * 1024^2)
  expect_error(
    .re_cov_output(result, posterior[1L, , drop = FALSE],
                   diagonal_only = NA),
    "'diagonal_only'",
    fixed = TRUE
  )
})

test_that("covariance defaults remain identical and diagonal singleton is stable", {

  df <- data.frame(id = factor("a"))
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, 2))

  default <- .re_cov_output(result, posterior)
  explicit <- .re_cov_output(
    result,
    posterior,
    diagonal_only = FALSE
  )
  diagonal <- .re_cov_output(
    result,
    posterior,
    diagonal_only = TRUE
  )

  expect_identical(explicit, default)
  expect_named(
    default$metadata,
    c(
      "parameter", "n_draws", "n_rows", "row_names", "row_order",
      "data_source", "dense", "dense_entries", "estimated_size_bytes",
      "potentially_expensive", "included_blocks", "skipped_blocks",
      "structures", "blocks"
    ),
    ignore.order = FALSE
  )
  expect_named(
    default$metadata$blocks$id,
    c(
      "block_name", "grouping", "structure", "compile_mode", "n_groups",
      "fitted_n_groups", "n_columns", "n_rows", "row_names", "row_order",
      "group_levels", "group_map", "new_levels", "new_group_levels",
      "new_level_rows", "column_names", "model_matrix",
      "sd_parameter_names", "row_varying_sd", "correlation_type",
      "group_covariance", "included", "skipped", "dense", "dense_entries"
    ),
    ignore.order = FALSE
  )
  expect_equal(dim(diagonal$samples), c(1L, 1L))
  expect_equal(unname(diagonal$samples), matrix(4, 1L, 1L))
  expect_identical(diagonal$metadata$sample_dim, c(1L, 1L))
  expect_equal(diagonal$metadata$equivalent_dense_entries, 1)
})

test_that("missing SD samples fail with block and structure context", {

  df <- data.frame(
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  result <- .re_cov_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    data = df,
    prior_random = prior_random(
      id = random_block(sd = .re_cov_sd_prior())
    )
  )
  random_term <- .re_cov_term(result, "id")
  posterior <- .re_cov_posterior(c(unrelated = 1))

  expect_error(
    .re_cov_output(result, posterior),
    "block 'id' with structure 'diag' cannot resolve SD draws",
    fixed = TRUE
  )

  duplicate_posterior <- .re_cov_posterior(.re_cov_sd_values(random_term, 2))
  duplicate_posterior <- cbind(
    duplicate_posterior,
    duplicate_posterior[, 1L, drop = FALSE]
  )
  colnames(duplicate_posterior)[2L] <- colnames(duplicate_posterior)[1L]
  expect_error(
    .re_cov_output(result, duplicate_posterior),
    "column names must be unique",
    fixed = TRUE
  )
})
