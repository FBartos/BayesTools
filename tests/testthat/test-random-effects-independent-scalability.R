skip_if_not_test_profile("unit")

.independent_backend_sd_prior <- function(){
  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.independent_backend_fixture <- function(structure = "diag"){
  data <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula <- stats::as.formula(paste0("~ 1 + x + ", structure, "(1 + x | id)"))
  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = .independent_backend_sd_prior(),
        monitor = random_monitor(
          latent = TRUE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  )
  sd_names <- unique(random_term$sd_parameter_names)
  posterior <- matrix(
    0,
    nrow = 2L,
    ncol = length(sd_names) + length(z_names),
    dimnames = list(NULL, c(sd_names, as.vector(z_names)))
  )
  if(identical(structure, "id")){
    posterior[, sd_names] <- c(2, 4)
  }else{
    posterior[, sd_names] <- rbind(c(2, 3), c(4, 5))
  }
  posterior[, as.vector(z_names)] <- rbind(
    c(1, 2, 3, 4),
    c(-1, 0.5, 2, -2)
  )

  list(
    data = data,
    formula = formula,
    result = result,
    random_term = random_term,
    posterior = posterior,
    z_names = z_names
  )
}

.independent_backend_existing_oracle <- function(fixture){
  random_term <- fixture$random_term
  posterior <- fixture$posterior
  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = posterior,
    prior_list = fixture$result$prior_list
  )
  out <- matrix(0, nrow(random_term$model_matrix), nrow(posterior))
  for(draw in seq_len(nrow(posterior))){
    z <- matrix(
      posterior[draw, as.vector(fixture$z_names)],
      nrow = random_term$n_groups,
      ncol = random_term$n_columns
    )
    effects <- sweep(z, 2L, sd_draws[draw, ], "*")
    out[, draw] <- rowSums(
      random_term$model_matrix *
        effects[random_term$group_map, , drop = FALSE]
    )
  }
  out
}

.independent_backend_covariance_oracle <- function(model_matrix, group_map,
                                                    column_sd, row_sd = NULL){
  out <- matrix(0, nrow(model_matrix), nrow(model_matrix))
  for(rows in split(seq_len(nrow(model_matrix)), group_map)){
    Z <- sweep(model_matrix[rows, , drop = FALSE], 2L, column_sd, "*")
    if(!is.null(row_sd)){
      Z <- Z * row_sd[rows]
    }
    out[rows, rows] <- tcrossprod(Z)
  }
  out
}

.independent_backend_mock_dense_helpers <- function(){
  testthat::local_mocked_bindings(
    .bt_random_effect_cholesky_draws = function(...){
      stop("dense Cholesky helper called", call. = FALSE)
    },
    .bt_random_effect_marginal_covariance_correlation_draws = function(...){
      stop("dense correlation helper called", call. = FALSE)
    },
    .bt_random_effect_mvn_group_draws_from_factor = function(...){
      stop("dense factor helper called", call. = FALSE)
    },
    .bt_JAGS_marglik_random_effect_cholesky = function(...){
      stop("dense bridge Cholesky helper called", call. = FALSE)
    },
    .package = "BayesTools"
  )
}

test_that("ID and DIAG reconstruction streams exact independent columns", {
  .independent_backend_mock_dense_helpers()

  for(structure in c("diag", "id")){
    fixture <- .independent_backend_fixture(structure)
    expected <- .independent_backend_existing_oracle(fixture)
    actual <- .bt_try_random_effect_contribution_from_latent(
      random_term = fixture$random_term,
      model_matrix = fixture$random_term$model_matrix,
      group_map = fixture$random_term$group_map,
      posterior = fixture$posterior,
      prior_list = fixture$result$prior_list
    )
    expect_equal(unname(actual), unname(expected), tolerance = 1e-12,
                 info = structure)

    bridge_value <- .bt_JAGS_marglik_random_effect_value(
      samples = fixture$posterior[1L, ],
      random_term = fixture$random_term,
      prior_list = fixture$result$prior_list
    )
    expect_equal(bridge_value, expected[, 1L], tolerance = 1e-12,
                 info = structure)
  }
})

test_that("latent reconstruction rejects duplicate posterior coordinates", {

  fixture <- .independent_backend_fixture("diag")
  sd_name <- fixture$random_term$sd_parameter_names[[1L]]
  duplicate <- fixture$posterior[, sd_name, drop = FALSE] + 100
  posterior <- cbind(duplicate, fixture$posterior)
  colnames(posterior)[1L] <- sd_name

  expect_error(
    .bt_try_random_effect_contribution_from_latent(
      random_term = fixture$random_term,
      model_matrix = fixture$random_term$model_matrix,
      group_map = fixture$random_term$group_map,
      posterior = posterior,
      prior_list = fixture$result$prior_list
    ),
    "'posterior_samples' column names must be unique",
    fixed = TRUE
  )
})

test_that("point-prior random scales support zero-dimensional bridge rows", {

  data <- data.frame(
    id = factor(c("a", "a", "b", "b"))
  )
  point_result <- JAGS_formula(
    formula = ~ 1 + (1 | id),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("point", list(location = 0))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 0.2))
      )
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  random_term <- point_result$formula_design$random_effects[[1L]]
  posterior <- matrix(numeric(), nrow = 1L, ncol = 0L)

  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = posterior,
    prior_list = point_result$prior_list
  )
  expect_equal(sd_draws, matrix(0.2, nrow = 1L, ncol = 1L))

  context <- .bt_JAGS_bridge_context_random_block(
    samples = numeric(),
    random_term = random_term,
    prior_list = point_result$prior_list,
    formula_prior_parameters = list(),
    data = data,
    parameters = list()
  )
  expect_equal(unname(context$scale$column_sd), 0.2)

  sampled_prior_list <- point_result$prior_list
  sampled_prior_list[[random_term$sd_parameter_names[[1L]]]] <-
    .independent_backend_sd_prior()
  expect_null(.bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = posterior,
    prior_list = sampled_prior_list
  ))
  expect_error(
    .bt_random_effect_marginal_covariance_validate_posterior(posterior),
    "must have non-empty column names",
    fixed = TRUE
  )
})

test_that("independent new-level sampling avoids identity and eigen matrices", {
  .independent_backend_mock_dense_helpers()
  fixture <- .independent_backend_fixture("diag")
  random_term <- fixture$random_term
  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = fixture$posterior,
    prior_list = fixture$result$prior_list
  )
  group_map <- c(3L, 3L, 4L, 4L)

  set.seed(917)
  expected <- matrix(0, nrow(random_term$model_matrix), nrow(fixture$posterior))
  for(draw in seq_len(nrow(fixture$posterior))){
    effects <- matrix(
      stats::rnorm(2L * random_term$n_columns),
      nrow = 2L,
      ncol = random_term$n_columns
    )
    effects <- sweep(effects, 2L, sd_draws[draw, ], "*")
    expected[, draw] <- rowSums(
      random_term$model_matrix * effects[c(1L, 1L, 2L, 2L), , drop = FALSE]
    )
  }
  set.seed(917)
  actual <- .bt_random_effect_group_contribution_sample(
    random_term = random_term,
    model_matrix = random_term$model_matrix,
    group_map = group_map,
    posterior = fixture$posterior,
    prior_list = fixture$result$prior_list,
    source_data = fixture$data
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("independent marginal covariance and bridge context stay compact", {
  .independent_backend_mock_dense_helpers()
  fixture <- .independent_backend_fixture("diag")
  random_term <- fixture$random_term
  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = fixture$posterior,
    prior_list = fixture$result$prior_list
  )
  covariance <- random_effects_marginal_vcov(
    fixture$result$formula_design,
    posterior_samples = fixture$posterior,
    prior_list = fixture$result$prior_list
  )
  for(draw in seq_len(nrow(fixture$posterior))){
    expect_equal(
      unname(covariance$samples[draw, , ]),
      .independent_backend_covariance_oracle(
        model_matrix = random_term$model_matrix,
        group_map = random_term$group_map,
        column_sd = sd_draws[draw, ]
      ),
      tolerance = 1e-12
    )
  }

  context <- .bt_JAGS_bridge_context_random_block(
    samples = fixture$posterior[1L, ],
    random_term = random_term,
    prior_list = fixture$result$prior_list,
    formula_prior_parameters = list(),
    data = fixture$data,
    parameters = as.list(fixture$posterior[1L, ])
  )
  expect_null(context$correlation$cholesky)
  expect_null(context$correlation$matrix)
  expect_null(context$covariance)
  expect_identical(context$correlation$structure, "diag")
  expect_identical(context$correlation$column_coordinates, 1:2)
  expect_equal(unname(context$scale$column_sd), sd_draws[1L, ])
})

test_that("row-indexed independent reconstruction streams scalar and column allocations", {
  .independent_backend_mock_dense_helpers()
  data <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  scalar <- JAGS_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = fixed_priors,
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source("tau", shape = "row"),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE)
      )
    )
  )
  component <- JAGS_formula(
    formula = ~ 1 + random(1 + x | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(name = "allocation",
        terms = "id",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  tau <- c(2, 4, 6, 8)
  z <- c(1, 2, 3, 4)
  for(specification in list(scalar = scalar, component = component)){
    random_term <- specification$formula_design$random_effects[[1L]]
    z_names <- as.vector(.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = random_term$n_groups,
      n_columns = random_term$n_columns
    ))
    values <- c(stats::setNames(tau, paste0("tau[", 1:4, "]")),
                stats::setNames(z, z_names))
    if(identical(specification, component)){
      weight_name <- random_term$sd_binding$allocations[[1L]]$weight_name
      values <- c(values, stats::setNames(
        c(0.25, 0.75),
        paste0(weight_name, "[", 1:2, "]")
      ))
      column_scale <- sqrt(c(0.25, 0.75))
    }else{
      column_scale <- c(1, 1)
    }
    posterior <- matrix(values, nrow = 1L,
                        dimnames = list(NULL, names(values)))
    effects <- matrix(z, nrow = 2L, ncol = 2L) *
      matrix(column_scale, nrow = 2L, ncol = 2L, byrow = TRUE)
    expected <- tau * rowSums(
      random_term$model_matrix *
        effects[random_term$group_map, , drop = FALSE]
    )

    reconstructed <- .bt_random_effect_row_indexed_contribution_from_latent(
      random_term = random_term,
      model_matrix = random_term$model_matrix,
      group_map = random_term$group_map,
      posterior = posterior,
      prior_list = specification$prior_list
    )
    marglik <- .bt_JAGS_marglik_random_effect_value(
      samples = values,
      random_term = random_term,
      prior_list = specification$prior_list
    )
    expect_equal(unname(drop(reconstructed)), unname(expected),
                 tolerance = 1e-12)
    expect_equal(unname(marglik), unname(expected), tolerance = 1e-12)

    covariance <- .bt_random_effect_marginal_covariance_row_indexed_block(
      random_term = random_term,
      model_matrix = random_term$model_matrix,
      group_map = random_term$group_map,
      source_data = data,
      prediction_rows = seq_len(nrow(data)),
      posterior = posterior,
      prior_list = specification$prior_list
    )
    expected_covariance <- .independent_backend_covariance_oracle(
      model_matrix = random_term$model_matrix,
      group_map = random_term$group_map,
      column_sd = column_scale,
      row_sd = tau
    )
    expect_equal(unname(covariance$samples[1L, , ]), expected_covariance,
                 tolerance = 1e-12)

    new_group_map <- c(3L, 3L, 4L, 4L)
    set.seed(193)
    effects <- matrix(stats::rnorm(4L), nrow = 2L, ncol = 2L) *
      matrix(column_scale, nrow = 2L, ncol = 2L, byrow = TRUE)
    expected_new <- tau * rowSums(
      random_term$model_matrix *
        effects[c(1L, 1L, 2L, 2L), , drop = FALSE]
    )
    set.seed(193)
    sampled <- .bt_random_effect_group_contribution_sample_row_indexed(
      random_term = random_term,
      model_matrix = random_term$model_matrix,
      group_map = new_group_map,
      rows = seq_len(nrow(random_term$model_matrix)),
      posterior = posterior,
      prior_list = specification$prior_list,
      source_data = data
    )
    expect_equal(unname(drop(sampled)), unname(expected_new),
                 tolerance = 1e-12)

    missing_latent <- posterior[, !colnames(posterior) %in% z_names,
                                drop = FALSE]
    expect_error(
      .bt_random_effect_row_indexed_contribution_from_latent(
        random_term = random_term,
        model_matrix = random_term$model_matrix,
        group_map = random_term$group_map,
        posterior = missing_latent,
        prior_list = specification$prior_list
      ),
      "cannot be reconstructed from the posterior samples",
      fixed = TRUE
    )
  }
})

test_that("large independent blocks remain linear in coefficient count", {
  .independent_backend_mock_dense_helpers()
  n_columns <- 4096L
  n_groups  <- 2L
  n_draws   <- 2L
  model_matrix <- matrix(0, nrow = 4L, ncol = n_columns)
  model_matrix[cbind(1:4, c(1L, n_columns, 2L, n_columns - 1L))] <-
    c(1, -2, 3, 0.5)
  random_term <- list(
    structure = "diag",
    block_name = "large",
    parameter_stem = "mu__xREx__large",
    group_levels = c("a", "b"),
    group_map = c(1L, 1L, 2L, 2L),
    n_groups = n_groups,
    n_columns = n_columns,
    model_matrix = model_matrix,
    sd_parameter_names = paste0("sd[", seq_len(n_columns), "]")
  )
  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  posterior <- matrix(
    0,
    nrow = n_draws,
    ncol = n_columns + length(z_names),
    dimnames = list(NULL, c(random_term$sd_parameter_names,
                           as.vector(z_names)))
  )
  posterior[, random_term$sd_parameter_names] <- 1
  posterior[, as.vector(z_names)] <- 0.25

  contribution <- .bt_try_random_effect_contribution_from_latent(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = random_term$group_map,
    posterior = posterior,
    prior_list = list()
  )
  covariance <- .bt_random_effect_marginal_covariance_block_samples(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = random_term$group_map,
    source_data = NULL,
    prediction_rows = NULL,
    posterior = posterior,
    prior_list = list()
  )
  set.seed(44)
  sampled <- .bt_random_effect_group_contribution_sample(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = c(3L, 3L, 4L, 4L),
    posterior = posterior,
    prior_list = list(),
    source_data = NULL
  )
  marglik <- .bt_JAGS_marglik_random_effect_value(
    samples = posterior[1L, ],
    random_term = random_term,
    prior_list = list()
  )

  expect_identical(dim(contribution), c(4L, n_draws))
  expect_identical(dim(covariance$samples), c(n_draws, 4L, 4L))
  expect_identical(dim(sampled), c(4L, n_draws))
  expect_length(marglik, 4L)
  expect_true(all(is.finite(c(contribution, covariance$samples,
                              sampled, marglik))))
})

test_that("prediction samples directly from authoritative factors", {

  factor <- matrix(c(1, 0.25, 0, sqrt(1 - 0.25^2)), nrow = 2L)
  set.seed(42)
  draws <- BayesTools:::.bt_random_effect_mvn_group_draws_from_factor(
    factor,
    n_groups = 5L
  )
  set.seed(42)
  latent <- matrix(stats::rnorm(10L), nrow = 5L, ncol = 2L)

  expect_identical(dim(draws), c(5L, 2L))
  expect_equal(draws, latent %*% t(factor), tolerance = 0)
  expect_null(attr(draws, "covariance_eigen_correction", exact = TRUE))

  expect_error(
    BayesTools:::.bt_random_effect_mvn_group_draws_from_factor(
      matrix(c(1, 0, 0, NA_real_), nrow = 2L),
      n_groups = 1L
    ),
    "finite numeric square matrix",
    fixed = TRUE
  )
})
