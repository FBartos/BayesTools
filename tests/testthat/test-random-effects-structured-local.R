skip_if_not_test_profile("unit")

test_that("structured local layout keeps only observed group-column cells", {

  model_matrix <- matrix(0, nrow = 8L, ncol = 5L)
  row_column   <- c(3L, 1L, 3L, 5L, 2L, 5L, 4L, 4L)
  group_map    <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L)
  model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)] <- 1

  layout <- BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = group_map,
    structure = "CS",
    exact_indicator = TRUE,
    parameter_stem = "mu__xREx__study"
  )

  expect_s3_class(layout, "BayesTools_random_effect_structured_local_layout")
  expect_equal(layout$global_n_columns, 5L)
  expect_equal(layout$n_groups, 3L)
  expect_equal(layout$n_local, 5L)
  expect_equal(layout$group_columns, list(`1` = c(1L, 3L), `2` = c(2L, 5L), `3` = 4L))
  expect_equal(layout$row_column, row_column)
  expect_equal(layout$row_local, c(2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L))
  expect_equal(
    layout$node_names,
    c(
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,3]",
      "mu__xREx__study_xRE_Zx[2,2]",
      "mu__xREx__study_xRE_Zx[2,5]",
      "mu__xREx__study_xRE_Zx[3,4]"
    )
  )
})

test_that("structured local latent names require unique character metadata", {

  model_matrix <- diag(3)
  layout <- BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = c(1L, 1L, 2L),
    structure = "cs",
    exact_indicator = TRUE,
    parameter_stem = "mu__xREx__id"
  )
  random_term <- list(
    block_name = "id",
    parameter_stem = "mu__xREx__id",
    latent_layout = layout
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = layout$n_groups,
      n_columns = layout$global_n_columns
    ),
    layout$node_names
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = 1L,
      n_columns = 1L
    ),
    layout$node_names[1L]
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = 1L,
      n_columns = 2L
    ),
    layout$node_names[1:2]
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = 2L,
      n_columns = 2L
    ),
    layout$node_names[1:2]
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = 0L,
      n_columns = layout$global_n_columns
    ),
    character()
  )
  expect_identical(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = layout$n_groups,
      n_columns = 0L
    ),
    character()
  )

  malformed_names <- list(
    non_character = seq_along(layout$node_names),
    missing = replace(layout$node_names, 1L, NA_character_),
    empty = replace(layout$node_names, 1L, ""),
    duplicate = replace(layout$node_names, 2L, layout$node_names[[1L]])
  )
  for(node_names in malformed_names){
    malformed_term <- random_term
    malformed_term$latent_layout$node_names <- node_names
    expect_error(
      BayesTools:::.bt_random_effect_latent_names(
        random_term = malformed_term,
        n_groups = layout$n_groups,
        n_columns = layout$global_n_columns
      ),
      paste0(
        "Random-effect local latent metadata for block 'id' must contain one ",
        "unique, non-missing, non-empty character node name per local latent cell."
      ),
      fixed = TRUE
    )
  }

  malformed_indices <- list(
    group_length = list(
      local_group = layout$local_group[-1L],
      local_column = layout$local_column
    ),
    column_length = list(
      local_group = layout$local_group,
      local_column = layout$local_column[-1L]
    ),
    zero_group = list(
      local_group = replace(layout$local_group, 1L, 0L),
      local_column = layout$local_column
    ),
    oversized_column = list(
      local_group = layout$local_group,
      local_column = replace(
        layout$local_column,
        1L,
        layout$global_n_columns + 1L
      )
    )
  )
  for(indices in malformed_indices){
    malformed_term <- random_term
    malformed_term$latent_layout$local_group <- indices$local_group
    malformed_term$latent_layout$local_column <- indices$local_column
    expect_error(
      BayesTools:::.bt_random_effect_latent_names(
        random_term = malformed_term,
        n_groups = layout$n_groups,
        n_columns = layout$global_n_columns
      ),
      paste0(
        "Random-effect local latent metadata for block 'id' must contain one ",
        "positive group and column index per local latent cell within the ",
        "stored dimensions."
      ),
      fixed = TRUE
    )
  }
})

test_that("structured local layout rejects non-one-hot and malformed mappings", {

  model_matrix <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))

  expect_null(BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = c(1L, 1L, 2L),
    structure = "cs"
  ))

  skipped <- BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = c(1L, 3L, 3L),
    structure = "cs",
    exact_indicator = TRUE,
    n_groups = 3L
  )
  expect_false(skipped$all_groups_observed)
  expect_equal(skipped$group_columns[[2L]], integer())
  expect_error(
    BayesTools:::.bt_random_effect_structured_local_layout(
      model_matrix = rbind(c(1, 1, 0), model_matrix[-1L, ]),
      group_map = c(1L, 1L, 2L),
      structure = "cs",
      exact_indicator = TRUE
    ),
    "exactly one non-zero",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_local_layout(
      model_matrix = rbind(c(2, 0, 0), model_matrix[-1L, ]),
      group_map = c(1L, 1L, 2L),
      structure = "cs",
      exact_indicator = TRUE
    ),
    "non-binary",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_local_layout(
      model_matrix = rbind(
        c(1 + .Machine$double.eps, 0, 0),
        model_matrix[-1L, ]
      ),
      group_map = c(1L, 1L, 2L),
      structure = "cs",
      exact_indicator = TRUE
    ),
    "non-binary",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_local_layout(
      model_matrix = model_matrix,
      group_map = c(1L, 1L),
      structure = "cs",
      exact_indicator = TRUE
    ),
    "one positive integer per model-matrix row",
    fixed = TRUE
  )
})

test_that("local CS and HCS Cholesky blocks equal global principal blocks", {

  K       <- 8L
  columns <- c(1L, 4L, 8L)
  for(structure in c("cs", "hcs")){
    for(rho in c(-0.99 / (K - 1L), -0.01, 0, 0.35, 0.9)){
      L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
        structure = structure,
        columns = columns,
        rho = rho,
        global_n_columns = K
      )
      global_R <- diag(1 - rho, K) + matrix(rho, nrow = K, ncol = K)
      expect_equal(tcrossprod(L), global_R[columns, columns], tolerance = 1e-12)
    }
  }

  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "cs",
      columns = columns,
      rho = -0.2,
      global_n_columns = K
    ),
    "global CS support for K = 8",
    fixed = TRUE
  )
})

test_that("local AR1 and HAR recurrences retain global index gaps", {

  K       <- 7L
  columns <- c(1L, 3L, 6L, 7L)
  for(structure in c("ar", "ar1", "har")){
    for(rho in c(-0.8, -0.2, 0, 0.4, 0.9)){
      L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
        structure = structure,
        columns = columns,
        rho = rho,
        global_n_columns = K
      )
      global_R <- rho^abs(outer(seq_len(K), seq_len(K), "-"))
      expect_equal(tcrossprod(L), global_R[columns, columns], tolerance = 1e-12)
      expect_equal(L[2L, 1L], rho^2, tolerance = 1e-12)
      expect_equal(L[3L, 1L], rho^5, tolerance = 1e-12)
    }
  }

  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "ar1",
      columns = c(1L, 2L),
      rho = -0.5,
      global_n_columns = 3L,
      column_coordinates = c(0, 0.5, 2)
    ),
    "integer ordered-index positions",
    fixed = TRUE
  )
})

test_that("local CAR recurrence retains irregular numeric time gaps", {

  coordinates <- c(0, 0.5, 2, 5)
  columns     <- c(1L, 3L, 4L)
  expect_equal(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = 2L,
      rho = 0.5,
      global_n_columns = length(coordinates),
      column_coordinates = coordinates
    ),
    matrix(1, nrow = 1L, ncol = 1L)
  )
  for(rho in c(0, 0.2, 0.75, 0.95)){
    L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = columns,
      rho = rho,
      global_n_columns = length(coordinates),
      column_coordinates = coordinates
    )
    global_R <- rho^abs(outer(coordinates, coordinates, "-"))
    expect_equal(tcrossprod(L), global_R[columns, columns], tolerance = 1e-12)
    expect_equal(L[2L, 1L], rho^2, tolerance = 1e-12)
    expect_equal(L[3L, 1L], rho^5, tolerance = 1e-12)
  }

  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = columns,
      rho = 0.5,
      global_n_columns = length(coordinates)
    ),
    "requires one numeric coordinate",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = 2L,
      rho = 0.5,
      global_n_columns = length(coordinates)
    ),
    "requires one numeric coordinate",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = 2L,
      rho = 0.5,
      global_n_columns = length(coordinates),
      column_coordinates = c(0, 0.5, 0.5, 5)
    ),
    "must be finite, unique, and match the global column count",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = "car",
      columns = rev(columns),
      rho = 0.5,
      global_n_columns = length(coordinates),
      column_coordinates = coordinates
    ),
    "ordered by increasing index coordinate",
    fixed = TRUE
  )
})

test_that("CAR transitions preserve representable boundary innovations", {

  rho_upper <- 1 - .Machine$double.eps / 2
  upper_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = rho_upper,
    left_coordinate = 0,
    right_coordinate = 0.5
  )
  expect_identical(upper_transition$log_phi, 0.5 * log(rho_upper))
  expect_identical(
    upper_transition$innovation_variance,
    -expm1(2 * upper_transition$log_phi)
  )
  expect_gt(upper_transition$innovation_variance, 0)

  upper_L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
    structure = "car",
    columns = 1:2,
    rho = rho_upper,
    global_n_columns = 2L,
    column_coordinates = c(0, 0.5)
  )
  expect_identical(
    upper_L[2L, 2L],
    sqrt(upper_transition$innovation_variance)
  )
  expect_gt(upper_L[2L, 2L], 0)
  expect_identical(
    BayesTools:::.bt_random_effect_structured_subset_transform(
      structure = "car",
      columns = 1:2,
      latent = c(0, 1),
      rho = rho_upper,
      global_n_columns = 2L,
      column_coordinates = c(0, 0.5)
    ),
    c(0, sqrt(upper_transition$innovation_variance))
  )

  tiny_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0.5,
    left_coordinate = 0,
    right_coordinate = 1e-320
  )
  expect_gt(tiny_transition$innovation_variance, 0)
  expect_identical(
    tiny_transition$innovation_variance,
    -expm1(2 * tiny_transition$log_phi)
  )
  tiny_L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
    structure = "car",
    columns = 1:2,
    rho = 0.5,
    global_n_columns = 2L,
    column_coordinates = c(0, 1e-320)
  )
  expect_identical(
    tiny_L[2L, 2L],
    sqrt(tiny_transition$innovation_variance)
  )
  expect_gt(tiny_L[2L, 2L], 0)

  zero_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0,
    left_coordinate = -2,
    right_coordinate = 3
  )
  expect_identical(
    zero_transition,
    list(log_phi = -Inf, phi = 0, innovation_variance = 1)
  )
})

test_that("CAR transition operation order and failures are explicit", {

  smallest <- .Machine$double.xmin * .Machine$double.eps
  transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0.5,
    left_coordinate = 0,
    right_coordinate = smallest
  )
  expect_identical(transition$log_phi, smallest * log(0.5))
  expect_identical(
    transition$innovation_variance,
    -expm1(2 * transition$log_phi)
  )
  expect_gt(
    transition$innovation_variance,
    -expm1(2 * smallest * log(0.5))
  )

  error_message <- tryCatch(
    BayesTools:::.bt_random_effect_markov_transition(
      rho = 0.9,
      left_coordinate = 0,
      right_coordinate = smallest
    ),
    error = function(error) conditionMessage(error)
  )
  smallest_label <- format(
    smallest,
    digits = 17L,
    scientific = TRUE,
    trim = TRUE
  )
  rho_label <- format(
    0.9,
    digits = 17L,
    scientific = TRUE,
    trim = TRUE
  )
  expect_match(
    error_message,
    "from coordinate 0.0000000000000000e+00",
    fixed = TRUE
  )
  expect_match(error_message, paste0("to ", smallest_label), fixed = TRUE)
  expect_match(error_message, paste0("rho = ", rho_label), fixed = TRUE)
  expect_match(error_message, paste0("gap = ", smallest_label), fixed = TRUE)
  expect_match(
    error_message,
    "non-positive or non-finite innovation variance",
    fixed = TRUE
  )
  expect_match(
    error_message,
    "requested coordinate/time resolution is not representable",
    fixed = TRUE
  )
})

test_that("stable CAR recurrence is shared by reconstruction and prediction", {

  rho <- 1 - .Machine$double.eps / 2
  coordinates <- c(0, 0.5)
  data <- data.frame(
    time = coordinates,
    id = factor(c("old", "old"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  result <- JAGS_formula(
    formula = ~ 1 + car(time | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior("normal", list(mean = 0, sd = 0.5))
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  latent_names <- as.vector(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = 1L,
      n_columns = 2L
    )
  )
  posterior_values <- c(0, 1, rho, 0, 1)
  names(posterior_values) <- c(
    "mu_intercept",
    unique(random_term$sd_parameter_names),
    random_term$correlation$rho_name,
    latent_names
  )
  posterior <- matrix(
    posterior_values,
    nrow = 1L,
    dimnames = list(NULL, names(posterior_values))
  )
  transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = rho,
    left_coordinate = coordinates[1L],
    right_coordinate = coordinates[2L]
  )
  expect_identical(
    unname(BayesTools:::.bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior
    )),
    unname(rho)
  )

  reconstructed <- BayesTools:::.bt_random_effect_structured_dense_contribution(
    random_term = random_term,
    model_matrix = random_term$model_matrix,
    group_map = random_term$group_map,
    posterior = posterior,
    scale_draws = matrix(1, nrow = 1L, ncol = 2L)
  )
  expect_equal(
    unname(drop(reconstructed)),
    c(0, sqrt(transition$innovation_variance)),
    tolerance = 0
  )

  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  new_data <- data.frame(
    time = coordinates,
    id = factor(c("new", "new"), levels = c("old", "new"))
  )
  set.seed(742)
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    prior_list = result$prior_list,
    formula_target = "conditional",
    new_levels = "sample"
  )
  set.seed(742)
  latent <- stats::rnorm(2L)
  expected_prediction <- c(
    latent[1L],
    transition$phi * latent[1L] +
      sqrt(transition$innovation_variance) * latent[2L]
  )
  expect_equal(
    unname(drop(prediction)),
    expected_prediction,
    tolerance = 0
  )

  bridge <- BayesTools:::.bt_random_effect_markov_conditional_unit_draw(
    observed_columns = 1L,
    missing_columns = 2L,
    observed_unit = 0,
    rho = rho,
    column_coordinates = coordinates,
    innovations = 1
  )
  expect_equal(
    unname(bridge),
    sqrt(transition$innovation_variance),
    tolerance = 0
  )
})

test_that("layout Cholesky blocks follow group-specific principal subsets", {

  model_matrix <- matrix(0, nrow = 7L, ncol = 6L)
  row_column   <- c(1L, 4L, 4L, 2L, 5L, 6L, 2L)
  group_map    <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)
  model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)] <- 1
  layout <- BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = group_map,
    structure = "ar1",
    exact_indicator = TRUE
  )
  blocks <- BayesTools:::.bt_random_effect_structured_local_cholesky_blocks(
    layout = layout,
    rho = -0.4
  )
  global_R <- (-0.4)^abs(outer(seq_len(6L), seq_len(6L), "-"))

  expect_length(blocks, 3L)
  for(group in seq_along(blocks)){
    columns <- layout$group_columns[[group]]
    expect_equal(tcrossprod(blocks[[group]]), global_R[columns, columns, drop = FALSE],
                 tolerance = 1e-12)
  }
})

test_that("direct subset transforms equal explicit Cholesky products", {

  set.seed(4813)
  cases <- list(
    cs = list(columns = c(1L, 8L, 30L, 678L), rho = -0.001),
    hcs = list(columns = c(1L, 8L, 30L, 678L), rho = 0.6),
    ar1 = list(columns = c(1L, 8L, 30L, 678L), rho = -0.4),
    har = list(columns = c(1L, 8L, 30L, 678L), rho = 0.8),
    car = list(columns = c(1L, 8L, 30L, 678L), rho = 0.7)
  )
  coordinates <- cumsum(seq_len(678L) / 678L)

  for(structure in names(cases)){
    case <- cases[[structure]]
    latent <- stats::rnorm(length(case$columns))
    structure_coordinates <- if(identical(structure, "car")) {
      coordinates
    } else {
      NULL
    }
    actual <- BayesTools:::.bt_random_effect_structured_subset_transform(
      structure = structure,
      columns = case$columns,
      latent = latent,
      rho = case$rho,
      global_n_columns = 678L,
      column_coordinates = structure_coordinates
    )
    L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
      structure = structure,
      columns = case$columns,
      rho = case$rho,
      global_n_columns = 678L,
      column_coordinates = structure_coordinates
    )
    expect_equal(actual, as.vector(L %*% latent), tolerance = 1e-12,
                 info = structure)
  }
})

test_that("irrelevant formula scaling does not materialize structured matrices", {

  random_term <- list(
    block_name = "study",
    parameter_stem = "mu__xREx__study",
    structure = "ar1",
    n_columns = 678L,
    sd_parameter_names = rep("mu__xREx__study_sd", 678L)
  )
  posterior <- cbind(mu__xREx__study_rho = c(0.2, 0.4))
  formula_scale <- list(mu = list(
    mu_x = list(mean = 1, sd = 2)
  ))
  expect_identical(
    BayesTools:::.bt_random_effect_summary_complete_scaled_samples(
      random_term = random_term,
      model_samples = posterior,
      prior_list = list(),
      parameter = "mu",
      formula_scale = formula_scale
    ),
    posterior
  )
})

test_that("structured conditional recurrences match exact Gaussian conditioning", {

  K <- 9L
  observed_columns <- c(1L, 4L, 8L)
  missing_columns  <- c(7L, 2L, 9L, 5L, 3L, 6L)
  observed_unit    <- c(0.3, -1.1, 0.6)
  car_coordinates  <- c(0, 0.2, 1.1, 2, 4.5, 4.75, 7, 10, 12.5)
  cases <- list(
    cs = list(rho = -0.12, coordinates = seq_len(K)),
    hcs = list(rho = 0.7, coordinates = seq_len(K)),
    ar1 = list(rho = -0.65, coordinates = seq_len(K)),
    har = list(rho = -0.35, coordinates = seq_len(K)),
    car = list(rho = 0.72, coordinates = car_coordinates)
  )

  for(structure in names(cases)){
    rho <- cases[[structure]]$rho
    coordinates <- cases[[structure]]$coordinates
    correlation <- if(structure %in% c("cs", "hcs")){
      out <- matrix(rho, nrow = K, ncol = K)
      diag(out) <- 1
      out
    }else{
      rho^abs(outer(coordinates, coordinates, "-"))
    }
    R_oo <- correlation[observed_columns, observed_columns, drop = FALSE]
    R_mo <- correlation[missing_columns, observed_columns, drop = FALSE]
    R_mm <- correlation[missing_columns, missing_columns, drop = FALSE]
    expected_mean <- as.vector(R_mo %*% solve(R_oo, observed_unit))
    expected_covariance <- R_mm - R_mo %*% solve(R_oo, t(R_mo))

    conditional <- function(innovations){
      if(structure %in% c("cs", "hcs")){
        BayesTools:::.bt_random_effect_cs_conditional_unit_draw(
          observed_unit = observed_unit,
          rho = rho,
          n_missing = length(missing_columns),
          innovations = innovations
        )
      }else{
        BayesTools:::.bt_random_effect_markov_conditional_unit_draw(
          observed_columns = observed_columns,
          missing_columns = missing_columns,
          observed_unit = observed_unit,
          rho = rho,
          column_coordinates = coordinates,
          innovations = innovations
        )
      }
    }
    base <- conditional(rep(0, length(missing_columns)))
    transform <- vapply(seq_along(missing_columns), function(column){
      innovation <- numeric(length(missing_columns))
      innovation[column] <- 1
      conditional(innovation) - base
    }, numeric(length(missing_columns)))

    expect_equal(base, expected_mean, tolerance = 1e-11, info = structure)
    expect_equal(
      tcrossprod(transform),
      expected_covariance,
      tolerance = 1e-11,
      info = structure
    )
  }
})

test_that("CAR bridge variance is evaluated ratio-first", {

  coordinates <- c(0, 1e-320, 2e-320)
  left_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0.5,
    left_coordinate = coordinates[1L],
    right_coordinate = coordinates[2L]
  )
  right_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0.5,
    left_coordinate = coordinates[2L],
    right_coordinate = coordinates[3L]
  )
  span_transition <- BayesTools:::.bt_random_effect_markov_transition(
    rho = 0.5,
    left_coordinate = coordinates[1L],
    right_coordinate = coordinates[3L]
  )
  multiply_first <-
    left_transition$innovation_variance *
      right_transition$innovation_variance /
      span_transition$innovation_variance
  expected_variance <-
    left_transition$innovation_variance *
      (right_transition$innovation_variance /
         span_transition$innovation_variance)
  expect_identical(multiply_first, 0)
  expect_gt(expected_variance, 0)

  draw <- BayesTools:::.bt_random_effect_markov_conditional_unit_draw(
    observed_columns = c(1L, 3L),
    missing_columns = 2L,
    observed_unit = c(0, 0),
    rho = 0.5,
    column_coordinates = coordinates,
    innovations = 1
  )
  expect_equal(unname(draw), sqrt(expected_variance), tolerance = 0)
})

test_that("structured conditional wrapper preserves innovation and SD semantics", {

  observed_columns <- c(1L, 4L, 8L)
  missing_columns  <- c(7L, 2L, 9L, 5L, 3L, 6L)
  observed_unit    <- c(0.3, -1.1, 0.6)
  sd <- seq(0, 2, length.out = 9L)

  set.seed(7291)
  actual <- BayesTools:::.bt_random_effect_structured_local_conditional_draw(
    structure = "ar1",
    observed_columns = observed_columns,
    missing_columns = missing_columns,
    observed_unit = observed_unit,
    sd = sd,
    rho = -0.65,
    global_n_columns = 9L,
    column_coordinates = seq_len(9L)
  )
  set.seed(7291)
  innovations <- stats::rnorm(length(missing_columns))
  expected_unit <- BayesTools:::.bt_random_effect_markov_conditional_unit_draw(
    observed_columns = observed_columns,
    missing_columns = missing_columns,
    observed_unit = observed_unit,
    rho = -0.65,
    column_coordinates = seq_len(9L),
    innovations = innovations
  )

  expect_identical(names(actual), as.character(missing_columns))
  expect_equal(unname(actual), expected_unit * sd[missing_columns], tolerance = 0)
})

test_that("Markov conditional reconstruction scales with requested coordinates", {

  K <- 100000L
  missing_columns <- c(2L, 7L, 101L, 999L, 1001L, 9999L,
                       25001L, 49999L, 50001L, 75000L, 99998L, 99999L)
  draw <- BayesTools:::.bt_random_effect_structured_local_conditional_draw(
    structure = "ar1",
    observed_columns = c(1L, 50000L, K),
    missing_columns = missing_columns,
    observed_unit = c(-0.2, 0.5, 1.1),
    sd = rep(1, K),
    rho = 0.999,
    global_n_columns = K,
    column_coordinates = seq_len(K)
  )

  expect_length(draw, length(missing_columns))
  expect_true(all(is.finite(draw)))
})

test_that("zero-SD conditional coordinates remain exactly deterministic", {

  set.seed(914)
  draw <- BayesTools:::.bt_random_effect_structured_local_conditional_draw(
    structure = "ar1",
    observed_columns = c(1L, 3L),
    missing_columns = c(2L, 4L),
    observed_unit = c(0.3, -0.7),
    sd = c(1, 0, 1, 0),
    rho = 0.6,
    global_n_columns = 4L,
    column_coordinates = seq_len(4L)
  )
  expect_equal(unname(draw), c(0, 0), tolerance = 0)
})

test_that("group-local bridge context retains only compact descriptors", {

  K <- 113L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(17L)), length.out = K))
  )
  result <- JAGS_formula(
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        cor = prior("point", list(location = 0.2))
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]
  samples <- stats::setNames(rep(0, term$latent_layout$n_local),
                             term$latent_layout$node_names)
  context <- BayesTools:::.bt_JAGS_bridge_context_random_block(
    samples = samples,
    random_term = term,
    prior_list = result$prior_list,
    formula_prior_parameters = list(),
    data = data,
    parameters = list()
  )

  expect_null(context$correlation$matrix)
  expect_null(context$correlation$cholesky)
  expect_null(context$correlation$blocks)
  expect_null(context$covariance)
  expect_null(context$covariance_blocks)
  expect_identical(context$correlation$structure, "cs")
  expect_s3_class(context$latent, "BayesTools_group_local_latent")
  expect_length(context$latent$values, K)
})

test_that("unobserved grouping levels preserve the dense fitted-level contract", {

  K <- 40L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(
      rep(paste0("group_", seq_len(4L)), length.out = K),
      levels = paste0("group_", seq_len(5L))
    )
  )
  result <- JAGS_formula(
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]

  expect_null(term$latent_layout)
  expect_match(
    result$formula_syntax,
    "mu__xREx__id_xRE_Zx[i,j] ~ dnorm(0, 1)",
    fixed = TRUE
  )
  expect_equal(term$n_groups, 5L)
})

test_that("row-indexed structured SD uses compact group-local projection", {

  K <- 113L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(17L)), length.out = K))
  )
  result <- JAGS_formula(
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source("tau", shape = "row"),
        cor = prior("point", list(location = 0.2))
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]

  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_equal(term$latent_layout$n_local, K)
  expect_equal(
    names(result$data),
    c("N_mu", "mu__xREx__id_xRE_MAPx", "mu__xREx__id_xRE_COLx")
  )
  expect_false(grepl("mu__xREx__id_xRE_DATAx", result$formula_syntax,
                     fixed = TRUE))
  expect_false(grepl("inprod(", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__id_xRE_COEFx", result$formula_syntax,
                     fixed = TRUE))
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id[i] = tau[i] * 1 * ",
      "mu__xREx__id_xRE_UNIT_COEFx[mu__xREx__id_xRE_MAPx[i],",
      "mu__xREx__id_xRE_COLx[i]]"
    ),
    fixed = TRUE
  )
})

test_that("row-indexed structured column allocation uses compact scale lookup", {

  K <- 113L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(17L)), length.out = K))
  )
  result <- JAGS_formula(
    formula = ~ 1 + hcs(index | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "id",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = rep(1, K)))
      ),
      id = random_block(cor = prior("point", list(location = 0.2)))
    )
  )
  term <- result$formula_design$random_effects[[1L]]

  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_equal(term$sd_binding$application, "column")
  expect_equal(
    names(result$data),
    c("N_mu", "mu__xREx__id_xRE_MAPx", "mu__xREx__id_xRE_COLx")
  )
  expect_false(grepl("mu__xREx__id_xRE_DATAx", result$formula_syntax,
                     fixed = TRUE))
  expect_false(grepl("inprod(", result$formula_syntax, fixed = TRUE))
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_ROW_COL_SCALEx[1] <- ",
      "sqrt(mu__xRE_ALLOCx_allocation__weight[1])"
    ),
    fixed = TRUE
  )
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id[i] = tau[i] * ",
      "mu__xREx__id_xRE_ROW_COL_SCALEx[mu__xREx__id_xRE_COLx[i]] * ",
      "mu__xREx__id_xRE_UNIT_COEFx[mu__xREx__id_xRE_MAPx[i],",
      "mu__xREx__id_xRE_COLx[i]]"
    ),
    fixed = TRUE
  )
})

test_that("structured complexity accounts for centered quadratic transforms", {

  monitor <- random_monitor(
    latent = FALSE,
    coefficients = FALSE,
    correlation = FALSE
  )
  estimates <- BayesTools:::.bt_random_effect_dense_complexity(
    structure = "cs",
    n_groups = 7L,
    n_columns = 5L,
    n_rows = 11L,
    monitor_policy = monitor,
    centered = TRUE
  )

  expect_equal(unname(estimates["transform_products"]), 7 * 5^2 + 11)
})

test_that("group-local complexity guard fails before oversized syntax emission", {

  model_matrix <- matrix(0, nrow = 8L, ncol = 5L)
  row_column   <- c(3L, 1L, 3L, 5L, 2L, 5L, 4L, 4L)
  group_map    <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L)
  model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)] <- 1
  layout <- BayesTools:::.bt_random_effect_structured_local_layout(
    model_matrix = model_matrix,
    group_map = group_map,
    structure = "cs",
    exact_indicator = TRUE
  )
  monitor <- random_monitor()

  expect_equal(
    BayesTools:::.bt_random_effect_group_local_complexity(
      layout = layout,
      n_rows = 8L,
      monitor_policy = monitor
    ),
    c(syntax_nodes = 25, monitored_values = 6)
  )
  expect_equal(
    BayesTools:::.bt_random_effect_group_local_complexity(
      layout = layout,
      n_rows = 8L,
      monitor_policy = monitor,
      row_indexed_external_sd = TRUE,
      column_allocation = TRUE
    ),
    c(syntax_nodes = 25, monitored_values = 11)
  )

  old_multiplier <- getOption("BayesTools.random_effects_complexity_multiplier")
  on.exit(options(
    BayesTools.random_effects_complexity_multiplier = old_multiplier
  ), add = TRUE)
  options(BayesTools.random_effects_complexity_multiplier = 1e-4)
  expect_error(
    BayesTools:::.bt_random_effect_check_group_local_complexity(
      random_term = list(block_name = "id"),
      layout = layout,
      n_rows = 8L,
      monitor_policy = monitor
    ),
    "exceeds the current group-local JAGS compiler limits",
    fixed = TRUE
  )

  K <- 40L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(5L)), length.out = K))
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + cs(index | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd_source = random_sd_source("tau", shape = "row"),
          cor = prior("point", list(location = 0.2))
        )
      )
    ),
    "exceeds the current group-local JAGS compiler limits",
    fixed = TRUE
  )
})

test_that("high-K sparse row-indexed structure compiles without dense JAGS data", {

  K <- 1024L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(128L)), length.out = K))
  )
  result <- JAGS_formula(
    formula = ~ 1 + ar1(index | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source("tau", shape = "row"),
        cor = prior("point", list(location = 0.5))
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]

  expect_equal(term$latent_layout$n_local, K)
  expect_equal(length(result$data$mu__xREx__id_xRE_MAPx), K)
  expect_equal(length(result$data$mu__xREx__id_xRE_COLx), K)
  expect_null(result$data$mu__xREx__id_xRE_DATAx)
  expect_false(grepl("inprod(", result$formula_syntax, fixed = TRUE))
})

test_that("row-indexed HAR and CAR use the compact structured projection", {

  K <- 113L
  data <- data.frame(
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    time = seq_len(K),
    id = factor(rep(paste0("group_", seq_len(17L)), length.out = K))
  )
  formulas <- list(
    har = ~ 1 + har(index | id),
    car = ~ 1 + car(0 + time | id)
  )

  for(structure in names(formulas)){
    result <- JAGS_formula(
      formula = formulas[[structure]],
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd_source = random_sd_source("tau", shape = "row"),
          cor = prior("point", list(location = 0.5))
        )
      )
    )
    term <- result$formula_design$random_effects[[1L]]

    expect_equal(term$latent_layout$n_local, K, info = structure)
    expect_null(result$data$mu__xREx__id_xRE_DATAx, info = structure)
    expect_equal(
      names(result$data),
      c("N_mu", "mu__xREx__id_xRE_MAPx", "mu__xREx__id_xRE_COLx"),
      info = structure
    )
    expect_false(
      grepl("inprod(", result$formula_syntax, fixed = TRUE),
      info = structure
    )
    if(identical(structure, "car")){
      expect_match(
        result$formula_syntax,
        "mu__xREx__id_xRE_CAR_LOG_PHIX[",
        fixed = TRUE
      )
      expect_match(
        result$formula_syntax,
        " <- pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[",
        fixed = TRUE
      )
      expect_false(grepl(
        "pow(mu__xREx__id_rho",
        result$formula_syntax,
        fixed = TRUE
      ))
      expect_false(grepl(
        "1 - pow(mu__xREx__id_xRE_CAR_PHIX",
        result$formula_syntax,
        fixed = TRUE
      ))
    }
  }
})
