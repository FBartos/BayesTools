skip_if_not_test_profile("unit")

.correlation_draws_sd_prior <- function(){

  prior("point", list(location = 1))
}

.correlation_draws_term <- function(structure, rho = NULL,
                                    rho_scale = "fisher_z"){

  data <- data.frame(
    id = factor(rep(c("g1", "g2"), each = 3L)),
    index = factor(rep(c("i1", "i2", "i3"), 2L)),
    time = rep(c(0, 0.5, 2), 2L)
  )
  formula <- switch(
    structure,
    cs = ~ 1 + cs(index | id),
    hcs = ~ 1 + hcs(index | id),
    ar1 = ~ 1 + ar1(index | id),
    har = ~ 1 + har(index | id),
    car = ~ 1 + car(time | id)
  )
  if(is.null(rho)){
    rho <- prior("normal", list(0, 0.5))
  }
  block <- if(identical(rho_scale, "fisher_z")){
    random_block(
      sd = .correlation_draws_sd_prior(),
      rho = rho
    )
  }else{
    random_block(
      sd = .correlation_draws_sd_prior(),
      covariance = random_covariance(
        rho = rho,
        rho_scale = rho_scale
      )
    )
  }
  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = block
    )
  )

  result$formula_design$random_effects[[1L]]
}

test_that("scalar structured correlation draws use compiled metadata", {

  cases <- list(
    cs = c(-0.25, 0.6),
    hcs = c(-0.25, 0.6),
    ar1 = c(-0.5, 0.6),
    har = c(-0.5, 0.6),
    car = c(0.25, 0.81)
  )
  for(structure in names(cases)){
    term <- .correlation_draws_term(structure)
    rho <- cases[[structure]]
    posterior <- cbind(mu = seq_along(rho), rho)
    colnames(posterior)[[2L]] <- term$correlation$rho_name
    actual <- random_effects_correlation_draws(term, posterior)
    distance <- if(structure %in% c("cs", "hcs")){
      1 - diag(3L)
    }else if(identical(structure, "car")){
      abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
    }else{
      abs(outer(seq_len(3L), seq_len(3L), "-"))
    }

    expect_equal(dim(actual), c(2L, 3L, 3L), info = structure)
    for(draw in seq_along(rho)){
      expected <- if(identical(structure, "car")){
        BayesTools:::.bt_random_effect_structured_subset_correlation(
          structure = structure,
          columns = seq_len(3L),
          rho = rho[draw],
          global_n_columns = 3L,
          column_coordinates = c(0, 0.5, 2)
        )
      }else{
        rho[[draw]]^distance
      }
      expect_equal(
        unname(actual[draw, , ]),
        expected,
        tolerance = 0,
        info = paste(structure, "draw", draw)
      )
    }
    expect_equal(dimnames(actual)[[2L]], term$column_names, info = structure)
    expect_null(term$correlation$distance_matrix, info = structure)
  }
})

test_that("fixed transformed rho and one-column blocks are reconstructed", {

  term <- .correlation_draws_term(
    "cs",
    rho = prior("point", list(location = 0.5))
  )
  posterior <- cbind(mu = 1:2)
  actual <- random_effects_correlation_draws(term, posterior)
  expected <- matrix(tanh(0.5), nrow = 3L, ncol = 3L)
  diag(expected) <- 1
  expect_equal(unname(actual[1L, , ]), expected, tolerance = 1e-15)
  expect_equal(unname(actual[2L, , ]), expected, tolerance = 1e-15)

  term$n_columns <- 1L
  term$column_names <- "intercept"
  term$correlation <- NULL
  scalar <- random_effects_correlation_draws(term, posterior)
  expect_equal(dim(scalar), c(2L, 1L, 1L))
  expect_equal(as.vector(scalar), c(1, 1), tolerance = 0)

  zero_column_posterior <- matrix(
    numeric(),
    nrow = 2L,
    ncol = 0L,
    dimnames = list(c("draw_1", "draw_2"), NULL)
  )
  scalar <- random_effects_correlation_draws(term, zero_column_posterior)
  expect_equal(dim(scalar), c(2L, 1L, 1L))
  expect_equal(as.vector(scalar), c(1, 1), tolerance = 0)
  expect_equal(dimnames(scalar), list(
    draw = c("draw_1", "draw_2"),
    row = "intercept",
    column = "intercept"
  ))

  term$n_columns <- 2L
  term$column_names <- c("intercept", "slope")
  expect_error(
    random_effects_correlation_draws(term, zero_column_posterior),
    "'posterior_samples' must have non-empty column names",
    fixed = TRUE
  )
})

test_that("correlation reconstruction validates support and compact CAR metadata", {

  cs_term <- .correlation_draws_term("cs")
  invalid_rho <- cbind(mu = 1, rho = -0.75)
  colnames(invalid_rho)[[2L]] <- cs_term$correlation$rho_name
  expect_error(
    random_effects_correlation_draws(cs_term, invalid_rho),
    "out-of-support draw",
    fixed = TRUE
  )
  expect_error(
    random_effects_correlation_draws(cs_term, cbind(mu = 1)),
    "missing canonical scalar correlation coordinates",
    fixed = TRUE
  )

  car_term <- .correlation_draws_term("car")
  car_term$correlation$time_values <- NULL
  car_term$car$time_values <- c(0, 0, 2)
  car_posterior <- cbind(mu = 1, rho = 0.5)
  colnames(car_posterior)[[2L]] <- car_term$correlation$rho_name
  expect_error(
    random_effects_correlation_draws(car_term, car_posterior),
    "ordered CAR time coordinates",
    fixed = TRUE
  )

  car_term$car$time_values <- NULL
  car_term$car$distance_matrix <-
    abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
  expect_error(
    random_effects_correlation_draws(car_term, car_posterior),
    "ordered CAR time coordinates",
    fixed = TRUE
  )

  cs_term$n_columns <- 2.5
  expect_error(
    random_effects_correlation_draws(cs_term, invalid_rho),
    "positive integer 'random_term$n_columns'",
    fixed = TRUE
  )
  expect_error(
    random_effects_correlation_draws(cs_term, list(mu = 1)),
    "matrix, data frame, mcmc, or mcmc.list",
    fixed = TRUE
  )

  unsupported <- .correlation_draws_term("cs")
  unsupported$structure <- "us"
  expect_error(
    random_effects_correlation_draws(unsupported, invalid_rho),
    "supports only scalar CS, HCS, AR1, HAR, and CAR",
    fixed = TRUE
  )
})

test_that("CAR correlation factors report canonical failing transitions", {

  smallest <- .Machine$double.xmin * .Machine$double.eps
  coordinates <- c(smallest, 2 * smallest, 1)
  rho <- 0.9
  random_term <- .correlation_draws_term("car")
  random_term$correlation$time_values <- coordinates
  random_term$car$time_values <- coordinates
  posterior <- matrix(
    rho,
    nrow = 1L,
    dimnames = list(NULL, random_term$correlation$rho_name)
  )
  distance <- abs(outer(coordinates, coordinates, "-"))
  latent_names <- as.vector(
    BayesTools:::.bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = length(random_term$group_levels),
      n_columns = 3L
    )
  )
  reconstruction_values <- c(rho, rep(0, length(latent_names)))
  names(reconstruction_values) <- c(
    random_term$correlation$rho_name,
    latent_names
  )
  reconstruction_posterior <- matrix(
    reconstruction_values,
    nrow = 1L,
    dimnames = list(NULL, names(reconstruction_values))
  )

  error_messages <- c(
    exported = tryCatch(
      random_effects_correlation_draws(random_term, posterior),
      error = function(error) conditionMessage(error)
    ),
    cholesky = tryCatch(
      BayesTools:::.bt_random_effect_cholesky_draws(
        random_term = random_term,
        n_columns = 3L,
        posterior = posterior
      ),
      error = function(error) conditionMessage(error)
    ),
    dense_reconstruction = tryCatch(
      BayesTools:::.bt_random_effect_structured_dense_contribution(
        random_term = random_term,
        model_matrix = random_term$model_matrix,
        group_map = random_term$group_map,
        posterior = reconstruction_posterior,
        scale_draws = matrix(1, nrow = 1L, ncol = 3L)
      ),
      error = function(error) conditionMessage(error)
    ),
    internal = tryCatch(
      BayesTools:::.bt_random_effect_structured_correlation_matrix(
        structure = "car",
        K = 3L,
        rho = rho,
        distance_matrix = distance,
        column_coordinates = coordinates,
        context = "Internal CAR correlation reconstruction"
      ),
      error = function(error) conditionMessage(error)
    )
  )
  expected_context <- c(
    exported = paste0(
      "Random-effect correlation reconstruction for block 'id', ",
      "posterior draw 1"
    ),
    cholesky = paste0(
      "Random-effect Cholesky reconstruction for block 'id', ",
      "posterior draw 1"
    ),
    dense_reconstruction = paste0(
      "Random-effect posterior reconstruction for block 'id', ",
      "posterior draw 1, group 1"
    ),
    internal = "Internal CAR correlation reconstruction"
  )
  labels <- format(
    c(coordinates[1L], coordinates[2L], rho, smallest),
    digits = 17L,
    scientific = TRUE,
    trim = TRUE
  )

  for(path in names(error_messages)){
    error_message <- error_messages[[path]]
    expect_match(
      error_message,
      expected_context[[path]],
      fixed = TRUE,
      info = path
    )
    expect_match(
      error_message,
      paste0("from coordinate ", labels[1L], " to ", labels[2L]),
      fixed = TRUE,
      info = path
    )
    expect_match(
      error_message,
      paste0("rho = ", labels[3L]),
      fixed = TRUE,
      info = path
    )
    expect_match(
      error_message,
      paste0("gap = ", labels[4L]),
      fixed = TRUE,
      info = path
    )
    expect_false(
      grepl(
        "from coordinate 0.0000000000000000e+00",
        error_message,
        fixed = TRUE
      ),
      info = path
    )
    expect_match(
      error_message,
      "requested coordinate/time resolution is not representable",
      fixed = TRUE,
      info = path
    )
  }
})

test_that("dense CAR correlation APIs are reconstructed from stable factors", {

  coordinates <- c(4, 4.5, 6)
  distance <- abs(outer(coordinates, coordinates, "-"))
  rho <- 1 - .Machine$double.eps / 2
  random_term <- .correlation_draws_term("car")
  random_term$correlation$time_values <- coordinates
  random_term$car$time_values <- coordinates
  posterior <- matrix(
    rho,
    nrow = 1L,
    dimnames = list(NULL, random_term$correlation$rho_name)
  )

  L <- BayesTools:::.bt_random_effect_structured_subset_cholesky(
    structure = "car",
    columns = seq_len(3L),
    rho = rho,
    global_n_columns = 3L,
    column_coordinates = coordinates
  )
  expected <- tcrossprod(L)
  diag(expected) <- 1
  exported <- random_effects_correlation_draws(random_term, posterior)
  internal <- BayesTools:::.bt_random_effect_structured_correlation_matrix(
    structure = "car",
    K = 3L,
    rho = rho,
    distance_matrix = distance,
    column_coordinates = coordinates
  )

  expect_equal(unname(exported[1L, , ]), expected, tolerance = 0)
  expect_equal(internal, expected, tolerance = 0)
  expect_identical(exported[1L, 1L, 2L], 1)
  expect_gt(L[2L, 2L], 0)

  expect_error(
    BayesTools:::.bt_random_effect_structured_correlation_matrix(
      structure = "car",
      K = 3L,
      rho = rho,
      distance_matrix = distance
    ),
    "a distance matrix alone cannot preserve the original coordinate values",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_correlation_matrix(
      structure = "car",
      K = 3L,
      rho = rho,
      distance_matrix = distance + diag(3L),
      column_coordinates = coordinates
    ),
    "zero diagonal",
    fixed = TRUE
  )
})

test_that("correlation reconstruction accepts coda posterior containers", {

  term <- .correlation_draws_term("ar1")
  values <- cbind(mu = 1:2, rho = c(-0.2, 0.4))
  colnames(values)[[2L]] <- term$correlation$rho_name
  samples <- coda::mcmc.list(coda::mcmc(values[1L, , drop = FALSE]),
                             coda::mcmc(values[2L, , drop = FALSE]))
  actual <- random_effects_correlation_draws(term, samples)

  expect_equal(dim(actual), c(2L, 3L, 3L))
  expect_equal(actual[1L, 1L, 3L], (-0.2)^2, tolerance = 0)
  expect_equal(actual[2L, 1L, 3L], 0.4^2, tolerance = 0)
})

test_that("transformed rho coordinates stay in the representable interior", {

  cases <- list(
    ar1_fisher = list(
      term = .correlation_draws_term("ar1"),
      coordinates = c(-1e300, 1e300)
    ),
    car_fisher = list(
      term = .correlation_draws_term("car"),
      coordinates = c(0, 1e300)
    ),
    car_logit = list(
      term = .correlation_draws_term("car", rho_scale = "logit"),
      coordinates = c(-1e300, 1e300)
    )
  )

  for(name in names(cases)){
    random_term <- cases[[name]]$term
    correlation <- random_term$correlation
    interior <- BayesTools:::.bt_random_effect_representable_rho_bounds(
      correlation$bounds,
      random_term$structure
    )
    posterior <- matrix(
      cases[[name]]$coordinates,
      ncol = 1L,
      dimnames = list(NULL, correlation$sample_name)
    )
    rho <- BayesTools:::.bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      out_of_support = "error"
    )

    expect_equal(rho, unname(interior), tolerance = 0, info = name)
    if(identical(random_term$structure, "car")){
      expect_gte(rho[[1L]], correlation$bounds[["lower"]])
    }else{
      expect_gt(rho[[1L]], correlation$bounds[["lower"]])
    }
    expect_lt(rho[[2L]], correlation$bounds[["upper"]])
    expect_equal(
      vapply(seq_len(nrow(posterior)), function(i){
        BayesTools:::.bt_JAGS_marglik_random_effect_rho(
          posterior[i, ],
          random_term
        )
      }, numeric(1)),
      unname(interior),
      tolerance = 0,
      info = paste(name, "bridge")
    )
  }

  near_zero <- BayesTools:::.bt_random_effect_representable_rho_bounds(
    c(lower = -1e-20, upper = 1),
    "cs"
  )
  expect_gt(near_zero[["lower"]], -1e-20)
  expect_lt(near_zero[["lower"]] - (-1e-20), 1e-30)
})

test_that("scalar structured Cholesky reconstruction requires canonical rho", {

  random_term <- .correlation_draws_term("ar1")
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(random_term, 3L)
  legacy_L <- diag(3L)
  legacy <- matrix(
    as.vector(legacy_L),
    nrow = 1L,
    dimnames = list(NULL, as.vector(L_names))
  )
  expect_error(
    BayesTools:::.bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = 3L,
      posterior = legacy
    ),
    "missing canonical scalar correlation coordinates",
    fixed = TRUE
  )

  invalid_rho <- cbind(legacy, 1)
  colnames(invalid_rho)[ncol(invalid_rho)] <- random_term$correlation$rho_name
  expect_error(
    BayesTools:::.bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = 3L,
      posterior = invalid_rho
    ),
    "out-of-support draw",
    fixed = TRUE
  )

  invalid_sample <- cbind(legacy, Inf)
  colnames(invalid_sample)[ncol(invalid_sample)] <-
    random_term$correlation$sample_name
  expect_error(
    BayesTools:::.bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = 3L,
      posterior = invalid_sample
    ),
    "scalar correlation coordinates must be finite",
    fixed = TRUE
  )

  invalid_metadata <- random_term
  invalid_metadata$correlation$type <- "lkj"
  expect_error(
    BayesTools:::.bt_random_effect_cholesky_draws(
      random_term = invalid_metadata,
      n_columns = 3L,
      posterior = legacy
    ),
    "does not define canonical scalar correlation metadata",
    fixed = TRUE
  )
})

test_that("scalar structured Cholesky reconstruction matches dense factors", {

  cases <- list(
    cs = c(-0.25, 0.6),
    hcs = c(-0.25, 0.6),
    ar1 = c(-0.5, 0.6),
    har = c(-0.5, 0.6),
    car = c(0.25, 0.81)
  )
  for(structure in names(cases)){
    random_term <- .correlation_draws_term(structure)
    rho <- cases[[structure]]
    posterior <- matrix(
      rho,
      ncol = 1L,
      dimnames = list(NULL, random_term$correlation$rho_name)
    )
    actual <- BayesTools:::.bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = 3L,
      posterior = posterior
    )
    distance <- if(structure %in% c("cs", "hcs")){
      1 - diag(3L)
    }else if(identical(structure, "car")){
      abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
    }else{
      abs(outer(seq_len(3L), seq_len(3L), "-"))
    }

    expect_equal(dim(actual), c(2L, 3L, 3L), info = structure)
    expect_null(dimnames(actual), info = structure)
    for(draw in seq_along(rho)){
      expected <- t(chol(rho[[draw]]^distance))
      expect_equal(
        actual[draw, , ],
        expected,
        tolerance = 1e-14,
        info = paste(structure, "draw", draw)
      )
    }
  }
})

test_that("CS Cholesky reconstruction is stable near its global lower bound", {

  n_columns <- 100L
  bounds <- BayesTools:::.bt_random_effect_structured_rho_bounds(
    K = n_columns,
    structure = "cs"
  )
  rho <- BayesTools:::.bt_random_effect_representable_rho_bounds(
    bounds,
    "cs"
  )[["lower"]]
  random_term <- list(
    structure = "cs",
    n_columns = n_columns,
    block_name = "near-boundary",
    correlation = list(
      type = "rho",
      rho_name = "rho",
      sample_name = "rho",
      rho_scale = "rho",
      sample_fixed = NULL,
      bounds = bounds
    )
  )
  posterior <- matrix(
    rho,
    nrow = 1L,
    dimnames = list(NULL, "rho")
  )

  actual <- BayesTools:::.bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )[1L, , ]
  expected <- matrix(rho, nrow = n_columns, ncol = n_columns)
  diag(expected) <- 1

  expect_true(all(is.finite(actual)))
  expect_gt(min(diag(actual)), 0)
  expect_equal(actual[upper.tri(actual)], rep(0, sum(upper.tri(actual))),
               tolerance = 0)
  expect_equal(tcrossprod(actual), expected, tolerance = 1e-12)
})
