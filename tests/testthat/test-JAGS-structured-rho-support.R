skip_if_not_test_profile("unit")

.structured_rho_support_term <- function(structure, K = 3L,
                                         distance_matrix = NULL){

  bounds <- BayesTools:::.bt_random_effect_structured_rho_bounds(
    K = K,
    structure = structure
  )
  correlation <- list(
    type = "rho",
    structure = structure,
    rho_name = "rho",
    sample_name = "rho",
    prior_name = "rho",
    sample_fixed = NULL,
    rho_scale = "rho",
    bounds = bounds,
    distance_matrix = distance_matrix,
    time_values = if(identical(structure, "car")) seq_len(K) else NULL
  )

  list(
    block_name = "id",
    structure = structure,
    n_columns = K,
    correlation = correlation
  )
}

.structured_rho_support_values <- function(random_term, rho){

  structure <- random_term$structure
  samples   <- c(rho = rho)
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_random_effect_scalar_rho_support(
    random_term = random_term,
    structure = structure
  )

  unname(c(
    marglik = bayestools_reference_random_effect_scalar_rho_support(
      samples = samples,
      random_term = random_term
    ),
    bridge = compiled(samples)
  ))
}

test_that("scalar structured rho support is exact without dense matrices", {

  car_distance <- abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
  cases <- list(
    cs  = list(lower = -0.5, interior = -0.5 + 1e-8),
    hcs = list(lower = -0.5, interior = -0.5 + 1e-8),
    ar1 = list(lower = -1, interior = -1 + 1e-8),
    har = list(lower = -1, interior = -1 + 1e-8),
    car = list(lower = 0, interior = 0)
  )

  for(structure in names(cases)){
    random_term <- .structured_rho_support_term(
      structure = structure,
      distance_matrix = if(identical(structure, "car")) car_distance else NULL
    )
    lower_expected <- if(identical(structure, "car")) c(0, 0) else c(-Inf, -Inf)

    expect_equal(
      .structured_rho_support_values(random_term, cases[[structure]]$lower),
      lower_expected,
      info = structure
    )
    expect_equal(
      .structured_rho_support_values(random_term, cases[[structure]]$interior),
      c(0, 0),
      info = structure
    )
    expect_equal(
      .structured_rho_support_values(random_term, 1),
      c(-Inf, -Inf),
      info = structure
    )
  }

  large_ar <- .structured_rho_support_term("ar1", K = 100000L)
  expect_equal(
    .structured_rho_support_values(large_ar, 0.25),
    c(0, 0)
  )
})

test_that("omitted structured rho priors cover each complete admissible range", {

  sd_prior <- prior("gamma", list(shape = 2, rate = 2))
  block_prior <- BayesTools:::.bt_random_prior_for_block(
    prior_random(sd = sd_prior),
    "id"
  )
  expected <- list(
    cs  = c(lower = -0.5, upper = 1),
    hcs = c(lower = -0.5, upper = 1),
    ar1 = c(lower = -1, upper = 1),
    har = c(lower = -1, upper = 1),
    car = c(lower = 0, upper = 1)
  )

  for(structure in names(expected)){
    rho_info <- BayesTools:::.bt_random_effect_structured_rho_prior(
      prior_prefix = "mu__xREx__id",
      node_prefix = "mu__xREx__id",
      K = 3L,
      structure = structure,
      block_prior = block_prior
    )
    rho_prior <- rho_info$prior_list[[1L]]

    expect_equal(rho_info$rho_scale, "rho", info = structure)
    expect_equal(rho_info$bounds, expected[[structure]], info = structure)
    expect_equal(rho_prior$distribution, "uniform", info = structure)
    expect_equal(
      rho_prior$parameters,
      list(a = expected[[structure]][["lower"]], b = 1),
      info = structure
    )
  }

  transformed_without_prior <- BayesTools:::.bt_random_prior_for_block(
    prior_random(
      sd = sd_prior,
      covariance = random_covariance(cor_scale = "fisher_z")
    ),
    "id"
  )
  expect_error(
    BayesTools:::.bt_random_effect_structured_rho_prior(
      prior_prefix = "mu__xREx__id",
      node_prefix = "mu__xREx__id",
      K = 3L,
      structure = "hcs",
      block_prior = transformed_without_prior
    ),
    "requires an explicit 'cor' prior",
    fixed = TRUE
  )
})

test_that("scalar rho support validates canonical bounds and CAR metadata", {

  malformed_bounds <- .structured_rho_support_term("har")
  malformed_bounds$correlation$bounds <- c(lower = -2, upper = 2)
  expect_error(
    bayestools_reference_random_effect_scalar_rho_support(
      samples = c(rho = 0),
      random_term = malformed_bounds
    ),
    "do not match the exact HAR support",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_random_effect_scalar_rho_support(
      random_term = malformed_bounds,
      structure = "har"
    ),
    "do not match the exact HAR support",
    fixed = TRUE
  )

  missing_car_distance <- .structured_rho_support_term("car")
  missing_car_distance$correlation$time_values <- NULL
  expect_error(
    bayestools_reference_random_effect_scalar_rho_support(
      samples = c(rho = 0.5),
      random_term = missing_car_distance
    ),
    "missing canonical ordered CAR time coordinates",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_random_effect_scalar_rho_support(
      random_term = missing_car_distance,
      structure = "car"
    ),
    "missing canonical ordered CAR time coordinates",
    fixed = TRUE
  )

  invalid_car_distance <- .structured_rho_support_term(
    "car",
    distance_matrix = matrix(0, nrow = 2L, ncol = 2L)
  )
  invalid_car_distance$correlation$time_values <- c(1, 1, 2)
  expect_error(
    bayestools_reference_random_effect_scalar_rho_support(
      samples = c(rho = 0.5),
      random_term = invalid_car_distance
    ),
    "missing canonical ordered CAR time coordinates",
    fixed = TRUE
  )
})

test_that("structured random-effect JAGS literals are locale independent", {

  withr::local_options(list(OutDec = ","))

  expect_equal(
    BayesTools:::.bt_JAGS_numeric_literal(0.5),
    "0.5"
  )

  data <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(rep(c("a", "b"), each = 3L))
  )
  result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
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

  expect_true(grepl(
    "mu__xREx__id_rho <- tanh(mu__xREx__id_rho_z)",
    result$formula_syntax,
    fixed = TRUE
  ))
  expect_true(grepl(
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[2] <- ",
      "0.5 * log(mu__xREx__id_rho)"
    ),
    result$formula_syntax,
    fixed = TRUE
  ))
  expect_true(grepl(
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[3] <- ",
      "1.5 * log(mu__xREx__id_rho)"
    ),
    result$formula_syntax,
    fixed = TRUE
  ))
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_INNOV_VARx[2] <- ",
      "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "pow(mu__xREx__id_rho",
    result$formula_syntax,
    fixed = TRUE
  ))
})

test_that("CAR innovation validation follows rho support and parameterization", {

  data <- data.frame(
    time = rep(c(0, 1e-320), 2L),
    id = factor(rep(c("a", "b"), each = 2L))
  )
  compile <- function(rho, parameterization = "noncentered",
                      sd = prior("point", list(location = 1)),
                      time_gap = 1e-320){
    compile_data <- data
    compile_data$time <- rep(c(0, time_gap), 2L)
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = compile_data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = sd,
          covariance = random_covariance(
            cor = rho,
            cor_scale = "cor"
          ),
          parameterization = parameterization
        )
      )
    )
  }

  expect_error(
    compile(prior("uniform", list(0, 1))),
    "stable innovation variance is zero or non-finite",
    fixed = TRUE
  )
  expect_error(
    compile(prior_mixture(list(
      prior("uniform", list(0, 0.25)),
      prior("uniform", list(0.9, 1))
    ))),
    "stable innovation variance is zero or non-finite",
    fixed = TRUE
  )

  noncentered <- compile(prior("point", list(location = 0.5)))
  expect_match(
    noncentered$formula_syntax,
    "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)",
    fixed = TRUE
  )
  expect_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered"
    ),
    "emitted conditional-normal precision"
  )
  expect_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      sd = prior("point", list(location = 0.5))
    ),
    "emitted conditional-normal precision"
  )
  expect_no_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      sd = prior("point", list(location = 1e8))
    )
  )
  expect_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      sd = prior(
        "normal",
        list(0, 1),
        truncation = list(lower = 0, upper = Inf)
      )
    ),
    "lower centered SD support"
  )

  equality_q <- 1 / .Machine$double.xmax
  equality_gap <- stats::qexp(equality_q) / (-2 * log(0.5))
  equality_log_phi <- equality_gap * log(0.5)
  expect_identical(
    stats::pexp(-2 * equality_log_phi, rate = 1),
    equality_q
  )
  expect_true(is.infinite(1 / equality_q))
  expect_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      time_gap = equality_gap
    ),
    "emitted conditional-normal precision"
  )

  bounded_q <- 2 / .Machine$double.xmax
  bounded_gap <- stats::qexp(bounded_q) / (-2 * log(0.5))
  bounded_log_phi <- bounded_gap * log(0.5)
  expect_identical(
    stats::pexp(-2 * bounded_log_phi, rate = 1),
    bounded_q
  )
  expect_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      sd = prior("uniform", list(0.5, 1)),
      time_gap = bounded_gap
    ),
    "lower centered SD support"
  )

  representable_q <- 0.5 / .Machine$double.xmax
  representable_gap <- stats::qexp(representable_q) / (-2 * log(0.5))
  representable_log_phi <- representable_gap * log(0.5)
  expect_identical(
    stats::pexp(-2 * representable_log_phi, rate = 1),
    representable_q
  )
  expect_true(is.finite(2 ^ -2 / representable_q))
  expect_no_error(
    compile(
      prior("point", list(location = 0.5)),
      parameterization = "centered",
      sd = prior("point", list(location = 2)),
      time_gap = representable_gap
    )
  )

  expect_error(
    compile(
      prior("point", list(location = 0)),
      parameterization = "centered",
      sd = prior("point", list(location = 1e-200))
    ),
    "unrepresentable initial JAGS precision"
  )

  independent <- compile(
    prior("point", list(location = 0)),
    parameterization = "centered"
  )
  expect_match(
    independent$formula_syntax,
    "mu__xREx__id_xRE_CAR_LOG_PHIX[2] <- ",
    fixed = TRUE
  )
  expect_match(
    independent$formula_syntax,
    " * log(mu__xREx__id_rho)",
    fixed = TRUE
  )
  expect_false(grepl(
    "dmnorm.vcov",
    independent$formula_syntax,
    fixed = TRUE
  ))
  expect_false(
    "mu__xREx__id_rho" %in% independent$add_parameters
  )
})

test_that("single-coordinate centered CAR validates initial precision without rho", {

  data <- data.frame(
    time = c(0, 0),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  compile <- function(sd){
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = sd,
          parameterization = "centered"
        )
      )
    )
  }

  expect_error(
    compile(prior("point", list(location = 1e-200))),
    "unrepresentable initial JAGS precision",
    fixed = TRUE
  )
  expect_error(
    compile(prior(
      "normal",
      list(0, 1),
      truncation = list(lower = 0, upper = Inf)
    )),
    "lower centered SD support",
    fixed = TRUE
  )

  result <- NULL
  expect_no_error(
    result <- compile(prior("point", list(location = 1)))
  )
  expect_false(any(grepl(
    "_rho",
    names(result$prior_list),
    fixed = TRUE
  )))
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_COEFx[g,i] ~ dnorm(0, ",
      "pow(mu__xREx__id_xRE_STDx[i], -2))"
    ),
    fixed = TRUE
  )
})

test_that("CAR Fisher-z syntax matches reconstruction at small positive values", {

  z <- 2 ^ -54
  data <- data.frame(
    time = rep(c(0, 0.5, 2), 2L),
    id = factor(rep(c("a", "b"), each = 3L))
  )
  result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        cor = prior("point", list(location = z))
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]

  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__id_rho <- tanh(mu__xREx__id_rho_z)"
    ),
    fixed = TRUE
  )
  reconstructed <- BayesTools:::.bt_random_effect_transform_rho(
    z,
    random_term$correlation,
    random_term = random_term
  )
  support_upper <- BayesTools:::.bt_random_effect_car_rho_support_upper(
    rho_prior = result$prior_list$mu__xREx__id_rho_z,
    rho_scale = random_term$correlation$rho_scale,
    bounds = random_term$correlation$bounds
  )
  expect_equal(support_upper, reconstructed, tolerance = 0)
  expect_equal(reconstructed, tanh(z), tolerance = 0)
  expect_gt(reconstructed, 0)
})

test_that("transformed scalar rho rejects saturated point coordinates", {

  data <- data.frame(
    index = factor(rep(c("i1", "i2", "i3"), 2L)),
    id = factor(rep(c("g1", "g2"), each = 3L))
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + ar1(index | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 1)),
          cor = prior("point", list(location = 1e300))
        )
      )
    ),
    "transformed point correlation is unavailable",
    fixed = TRUE
  )
})
