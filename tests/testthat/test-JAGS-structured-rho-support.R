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
    marglik = BayesTools:::.bt_JAGS_marglik_random_effect_scalar_rho_support(
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

test_that("scalar rho support validates canonical bounds and CAR metadata", {

  malformed_bounds <- .structured_rho_support_term("har")
  malformed_bounds$correlation$bounds <- c(lower = -2, upper = 2)
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_scalar_rho_support(
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
    BayesTools:::.bt_JAGS_marglik_random_effect_scalar_rho_support(
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
    BayesTools:::.bt_JAGS_marglik_random_effect_scalar_rho_support(
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
        rho = prior("normal", list(0, 0.5))
      )
    )
  )

  expect_true(grepl(
    "max(0, min(0.99999999999999989,",
    result$formula_syntax,
    fixed = TRUE
  ))
  expect_true(grepl(
    "pow(mu__xREx__id_rho, 0.5)",
    result$formula_syntax,
    fixed = TRUE
  ))
  expect_true(grepl(
    "pow(mu__xREx__id_rho, 1.5)",
    result$formula_syntax,
    fixed = TRUE
  ))
})
