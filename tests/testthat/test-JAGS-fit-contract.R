skip_if_not_test_profile("unit")

test_that("structured JAGS parameter encoding is injective and reversible", {

  cases <- list(
    list(kind = "fixed", formula_parameter = "mu", term = "a_b:c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu_a", term = "b c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu", term = "a__xXx__b", role = "coefficient"),
    list(kind = "random", formula_parameter = "mu", term = "group[one]", role = "sd"),
    list(kind = "fixed", formula_parameter = "mu", term = intToUtf8(c(946L, 233L)), role = "coefficient")
  )

  encoded <- vapply(cases, .bt_parameter_encode, character(1))
  expect_length(unique(encoded), length(cases))
  expect_true(all(grepl("^[A-Za-z][A-Za-z0-9_.]*$", encoded)))
  for(i in seq_along(cases)){
    decoded <- .bt_parameter_decode(encoded[i])
    expect_identical(decoded[names(cases[[i]])], cases[[i]])
    expect_identical(decoded$encoding_version, 1L)
  }

  left <- .bt_parameter_encode(list(
    kind = "fixed", formula_parameter = "a_b", term = "c", role = "coefficient"
  ))
  right <- .bt_parameter_encode(list(
    kind = "fixed", formula_parameter = "a", term = "b_c", role = "coefficient"
  ))
  expect_false(identical(left, right))
  expect_error(.bt_parameter_decode("BT2_00_00_00_00"), "unsupported")
  expect_error(.bt_parameter_decode("BT1_0_00_00_00"), "malformed")
  for(invalid_utf8 in c("FF", "C0AF", "EDA080", "C3")){
    expect_error(
      .bt_parameter_decode(paste0("BT1_", invalid_utf8, "_6D75__636F6566")),
      "invalid UTF-8",
      fixed = TRUE
    )
  }
})

test_that("formula designs persist a validated semantic name map", {

  result <- JAGS_formula(
    formula = ~ x + f,
    parameter = "mu",
    data = data.frame(
      x = c(-1, 0, 1, 2),
      f = factor(c("a", "b", "a", "b"))
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  map <- result$formula_design$name_map

  expect_identical(result$formula_design$schema_version, 4L)
  expect_s3_class(map, "BayesTools_formula_name_map")
  expect_identical(attr(map, "schema_version"), 1L)
  expect_setequal(map$jags_name, c("mu", "mu_intercept", "mu_x", "mu_f"))
  expect_identical(map$term[map$jags_name == "mu_f"], "f")

  auxiliary <- JAGS_formula(
    formula = ~ x,
    parameter = "mu",
    data = data.frame(x = c(-1, 0, 1, 2)),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior_spike_and_slab(prior("normal", list(0, 1)))
    )
  )$formula_design$name_map
  indicator <- auxiliary[auxiliary$jags_name == "mu_x_indicator", ]
  expect_identical(indicator$kind, "fixed_auxiliary")
  expect_identical(indicator$term, "x")
  expect_identical(indicator$role, "_indicator")

  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit <- .bt_attach_fit_contract(fit)
  expect_identical(JAGS_formula_name_map(fit, "mu"), map)
  expect_identical(
    JAGS_fit_contract(fit)$formula_design_version,
    4L
  )
  expect_silent(JAGS_validate_fit_contract(
    fit,
    requires = c("name_encoding", "formula_name_map", "formula_design")
  ))

  broken <- fit
  attr(broken, "fit_contract")$formula_name_map_version <- 2L
  expect_error(
    JAGS_formula_name_map(broken, "mu"),
    "Refit the model"
  )

  missing_map <- list(mu = result$formula_design)
  missing_map$mu$name_map <- NULL
  attr(broken, "formula_design") <- missing_map
  broken <- .bt_attach_fit_contract(broken)
  expect_error(JAGS_formula_name_map(broken), "name-map metadata are missing")
  expect_error(
    .bt_build_parameter_coordinates(
      columns = "mu_intercept", formula_design = missing_map
    ),
    "name-map metadata are missing"
  )
})

test_that("functions reading fitted metadata refuse fits without the current contract", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{}",
    prior_list = list(
      mu    = prior("normal", list(0, 1)),
      sigma = prior("normal", list(0, 1), list(0, Inf))
    ),
    chains = 2, adapt = 50, burnin = 50, sample = 200, seed = 1
  ))
  models <- function(object){
    list(list(fit = object, marglik = bridgesampling_object(0), prior_weights = 1))
  }
  calls <- list(
    JAGS_check_convergence = function(object){
      JAGS_check_convergence(object)
    },
    JAGS_diagnostics = function(object){
      JAGS_diagnostics(object, "mu", type = "trace", plot_type = "ggplot")
    },
    as_mixed_posteriors = function(object){
      as_mixed_posteriors(object, "mu")
    },
    mix_posteriors = function(object){
      mix_posteriors(models(object), "mu", list(mu = FALSE))
    },
    JAGS_extend = function(object){
      JAGS_extend(object, autofit_control = list(max_extend = 1, sample_extend = 10))
    }
  )

  # The metadata state of a BayesTools 0.3.0 fit: no parameter map, contract,
  # or draw geometry. The message is that of the summary tables.
  stripped <- fit
  for(name in c("parameter_map", "fit_contract", "draw_geometry")){
    attr(stripped, name) <- NULL
  }
  missing_map <- paste0(
    "The fitted object does not contain parameter-map metadata. ",
    "Refit the model with the current BayesTools version."
  )
  expect_error(runjags_estimates_table(stripped), missing_map, fixed = TRUE)

  without_contract <- fit
  attr(without_contract, "fit_contract") <- NULL
  missing_contract <- paste0(
    "The fitted object does not contain a supported schema contract. ",
    "Refit the model with this version of BayesTools."
  )

  previous_map <- fit
  map <- attr(previous_map, "parameter_map", exact = TRUE)
  map$schema_version <- map$schema_version - 1L
  attr(previous_map, "parameter_map") <- map
  unsupported_map <- paste0(
    "Parameter-map metadata are missing, malformed, or unsupported. ",
    "Refit the model with the current BayesTools version."
  )

  for(name in names(calls)){
    expect_error(calls[[name]](stripped), missing_map, fixed = TRUE, info = name)
    expect_error(calls[[name]](without_contract), missing_contract, fixed = TRUE, info = name)
    expect_error(calls[[name]](previous_map), unsupported_map, fixed = TRUE, info = name)
  }

  # Plain runjags objects carry no fitted metadata at all.
  plain <- fit
  class(plain) <- "runjags"
  expect_error(
    JAGS_check_convergence(plain),
    "'fit' must be a 'BayesTools_fit' created by JAGS_fit(). Refit the model with this version of BayesTools.",
    fixed = TRUE
  )
  expect_error(
    mix_posteriors(models(plain), "mu", list(mu = FALSE)),
    "'model_list:fit' must be a 'BayesTools_fit' created by JAGS_fit(). Refit the model with this version of BayesTools.",
    fixed = TRUE
  )

  # Fits of this version pass.
  expect_true(JAGS_check_convergence(
    fit, max_Rhat = 2, min_ESS = 1, max_error = NULL, max_SD_error = NULL
  ))
  expect_s3_class(as_mixed_posteriors(fit, "mu")$mu, "mixed_posteriors")
})
