skip_if_not_test_profile("unit")

test_that("structured JAGS parameter encoding is injective and reversible", {

  cases <- list(
    list(kind = "fixed", formula_parameter = "mu", term = "a_b:c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu_a", term = "b c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu", term = "a__xXx__b", role = "coefficient"),
    list(kind = "random", formula_parameter = "mu", term = "group[one]", role = "sd"),
    list(kind = "fixed", formula_parameter = "mu", term = intToUtf8(c(946L, 233L)), role = "coefficient")
  )

  encoded <- vapply(cases, JAGS_parameter_encode, character(1))
  expect_length(unique(encoded), length(cases))
  expect_true(all(grepl("^[A-Za-z][A-Za-z0-9_.]*$", encoded)))
  for(i in seq_along(cases)){
    decoded <- JAGS_parameter_decode(encoded[i])
    expect_identical(decoded[names(cases[[i]])], cases[[i]])
    expect_identical(decoded$encoding_version, 1L)
  }

  left <- JAGS_parameter_encode(list(
    kind = "fixed", formula_parameter = "a_b", term = "c", role = "coefficient"
  ))
  right <- JAGS_parameter_encode(list(
    kind = "fixed", formula_parameter = "a", term = "b_c", role = "coefficient"
  ))
  expect_false(identical(left, right))
  expect_error(JAGS_parameter_decode("BT2_00_00_00_00"), "unsupported")
  expect_error(JAGS_parameter_decode("BT1_0_00_00_00"), "malformed")
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
})
