skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

# Public entry points to BayesTools helpers that downstream packages (RoBMA)
# previously reached through the unexported namespace. Each export is a
# validated wrapper: it must return what the internal returns for valid input
# and refuse invalid input with its own message.

test_that("format_BF_from_log() returns the internal Bayes factor carrier", {

  logs <- c(a = 2, b = 800, c = -800, d = NA, e = 0)
  bounds <- c(">", NA, "<", NA, ">")
  diagnostics <- list(list(x = 1), NULL, "note", NULL, list())
  for(logBF in c(FALSE, TRUE)){
    for(BF01 in c(FALSE, TRUE)){
      for(inclusion in c(FALSE, TRUE)){
        expect_identical(
          format_BF_from_log(logs, logBF, BF01, inclusion, bounds, diagnostics),
          BayesTools:::.format_BF_from_log(logs, logBF, BF01, inclusion, bounds,
                                           diagnostics = diagnostics)
        )
      }
    }
  }
  expect_identical(
    format_BF_from_log(c(1, NaN, -Inf, Inf)),
    BayesTools:::.format_BF_from_log(c(1, NaN, -Inf, Inf))
  )

  # the carrier keeps the logarithm of Bayes factors beyond the double range
  # and is what format_BF() returns for a carrier; a plain vector is not
  carrier <- format_BF_from_log(c(a = 2, b = 800))
  expect_s3_class(carrier, "BayesTools_BF")
  expect_identical(names(carrier), c("a", "b"))
  expect_identical(as.numeric(carrier)[2L], Inf)
  expect_identical(
    format_BF(carrier, logBF = TRUE),
    format_BF_from_log(c(a = 2, b = 800), logBF = TRUE)
  )
  expect_identical(as.numeric(format_BF(carrier, logBF = TRUE)), c(2, 800))
  expect_identical(as.numeric(format_BF(carrier, BF01 = TRUE, logBF = TRUE)), c(-2, -800))
  plain <- format_BF(exp(c(a = 2, b = 800)), logBF = TRUE)
  expect_false(inherits(plain, "BayesTools_BF"))
  expect_identical(as.numeric(plain)[2L], Inf)

  # bound operators refer to the BF10 values and follow the output orientation
  expect_identical(
    attr(format_BF_from_log(c(1, 2), bound_operator = c("<", NA)), "bound_operator"),
    c("<", NA)
  )
  expect_identical(
    attr(format_BF_from_log(c(1, 2), BF01 = TRUE, bound_operator = c("<", NA)), "bound_operator"),
    c(">", NA)
  )
  expect_identical(
    attr(format_BF_from_log(1, inclusion = TRUE, BF01 = TRUE), "name"),
    "Exclusion BF"
  )
})

test_that("format_BF_from_log() validates its arguments", {

  expect_error(format_BF_from_log("a"), "The 'log_BF' argument must be a numeric vector.", fixed = TRUE)
  expect_error(format_BF_from_log(NULL), "The 'log_BF' argument cannot be NULL.", fixed = TRUE)
  expect_error(format_BF_from_log(numeric()), "The 'log_BF' argument cannot be NULL.", fixed = TRUE)
  expect_error(format_BF_from_log(NA), "The 'log_BF' argument must be a numeric vector.", fixed = TRUE)
  expect_error(format_BF_from_log(1, logBF = NA), "The 'logBF' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(format_BF_from_log(1, logBF = c(TRUE, FALSE)), "The 'logBF' argument must have length '1'.", fixed = TRUE)
  expect_error(format_BF_from_log(1, BF01 = "a"), "The 'BF01' argument must be a logical vector.", fixed = TRUE)
  expect_error(format_BF_from_log(1, inclusion = NA), "The 'inclusion' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(format_BF_from_log(1, inclusion = 1), "The 'inclusion' argument must be a logical vector.", fixed = TRUE)
  expect_error(format_BF_from_log(c(1, 2), bound_operator = "="), "The '=' values are not recognized by the 'bound_operator' argument.", fixed = TRUE)
  expect_error(format_BF_from_log(c(1, 2), bound_operator = 1), "The 'bound_operator' argument must be a character vector", fixed = TRUE)
  expect_error(format_BF_from_log(c(1, 2), bound_operator = "<"), "The 'bound_operator' argument must have one element per value of 'log_BF'.", fixed = TRUE)
  expect_error(format_BF_from_log(c(1, 2), diagnostics = 1), "The 'diagnostics' argument must be NULL or a list", fixed = TRUE)
  expect_error(format_BF_from_log(c(1, 2), diagnostics = list(1)), "The 'diagnostics' argument must be NULL or a list", fixed = TRUE)
})

test_that("BF_carrier_drop() removes the private carrier attributes", {

  carrier <- format_BF_from_log(c(a = 2, b = 800), inclusion = TRUE, diagnostics = list(1, 2))
  dropped <- BF_carrier_drop(carrier)
  expect_identical(dropped, BayesTools:::.BF_carrier_drop(carrier))
  expect_false(inherits(dropped, "BayesTools_BF"))
  expect_setequal(names(attributes(dropped)), c("names", "name", "logBF", "BF01", "bound_operator"))
  expect_identical(attr(dropped, "name"), "Inclusion BF")
  expect_identical(names(dropped), c("a", "b"))

  # the class stays while a bound operator is present
  bounded <- format_BF_from_log(c(1, 2), bound_operator = c("<", NA))
  expect_identical(BF_carrier_drop(bounded), BayesTools:::.BF_carrier_drop(bounded))
  expect_s3_class(BF_carrier_drop(bounded), "BayesTools_BF")
  expect_setequal(
    names(attributes(BF_carrier_drop(bounded))),
    c("class", "name", "logBF", "BF01", "bound_operator")
  )

  # plain numbers pass unchanged
  expect_identical(BF_carrier_drop(c(a = 1, b = 2)), c(a = 1, b = 2))
  expect_identical(BF_carrier_drop(1:3), 1:3)

  expect_error(BF_carrier_drop("a"), "The 'x' argument must be a numeric vector.", fixed = TRUE)
  expect_error(BF_carrier_drop(NULL), "The 'x' argument must be a numeric vector.", fixed = TRUE)
  expect_error(BF_carrier_drop(NA), "The 'x' argument must be a numeric vector.", fixed = TRUE)
})

test_that("invert_BF_bound_operator() swaps the operators and rejects other values", {

  operators <- c(a = "<", b = ">", c = NA, d = "<")
  expect_identical(invert_BF_bound_operator(operators), c(">", "<", NA, ">"))
  expect_identical(
    invert_BF_bound_operator(operators),
    BayesTools:::.invert_BF_bound_operator(operators)
  )
  expect_identical(invert_BF_bound_operator(c(NA, NA)), c(NA_character_, NA_character_))
  expect_identical(invert_BF_bound_operator(NA), NA_character_)
  expect_identical(invert_BF_bound_operator(character()), character())
  expect_identical(invert_BF_bound_operator(invert_BF_bound_operator(operators)), unname(operators))

  # the internal passes other strings through; the export names them
  expect_identical(BayesTools:::.invert_BF_bound_operator(c("<", "x")), c(">", "x"))
  expect_error(
    invert_BF_bound_operator(c("<", "x", "y", "x")),
    "The 'x', 'y' values are not recognized by the 'operator' argument.",
    fixed = TRUE
  )
  expect_error(invert_BF_bound_operator(""), "The '' values are not recognized", fixed = TRUE)
  expect_error(invert_BF_bound_operator("<="), "The '<=' values are not recognized", fixed = TRUE)
  expect_error(invert_BF_bound_operator(1), "The 'operator' argument must be a character vector", fixed = TRUE)
  expect_error(invert_BF_bound_operator(NULL), "The 'operator' argument must be a character vector", fixed = TRUE)
  expect_error(invert_BF_bound_operator(c(NA, 1)), "The 'operator' argument must be a character vector", fixed = TRUE)
  expect_error(invert_BF_bound_operator(TRUE), "The 'operator' argument must be a character vector", fixed = TRUE)
})

test_that("posterior_log_ordinate_attribute() equals the internal producer", {

  value <- c(0, .5)
  log_ordinate <- c(-1000, -2)
  diagnostics <- list(BF_error_percent = 1.5)
  attribute <- posterior_log_ordinate_attribute(
    value, log_ordinate, "unit point part", "precomputed",
    diagnostics = diagnostics, parameter = "theta", conditional = "bias"
  )
  expect_identical(
    attribute,
    BayesTools:::.posterior_log_ordinate_attribute(
      value, log_ordinate, "unit point part", "precomputed",
      diagnostics = diagnostics, parameter = "theta", conditional = "bias"
    )
  )
  expect_identical(
    class(attribute),
    c("BayesTools_posterior_log_ordinate", "BayesTools_posterior_ordinate", "list")
  )
  expect_identical(
    names(attribute),
    c("status", "value", "ordinate", "log_ordinate", "method", "density_method",
      "diagnostics", "parameter", "conditional")
  )
  expect_identical(attribute$status, "ok")
  expect_identical(attribute$ordinate, exp(log_ordinate))
  expect_identical(attribute$ordinate[1L], 0)
  expect_identical(attribute$log_ordinate, log_ordinate)
  expect_identical(
    posterior_log_ordinate_attribute(value, log_ordinate, "m", "KDE"),
    BayesTools:::.posterior_log_ordinate_attribute(value, log_ordinate, "m", "KDE")
  )

  # the natural-height sibling refuses the underflowing height
  expect_error(
    posterior_ordinate_attribute(0, exp(-1000), "m", "precomputed"),
    "Posterior ordinates must be finite and positive.", fixed = TRUE
  )

  # the attribute is accepted as draw metadata and keeps its log heights
  draws <- c(.1, .2, .3)
  attr(draws, "parameter") <- "theta"
  posterior_metadata(draws, "posterior_ordinate") <- attribute
  expect_identical(posterior_metadata(draws, "posterior_ordinate")$log_ordinate, log_ordinate)
})

test_that("posterior_log_ordinate_attribute() validates like its natural-height sibling", {

  expect_error(posterior_log_ordinate_attribute(0, -1, 1, "precomputed"),
               "The 'method' argument must be a character vector.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, c("a", "b"), "precomputed"),
               "The 'method' argument must have length '1'.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, NA_character_, "precomputed"),
               "The 'method' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, "m", NULL),
               "The 'density_method' argument cannot be NULL.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, "m", NA_character_),
               "The 'density_method' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, "m", "precomputed", NULL, "unnamed"),
               "must have unique, nonmissing names", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, "m", "precomputed", NULL, a = 1, a = 2),
               "must have unique, nonmissing names", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -1, "m", "precomputed", NULL, a = 1, 2),
               "must have unique, nonmissing names", fixed = TRUE)
  for(reserved in c("status", "ordinate", "x", "y", "null_hypothesis", "height", "posterior_height",
                    "ordinates", "posterior_ordinate", "posterior_ordinates", "estimator")){
    arguments <- list(0, -1, "m", "precomputed")
    arguments[[reserved]] <- 1
    expect_error(do.call(posterior_log_ordinate_attribute, arguments),
                 "cannot replace reserved fields", fixed = TRUE)
  }

  # the checks of the internal producer remain
  expect_error(posterior_log_ordinate_attribute(c(0, 0), c(-1, -2), "m", "precomputed"),
               "Log posterior ordinates require unique finite values and finite log heights.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, c(-1, -2), "m", "precomputed"),
               "Log posterior ordinates require unique finite values and finite log heights.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(0, -Inf, "m", "precomputed"),
               "Log posterior ordinates require unique finite values and finite log heights.", fixed = TRUE)
  expect_error(posterior_log_ordinate_attribute(NA_real_, -1, "m", "precomputed"),
               "Log posterior ordinates require unique finite values and finite log heights.", fixed = TRUE)
})

test_that("posterior_measure_check() raises the declared typed refusals", {

  declared <- data.frame(
    column  = c("a", "a", "b", "b"),
    measure = c("atoms", "prior_density", "support", "prior_density"),
    reason  = c("atoms refused", "density refused", "support refused", "map refused"),
    cause   = c(NA, NA, NA, "nonlinear_map"),
    stringsAsFactors = FALSE
  )
  draws <- cbind(a = c(.1, .2, .3), b = c(1, 2, 3))
  posterior_metadata(draws, "measure_unavailable") <- declared

  # nothing declared, another measure, or another column pass
  checked <- withVisible(posterior_measure_check(draws, "atoms", "b"))
  expect_false(checked$visible)
  expect_true(checked$value)
  expect_true(posterior_measure_check(draws, "support", "a"))
  expect_true(posterior_measure_check(draws, "atoms", "c"))
  expect_true(posterior_measure_check(draws, "atoms", character()))
  expect_true(posterior_measure_check(c(.1, .2), "atoms"))
  expect_identical(
    posterior_measure_check(draws, "support", "a"),
    BayesTools:::.bt_formula_measure_check(draws, "support", "a")
  )

  refusal <- function(x, measure, column = NULL) tryCatch(posterior_measure_check(x, measure, column), error = identity)
  atoms <- refusal(draws, "atoms", "a")
  expect_identical(class(atoms), c("BayesTools_formula_atoms_unavailable", "error", "condition"))
  expect_identical(conditionMessage(atoms), "Formula atoms for 'a' are unavailable: atoms refused.")
  expect_null(conditionCall(atoms))
  expect_identical(atoms$target, "a")
  expect_identical(atoms$reason, "atoms refused")
  expect_identical(atoms$detail, "atoms refused")
  expect_null(atoms$diagnostics)
  expect_identical(
    atoms,
    tryCatch(BayesTools:::.bt_formula_measure_check(draws, "atoms", "a"), error = identity)
  )

  density <- refusal(draws, "prior_density", "a")
  expect_identical(
    class(density),
    c("BayesTools_formula_prior_density_unavailable", "BayesTools_formula_transform_unavailable", "error", "condition")
  )
  expect_identical(density$reason, "density refused")
  mapped <- refusal(draws, "prior_density", "b")
  expect_identical(
    class(mapped),
    c("BayesTools_formula_prior_density_unavailable", "BayesTools_formula_measure_unavailable",
      "BayesTools_formula_transform_unavailable", "error", "condition")
  )
  expect_identical(mapped$target, "b")
  expect_identical(mapped$reason, "nonlinear_map")
  expect_identical(mapped$detail, "map refused")

  support <- refusal(draws, "support", "b")
  expect_identical(class(support), c("BayesTools_formula_support_unavailable", "error", "condition"))
  expect_identical(support$target, "b")

  # without a column, the column names of matrix draws are checked
  expect_identical(class(refusal(draws, "atoms")), class(atoms))
  expect_identical(refusal(draws, "support")$target, "b")

  # vector draws check their declared parameter
  vector_draws <- c(.1, .2)
  attr(vector_draws, "parameter") <- "theta"
  posterior_metadata(vector_draws, "measure_unavailable") <- data.frame(
    column = "theta", measure = "support", reason = "not certified"
  )
  expect_identical(refusal(vector_draws, "support")$target, "theta")
  expect_true(posterior_measure_check(vector_draws, "atoms"))
})

test_that("posterior_measure_check() messages end in one full stop and agree in number", {

  # a stored reason that is already a sentence adds no second full stop; 'detail' keeps the raw reason
  draws <- c(.1, .2, .3)
  attr(draws, "parameter") <- "a"
  posterior_metadata(draws, "measure_unavailable") <- data.frame(
    column  = "a",
    measure = c("atoms", "support", "prior_density"),
    reason  = c("atoms refused.", "support refused.", "density refused."),
    stringsAsFactors = FALSE
  )
  refusal <- function(measure) tryCatch(posterior_measure_check(draws, measure), error = identity)

  atoms <- refusal("atoms")
  expect_identical(conditionMessage(atoms), "Formula atoms for 'a' are unavailable: atoms refused.")
  expect_identical(atoms$detail, "atoms refused.")
  support <- refusal("support")
  expect_identical(conditionMessage(support), "Formula support for 'a' is unavailable: support refused.")
  expect_identical(support$detail, "support refused.")
  density <- refusal("prior_density")
  expect_identical(conditionMessage(density), "Prior density for 'a' is unavailable: density refused.")
  expect_identical(density$detail, "density refused.")
  expect_identical(density$reason, "density refused.")
})

test_that("posterior_measure_check() validates its arguments", {

  draws <- c(.1, .2)
  expect_error(posterior_measure_check(draws, "density"),
               "The 'density' values are not recognized by the 'measure' argument.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, c("atoms", "support")),
               "The 'measure' argument must have length '1'.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, NA_character_),
               "The 'measure' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, NULL),
               "The 'measure' argument cannot be NULL.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, 1),
               "The 'measure' argument must be a character vector.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, "atoms", 1),
               "The 'column' argument must be a character vector.", fixed = TRUE)
  expect_error(posterior_measure_check(draws, "atoms", NA_character_),
               "The 'column' argument cannot contain NA/NaN values.", fixed = TRUE)
})

test_that("JAGS_formula_point_finalize() equals the internal point finalization", {

  data <- data.frame(x = c(-2, -1, 1, 2))
  priors <- list(
    intercept = prior("point", list(0)),
    x         = prior("point", list(location = expression(2 * theta + model_c[2])))
  )
  output <- JAGS_formula(~ x, "mu", data, priors)
  owning <- c(output$prior_list, list(theta = prior("normal", list(0, 1))))
  model_data <- list(model_c = c(10, 20))
  parameters <- c("theta", "mu_x")

  finalized <- JAGS_formula_point_finalize(
    output$formula_design, owning, data, model_data, parameters
  )
  expect_identical(
    finalized,
    BayesTools:::.bt_formula_point_finalize(
      output$formula_design, owning, data, model_data, parameters
    )
  )
  owner <- finalized$point_expression_owner
  expect_identical(owner$targets, "mu_x")
  expect_setequal(names(owner$prior_list), c("mu_x", "theta"))
  expect_identical(owner$extra_parents, character())
  expect_identical(names(owner$data), "model_c")
  expect_identical(names(owner$points), "mu_x")
  expect_null(output$formula_design$point_expression_owner)

  # the defaults of the optional arguments
  expect_identical(
    JAGS_formula_point_finalize(output$formula_design, owning, data),
    BayesTools:::.bt_formula_point_finalize(output$formula_design, owning, data)
  )
  expect_identical(
    JAGS_formula_point_finalize(output$formula_design, owning, NULL, NULL, NULL),
    BayesTools:::.bt_formula_point_finalize(output$formula_design, owning, NULL, NULL, character())
  )

  # a design without expression-point coefficients is returned unchanged
  plain <- JAGS_formula(~ x, "mu", data, list(
    intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))
  ))
  expect_identical(
    JAGS_formula_point_finalize(plain$formula_design, plain$prior_list, data),
    plain$formula_design
  )

  # the typed refusals of the internal finalization remain
  expect_error(
    JAGS_formula_point_finalize(output$formula_design, unname(owning), data, model_data, parameters),
    class = "BayesTools_formula_point_unavailable"
  )
  expect_error(
    JAGS_formula_point_finalize(output$formula_design, owning, list(model_c = c(1, 2)), model_data, parameters),
    class = "BayesTools_formula_point_unavailable"
  )
})

test_that("JAGS_formula_point_finalize() validates its arguments", {

  data <- data.frame(x = c(-2, -1, 1, 2))
  output <- JAGS_formula(~ x, "mu", data, list(
    intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))
  ))
  design <- output$formula_design
  finalize <- function(...) JAGS_formula_point_finalize(...)

  for(invalid in list(NULL, list(), unclass(design), "design", data)){
    expect_error(finalize(invalid, output$prior_list, data),
                 "The 'design' argument must be a formula design returned by 'JAGS_formula()'.", fixed = TRUE)
  }
  expect_error(finalize(design, NULL, data), "The 'prior_list' argument cannot be NULL.", fixed = TRUE)
  expect_error(finalize(design, list(), data), "The 'prior_list' argument cannot be NULL.", fixed = TRUE)
  expect_error(finalize(design, 1, data), "The 'prior_list' argument must be a list.", fixed = TRUE)
  expect_error(finalize(design, output$prior_list, 1),
               "The 'formula_data' argument must be NULL, a data frame, or a named list.", fixed = TRUE)
  expect_error(finalize(design, output$prior_list, "data"),
               "The 'formula_data' argument must be NULL, a data frame, or a named list.", fixed = TRUE)
  expect_error(finalize(design, output$prior_list, data, model_data = 1),
               "The 'model_data' argument must be a list.", fixed = TRUE)
  expect_error(finalize(design, output$prior_list, data, parameter_names = 1),
               "The 'parameter_names' argument must be a character vector.", fixed = TRUE)
  expect_error(finalize(design, output$prior_list, data, parameter_names = c("a", NA)),
               "The 'parameter_names' argument cannot contain NA/NaN values.", fixed = TRUE)
})

test_that("JAGS_formula_scale_finalize() equals the internal scale finalization", {

  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "sigma"
  data <- data.frame(x = c(10, 20, 30, 45))
  output <- JAGS_formula(
    ~ x, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = slope),
    formula_scale = TRUE
  )
  design <- JAGS_formula_point_finalize(output$formula_design, output$prior_list, data)
  expect_false(isTRUE(attr(output$formula_scale, "unscale_design")$complete))

  scale <- JAGS_formula_scale_finalize(
    output$formula_scale, design, output$prior_list, list(sigma = 3)
  )
  expect_identical(
    scale,
    BayesTools:::.bt_formula_scale_finalize(
      output$formula_scale, design, prior_list = output$prior_list,
      model_data = list(sigma = 3), owner_scope = "fit"
    )
  )
  specification <- attr(scale, "unscale_design")
  expect_identical(specification$state_constants, c(sigma = 3))
  expect_identical(specification$owner_scope, "fit")
  expect_true(specification$complete)

  prior_only <- JAGS_formula_scale_finalize(
    output$formula_scale, design, output$prior_list, list(sigma = 3), owner_scope = "prior_only"
  )
  expect_identical(attr(prior_only, "unscale_design")$owner_scope, "prior_only")
  expect_identical(
    prior_only,
    BayesTools:::.bt_formula_scale_finalize(
      output$formula_scale, design, output$prior_list, list(sigma = 3), "prior_only"
    )
  )

  # the scalar state also comes from a point prior of the model
  from_prior <- JAGS_formula_scale_finalize(
    output$formula_scale, design, c(output$prior_list, list(sigma = prior("point", list(3))))
  )
  expect_identical(attr(from_prior, "unscale_design")$state_constants, c(sigma = 3))

  # the defaults and the empty scale
  expect_identical(
    JAGS_formula_scale_finalize(output$formula_scale),
    BayesTools:::.bt_formula_scale_finalize(output$formula_scale)
  )
  expect_null(JAGS_formula_scale_finalize(NULL))
  unscaled <- JAGS_formula(~ x, "mu", data, list(
    intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))
  ))
  expect_identical(
    JAGS_formula_scale_finalize(unscaled$formula_scale, unscaled$formula_design, unscaled$prior_list),
    BayesTools:::.bt_formula_scale_finalize(unscaled$formula_scale, unscaled$formula_design, unscaled$prior_list)
  )

  # contradictory owners of one state keep the typed refusal
  contradictory <- tryCatch(
    JAGS_formula_scale_finalize(
      output$formula_scale, design,
      c(output$prior_list, list(sigma = prior("point", list(2)))), list(sigma = 3)
    ),
    error = identity
  )
  expect_s3_class(contradictory, "BayesTools_formula_transform_unavailable")
  expect_identical(contradictory$reason, "contradictory_state_constants")
})

test_that("JAGS_formula_scale_finalize() validates its arguments", {

  data <- data.frame(x = c(10, 20, 30, 45))
  output <- JAGS_formula(~ x, "mu", data, list(
    intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))
  ), formula_scale = TRUE)
  design <- output$formula_design
  scale <- output$formula_scale
  finalize <- function(...) JAGS_formula_scale_finalize(...)

  expect_error(finalize(1), "The 'formula_scale' argument must be a list.", fixed = TRUE)
  expect_error(finalize("scale"), "The 'formula_scale' argument must be a list.", fixed = TRUE)
  for(invalid in list(list(), unclass(design), "design", data)){
    expect_error(finalize(scale, invalid),
                 "The 'design' argument must be a formula design returned by 'JAGS_formula()' or NULL.", fixed = TRUE)
  }
  expect_error(finalize(scale, design, 1), "The 'prior_list' argument must be a list.", fixed = TRUE)
  expect_error(finalize(scale, design, list(), 1), "The 'model_data' argument must be a list.", fixed = TRUE)
  expect_error(finalize(scale, design, owner_scope = "compiler"),
               "The 'compiler' values are not recognized by the 'owner_scope' argument.", fixed = TRUE)
  expect_error(finalize(scale, design, owner_scope = c("fit", "prior_only")),
               "The 'owner_scope' argument must have length '1'.", fixed = TRUE)
  expect_error(finalize(scale, design, owner_scope = NA_character_),
               "The 'owner_scope' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(finalize(scale, design, owner_scope = TRUE),
               "The 'owner_scope' argument must be a character vector.", fixed = TRUE)
})

.robma_facing_random_term <- function(formula, data, ...){

  JAGS_formula(
    formula, "mu", data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(...)
  )$formula_design$random_effects[[1L]]
}

test_that("random group helpers equal the internal chain on compiled terms", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  fitted <- data.frame(g = c("a", "b", "c", "a"))
  term <- .robma_facing_random_term(~ 1 + diag(1 | g), fitted, sd = sd_prior)

  expect_identical(random_group_component_names(term), "g")
  expect_identical(random_group_component_names(term), BayesTools:::.bt_random_group_component_names(term))

  new_data <- data.frame(g = c("b", "z", "a", "b", "z"))
  keys <- random_group_prediction_keys(term, new_data)
  expect_identical(
    keys,
    BayesTools:::.bt_random_group_prediction_tuple_keys(
      term, BayesTools:::.bt_random_group_observations(term, new_data)
    )
  )
  expect_type(keys, "character")
  expect_length(keys, nrow(new_data))
  # fitted groups carry their fitted key, a new group a key of its own
  fitted_keys <- term$group_tuple_keys[match(c("b", "a"), term$group_levels)]
  expect_identical(keys[c(1L, 3L)], fitted_keys)
  expect_identical(keys[4L], keys[1L])
  expect_identical(keys[5L], keys[2L])
  expect_false(keys[2L] %in% term$group_tuple_keys)
  expect_identical(
    random_group_prediction_keys(term, fitted),
    term$group_tuple_keys[match(fitted$g, term$group_levels)]
  )
  # a factor column and a data frame with other columns give the same keys
  expect_identical(
    random_group_prediction_keys(term, data.frame(other = 1:3, g = factor(c("c", "a", "b")))),
    term$group_tuple_keys[match(c("c", "a", "b"), term$group_levels)]
  )

  # the stored design of draws carries the same term
  draws <- JAGS_formula_draws(
    cbind(mu_intercept = c(1, 2)), ~ 1 + diag(1 | g), "mu", fitted,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = sd_prior)
  )
  stored <- JAGS_formula_design(draws, "mu")$random_effects[[1L]]
  expect_identical(random_group_component_names(stored), "g")
  expect_identical(random_group_prediction_keys(stored, new_data), keys)

  # grouping by several variables
  tuples <- data.frame(
    g1 = factor(rep(c("b", "a", "b", "a"), 2), levels = c("a", "b")),
    g2 = factor(rep(c("y", "y", "x", "x"), 2), levels = c("x", "y"))
  )
  prediction <- data.frame(
    g1 = as.character(tuples$g1[c(4L, 1L, 2L)]),
    g2 = c("x", "y", "w")
  )
  check_tuple_term <- function(block){

    expect_identical(random_group_component_names(block), c("g1", "g2"))
    expect_identical(random_group_component_names(block), BayesTools:::.bt_random_group_component_names(block))
    block_keys <- random_group_prediction_keys(block, prediction)
    expect_identical(
      block_keys,
      BayesTools:::.bt_random_group_prediction_tuple_keys(
        block, BayesTools:::.bt_random_group_observations(block, prediction)
      )
    )
    expect_identical(
      block_keys[1:2],
      block$group_tuple_keys[match(c("a:x", "b:y"), block$group_levels)]
    )
    expect_false(block_keys[3L] %in% block$group_tuple_keys)
  }
  check_tuple_term(.robma_facing_random_term(~ 1 + diag(1 | g1:g2), tuples, g1_g2 = random_block(sd = sd_prior)))
  check_tuple_term(.robma_facing_random_term(~ 1 + diag(1 | g1 / g2), tuples, sd = sd_prior))
})

test_that("random_group_prediction_keys() keeps the matching rules of the fitted levels", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  adjacent <- c(.3, .3 + .Machine$double.eps)
  term <- .robma_facing_random_term(~ 1 + diag(1 | g), data.frame(g = adjacent), sd = sd_prior)
  expect_identical(
    random_group_prediction_keys(term, data.frame(g = adjacent[2:1])),
    term$group_tuple_keys[2:1]
  )
  expect_identical(
    random_group_prediction_keys(term, data.frame(g = term$group_levels[2L])),
    term$group_tuple_keys[2L]
  )
  expect_error(
    random_group_prediction_keys(term, data.frame(g = "0.3")),
    "Random-effect grouping value '0.3' for component 'g' is ambiguous", fixed = TRUE
  )
  new_group <- random_group_prediction_keys(term, data.frame(g = .3 + 2 * .Machine$double.eps))
  expect_false(new_group %in% term$group_tuple_keys)
})

test_that("random group helpers validate the term and the data", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  term <- .robma_facing_random_term(~ 1 + diag(1 | g), data.frame(g = c("a", "b")), sd = sd_prior)

  expect_error(random_group_component_names(NULL), "The 'term' argument cannot be NULL.", fixed = TRUE)
  expect_error(random_group_component_names(list()), "The 'term' argument cannot be NULL.", fixed = TRUE)
  expect_error(random_group_component_names("g"), "The 'term' argument must be a list.", fixed = TRUE)
  expect_error(random_group_component_names(1), "The 'term' argument must be a list.", fixed = TRUE)
  expect_error(random_group_component_names(list(group_label = "g")),
               "The 'term' argument must be a compiled random-effect term with a grouping expression.", fixed = TRUE)
  expect_error(random_group_prediction_keys(NULL, data.frame(g = "a")),
               "The 'term' argument cannot be NULL.", fixed = TRUE)
  expect_error(random_group_prediction_keys("g", data.frame(g = "a")),
               "The 'term' argument must be a list.", fixed = TRUE)

  # the ownership validator runs first
  expect_error(random_group_prediction_keys(list(group_expr = quote(g)), "not data"),
               class = "BayesTools_refit_required")
  stale <- term
  stale$group_component_identity_keys <- NULL
  expect_error(random_group_prediction_keys(stale, data.frame(g = "a")),
               class = "BayesTools_refit_required")
  stale <- term
  stale$group_components <- "h"
  expect_error(random_group_prediction_keys(stale, data.frame(g = "a")),
               class = "BayesTools_refit_required")

  expect_error(random_group_prediction_keys(term, list(g = "a")),
               "The 'data' argument must be a data frame.", fixed = TRUE)
  expect_error(random_group_prediction_keys(term, matrix("a", dimnames = list(NULL, "g"))),
               "The 'data' argument must be a data frame.", fixed = TRUE)
  expect_error(random_group_prediction_keys(term, NULL),
               "The 'data' argument must be a data frame.", fixed = TRUE)

  # the checks of the internal observations remain
  expect_error(random_group_prediction_keys(term, data.frame(h = "a")),
               "The 'g' random-effect grouping variable is missing in the data set.", fixed = TRUE)
  expect_error(random_group_prediction_keys(term, data.frame(g = c("a", NA))),
               "Random-effect grouping expression 'g' must not contain missing values.", fixed = TRUE)
})

test_that("hypothesis_linear_target_side() equals the internal side reading", {

  statement <- function(text) hypothesis_parse(text)$statements[[1L]]

  region <- statement("f[A] - 2*f[B] > 0.5")
  reading <- hypothesis_linear_target_side(region$left, "f")
  expect_identical(reading, BayesTools:::.hypothesis_linear_target_side(region$left, "f"))
  expect_identical(names(reading), c("operator", "coefficients", "numerator_coefficients", "value", "divisor"))
  expect_identical(reading$operator, ">")
  expect_identical(reading$coefficients, c(A = 1, B = -2))
  expect_identical(reading$numerator_coefficients, c(A = 1, B = -2))
  expect_identical(reading$value, 0.5)
  expect_identical(reading$divisor, 1)
  # the implicit complement side reads the same combination
  complement <- hypothesis_linear_target_side(region$right, "f")
  expect_identical(complement, BayesTools:::.hypothesis_linear_target_side(region$right, "f"))
  expect_identical(complement$operator, "<=")
  expect_identical(complement$coefficients, reading$coefficients)
  expect_identical(complement$value, reading$value)

  # the divisor scales the coefficients, not the numerator
  point <- statement("(f[A] + f[B])/2 = 0")
  halves <- hypothesis_linear_target_side(point$left, "f")
  expect_identical(halves, BayesTools:::.hypothesis_linear_target_side(point$left, "f"))
  expect_identical(halves$operator, "=")
  expect_identical(halves$coefficients, c(A = .5, B = .5))
  expect_identical(halves$numerator_coefficients, c(A = 1, B = 1))
  expect_identical(halves$value, 0)
  expect_identical(halves$divisor, 2)
  expect_identical(hypothesis_linear_target_side(point$right, "f")$operator, "!=")

  # constants move to the value, levels sort by name and cancel to zero
  moved <- hypothesis_linear_target_side(statement("f[B] + 1 < 2*f[A] + 3 - f[B]")$left, "f")
  expect_identical(moved, BayesTools:::.hypothesis_linear_target_side(
    statement("f[B] + 1 < 2*f[A] + 3 - f[B]")$left, "f"
  ))
  expect_identical(moved$operator, "<")
  expect_identical(moved$coefficients, c(A = -2, B = 2))
  expect_identical(moved$value, 2)
  expect_identical(
    hypothesis_linear_target_side(statement("f[A] + f[B] - f[B] != 1")$left, "f")$coefficients,
    c(A = 1)
  )

  # a side of another parameter is not a target of this parameter
  expect_error(
    hypothesis_linear_target_side(statement("g[A] > 0")$left, "f"),
    "A linear target may reference levels of only 'f'.", fixed = TRUE
  )
  expect_error(
    hypothesis_linear_target_side(statement("f[A]^2 > 0")$left, "f"),
    "A linear target must be a linear combination of levels of 'f' and numbers.", fixed = TRUE
  )
  expect_error(
    hypothesis_linear_target_side(statement("f[A] > 0 & f[B] > 0")$left, "f"),
    "Linear targets support only point and simple region statements.", fixed = TRUE
  )
})

test_that("hypothesis_linear_target_side() validates its arguments", {

  side <- hypothesis_parse("f[A] > 0")$statements[[1L]]$left

  expect_error(hypothesis_linear_target_side(side, ""),
               "The 'parameter' argument cannot be an empty string.", fixed = TRUE)
  expect_error(hypothesis_linear_target_side(side, NULL),
               "The 'parameter' argument cannot be NULL.", fixed = TRUE)
  expect_error(hypothesis_linear_target_side(side, c("f", "g")),
               "The 'parameter' argument must have length '1'.", fixed = TRUE)
  expect_error(hypothesis_linear_target_side(side, NA_character_),
               "The 'parameter' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(hypothesis_linear_target_side(side, 1),
               "The 'parameter' argument must be a character vector.", fixed = TRUE)

  malformed <- "Hypothesis AST contains a malformed side. Parse the hypothesis again."
  expect_error(hypothesis_linear_target_side(NULL, "f"), malformed, fixed = TRUE)
  expect_error(hypothesis_linear_target_side("f[A] > 0", "f"), malformed, fixed = TRUE)
  expect_error(hypothesis_linear_target_side(side$expression, "f"), malformed, fixed = TRUE)
  expect_error(hypothesis_linear_target_side(unclass(side)[c("type", "expression")], "f"), malformed, fixed = TRUE)
  altered <- side
  altered$type <- "other"
  expect_error(hypothesis_linear_target_side(altered, "f"), malformed, fixed = TRUE)
  stale <- side
  stale$source <- "f[B] > 0"
  expect_error(hypothesis_linear_target_side(stale, "f"),
               "Hypothesis AST side source is inconsistent with its expression.", fixed = TRUE)
})

test_that("prior_ordered_total_name() names the ordered total node", {

  expect_identical(prior_ordered_total_name("mu_f"), "mu_f_ordered_total")
  names <- c("mu_f", "g", "a b")
  expect_identical(prior_ordered_total_name(names), paste0(names, "_ordered_total"))
  expect_identical(prior_ordered_total_name(names), BayesTools:::.prior_ordered_total_name(names))

  # the name is the one of the total node of the generated model; the
  # inclusion indicator of a mixture total appends '_indicator'
  data <- data.frame(f = ordered(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")))
  compile <- function(total){
    JAGS_formula(~ f, "mu", data, list(
      intercept = prior("normal", list(0, 1)),
      f         = prior_ordered(total)
    ))$prior_list
  }
  simple <- compile(prior("normal", list(1, 1)))
  expect_true(prior_ordered_total_name("mu_f") %in% JAGS_to_monitor(simple))
  expect_true(grepl(prior_ordered_total_name("mu_f"), JAGS_add_priors("model{}", simple), fixed = TRUE))
  mixture <- compile(prior_mixture(
    list(prior("spike", list(0)), prior("normal", list(1, 1))), is_null = c(TRUE, FALSE)
  ))
  expect_true(paste0(prior_ordered_total_name("mu_f"), "_indicator") %in% JAGS_to_monitor(mixture))

  expect_error(prior_ordered_total_name(NULL), "The 'parameter_name' argument cannot be NULL.", fixed = TRUE)
  expect_error(prior_ordered_total_name(character()), "The 'parameter_name' argument cannot be NULL.", fixed = TRUE)
  expect_error(prior_ordered_total_name(1), "The 'parameter_name' argument must be a character vector.", fixed = TRUE)
  expect_error(prior_ordered_total_name(NA_character_), "The 'parameter_name' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(prior_ordered_total_name(c("a", NA)), "The 'parameter_name' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(prior_ordered_total_name(""), "The 'parameter_name' argument cannot contain empty strings.", fixed = TRUE)
  expect_error(prior_ordered_total_name(c("a", "")), "The 'parameter_name' argument cannot contain empty strings.", fixed = TRUE)
})

test_that("simplify_as_mixed_posterior_bias() extracts the scalar bias components", {

  prior_list <- list(bias = prior_mixture(list(
    prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5))),
    prior_PET("normal", list(0, 1)), prior_PEESE("normal", list(0, 1))
  )))
  posterior <- cbind(
    bias_indicator = rep(1:3, each = 6),
    "omega[1]" = 1,
    "omega[2]" = c(rep(.5, 6), rep(1, 12)),
    PET = c(rep(0, 6), seq(.1, 1, length.out = 6), rep(0, 6)),
    PEESE = c(rep(0, 12), seq(.2, 2, length.out = 6))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  fit <- attach_test_parameter_map(fit)
  samples <- as_mixed_posteriors(fit, "bias")
  samples$other <- c(10, 20)

  for(parameter in c("PET", "PEESE", "omega")){
    scalar <- simplify_as_mixed_posterior_bias(samples, parameter)
    expect_equal(scalar, BayesTools:::.simplify_as_mixed_posterior_bias(samples, parameter))
    expect_setequal(names(scalar), c(parameter, "other"))
    expect_identical(scalar$other, c(10, 20))
    expect_false(inherits(scalar[[parameter]], "mixed_posteriors.bias"))
    expect_identical(nrow(scalar[[parameter]]), nrow(samples$bias))
    expect_identical(length(attr(scalar[[parameter]], "prior_list")), 3L)
  }
  for(parameter in c("PET", "PEESE")){
    scalar <- simplify_as_mixed_posterior_bias(samples, parameter)[[parameter]]
    expect_identical(colnames(scalar), parameter)
    expect_equal(as.numeric(scalar), unname(posterior[, parameter]))
    expect_identical(colnames(BayesTools:::.posterior_atoms_get(scalar)$locations), parameter)
    # the branches without the component are null values with their weights
    expect_true(all(vapply(attr(scalar, "prior_list")[c(1L, if(parameter == "PET") 3L else 2L)],
                           is.prior.point, logical(1))))
  }
  omega <- simplify_as_mixed_posterior_bias(samples, "omega")$omega
  omega_columns <- grep("omega", colnames(samples$bias), value = TRUE)
  expect_length(omega_columns, 2L)
  expect_identical(colnames(omega), omega_columns)
  expect_s3_class(omega, "mixed_posteriors.weightfunction")
  expect_equal(as.numeric(omega[, 1L]), unname(posterior[, "omega[1]"]))
  expect_equal(as.numeric(omega[, 2L]), unname(posterior[, "omega[2]"]))

  # the typed refusal of the internal extraction remains
  incomplete <- samples
  incomplete$bias <- incomplete$bias[, colnames(incomplete$bias) != "PET", drop = FALSE]
  for(attribute in setdiff(names(attributes(samples$bias)), c("dim", "dimnames"))){
    attr(incomplete$bias, attribute) <- attributes(samples$bias)[[attribute]]
  }
  incomplete$bias <- BayesTools:::.bt_meta_refresh(incomplete$bias)
  expect_error(
    simplify_as_mixed_posterior_bias(incomplete, "PET"),
    "Posterior samples for 'PET' are unavailable because an active bias branch is not a point prior.",
    fixed = TRUE
  )

  # the arguments are validated before the extraction
  expect_error(simplify_as_mixed_posterior_bias(samples, "PETPEESE"),
               "The 'PETPEESE' values are not recognized by the 'parameter' argument.", fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples, "bias"),
               "The 'bias' values are not recognized by the 'parameter' argument.", fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples, c("PET", "PEESE")),
               "The 'parameter' argument must have length '1'.", fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples, NA_character_),
               "The 'parameter' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples, NULL),
               "The 'parameter' argument cannot be NULL.", fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples, 1),
               "The 'parameter' argument must be a character vector.", fixed = TRUE)
  bias_message <- "The 'samples' argument must be a list with a 'bias' element of class 'mixed_posteriors.bias'"
  expect_error(simplify_as_mixed_posterior_bias(NULL, "PET"), bias_message, fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(samples$bias, "PET"), bias_message, fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(list(), "PET"), bias_message, fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(list(bias = 1:3), "PET"), bias_message, fixed = TRUE)
  expect_error(simplify_as_mixed_posterior_bias(list(PET = samples$bias), "PET"), bias_message, fixed = TRUE)
})

.robma_facing_prior_density <- function(){

  BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    )),
    weights = c(beta = 1),
    n_grid  = 256
  )
}

test_that("plot() of a prior linear density draws its plot data", {

  density <- .robma_facing_prior_density()
  expect_s3_class(density, "prior_linear_density")
  expect_false(inherits(density, "prior"))
  expect_identical(
    utils::getS3method("plot", "prior_linear_density"),
    BayesTools:::plot.prior_linear_density
  )
  # plot() of a prior keeps its own method
  expect_identical(
    utils::getS3method("plot", "prior"),
    BayesTools:::plot.prior
  )
  expect_null(utils::getS3method("plot", "prior_linear_density_deferred", optional = TRUE))

  # ggplot: the layers of the internal renderer on the same plot data
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(density, n_points = 64, x_range = c(-3, 3))
  direct <- BayesTools:::.plot_prior_list.both(plot_data, plot_type = "ggplot", par_name = "beta")
  plotted <- plot(density, plot_type = "ggplot", n_points = 64, xlim = c(-3, 3), par_name = "beta")
  expect_s3_class(plotted, "ggplot")
  geoms <- function(plot) unname(vapply(plot$layers, function(layer) class(layer$geom)[[1L]], character(1)))
  expect_identical(geoms(plotted), c("GeomLine", "GeomSegment"))
  expect_identical(geoms(plotted), geoms(direct))
  expect_equal(lapply(plotted$layers, `[[`, "data"), lapply(direct$layers, `[[`, "data"))
  curve <- plotted$layers[[1L]]$data
  expect_equal(curve$y, .5 * stats::dnorm(curve$x), tolerance = 1e-3)
  expect_equal(range(curve$x), c(-3, 3))
  atom <- plotted$layers[[2L]]$data
  expect_equal(atom$x, 0)
  expect_gt(atom$yend, 0)

  # the transformation arguments and the graphical arguments reach the renderer
  transformed <- plot(density, plot_type = "ggplot", n_points = 64, xlim = c(0, 4),
                      transformation = "exp", transformation_settings = TRUE, ylim = c(0, 1))
  transformed_data <- BayesTools:::.prior_linear_density_to_plot_data(
    density, n_points = 64, x_range = c(0, 4), transformation = "exp", transformation_settings = TRUE
  )
  expect_equal(
    lapply(transformed$layers, `[[`, "data"),
    lapply(BayesTools:::.plot_prior_list.both(transformed_data, "ggplot", ylim = c(0, 1))$layers, `[[`, "data")
  )
  expect_s3_class(
    plot(density, plot_type = "ggplot", n_points = 64, xlim = c(-3, 3), add = TRUE),
    "BayesTools_prior_overlay"
  )
  expect_s3_class(plot(BayesTools:::.prior_linear_density_point(2), plot_type = "ggplot"), "ggplot")

  # base graphics: invisible NULL, one device, and overlays keep the axis scale
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- withVisible(plot(density, n_points = 64, xlim = c(-3, 3)))
  expect_false(result$visible)
  expect_null(result$value)
  expect_silent(plot(density, plot_type = "base", n_points = 64, xlim = c(-3, 3), add = TRUE, col = "red"))
  expect_false(withVisible(plot(density, n_points = 64, par_name = "beta", lwd = 2))$visible)
})

test_that("plot() of a prior linear density validates its arguments", {

  density <- .robma_facing_prior_density()

  expect_error(plot(density, plot_type = "plotly"),
               "The 'plotly' values are not recognized by the 'plot_type' argument.", fixed = TRUE)
  expect_error(plot(density, plot_type = c("base", "ggplot")),
               "The 'plot_type' argument must have length '1'.", fixed = TRUE)
  expect_error(plot(density, n_points = 1), "The 'n_points' must be equal or higher than 2.", fixed = TRUE)
  expect_error(plot(density, n_points = 10.5), "The 'n_points' argument must be an integer vector.", fixed = TRUE)
  expect_error(plot(density, n_points = c(10, 20)), "The 'n_points' argument must have length '1'.", fixed = TRUE)
  expect_error(plot(density, xlim = 1), "The 'xlim' argument must have length '2'.", fixed = TRUE)
  expect_error(plot(density, xlim = c("a", "b")), "The 'xlim' argument must be a numeric vector.", fixed = TRUE)
  expect_error(plot(density, transformation = 3), "Uknown format of the 'transformation' argument.", fixed = TRUE)

  # a density without a curve and without a point mass has nothing to draw
  testthat::local_mocked_bindings(
    .prior_linear_density_to_plot_data = function(...) list(),
    .package = "BayesTools"
  )
  expect_error(plot(density, plot_type = "ggplot"),
               "The prior density has neither a continuous curve nor a point mass to plot.", fixed = TRUE)
})
