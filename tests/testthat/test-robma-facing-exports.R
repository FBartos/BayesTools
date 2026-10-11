skip_if_not_test_profile("unit")

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
