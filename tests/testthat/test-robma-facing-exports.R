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
