# ============================================================================ #
# TEST FILE: Input Validation Tests for Tools
# ============================================================================ #
#
# PURPOSE:
#   Tests for input validation functions (check_bool, check_char, check_real,
#   check_int, check_list) in R/tools.R
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#   - Tests the check_* functions exported from BayesTools
#
# SKIP CONDITIONS:
#   - None (fast, pure R tests)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @input-validation, @fast
# ============================================================================ #


test_that("check_bool validates logical inputs", {


  # Valid inputs

  expect_null(check_bool(NULL,  "", allow_NULL = TRUE))
  expect_null(check_bool(TRUE,  ""))
  expect_null(check_bool(FALSE, ""))
  expect_null(check_bool(as.logical(stats::rbinom(5, 1, .5)), "", check_length = 0))
  expect_null(check_bool(c(FALSE, FALSE), "", check_length = 2))
  expect_null(check_bool(NA,  ""))

  # Invalid type: matrix
  expect_error(
    check_bool(as.matrix(as.logical(rbinom(5, 1, .5))), "test object"),
    "The 'test object' argument must be a logical vector."
  )
  # Invalid type: string
  expect_error(
    check_bool("string", "test object"),
    "The 'test object' argument must be a logical vector."
  )
  # Invalid type: numeric
  expect_error(
    check_bool(1, "test object"),
    "The 'test object' argument must be a logical vector."
  )
  # Invalid type: list
  expect_error(
    check_bool(list(TRUE), "test object"),
    "The 'test object' argument must be a logical vector."
  )
  # Invalid length
  expect_error(
    check_bool(TRUE, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  # NULL not allowed
  expect_error(
    check_bool(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
  # NA not allowed
  expect_error(
    check_bool(NA, "test object", allow_NA = FALSE),
    "The 'test object' argument cannot contain NA/NaN values."
  )
})


test_that("check_char validates character inputs", {

  # Valid inputs
  expect_null(check_char(NULL,  "", allow_NULL = TRUE))
  expect_null(check_char("string",  ""))
  expect_null(check_char(c("string", "string1"),  "", check_length = 0))
  expect_null(check_char(as.character(stats::rbinom(5, 1, .5)), "", check_length = 5))
  expect_null(check_char(c("string", "string1"),  "", check_length = 0, allow_values = c("string", "string1")))
  expect_null(check_char(c(NA, ""),  "", check_length = 2))

  # Invalid type: matrix
  expect_error(
    check_char(as.matrix(as.logical(as.character(5, 1, .5))), "test object"),
    "The 'test object' argument must be a character vector."
  )
  # Invalid type: logical
  expect_error(
    check_char(TRUE, "test object"),
    "The 'test object' argument must be a character vector."
  )
  # Invalid type: numeric
  expect_error(
    check_char(1, "test object"),
    "The 'test object' argument must be a character vector."
  )
  # Invalid type: list
  expect_error(
    check_char(list("string"), "test object"),
    "The 'test object' argument must be a character vector."
  )
  # Invalid length
  expect_error(
    check_char("string", "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  # Invalid allowed values
  expect_error(
    check_char(c("string", "string1"),  "test object", check_length = 0, allow_values = c("string")),
    "The 'string1' values are not recognized by the 'test object' argument."
  )
  # NULL not allowed
  expect_error(
    check_char(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
  # NA not allowed
  expect_error(
    check_char(c("a", NA), "test object", allow_NA = FALSE, check_length = FALSE),
    "The 'test object' argument cannot contain NA/NaN values."
  )
})


test_that("check_real validates numeric inputs", {

  # Valid inputs
  expect_null(check_real(NULL,  "", allow_NULL = TRUE))
  expect_null(check_real(pi,  ""))
  expect_null(check_real(c(pi, 2),  "", check_length = 0))
  expect_null(check_real(stats::rnorm(4, 1, .5), "", check_length = 4))
  expect_null(check_real(stats::rgamma(1, 1, 1), "", lower = 0))
  expect_null(check_real(stats::rbeta(1, 1, 1), "",  upper = 1))
  expect_null(check_real(c(0, 1), "", lower = 0, upper = 1, check_length = 2))
  expect_null(check_real(c(NA, NaN),  "", check_length = 2))

  # Invalid type: matrix
  expect_error(
    check_real(as.matrix(stats::rnorm(4, 1, .5)), "test object", check_length = FALSE),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: logical
  expect_error(
    check_real(TRUE, "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: string
  expect_error(
    check_real("string", "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: list
  expect_error(
    check_real(list(3.2), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid length
  expect_error(
    check_real(1, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  # Upper bound violation
  expect_error(
    check_real(stats::rgamma(1, 1, 1), "test object", upper = 0),
    "The 'test object' must be equal or lower than 0."
  )
  # Lower bound violation
  expect_error(
    check_real(stats::rbeta(1, 1, 1), "test object",  lower = 1),
    "The 'test object' must be equal or higher than 1."
  )
  # Boundary not allowed (lower)
  expect_error(
    check_real(0, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be higher than 0."
  )
  # Boundary not allowed (upper)
  expect_error(
    check_real(1, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be lower than 1."
  )
  # NULL not allowed
  expect_error(
    check_real(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
  # NA not allowed
  expect_error(
    check_real(NaN, "test object", allow_NA = FALSE),
    "The 'test object' argument cannot contain NA/NaN values."
  )
})


test_that("check_int validates integer inputs", {

  # Valid inputs
  expect_null(check_int(NULL,  "", allow_NULL = TRUE))
  expect_null(check_int(0,  ""))
  expect_null(check_int(c(-1, 2),  "", check_length = 0))
  expect_null(check_int(stats::rpois(4, 1), "", check_length = 4))
  expect_null(check_int(c(0, 5, 10), "", lower = 0, check_length = 3))
  expect_null(check_int(c(0, -2, 2), "", upper = 2, check_length = 3))
  expect_null(check_int(c(-3, -1), "", lower = -3, upper = -1, check_length = 2))
  expect_null(check_int(c(NA, NaN),  "", check_length = 2))

  # Invalid type: matrix
  expect_error(
    check_int(as.matrix(stats::rpois(4, 1)), "test object", check_length = FALSE),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: logical
  expect_error(
    check_int(TRUE, "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: string
  expect_error(
    check_int("string", "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid type: list
  expect_error(
    check_int(list(3.2), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  # Invalid length
  expect_error(
    check_int(1, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  # Upper bound violation
  expect_error(
    check_int(1, "test object", upper = 0),
    "The 'test object' must be equal or lower than 0."
  )
  # Lower bound violation
  expect_error(
    check_int(0, "test object",  lower = 1),
    "The 'test object' must be equal or higher than 1."
  )
  # Boundary not allowed (lower)
  expect_error(
    check_int(0, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be higher than 0."
  )
  # Boundary not allowed (upper)
  expect_error(
    check_int(1, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be lower than 1."
  )
  # NULL not allowed
  expect_error(
    check_int(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
  # NA not allowed
  expect_error(
    check_int(c(1, NA), "test object", allow_NA = FALSE, check_length = FALSE),
    "The 'test object' argument cannot contain NA/NaN values."
  )
  # Non-integer values rejected
  expect_error(check_int(1.5, "test"), "must be an integer")
  expect_error(check_int(c(1, 2.5, 3), "test", check_length = 3), "must be an integer")
})


test_that("check_list validates list inputs", {

  # Valid inputs
  expect_null(check_list(NULL,  "", allow_NULL = TRUE))
  expect_null(check_list(list(), "", allow_NULL = TRUE))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), ""))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_names = c("a", "b")))
  expect_null(check_list(list("a" = c("a", "b")), "", check_names = c("a", "b")))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_names = c("a", "b"), all_objects = TRUE))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_length = 2))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1, "c" = c("a", "b")), "", check_names = c("a", "b"), all_objects = TRUE, allow_other = TRUE))

  # Invalid type: string
  expect_error(
    check_list("string", "test object"),
    "The 'test object' argument must be a list."
  )
  # Invalid type: numeric
  expect_error(
    check_list(1, "test object"),
    "The 'test object' argument must be a list."
  )
  # Invalid type: logical
  expect_error(
    check_list(TRUE, "test object"),
    "The 'test object' argument must be a list."
  )
  # Empty list with length requirement
  expect_error(
    check_list(list(), "test object", check_length = 2),
    "The 'test object' argument cannot be NULL."
  )
  # Unrecognized names
  expect_error(
    check_list(list("c" = c("a", "b")), "test object", check_names = c("a", "b")),
    "The 'c' objects are not recognized by the 'test object' argument."
  )
  # Missing required names
  expect_error(
    check_list(list("a" = c("a", "b")), "test object", check_names = c("a", "b"), all_objects = TRUE),
    "The 'b' objects are missing in the 'test object' argument."
  )
  # NULL not allowed
  expect_error(
    check_list(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})


test_that("check functions handle empty vectors correctly", {

  # Empty vectors treated like NULL
  expect_error(check_bool(logical(0), "test"))
  expect_error(check_char(character(0), "test"))
  expect_error(check_real(numeric(0), "test"))
  expect_error(check_int(integer(0), "test"))

  # Empty vectors allowed with allow_NULL
  expect_null(check_bool(logical(0), "test", allow_NULL = TRUE))
  expect_null(check_char(character(0), "test", allow_NULL = TRUE))
  expect_null(check_real(numeric(0), "test", allow_NULL = TRUE))
  expect_null(check_int(integer(0), "test", allow_NULL = TRUE))
})


test_that("check functions support custom error prefix", {

  expect_error(
    check_bool("string", "test", call = "[custom] "),
    "\\[custom\\] The 'test' argument must be a logical"
  )
  expect_error(
    check_char(1, "test", call = "[custom] "),
    "\\[custom\\] The 'test' argument must be a character"
  )
  expect_error(
    check_real("a", "test", call = "[custom] "),
    "\\[custom\\] The 'test' argument must be a numeric"
  )
  expect_error(
    check_int("a", "test", call = "[custom] "),
    "\\[custom\\] The 'test' argument must be a numeric"
  )
  expect_error(
    check_list("a", "test", call = "[custom] "),
    "\\[custom\\] The 'test' argument must be a list"
  )
})
