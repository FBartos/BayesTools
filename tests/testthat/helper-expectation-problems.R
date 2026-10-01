# The failed expectations of 'code', as a character vector of their messages.
#
# Every expectation the code makes is still checked with its own tolerance and
# message, but none is recorded by the reporter: a test that checks many
# elements (one per row, level, or grid point) asserts the returned vector, which
# is empty when everything held, once per unit it checks. That keeps a failure as
# specific as the per-element expectation was (the message names the element
# and the values), without one recorded result per element. An error, a skip or a
# warning of the code is not an expectation and passes through.
expectation_problems <- function(code) {
  problems <- character()
  withCallingHandlers(
    code,
    expectation_success = function(condition) {
      invokeRestart("continue_test")
    },
    expectation_failure = function(condition) {
      problems <<- c(problems, conditionMessage(condition))
      invokeRestart("continue_test")
    }
  )

  problems
}

# 'expect_equal()' for every element of two numeric vectors in one expectation,
# with the criterion of 'expect_equal()' on a single number: the difference of
# an element is relative to its expected value where that exceeds the
# tolerance, and absolute where it does not. 'expect_equal()' of the vectors
# compares the mean difference of the elements that differ with their mean
# expected value, so an element with a large expected value hides a large
# difference of another; here every element has to meet the tolerance.
# Elements that are both NA, or equal (including equal infinities), agree. A
# failure lists the elements that differ with their values.
expect_equal_each <- function(object, expected,
                              tolerance = testthat::testthat_tolerance(),
                              info = NULL) {
  object   <- as.numeric(object)
  expected <- as.numeric(expected)
  if (length(object) != length(expected)) {
    return(testthat::expect(
      FALSE,
      sprintf("The numbers of elements differ: %d and %d.", length(object), length(expected)),
      info = info
    ))
  }
  both_na   <- is.na(object) & is.na(expected)
  one_na    <- xor(is.na(object), is.na(expected))
  difference <- abs(object - expected)
  scale     <- ifelse(is.finite(expected) & abs(expected) > tolerance, abs(expected), 1)
  agree     <- both_na | (!one_na & (object == expected | difference / scale < tolerance))
  agree[is.na(agree)] <- FALSE
  failed    <- which(!agree)
  testthat::expect(
    length(failed) == 0L,
    sprintf(
      "%d of %d elements differ by more than %s: %s",
      length(failed), length(object), format(tolerance),
      paste(
        sprintf(
          "[%d] %s vs %s", utils::head(failed, 5L),
          format(object[utils::head(failed, 5L)], digits = 15L),
          format(expected[utils::head(failed, 5L)], digits = 15L)
        ),
        collapse = "; "
      )
    ),
    info = info
  )
  invisible(object)
}

# Passing comparisons skip testthat's per-call output setup. Every
# 'expect_identical()', 'expect_equal()', 'expect_true()' and 'expect_false()'
# reaches 'waldo_compare()', which sets up reproducible output (language,
# options, environment variables, and a gettext cache reset that creates and
# removes a temporary directory) before it compares anything. That setup costs
# about 2 ms per call on Windows and was the larger part of the cost of the
# passing expectations of the unit profile. The versions below record a success
# only where testthat's comparison cannot fail: the values are identical() (which
# implies equality at any tolerance) and the call passes no comparison option
# (one can be stricter than identical(), e.g. 'ignore_encoding = FALSE', or
# deprecated, which testthat warns about) and a tolerance testthat accepts; or
# the object of 'expect_true()' / 'expect_false()' is a logical vector without a
# class whose value is the constant (testthat ignores its other attributes).
# Otherwise they make testthat's own expectation of the same arguments, labelled
# as testthat labels the original call, so a failure reads and is located as
# before. They mask the testthat functions of the same names for the test files
# only.
.expectation_label <- utils::getFromNamespace("expr_label", "testthat")

# Whether testthat's comparison of identical values passes for certain: no
# comparison options, and no 'waldo_opts' attribute that sets one.
.expectation_identical <- function(object, expected, n_options) {
  n_options == 0L && identical(object, expected) &&
    is.null(attr(object, "waldo_opts", exact = TRUE))
}

# Whether testthat accepts the tolerance (a number of at least zero, or NULL).
.expectation_tolerance_valid <- function(tolerance) {
  is.null(tolerance) ||
    (is.numeric(tolerance) && !is.object(tolerance) && length(tolerance) == 1L &&
       !is.na(tolerance) && tolerance >= 0)
}

# Whether 'object' is a logical vector without a class whose value is 'constant'.
.expectation_constant <- function(object, constant) {
  is.logical(object) && !is.object(object) && identical(as.vector(object), constant)
}

expect_identical <- function(object, expected, info = NULL, label = NULL,
                             expected.label = NULL, ...) {
  if (.expectation_identical(object, expected, ...length())) {
    testthat::succeed()
    return(invisible(object))
  }
  testthat::expect_identical(
    object, expected, info = info,
    label = if (is.null(label)) .expectation_label(substitute(object)) else label,
    expected.label = if (is.null(expected.label)) .expectation_label(substitute(expected)) else expected.label,
    ...
  )
}

expect_equal <- function(object, expected, ..., tolerance, info = NULL,
                         label = NULL, expected.label = NULL) {
  if (.expectation_identical(object, expected, ...length()) &&
      (missing(tolerance) || .expectation_tolerance_valid(tolerance))) {
    testthat::succeed()
    return(invisible(object))
  }
  label <- if (is.null(label)) .expectation_label(substitute(object)) else label
  expected.label <- if (is.null(expected.label)) .expectation_label(substitute(expected)) else expected.label
  if (missing(tolerance)) {
    testthat::expect_equal(
      object, expected, ..., info = info, label = label, expected.label = expected.label
    )
  } else {
    testthat::expect_equal(
      object, expected, ..., tolerance = tolerance, info = info,
      label = label, expected.label = expected.label
    )
  }
}

expect_true <- function(object, info = NULL, label = NULL) {
  if (.expectation_constant(object, TRUE)) {
    testthat::succeed()
    return(invisible(object))
  }
  testthat::expect_true(
    object, info = info,
    label = if (is.null(label)) .expectation_label(substitute(object)) else label
  )
}

expect_false <- function(object, info = NULL, label = NULL) {
  if (.expectation_constant(object, FALSE)) {
    testthat::succeed()
    return(invisible(object))
  }
  testthat::expect_false(
    object, info = info,
    label = if (is.null(label)) .expectation_label(substitute(object)) else label
  )
}
