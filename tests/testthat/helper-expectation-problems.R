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
