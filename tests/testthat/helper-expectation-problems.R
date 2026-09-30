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
