# ============================================================================ #
# TEST FILE: Distribution Tools Helpers
# ============================================================================ #
#
# PURPOSE:
#   Tests for internal distribution helper functions like .check_log,
#   .check_log.p, and .check_lower.tail.
#
# DEPENDENCIES:
#   - None (pure R)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# TAGS: @evaluation, @distributions, @tools
# ============================================================================ #


test_that(".check_log works", {

  expect_null(BayesTools:::.check_log(TRUE))
  expect_null(BayesTools:::.check_log(FALSE))

  expect_error(BayesTools:::.check_log("TRUE"), "must be a logical")
  expect_error(BayesTools:::.check_log(1), "must be a logical")
  expect_error(BayesTools:::.check_log(NULL), "cannot be NULL")

})


test_that(".check_log.p works", {

  expect_null(BayesTools:::.check_log.p(TRUE))
  expect_null(BayesTools:::.check_log.p(FALSE))

  expect_error(BayesTools:::.check_log.p("TRUE"), "must be a logical")
  expect_error(BayesTools:::.check_log.p(1), "must be a logical")

})


test_that(".check_lower.tail works", {

  expect_null(BayesTools:::.check_lower.tail(TRUE))
  expect_null(BayesTools:::.check_lower.tail(FALSE))

  expect_error(BayesTools:::.check_lower.tail("TRUE"), "must be a logical")
  expect_error(BayesTools:::.check_lower.tail(1), "must be a logical")

})


test_that(".check_x works", {

  expect_null(BayesTools:::.check_x(0.5))
  expect_null(BayesTools:::.check_x(c(0, 0.5, 1)))
  expect_null(BayesTools:::.check_x(0.5, lower = 0, upper = 1))

  expect_error(BayesTools:::.check_x(-1, lower = 0), "must be equal or higher than 0")
  expect_error(BayesTools:::.check_x(2, upper = 1), "must be equal or lower than 1")
  expect_error(BayesTools:::.check_x("a"), "must be a numeric")

})


test_that(".check_n works", {

  expect_null(BayesTools:::.check_n(1))
  expect_null(BayesTools:::.check_n(100))

  expect_error(BayesTools:::.check_n(0), "must be equal or higher than 1")
  expect_error(BayesTools:::.check_n(-1), "must be equal or higher than 1")
  expect_error(BayesTools:::.check_n(c(1, 2)), "must have length '1'")
  expect_error(BayesTools:::.check_n("a"), "must be a numeric")

})


test_that(".check_q works", {

  expect_null(BayesTools:::.check_q(0.5))
  expect_null(BayesTools:::.check_q(c(-1, 0, 1)))
  expect_null(BayesTools:::.check_q(0.5, lower = 0, upper = 1))

  expect_error(BayesTools:::.check_q(-1, lower = 0), "must be equal or higher than 0")
  expect_error(BayesTools:::.check_q(2, upper = 1), "must be equal or lower than 1")
  expect_error(BayesTools:::.check_q("a"), "must be a numeric")

})


test_that(".check_p works", {

  # Standard probability checks (log.p = FALSE)
  expect_null(BayesTools:::.check_p(0.5, FALSE))
  expect_null(BayesTools:::.check_p(c(0, 0.5, 1), FALSE))

  expect_error(BayesTools:::.check_p(-0.1, FALSE), "must be equal or higher than 0")
  expect_error(BayesTools:::.check_p(1.1, FALSE), "must be equal or lower than 1")
  expect_error(BayesTools:::.check_p("a", FALSE), "must be a numeric")

  # Log probability checks (log.p = TRUE)
  expect_null(BayesTools:::.check_p(-1, TRUE))
  expect_null(BayesTools:::.check_p(c(-2, -1, 0), TRUE))

  expect_error(BayesTools:::.check_p(0.1, TRUE), "must be equal or lower than 0")
  expect_error(BayesTools:::.check_p("a", TRUE), "must be a numeric")

})
