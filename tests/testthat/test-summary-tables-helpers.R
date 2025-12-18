# ============================================================================ #
# TEST FILE: Summary Tables Helper Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for format_BF, format_estimates, and other summary table
#   formatting utilities.
#
# DEPENDENCIES:
#   - common-functions.R: test_reference_table, REFERENCE_DIR
#
# SKIP CONDITIONS:
#   - None (can run on CRAN - pure R with reference file testing)
#
# TAGS: @evaluation, @summary-tables, @formatting
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "summary-tables-helpers")
source(testthat::test_path("common-functions.R"))


test_that("format_BF works correctly", {

  # Basic usage
  BF <- format_BF(3.5)
  expect_equal(as.numeric(BF), 3.5)
  expect_equal(attr(BF, "name"), "BF")
  expect_false(attr(BF, "logBF"))
  expect_false(attr(BF, "BF01"))

  # With BF01 = TRUE (inverted)
  BF_01 <- format_BF(2, BF01 = TRUE)
  expect_equal(as.numeric(BF_01), 0.5)
  expect_equal(attr(BF_01, "name"), "1/BF")
  expect_true(attr(BF_01, "BF01"))

  # With logBF = TRUE
  BF_log <- format_BF(exp(2), logBF = TRUE)
  expect_equal(as.numeric(BF_log), 2, tolerance = 1e-10)
  expect_match(attr(BF_log, "name"), "log\\(BF\\)")
  expect_true(attr(BF_log, "logBF"))

  # With inclusion = TRUE
  BF_incl <- format_BF(5, inclusion = TRUE)
  expect_equal(attr(BF_incl, "name"), "Inclusion BF")

  # With BF01 = TRUE and inclusion = TRUE
  BF_excl <- format_BF(5, BF01 = TRUE, inclusion = TRUE)
  expect_equal(attr(BF_excl, "name"), "Exclusion BF")

  # Combined logBF and BF01
  BF_both <- format_BF(10, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(BF_both), log(0.1), tolerance = 1e-10)
  expect_match(attr(BF_both, "name"), "log\\(1/BF\\)")

  # Vector input with NA - NA_real_ must be part of numeric vector
  BF_vec_na <- format_BF(c(1, 2, NA_real_))
  expect_equal(as.numeric(BF_vec_na)[1:2], c(1, 2))
  expect_true(is.na(as.numeric(BF_vec_na)[3]))

  # Vector input
  BF_vec <- format_BF(c(1, 2, 3))
  expect_equal(as.numeric(BF_vec), c(1, 2, 3))

})


test_that("format_BF input validation works", {

  expect_error(format_BF(-1), "must be equal or higher than 0")
  expect_error(format_BF("3"), "must be a numeric")
  expect_error(format_BF(3, logBF = "TRUE"), "must be a logical")
  expect_error(format_BF(3, BF01 = "TRUE"), "must be a logical")

})


test_that("add_column works correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Create a proper BayesTools_table with 3 columns (needed for middle position test)
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Add column at end (default)
  result1 <- add_column(test_data, "CI_lower", c(-0.5, 0.8))
  expect_equal(ncol(result1), 4)
  expect_true("CI_lower" %in% names(result1))
  expect_s3_class(result1, "BayesTools_table")
  expect_equal(attr(result1, "type"), c("estimate", "estimate", "estimate", "estimate"))
  test_reference_table(result1, "add_column_end.txt")

  # Add column at specific position
  result2 <- add_column(test_data, "CI_lower", c(-0.5, 0.8), column_position = 2)
  expect_equal(ncol(result2), 4)
  expect_equal(names(result2)[2], "CI_lower")
  test_reference_table(result2, "add_column_position2.txt")

  # Add column at position 1
  result3 <- add_column(test_data, "ID", c(1, 2), column_position = 1)
  expect_equal(ncol(result3), 4)
  expect_equal(names(result3)[1], "ID")
  expect_equal(attr(result3, "type")[1], "integer")
  test_reference_table(result3, "add_column_position1.txt")

  # With specified column type
  result4 <- add_column(test_data, "Prob", c(0.5, 0.8), column_type = "probability")
  expect_equal(ncol(result4), 4)
  expect_equal(attr(result4, "type")[4], "probability")
  test_reference_table(result4, "add_column_probability.txt")

  # Add column with string values (must specify column_type)
  result5 <- add_column(test_data, "Category", c("A", "B"), column_type = "string")
  expect_equal(ncol(result5), 4)
  expect_true("Category" %in% names(result5))
  expect_equal(attr(result5, "type")[4], "string")
  test_reference_table(result5, "add_column_string.txt")

})


test_that("add_column input validation works", {

  # Create a proper BayesTools_table with data (3 columns)
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Wrong table class
  expect_error(add_column(data.frame(a = 1), "b", 2), "must be of class 'BayesTools_table'")

  # Wrong column_values length
  expect_error(add_column(test_data, "new", c(1, 2, 3)), "must be a vector of the same length")

})


test_that("remove_column works correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Create a proper BayesTools_table with data
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Remove last column (default)
  result0 <- remove_column(test_data)
  expect_equal(ncol(result0), 2)
  expect_false("SD" %in% names(result0))
  expect_s3_class(result0, "BayesTools_table")
  expect_equal(attr(result0, "type"), c("estimate", "estimate"))
  test_reference_table(result0, "remove_column_last.txt")

  # Remove by position
  result2 <- remove_column(test_data, column_position = 2)
  expect_equal(ncol(result2), 2)
  expect_false("Median" %in% names(result2))
  test_reference_table(result2, "remove_column_position2.txt")

})


test_that("remove_column input validation works", {

  # Create a proper BayesTools_table with data
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Wrong table class
  expect_error(remove_column(data.frame(a = 1)), "must be of class 'BayesTools_table'")

  # Invalid column position
  expect_error(remove_column(test_data, column_position = 10), "'column_position'")

})


test_that("ensemble_estimates_empty_table works correctly", {

  empty_table <- ensemble_estimates_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_estimates_empty.txt")

})


test_that("ensemble_inference_empty_table works correctly", {

  empty_table <- ensemble_inference_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_inference_empty.txt")

})


test_that("ensemble_summary_empty_table works correctly", {

  empty_table <- ensemble_summary_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_summary_empty.txt")

})


test_that("ensemble_diagnostics_empty_table works correctly", {

  empty_table <- ensemble_diagnostics_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_diagnostics_empty.txt")

})
