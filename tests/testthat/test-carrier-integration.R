skip_if_not_test_profile("unit")

test_that("public formatting retains a live carrier's canonical source logs", {
  producer <- hypothesis_BF(c(1, 2, 3), prior("normal", list(0, 1)),
    "a = 356 vs a != 356", parameter = "a", density_method = "normal", columns = "all")
  expect_equal(as.numeric(format_BF(producer$BF, logBF = TRUE)), 710, tolerance = 1e-12)
  inverse <- as.numeric(format_BF(producer$BF, BF01 = TRUE))
  expect_gt(inverse, 0)
  expect_equal(inverse / exp(-710), 1, tolerance = 1e-12)
  logs <- c(overflow = 800, tiny = -800, ordinary = -2, missing = NA_real_)
  diagnostics <- list(list(marker = 1), list(marker = 2), NULL, NULL)
  source_bounds <- c(">", "<", NA_character_, NA_character_)
  for(start_log in c(FALSE, TRUE)) for(start_inverse in c(FALSE, TRUE)){
    carrier <- .format_BF_from_log(logs, start_log, start_inverse,
      bound_operator = source_bounds, diagnostics = diagnostics)
    for(want_log in c(FALSE, TRUE)) for(want_inverse in c(FALSE, TRUE)){
      actual <- format_BF(carrier, want_log, want_inverse)
      expected <- .format_BF_from_log(logs, want_log, want_inverse,
        bound_operator = source_bounds, diagnostics = diagnostics)
      expect_identical(actual, expected)
      expect_identical(names(actual), names(logs))
      expect_identical(.BF_carrier_log(actual), unname(logs))
      expect_identical(attr(actual, "numerical_diagnostics"), diagnostics)
    }
    selected <- carrier[c(3L, 1L, 1L, 4L)]
    actual <- format_BF(selected, logBF = TRUE)
    expect_identical(as.numeric(actual), unname(logs[c(3L, 1L, 1L, 4L)]))
    expect_identical(attr(actual, "numerical_diagnostics"), diagnostics[c(3L, 1L, 1L, 4L)])
  }
  for(value in list(c(0, Inf, NA, 2), unclass(producer$BF))){
    plain <- as.numeric(value)
    expect_identical(as.numeric(format_BF(plain, logBF = TRUE)), log(plain))
  }
  for(invalidated in list(producer$BF + 0, exp(producer$BF), {x <- producer$BF; x[1L] <- 2; x})){
    expect_null(.BF_carrier_log(invalidated))
    expect_identical(as.numeric(format_BF(invalidated, logBF = TRUE)), log(as.numeric(invalidated)))
  }
})
