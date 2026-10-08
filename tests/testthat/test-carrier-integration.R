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

test_that("removing a public linear recipe clears its companion space", {
  x <- c(1, 2, 3)
  posterior_metadata(x, "linear_weights") <- c(a = 1)
  expect_identical(posterior_metadata(x, "linear_weight_space"), "coefficient")
  posterior_metadata(x, "linear_weights") <- NULL
  expect_null(posterior_metadata(x, "linear_weights"))
  expect_null(posterior_metadata(x, "linear_weight_space"))
})

test_that("public linear recipes reject nonfinite weights without mutation", {
  x <- c(1, 2, 3)
  for(value in c(NA_real_, NaN, Inf, -Inf)){
    for(weights in list(c(a = value), matrix(value, 1L, 1L, dimnames = list("target", "a")))){
      original <- serialize(x, NULL)
      condition <- tryCatch({posterior_metadata(x, "linear_weights") <- weights; NULL}, error = identity)
      expect_s3_class(condition, "error")
      if(inherits(condition, "condition")){
        expect_match(conditionMessage(condition), "finite numeric vector or matrix", fixed = TRUE)
        expect_null(conditionCall(condition))
      }
      expect_identical(serialize(x, NULL), original)
    }
  }
  finite <- matrix(c(1, 0, -2, 3), 2L, dimnames = list(c("A", "B"), c("a", "b")))
  posterior_metadata(x, "linear_weights") <- finite
  expect_identical(posterior_metadata(x, "linear_weights"), finite)
})

test_that("list measure refusals emit exactly once at the actual level", {
  good <- structure(seq(-1, 1, length.out = 100L), class = c("marginal_posterior.simple", "marginal_posterior"))
  attr(good, "parameter") <- "theta"
  posterior_metadata(good, "atoms") <- posterior_atom_attribute()
  posterior_metadata(good, "prior_density") <- prior("normal", list(0, 1))
  bad <- good
  posterior_metadata(bad, "measure_unavailable") <- data.frame(column = "theta", measure = "atoms",
    reason = "Declared contribution measure is unavailable", cause = "unsupported_contribution_measure")
  parent <- structure(list(bad = bad, good = good), class = c("marginal_posterior.factor", "marginal_posterior", "list"), parameter = "theta")
  messages <- character()
  actual <- withCallingHandlers(Savage_Dickey_BF(parent, silent = FALSE), warning = function(w){messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")})
  expect_length(messages, 1L)
  expect_match(messages, "theta[bad]:", fixed = TRUE)
  expect_true(is.na(actual$bad))
  expect_true(is.finite(actual$good))
  expect_s3_class(attr(actual$bad, "numerical_diagnostics"), "BayesTools_formula_measure_unavailable")
  expect_silent(Savage_Dickey_BF(parent, silent = TRUE))
  expect_error(Savage_Dickey_BF(bad, silent = TRUE), class = "BayesTools_formula_measure_unavailable")
})
