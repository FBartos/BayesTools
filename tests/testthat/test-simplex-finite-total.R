skip_if_not_test_profile("unit")

test_that("R133 finite simplex entries cannot normalize a nonfinite total", {
  huge <- c(1e308, 1e308)
  expect_error(.canonicalize_simplex(huge, "weights"),
    "The 'weights' simplex values must sum to one; row 1 has a non-finite total or roundoff bound.", fixed = TRUE)
  rows <- rbind(valid = c(.5, .5), invalid = huge)
  expect_error(.canonicalize_simplex(rows, "weights"),
    "The 'weights' simplex values must sum to one; row 2 has a non-finite total or roundoff bound.", fixed = TRUE)
  expect_identical(lpdf(prior("dirichlet", list(alpha = c(1, 1))), huge), -Inf)
  expect_error(prior_ordered(prior("normal", list(0, 1)), allocation = huge),
    "The 'allocation fixed allocation' simplex values must sum to one; row 1 has a non-finite total or roundoff bound.", fixed = TRUE)
  expect_error(.canonicalize_simplex(c(1, 1)), "exceeding the roundoff bound")
  expect_error(.canonicalize_simplex(c(Inf, 0)), "must be finite numeric values")
})

test_that("R133 near-one simplex values retain normalization diagnostics and labels", {
  values <- c(first = .2, second = .3, third = .5 + .Machine$double.eps)
  normalized <- .canonicalize_simplex(values, diagnostics = TRUE)
  expect_identical(normalized$values, values / sum(values))
  expect_identical(normalized$diagnostics$original_sum, sum(values))
  expect_identical(normalized$diagnostics$roundoff_bound, .simplex_roundoff_bound(values))
  expect_identical(normalized$diagnostics$max_correction, max(abs(values / sum(values) - values)))
  rows <- rbind(a = values, b = c(first = .2, second = .3, third = .5))
  expect_identical(.canonicalize_simplex(rows), rows / rowSums(rows))
})
