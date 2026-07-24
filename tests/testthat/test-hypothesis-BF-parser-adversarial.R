skip_if_not_test_profile("unit")


test_that("level-reference normalization preserves explicit escaping", {

  references <- c(
    "`mu alloc[level A]`",
    "`mu alloc[a & b]`",
    "`mu alloc[a vs b]`"
  )

  expect_equal(hypothesis_normalize_level_references(references), references)
  expect_equal(
    hypothesis_normalize_level_references(
      hypothesis_normalize_level_references("theta[level A]")
    ),
    "`theta[level A]`"
  )

  parsed <- hypothesis_parse_level_reference(references)
  expect_true(all(parsed[["direct"]]))
  expect_equal(parsed[["symbol"]], c(
    "mu alloc[level A]",
    "mu alloc[a & b]",
    "mu alloc[a vs b]"
  ))
  expect_equal(parsed[["parameter"]], rep("mu alloc", 3L))
  expect_equal(parsed[["level"]], c("level A", "a & b", "a vs b"))
})


test_that("hypothesis parser rejects malformed call nodes cleanly", {

  expect_error(
    hypothesis_parse_point_reference("abs() > 0"),
    "Hypothesis expression call 'abs' requires at least one argument.",
    fixed = TRUE
  )
  expect_error(
    hypothesis_parse_point_reference("base::log(theta) > 0"),
    "Unsupported hypothesis expression call 'base::log'.",
    fixed = TRUE
  )
})
