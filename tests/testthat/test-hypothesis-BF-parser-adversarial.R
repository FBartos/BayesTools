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


test_that("hypothesis parser implements the documented region grammar", {

  expect_s3_class(
    BayesTools:::.parse_hypothesis_BF("(theta > 0)"),
    "BayesTools_hypothesis_BF_parsed"
  )
  expect_s3_class(
    BayesTools:::.parse_hypothesis_BF("!(theta > 0)"),
    "BayesTools_hypothesis_BF_parsed"
  )
  expect_s3_class(
    BayesTools:::.parse_hypothesis_BF(
      "(theta > 0) & (!(phi >= 1) | abs(eta) < 2)"
    ),
    "BayesTools_hypothesis_BF_parsed"
  )
  expect_s3_class(
    BayesTools:::.parse_hypothesis_BF(
      "theta > 0 vs (phi <= 1 | eta >= 2)"
    ),
    "BayesTools_hypothesis_BF_parsed"
  )

  draws <- data.frame(
    theta = c(-1, 1, 2),
    phi   = c(0, 2, 0),
    eta   = c(0, 3, 0)
  )
  expect_identical(
    BayesTools:::.hypothesis_eval_condition(
      "(theta > 0) & (!(phi >= 1) | abs(eta) < 2)",
      draws
    ),
    c(FALSE, FALSE, TRUE)
  )

  functions <- c("abs", "exp", "log", "sqrt", "plogis", "qlogis")
  function_draws <- data.frame(theta = c(0.25, 0.75))
  for(fun in functions){
    condition <- paste0(fun, "(theta) > -1e300")
    expect_true(all(BayesTools:::.hypothesis_eval_condition(
      condition,
      function_draws
    )))
  }
})


test_that("point hypotheses accept literals or symbolic right-hand sides", {

  parsed <- lapply(
    c("theta = -0.5", "theta == +2", "theta != 3"),
    BayesTools:::.parse_hypothesis_BF
  )
  expect_equal(
    vapply(parsed, function(x) x[["left"]][["value"]], numeric(1)),
    c(-0.5, 2, 3)
  )

  symbolic <- hypothesis_parse("theta = other")
  expect_identical(
    hypothesis_render(symbolic),
    "theta - other = 0"
  )
  expect_identical(symbolic$statements[[1L]]$left$label, "theta = other")

  for(hypothesis in c(
    "theta = 1 + 1",
    "theta = --1",
    "theta = Inf"
  )){
    expect_error(
      BayesTools:::.parse_hypothesis_BF(hypothesis),
      "numeric value written as one finite literal",
      fixed = TRUE
    )
  }
})


test_that("constant-left relations use the canonical scalar target", {

  reversed <- hypothesis_parse("0 > theta vs 0 = theta")
  expect_identical(
    hypothesis_render(reversed),
    "theta < 0 vs theta = 0"
  )
  expect_identical(
    reversed$statements[[1L]]$left$label,
    "0 > theta"
  )
  expect_identical(
    reversed$statements[[1L]]$right$label,
    "0 = theta"
  )

  point_reference <- hypothesis_parse_point_reference(
    "0 > theta vs 0 = theta"
  )
  expect_true(point_reference[["direct"]])
  expect_identical(point_reference[["symbol"]], "theta")

  set.seed(10)
  prior     <- stats::rnorm(2000)
  posterior <- stats::rnorm(2000, mean = 0.3)
  reversed_BF <- hypothesis_BF(
    posterior      = posterior,
    prior          = prior,
    hypothesis     = "0 > theta vs 0 = theta",
    parameter      = "theta",
    density_method = "normal"
  )
  canonical_BF <- hypothesis_BF(
    posterior      = posterior,
    prior          = prior,
    hypothesis     = "theta < 0 vs theta = 0",
    parameter      = "theta",
    density_method = "normal"
  )
  expect_equal(
    attr(reversed_BF, "raw_BF"),
    attr(canonical_BF, "raw_BF")
  )
})


test_that("escaped reserved identifiers take precedence over R constants", {

  reserved <- c("Inf", "NaN", "NA", "TRUE", "FALSE")
  for(name in reserved){
    draws <- data.frame(c(1, -1), check.names = FALSE)
    names(draws) <- name
    condition <- paste0("`", name, "` > 0")

    expect_identical(
      BayesTools:::.hypothesis_expression_symbols(
        BayesTools:::.hypothesis_parse_expression(condition)
      ),
      name
    )
    expect_identical(
      BayesTools:::.hypothesis_eval_condition(condition, draws),
      c(TRUE, FALSE)
    )
  }

  for(hypothesis in c(
    "Inf > theta",
    "NaN > theta",
    "NA > theta",
    "TRUE > theta",
    "FALSE > theta"
  )){
    expect_error(
      BayesTools:::.parse_hypothesis_BF(hypothesis),
      "finite numeric literals|Unescaped reserved literal"
    )
  }
})


test_that("hypothesis parser rejects constructs outside the grammar", {

  expect_error(
    BayesTools:::.parse_hypothesis_BF("!(theta == 0)"),
    "Point equalities cannot be negated",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("(theta == 0)"),
    "point equalities cannot be parenthesized",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("theta > 0 & abs(phi)"),
    "Unsupported region operator 'abs'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("theta > 0 && phi < 1"),
    "Unsupported region operator '&&'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("theta > 1 / 0"),
    "Constant hypothesis arithmetic must evaluate to one finite numeric value.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("abs(theta, phi) > 0"),
    "unsupported number of arguments",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.parse_hypothesis_BF("sin(theta) > 0"),
    "Unsupported hypothesis expression operator or function 'sin'.",
    fixed = TRUE
  )
})


test_that("parentheses do not change comparison compatibility", {

  point <- BayesTools:::.parse_hypothesis_BF("theta = 0")[["left"]]
  explicit_not_point <- BayesTools:::.parse_hypothesis_BF(
    "(theta) != 0"
  )[["left"]]
  expect_true(BayesTools:::.hypothesis_sides_point_complement(
    point,
    explicit_not_point
  ))

  simple_region <- BayesTools:::.parse_hypothesis_BF(
    "(theta) > 1"
  )[["left"]]
  compound_region <- BayesTools:::.parse_hypothesis_BF(
    "(theta > 0) & (theta < 2)"
  )[["left"]]
  expect_true(BayesTools:::.hypothesis_point_region_compatible(
    point,
    simple_region
  ))
  expect_true(BayesTools:::.hypothesis_point_region_compatible(
    point,
    compound_region
  ))

  posterior <- c(-2, -1, 1, 2)
  prior <- c(-3, -1, 1, 3)
  expect_s3_class(
    hypothesis_BF(
      posterior,
      prior,
      "theta = 0 vs (theta) != 0",
      parameter = "theta"
    ),
    "BayesTools_hypothesis_BF"
  )
})


test_that("whitelisted function names remain valid quantity identifiers", {

  functions <- c("abs", "exp", "log", "sqrt", "plogis", "qlogis")
  for(name in functions){
    expect_identical(
      BayesTools:::.hypothesis_expression_symbols(
        BayesTools:::.hypothesis_parse_expression(paste(name, "> 0"))
      ),
      name
    )
  }
  expect_identical(
    BayesTools:::.hypothesis_expression_symbols(
      BayesTools:::.hypothesis_parse_expression("exp(abs) > 1")
    ),
    "abs"
  )

  posterior <- data.frame(abs = c(-2, -1, 1, 2))
  prior <- data.frame(abs = c(-3, -1, 1, 3))
  expect_s3_class(
    hypothesis_BF(posterior, prior, "exp(abs) > 1"),
    "BayesTools_hypothesis_BF"
  )
})
