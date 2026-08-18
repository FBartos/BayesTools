skip_if_not_test_profile("unit")

test_that("hypothesis AST preserves structure and quoted symbols", {

  hypothesis <- c(
    "(theta > 0) & (!(phi >= 1) | abs(eta) < 2)",
    "`a vs b` = 0 vs theta != 0",
    "alpha.beta + alpha_beta + `a b` + `TRUE` + `a>=b` > `mu alloc[level A]`"
  )
  ast <- hypothesis_parse(hypothesis)

  expect_s3_class(ast, "BayesTools_hypothesis_ast")
  expect_identical(ast$schema_version, 1L)
  expect_identical(hypothesis_render(hypothesis_parse(hypothesis_render(ast))),
                   hypothesis_render(ast))
  expect_identical(
    hypothesis_symbols(ast),
    c("theta", "phi", "eta", "a vs b", "alpha.beta", "alpha_beta",
      "a b", "TRUE", "a>=b", "mu alloc")
  )

  occurrences <- hypothesis_symbols(ast, occurrences = TRUE)
  level <- occurrences[occurrences$parameter == "mu alloc", , drop = FALSE]
  expect_true(nrow(level) > 0L)
  expect_true(all(level$level == "level A"))
  expect_true(all(level$symbol == "mu alloc[level A]"))
  expect_true(all(c("ast", "statement", "side", "node", "resolution") %in%
                    hypothesis_ast_schema()$object))
})

test_that("hypothesis parsing recognizes exact non-syntactic catalog aliases", {

  coordinates <- .bt_build_parameter_coordinates(columns = "theta")
  catalog     <- .bt_build_parameter_catalog(coordinates)
  quantity <- .bt_parameter_catalog_quantity(
    canonical_name = "(mu) random_total: var_prop(study)",
    namespace      = "mu",
    role           = "random_var_prop",
    component      = "random",
    owner_type     = "variance_allocation",
    owner_name     = "random_total",
    quantity       = "var_prop",
    arguments      = "study",
    source_type    = "composite",
    extraction_key = list(type = "test", dependencies = character())
  )
  quantity$provider    <- "RoBMA"
  quantity$quantity_id <- "RoBMA::random_fraction"
  aliases <- data.frame(
    alias       = "random_total: var_prop(study)",
    quantity_id = quantity$quantity_id,
    namespace   = quantity$namespace,
    component   = quantity$component,
    stringsAsFactors = FALSE
  )
  catalog <- parameter_catalog_extend(
    catalog,
    quantities = quantity,
    aliases    = aliases,
    provider   = "RoBMA"
  )
  hypothesis <- c(
    "random_total: var_prop(study) != 0 vs random_total: var_prop(study) = 0",
    "random_total: var_prop(study) != 1 vs random_total: var_prop(study) = 1"
  )

  ast <- hypothesis_parse(
    hypothesis,
    catalog   = catalog,
    component = "random"
  )

  expect_identical(
    hypothesis_symbols(ast),
    "random_total: var_prop(study)"
  )
  expect_identical(
    hypothesis_render(ast),
    gsub(
      "random_total: var_prop(study)",
      "`random_total: var_prop(study)`",
      hypothesis,
      fixed = TRUE
    )
  )
  expect_identical(
    hypothesis_render(hypothesis_parse(
      hypothesis_render(ast),
      catalog   = catalog,
      component = "random"
    )),
    hypothesis_render(ast)
  )
  expect_error(
    hypothesis_parse(hypothesis, component = "random"),
    "require 'catalog'",
    fixed = TRUE
  )
})

test_that("hypothesis rewriting edits exact symbol roots only", {

  ast <- hypothesis_parse(
    "a + aa + `a b` + abs(theta) > `a[level a]`"
  )
  rewritten <- hypothesis_rewrite(ast, c(a = "new root"))
  rendered <- hypothesis_render(rewritten)

  expect_match(rendered, "`new root` + aa + `a b` + abs(theta)",
               fixed = TRUE)
  expect_match(rendered, "`new root[level a]`", fixed = TRUE)
  expect_identical(
    hypothesis_symbols(rewritten),
    c("new root", "aa", "a b", "theta")
  )
  expect_error(
    hypothesis_rewrite(ast, c(a = "aa")),
    "duplicate or colliding"
  )
  expect_error(
    hypothesis_rewrite(ast, c(missing = "new")),
    "unknown hypothesis symbol"
  )
  expect_error(
    hypothesis_rewrite(ast, c(a = "new[level]")),
    "preserve their exact identity"
  )

  reserved <- hypothesis_rewrite(
    hypothesis_parse("theta > 0"),
    c(theta = "TRUE")
  )
  expect_identical(hypothesis_render(reserved), "`TRUE` > 0")
})

test_that("hypothesis resolution delegates ambiguity to the catalog", {

  coordinates <- .bt_build_parameter_coordinates(columns = "theta")
  catalog <- .bt_build_parameter_catalog(coordinates)
  location <- .bt_parameter_catalog_quantity(
    canonical_name = "effect_location",
    namespace = "location",
    role = "coefficient",
    extraction_key = list(type = "test", dependencies = character())
  )
  scale <- .bt_parameter_catalog_quantity(
    canonical_name = "effect_scale",
    namespace = "scale",
    role = "coefficient",
    extraction_key = list(type = "test", dependencies = character())
  )
  quantities <- rbind(location, scale)
  quantities$provider <- "RoBMA"
  quantities$quantity_id <- c("RoBMA::location", "RoBMA::scale")
  aliases <- data.frame(
    alias = c("effect", "effect"),
    quantity_id = quantities$quantity_id,
    namespace = quantities$namespace,
    component = c("", ""),
    stringsAsFactors = FALSE
  )
  catalog <- parameter_catalog_extend(
    catalog,
    quantities = quantities,
    aliases = aliases,
    provider = "RoBMA"
  )
  ast <- hypothesis_parse("effect > 0")

  condition <- tryCatch(
    hypothesis_resolve(ast, catalog),
    error = identity
  )
  expect_s3_class(condition, "BayesTools_parameter_ambiguous")
  expect_identical(condition$alias, "effect")
  expect_setequal(
    condition$candidates$quantity_id,
    c("RoBMA::location", "RoBMA::scale")
  )

  resolved <- hypothesis_resolve(ast, catalog, namespace = "scale")
  effect <- resolved$occurrences[
    resolved$occurrences$parameter == "effect",
    ,
    drop = FALSE
  ]
  direct <- parameter_catalog_resolve(catalog, "effect", namespace = "scale")
  expect_identical(resolved$schema_version, 1L)
  expect_true(all(effect$quantity_id == direct$quantity_id))
  expect_true(all(effect$canonical_name == "effect_scale"))
  expect_error(
    hypothesis_resolve(hypothesis_parse("missing > 0"), catalog),
    class = "BayesTools_parameter_not_found"
  )
})

test_that("public reference helpers agree with AST nodes", {

  hypothesis <- c(
    "theta[level A] = 0",
    "theta + phi != -1",
    "theta = 0 vs phi != 1"
  )
  references <- hypothesis_parse_point_reference(hypothesis)
  ast <- hypothesis_parse(hypothesis)

  expect_identical(references$direct, c(TRUE, FALSE, TRUE, TRUE))
  expect_identical(references$parameter,
                   c("theta", NA_character_, "theta", "phi"))
  expect_identical(references$level,
                   c("level A", NA_character_, NA_character_, NA_character_))
  expect_identical(
    references$operator,
    vapply(
      list(
        ast$statements[[1L]]$left,
        ast$statements[[2L]]$left,
        ast$statements[[3L]]$left,
        ast$statements[[3L]]$right
      ),
      function(side) if(side$type == "point") "=" else "!=",
      character(1)
    )
  )
})

test_that("hypothesis BF consumes validated ASTs without changing results", {

  prior <- seq(-2, 2, length.out = 201)
  posterior <- seq(-1.8, 2.2, length.out = 201)
  text <- "theta = 0"
  ast <- hypothesis_parse(text)
  arguments <- list(
    posterior = posterior,
    prior = prior,
    parameter = "theta",
    density_method = "normal",
    columns = "all"
  )

  from_text <- do.call(hypothesis_BF, c(arguments, list(hypothesis = text)))
  from_ast <- do.call(hypothesis_BF, c(arguments, list(hypothesis = ast)))
  expect_identical(from_ast, from_text)

  expect_identical(unserialize(serialize(ast, NULL)), ast)
  stale <- ast
  stale$schema_version <- 2L
  expect_error(
    hypothesis_render(stale),
    "missing or unsupported"
  )
  malformed <- unclass(ast)
  expect_error(
    hypothesis_BF(posterior, prior, malformed, parameter = "theta"),
    "'hypothesis'"
  )
})
