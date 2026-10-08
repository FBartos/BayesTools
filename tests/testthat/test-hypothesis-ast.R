skip_if_not_test_profile("unit")

test_that("escaped literal roots preserve statement grammar on roundtrip", {
  roots <- c("a`b", "a` vs b", "a`=b", "a`>=b", "a`&b", "a`|b",
    "a`<-b", "a\\b", "a`\\ vs b", "a\\`b", "a\\\\`b", "a vs b", "TRUE")
  for(root in roots){
    for(text in c("theta = 0", "theta > 0", "theta > 0 vs theta <= 0")){
      rewritten <- hypothesis_rewrite(hypothesis_parse(text), c(theta = root))
      reparsed <- tryCatch(hypothesis_parse(hypothesis_render(rewritten)), error = identity)
      expect_false(inherits(reparsed, "error"), info = root)
      if(inherits(reparsed, "error")) next
      expect_identical(hypothesis_symbols(reparsed), root)
      expect_identical(hypothesis_render(reparsed), hypothesis_render(rewritten))
      expect_identical(lapply(reparsed$statements, function(x) x$left$expression),
        lapply(rewritten$statements, function(x) x$left$expression))
    }
  }
  expect_identical(hypothesis_symbols(hypothesis_parse("theta[level A] > 0")), "theta")
  for(text in c("theta <- 0", "theta > 0 vs phi <- 1", "sin(theta) > 0",
                "theta + TRUE > 0", "theta + Inf > 0")) expect_error(hypothesis_parse(text))
})

test_that("integer hypothesis literals share double semantic rendering", {
  integer_text <- c("theta = 1L", "theta > 1L", "2L * theta = 0",
    "theta > -1L", "theta / 2L > 0", "theta^2L > 1L")
  double_text <- gsub("L", "", integer_text, fixed = TRUE)
  for(i in seq_along(integer_text)){
    integer_ast <- tryCatch(hypothesis_parse(integer_text[[i]]), error = identity)
    expect_false(inherits(integer_ast, "error"), info = integer_text[[i]])
    if(inherits(integer_ast, "error")) next
    double_ast <- hypothesis_parse(double_text[[i]])
    expect_identical(integer_ast$statements[[1L]]$left$label, integer_text[[i]])
    expect_identical(integer_ast$statements[[1L]]$left$expression,
      double_ast$statements[[1L]]$left$expression)
    expect_identical(hypothesis_render(integer_ast), hypothesis_render(double_ast))
    reparsed <- hypothesis_parse(hypothesis_render(integer_ast))
    expect_identical(reparsed$statements[[1L]]$left$expression,
      double_ast$statements[[1L]]$left$expression)
    expect_identical(reparsed$statements[[1L]]$left$value,
      double_ast$statements[[1L]]$left$value)
  }
})

test_that("hypothesis rendering preserves former codec and placeholder names", {
  hypotheses <- c(".BayesToolsHypothesisLiteral1. + 1 = 0",
    "prefix.BayesToolsHypothesisLiteral1.suffix + 1 = 0",
    ".BayesTools_escaped_constant_Inf = 0",
    "`Inf` + .BayesTools_escaped_constant_Inf + 0.30000000000000004 > `TRUE`",
    "`NA` + `NaN` + `FALSE` > 0")
  for(text in hypotheses){
    ast <- hypothesis_parse(text)
    reparsed <- hypothesis_parse(hypothesis_render(ast))
    expect_identical(hypothesis_symbols(reparsed), hypothesis_symbols(ast))
    expect_identical(reparsed$statements[[1L]]$left$expression, ast$statements[[1L]]$left$expression)
    expect_identical(hypothesis_render(reparsed), hypothesis_render(ast))
  }
  expect_identical(hypothesis_symbols(hypothesis_parse(hypotheses[3L])), ".BayesTools_escaped_constant_Inf")
  expression <- quote(theta[".BayesToolsHypothesisLiteral1."] + 1)
  expect_identical(parse(text = .hypothesis_expression_text(expression))[[1L]], expression)
  expression <- quote(.BayesToolsHypothesisLiteral1.(theta) + 1)
  expect_identical(parse(text = .hypothesis_expression_text(expression))[[1L]], expression)
  expression <- quote(`Inf` + .BayesTools_escaped_constant_Inf + 1)
  expect_identical(.hypothesis_parse_expression(.hypothesis_expression_text(expression)), expression)
  for(text in c("theta + TRUE > 0", "theta + Inf > 0", "theta <- 0", "sin(theta) > 0")) expect_error(hypothesis_parse(text))
})

test_that("quoted reserved names and former aliases retain separate draw values", {
  draws <- data.frame(list("Inf" = c(1, 2), .BayesTools_escaped_constant_Inf = c(3, 4)), check.names = FALSE)
  expression <- .hypothesis_parse_expression(".BayesTools_escaped_constant_Inf - `Inf`")
  expect_identical(.hypothesis_eval_expression(expression, draws), c(2, 2))
  expect_identical(.hypothesis_affine_value(.hypothesis_affine_read(expression, names(draws)), draws), c(2, 2))
  ast <- hypothesis_parse("theta > 0.2")
  from_ast <- hypothesis_BF(c(-1, 1), prior("normal", list(0, 1)), ast, parameter = "theta", seed = 17, columns = "all")
  from_text <- hypothesis_BF(c(-1, 1), prior("normal", list(0, 1)), hypothesis_render(ast), parameter = "theta", seed = 17, columns = "all")
  expect_identical(attr(from_ast, "raw_BF"), attr(from_text, "raw_BF"))
  expect_identical(attr(from_ast, "raw_log_BF"), attr(from_text, "raw_log_BF"))
})

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

test_that("hypothesis rendering preserves numeric values exactly", {

  values <- c(0.123456789, 123456789.123, pi, 1 + .Machine$double.eps,
              .Machine$double.xmin, .Machine$double.xmax)
  for(value in values){
    literal <- sprintf("%.17g", value)
    for(hypothesis in c(paste("theta =", literal),
                        paste("theta >", literal),
                        paste("theta +", literal, "> 0"))){
      ast <- hypothesis_parse(hypothesis)
      reparsed <- hypothesis_parse(hypothesis_render(ast))
      expect_identical(reparsed$statements[[1L]]$left$value,
                       ast$statements[[1L]]$left$value)
      expect_identical(reparsed$statements[[1L]]$left$expression,
                       ast$statements[[1L]]$left$expression)
      expect_identical(hypothesis_render(reparsed), hypothesis_render(ast))
    }
  }

  withr::local_options(OutDec = ",")
  ast <- hypothesis_parse("theta = 0.123456789")
  expect_identical(hypothesis_render(ast), "theta = 0.123456789")
  expect_identical(hypothesis_parse(hypothesis_render(ast)), ast)
})

test_that("hypothesis labels render literals with shortest round-trip text", {

  hypothesis <- c(
    "theta > 0.2",
    "theta - 0.1 > 0.3 vs theta - 0.1 = 0.2",
    "exp(theta) > 1.1",
    "-theta^2 > -0.04",
    "theta > 1e-20",
    "theta > 0.30000000000000004"
  )
  ast <- hypothesis_parse(hypothesis)
  expect_identical(hypothesis_render(ast), hypothesis)
  expect_identical(
    vapply(ast$statements, function(statement) statement$right$label, ""),
    c("theta <= 0.2", "theta - 0.1 = 0.2", "exp(theta) <= 1.1",
      "-theta^2 <= -0.04", "theta <= 1e-20", "theta <= 0.30000000000000004")
  )
  expect_identical(hypothesis_parse(hypothesis_render(ast)), ast)

  # Exact values are retained where labels are shortened or not.
  precise <- hypothesis_parse("theta > 0.30000000000000004")$statements[[1L]]
  expect_identical(precise$left$expression$right$value, 0.1 + 0.2)
  expect_identical(
    .hypothesis_simple_parameter_comparison(precise$left, "theta")$value,
    0.1 + 0.2
  )

  rewritten <- hypothesis_rewrite(
    hypothesis_parse(c("theta > 0.2", "theta = 0.1 vs theta > 0.3")),
    c(theta = "phi")
  )
  expect_identical(
    lapply(rewritten$statements, function(statement){
      c(statement$left$label, statement$right$label)
    }),
    list(c("phi > 0.2", "phi <= 0.2"), c("phi = 0.1", "phi > 0.3"))
  )

  set.seed(1)
  out <- hypothesis_BF(
    stats::rnorm(2000, 0.3, 0.1),
    prior("normal", list(0, 1)),
    hypothesis = "theta > 0.2",
    parameter  = "theta"
  )
  expect_identical(out[["Null"]], "theta <= 0.2")
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
    simplified  = FALSE,
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
  expect_error(
    hypothesis_parse("theta > 0", simplify_names = TRUE),
    "'simplify_names = TRUE' requires 'catalog'.",
    fixed = TRUE
  )
})

test_that("catalog aliases containing brackets resolve as whole symbols", {

  data <- data.frame(
    f = factor(c("a", "b", "c", "a", "b", "c")),
    g = factor(c("u", "v", "u", "v", "u", "v")),
    x = c(1, 2, 3, 4, 5, 6)
  )
  formula_result <- JAGS_formula(~ f * g + x, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"),
    x = prior("normal", list(0, 1)),
    "f:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
  ))
  coordinates <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(is.prior.factor(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  coordinates <- .bt_build_parameter_coordinates(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )
  resolved_names <- function(hypothesis, ...){
    ast <- hypothesis_parse(hypothesis, catalog = catalog, namespace = "mu")
    unique(hypothesis_resolve(ast, catalog, namespace = "mu", ...)$
             occurrences$canonical_name)
  }

  # Exact aliases with brackets (canonical names and display labels). The
  # bracket holds the level label: `mu_f[b]` is level "b", and the backend
  # coordinate `mu_f[1]` (the same level cell) is not a public selector.
  expect_identical(resolved_names("mu_f[b] > 0"), "mu_f[b]")
  expect_error(resolved_names("mu_f[1] > 0"),
               class = "BayesTools_parameter_not_found")
  expect_identical(resolved_names("(mu) f[b] > 0"), "mu_f[b]")
  expect_identical(resolved_names("mu_f[b] > mu_f[c]"),
                   c("mu_f[b]", "mu_f[c]"))
  expect_identical(resolved_names("(mu) f[b] - (mu) f[c] = 0"),
                   c("mu_f[b]", "mu_f[c]"))
  expect_identical(resolved_names("mu_f[b] > 0", component = "b"), "mu_f[b]")
  # Parameter-level splitting remains available.
  expect_identical(resolved_names("f[c] > 0"), "mu_f[c]")

  # A non-syntactic alias followed by a level is not quoted on its own.
  for(hypothesis in c("f:g[f=b, g=v] > 0", "f:g[ f=b, g=v ] > 0")){
    ast <- hypothesis_parse(hypothesis, catalog = catalog, namespace = "mu")
    expect_identical(hypothesis_render(ast), "`f:g[f=b, g=v]` > 0",
                     info = hypothesis)
    expect_identical(resolved_names(hypothesis), "mu_f__xXx__g[f=b, g=v]",
                     info = hypothesis)
  }
  expect_identical(
    hypothesis_render(hypothesis_parse(
      "f:g[f=b,g=v] > 0", catalog = catalog, namespace = "mu"
    )),
    hypothesis_render(hypothesis_parse("f:g[f=b,g=v] > 0"))
  )
})

test_that("bracketed level names parse without splitting interaction symbols", {

  # Level names containing brackets (cut() intervals) are level references.
  levels <- hypothesis_parse_level_reference(
    c("`mu[(0,1]]`", "`mu[[0,1)]`", "`mu[[0,1]]`")
  )
  expect_identical(levels$parameter, rep("mu", 3L))
  expect_identical(levels$level, c("(0,1]", "[0,1)", "[0,1]"))
  expect_true(all(levels$direct))

  # Interaction names with several bracket groups remain single symbols.
  interactions <- c("mu_f[dif: a]__xXx__g[dif: u]",
                    "mu_a[a1]__xXx__year__xXx__b[b1]")
  ast <- hypothesis_parse(paste0("`", interactions[1L], "` > `",
                                 interactions[2L], "`"))
  expect_identical(hypothesis_symbols(ast), interactions)
  occurrences <- hypothesis_symbols(ast, occurrences = TRUE)
  expect_true(all(is.na(occurrences$level)))
  expect_identical(unique(occurrences$parameter), interactions)
  expect_identical(
    hypothesis_render(hypothesis_rewrite(ast, setNames("z", interactions[1L]))),
    paste0("z > `", interactions[2L], "`")
  )
  expect_false(any(hypothesis_parse_level_reference(
    paste0("`", interactions, "`")
  )$direct))
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

  qualified <- hypothesis_rewrite(
    hypothesis_parse("mu - fac[A] > 0"),
    c(fac = "mu")
  )
  expect_identical(hypothesis_render(qualified), "mu - `mu[A]` > 0")
  expect_identical(
    unique(hypothesis_symbols(qualified, occurrences = TRUE)$symbol),
    c("mu", "mu[A]")
  )
  expect_error(
    hypothesis_rewrite(hypothesis_parse("mu[A] - fac[A] > 0"),
                       c(fac = "mu")),
    "Rewrite mapping creates duplicate or colliding hypothesis symbols.",
    fixed = TRUE
  )
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
    simplified = c(FALSE, FALSE),
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

test_that("hypothesis resolution classes statements without parameter symbols", {

  catalog <- .bt_build_parameter_catalog(
    .bt_build_parameter_coordinates(columns = "theta")
  )
  for(hypothesis in c("1 > 0", "0 = 0", "2 > 1 & 1 > 0", "1 = 0 vs 2 > 1")){
    condition <- tryCatch(
      hypothesis_resolve(hypothesis_parse(hypothesis), catalog),
      error = identity
    )
    expect_identical(
      class(condition),
      c("BayesTools_hypothesis_no_parameters",
        "BayesTools_parameter_resolution_error", "error", "condition"),
      info = hypothesis
    )
    expect_identical(
      conditionMessage(condition),
      "The hypothesis contains no parameter symbols to resolve.",
      info = hypothesis
    )
  }
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

  precise <- hypothesis_parse("theta = 0.701406683025")
  precise_reference <- hypothesis_parse_point_reference(precise)
  expect_identical(
    precise_reference$value,
    precise$statements[[1L]]$left$value
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

  # raw prior draws: the prior ordinate is an inexact normal estimate
  expect_warning(
    from_text <- do.call(hypothesis_BF, c(arguments, list(hypothesis = text))),
    class = "BayesTools_inexact_ordinate"
  )
  expect_warning(
    from_ast <- do.call(hypothesis_BF, c(arguments, list(hypothesis = ast))),
    class = "BayesTools_inexact_ordinate"
  )
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

# Reference: the validator without its memo, which is the pre-memo behaviour.
test_that("an AST is validated once per content and modifications are rechecked", {

  .BayesTools_private$content_memo <- NULL
  withr::defer(.BayesTools_private$content_memo <- NULL)
  original <- .bt_validate_hypothesis_ast_uncached
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_validate_hypothesis_ast_uncached = function(ast){
      calls <<- calls + 1L
      original(ast)
    },
    .package = "BayesTools"
  )

  # the constructor validates, so the object is known to every later call
  ast <- hypothesis_parse(c("theta > 0", "abs(phi) + 1 = 2 vs phi != 1"))
  expect_identical(calls, 1L)
  for(i in 1:3){
    hypothesis_render(ast)
    hypothesis_symbols(ast, occurrences = TRUE)
  }
  expect_identical(calls, 1L)
  # an equal copy (as after saving and loading) is recognised too
  expect_identical(hypothesis_render(unserialize(serialize(ast, NULL))), hypothesis_render(ast))
  expect_identical(calls, 1L)

  # a modified AST is never covered by the original: invalid ones fail on every
  # call, and the memo and the unmemoized validator agree on the refusal
  wrong_source <- ast
  wrong_source$statements[[1L]]$source <- "phi > 0"
  wrong_side <- ast
  wrong_side$statements[[2L]]$left$label <- ""
  wrong_node <- ast
  wrong_node$statements[[1L]]$left$expression$type <- "unknown"
  wrong_version <- ast
  wrong_version$schema_version <- 2L
  tampered <- list(
    source = wrong_source, side = wrong_side, node = wrong_node,
    version = wrong_version
  )
  for(name in names(tampered)){
    for(i in 1:2){
      memoized <- tryCatch(hypothesis_render(tampered[[name]]), error = identity)
      unmemoized <- tryCatch(original(tampered[[name]]), error = identity)
      expect_s3_class(memoized, "error")
      expect_identical(conditionMessage(memoized), conditionMessage(unmemoized))
    }
  }
  # one check per refusal and repetition (the reference validator is not counted)
  expect_identical(calls, 1L + 2L * length(tampered))
  # the original is still recognised after the refused modifications
  before <- calls
  hypothesis_render(ast)
  expect_identical(calls, before)

  # a valid modification is checked (once) when it is made
  rewritten <- hypothesis_rewrite(ast, c(theta = "eta"))
  before <- calls
  expect_identical(
    hypothesis_render(rewritten),
    c("eta > 0", "abs(phi) + 1 = 2 vs phi != 1")
  )
  expect_identical(calls, before)
})
