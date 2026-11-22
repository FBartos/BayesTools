skip_if_not_test_profile("unit")

test_that("JAGS_fit binds package sampler declarations before creating the backend", {

  skip_if_not_installed("rjags")
  called <- 0L
  supplied <- NULL
  backend <- NULL
  testthat::local_mocked_bindings(
    .JAGS_load_modules = function(...) invisible(TRUE),
    .bt_rjags_fit = function(...){
      backend <<- list(...)
      simpleError("Captured the backend input")
    },
    .package = "BayesTools"
  )
  declare <- function(data, formula_design, prior_list){
    called <<- called + 1L
    supplied <<- list(priors = prior_list, design = formula_design, data = data)
    list(markers = list(list(name = "Package_joint_marker", distribution = "dpackage_marker",
      arguments = list(c("mu_intercept", "mu_x"), "Package_joint_data", c(1, 2)))),
      data = list(Package_joint_data = matrix(1:6, 2)))
  }
  observed <- c(1, 2, 4)
  normal <- prior("normal", list(0, 1))
  fit <- JAGS_fit(
    "model { for(i in 1:N_mu) { y[i] ~ dnorm(mu[i], tau) } }",
    data = list(y = observed), prior_list = list(tau = prior("gamma", list(1, 1))),
    formula_list = list(mu = ~ x), formula_data_list = list(mu = data.frame(x = c(0, 1, 3))),
    formula_prior_list = list(mu = list(intercept = normal, x = normal)),
    sampler_joint = declare, jags_modules = "Package", seed = 1L
  )
  expect_identical(called, 1L)
  expect_named(supplied$priors, c("mu_intercept", "mu_x", "tau"))
  expect_named(supplied$design, "mu")
  expect_identical(supplied$data$y, observed)
  expect_true("N_mu" %in% names(supplied$data))
  expect_identical(backend$data$Package_joint_marker, 0)
  expect_identical(backend$data$Package_joint_data, matrix(1:6, 2))
  expect_match(backend$model,
    "Package_joint_marker ~ dpackage_marker(c(mu_intercept, mu_x), Package_joint_data, c(1, 2))",
    fixed = TRUE)
  expect_false(grepl("Package_joint_marker", attr(fit, "model_syntax"), fixed = TRUE))
  declarations <- attr(fit, "sampler_declarations", exact = TRUE)
  expect_identical(declarations$joint$data, "Package_joint_data")
  expect_identical(declarations$joint$markers[[1L]]$name, "Package_joint_marker")
})

test_that("package markers leave the scientific gate scan intact", {

  priors <- list(bias = prior_mixture(
    list(prior_none(), prior_PET("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  ))
  data <- list(y = 0)
  syntax <- JAGS_add_priors("model { y ~ dnorm(PET, 1) }", priors)
  legacy <- .bt_sampler_declarations(syntax, data, priors, NULL, sampler_gates = TRUE)
  plan <- .bt_sampler_apply_plan(legacy, NULL, priors, NULL, data,
    sampler_joint = function(prior_list, formula_design, data){
      list(markers = list(list(name = "Package_joint_marker", distribution = "dbt_block",
        arguments = list("PET_1", 1))))
    },
    sampler_probabilities = "bias_indicator")
  expect_identical(plan$declarations$gates, legacy$declarations$gates)
  expect_identical(plan$declarations$gates$node, "PET_1")
  expect_match(plan$model, "BayesTools_finite_1 <- bt_finite_tag(bias_indicator)", fixed = TRUE)
  expect_identical(plan$declarations$probabilities$readout, "BayesTools_finite_readout_1")
})

test_that("invalid sampler declarations fail before JAGS compilation", {

  expect_error(.bt_sampler_check_plan_arguments(1, NULL), class = "BayesTools_sampler_declaration_error")
  for(nodes in list(c("g", "g "), "g[1:2]", "g[3 2]", "g\n", NA_character_)){
    expect_error(.bt_sampler_check_plan_arguments(NULL, nodes), class = "BayesTools_sampler_declaration_error")
  }
  expect_identical(.bt_sampler_check_plan_arguments(NULL, c("a", "b [ 2, 1 ]")), c("a", "b[2,1]"))

  marker <- list(name = "package_marker", distribution = "dpackage_marker", arguments = list("x", 1))
  check <- function(fields){
    marker[names(fields)] <- fields
    .bt_sampler_joint_markers(list(markers = list(marker)), c("x", "y"), "Package")
  }
  for(bad in list("y", "BayesTools_finite_readout_1", "mark\n", "x <- 1")){
    expect_error(check(list(name = bad)), class = "BayesTools_sampler_declaration_error")
  }
  for(bad in list("equals(x, 1)", "unknown", x ~ y, Inf, list(x = 2))){
    expect_error(check(list(arguments = list(bad))), class = "BayesTools_sampler_declaration_error")
  }
  expect_error(check(list(distribution = "dnorm")),
    "must be a BayesTools marker or use the prefix", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker, marker)), c("x", "y"), "Package"),
    "used more than once", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker), data = list(x = 0)), c("x", "y"), "Package"),
    "already used or reserved", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker)), c("x", "y"), "Package", retired = "x"),
    "use 'sampler_nested_weights' = FALSE", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker), data = list(spare = 1)), c("x", "y"), "Package"),
    "data 'spare' is not referenced by any marker", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
})

test_that("sampler references come from model relations and inputs", {

  syntax <- paste(
    "data { n ~ dpois(4) }", "model {",
    "for(i in 1:2) { logit(prob[i]) <- theta }",
    "y ~ dbin(prob[1], n)", "theta ~ dnorm(0, 1)", "}", sep = "\n"
  )
  data <- list(y = 1, known_constant = 9)
  prior_list <- list(bias = prior("normal", list(0, 1)))
  legacy <- list(model = syntax, data = data)
  declare <- function(node){
    function(prior_list, formula_design, data){
      list(markers = list(list(name = "trace_marker", distribution = "dpackage_marker",
        arguments = list(node, "n", "known_constant"))))
    }
  }
  plan <- function(nodes = NULL, joint = NULL){
    .bt_sampler_apply_plan(legacy, "undeclared", prior_list, NULL, data,
      sampler_joint = joint, sampler_probabilities = nodes, jags_modules = "Package")
  }
  for(node in c("bias", "dnorm", "i", "logit", "undeclared")){
    expect_error(plan(nodes = node), "unavailable in the model", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
    expect_error(plan(joint = declare(node)), paste0("unknown node '", node, "'"), fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }
  accepted <- plan(nodes = c("prob[1]", "n", "theta"), joint = declare("prob[1:2]"))
  expect_identical(accepted$declarations$probabilities$node, c("prob[1]", "n", "theta"))
  expect_match(accepted$model, "trace_marker ~ dpackage_marker(prob[1:2], n, known_constant)", fixed = TRUE)
})

test_that("callback and ordinary FactorSlice block markers cannot share one fit", {

  priors <- stats::setNames(rep(list(prior("normal", list(0, 1))), 3L), c("x", "w", "z"))
  data <- list(y = 1)
  syntax <- JAGS_add_priors("model { alias <- x\n y ~ dnorm(alias + w + z, 1) }", priors)
  legacy <- .bt_sampler_declarations(syntax, data, priors, NULL, sampler_blocks = list(c("x", "w")))
  declare <- function(target, id, distribution = "dbt_block"){
    function(prior_list, formula_design, data){
      list(markers = list(list(name = "package_marker", distribution = distribution,
        arguments = list(target, id))))
    }
  }
  plan <- function(joint) .bt_sampler_apply_plan(legacy, NULL, priors, NULL, data,
    sampler_joint = joint, sampler_probabilities = NULL, jags_modules = "Package")
  expect_error(plan(declare("alias", 2L)), "blocks from 'sampler_blocks'", fixed = TRUE,
    class = "BayesTools_sampler_declaration_error")
  expect_error(plan(declare("z", 1L)), "blocks from 'sampler_blocks'", fixed = TRUE,
    class = "BayesTools_sampler_declaration_error")
  expect_true(plan(declare("z", 1L, "dpackage_marker"))$uses_module)
})

test_that("typed callbacks cannot bypass the scientific gate check", {

  callback <- function(prior_list, formula_design, data){
    list(markers = list(list(name = "manual_gate", distribution = "dbt_gate",
      arguments = list("x", "open"))))
  }
  expect_error(JAGS_fit(
    "model { open ~ dbern(0.99)\n y ~ dnorm(x, 100) }", data = list(y = 10),
    prior_list = list(x = prior("normal", list(0, 1))),
    sampler_joint = callback, jags_modules = "bt",
    runtime_setup = function(context) stop("Unexpected backend setup")
  ), "distribution 'dbt_gate' is reserved for gates verified from the model. Use 'sampler_gates' to enable those gates.",
    fixed = TRUE, class = "BayesTools_sampler_declaration_error")
})

test_that("package tags are terminal declarations with no observed marker data", {

  tag <- list(name = "package_terminal", tag = "package_zero", arguments = list("x"))
  marker <- list(name = "package_observed", distribution = "dpackage_marker",
    arguments = list("x", "geometry"))
  legacy <- list(model = "model { x ~ dnorm(0, 1) }", data = NULL)
  plan <- function(markers, payload = NULL, model = legacy$model, monitor = NULL){
    legacy$model <- model
    .bt_sampler_apply_plan(legacy, monitor, NULL, NULL, NULL,
      sampler_joint = function(prior_list, formula_design, data){
        list(markers = markers, data = payload)
      },
      sampler_probabilities = NULL, jags_modules = "Package")
  }

  tag_only <- plan(list(tag))
  expect_identical(tag_only$data, legacy$data)
  expect_match(tag_only$model, "package_terminal <- package_zero(x)", fixed = TRUE)
  expect_identical(tag_only$declarations$joint$markers, list(tag))
  combined <- plan(list(tag, marker), list(geometry = 2))
  expect_identical(combined$data, list(geometry = 2, package_observed = 0))
  expect_match(combined$model, "package_observed ~ dpackage_marker(x, geometry)", fixed = TRUE)
  expect_identical(combined$declarations$joint$markers, list(tag, marker))

  for(symbol in c("other_zero", "Package_zero", "exp")){
    wrong_provider <- tag
    wrong_provider$tag <- symbol
    expect_error(plan(list(wrong_provider)), "must use the lowercase prefix of a provider",
      fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  }
  expect_error(plan(list(c(tag, list(distribution = "dbt_gate")))), "exactly one of 'distribution' and 'tag'",
    fixed = TRUE, class = "BayesTools_sampler_declaration_error")

  marker$arguments <- list("package_terminal")
  for(markers in list(list(tag, marker), list(marker, tag))){
    expect_error(plan(markers), "unknown node 'package_terminal'", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }
  expect_error(plan(list(tag, tag)), "is used more than once", fixed = TRUE,
    class = "BayesTools_sampler_declaration_error")
  expect_error(plan(list(tag), payload = list(package_terminal = 0)), "is used more than once",
    fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(plan(list(tag), model = "model { x ~ dnorm(package_terminal, 1) }"),
    "is already used or reserved", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(plan(list(tag), monitor = " package_terminal [1] "),
    "is already used or reserved", fixed = TRUE, class = "BayesTools_sampler_declaration_error")

  for(line_end in c("\n", "\r\n", "\r")){
    calling <- paste0("model { x ~ dnorm(0, 1)\n# earlier statement", line_end,
      " y <- package_zero /* bring the argument next */ (x)\n}")
    expect_error(plan(list(tag), model = calling),
      "tag function 'package_zero' is already called by the model.", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }
  separate <- paste(
    "model { x ~ dnorm(0, 1)",
    "# package_zero(x)",
    "/* package_zero(x) */ y <- package_other_zero(x)",
    "}", "/* package_zero(x) }", sep = "\n"
  )
  expect_match(plan(list(tag), model = separate)$model,
    "package_terminal <- package_zero(x)", fixed = TRUE)
  for(symbol in c("bt_finite_tag", "bt_finite_witness_tag")){
    reserved <- tag
    reserved$tag <- symbol
    expect_error(.bt_sampler_joint_markers(list(markers = list(reserved)), "x", "bt"),
      paste0("tag function '", symbol, "' is reserved by BayesTools"), fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }

  for(keyword in c("model", "data", "var", "in", "for")){
    reserved <- tag
    reserved$name <- keyword
    expect_error(plan(list(reserved)), class = "BayesTools_sampler_declaration_error")
    expect_error(plan(list(tag), payload = stats::setNames(list(0), keyword)),
      class = "BayesTools_sampler_declaration_error")
    expect_error(.bt_sampler_check_plan_arguments(NULL, paste0(keyword, "[1]")),
      class = "BayesTools_sampler_declaration_error")
  }
})

test_that("comment masking preserves native line breaks and source offsets", {

  heading <- "# comment with a closing token }"
  for(line_end in c("\n", "\r\n", "\r")){
    suffix <- paste0(line_end, "z ~ dbern(0.5)")
    expect_identical(.JAGS_mask_comments(paste0(heading, suffix)),
      paste0(strrep(" ", nchar(heading)), suffix))
    model <- "model { z ~ dbern(0.5) }"
    trailing <- paste0("/* Neither this brace } nor", line_end, "ghost ~ dbern(0.5) is code")
    expect_identical(.JAGS_mask_comments(paste0(model, trailing)),
      paste0(model, strrep(" ", nchar("/* Neither this brace } nor")), line_end,
        strrep(" ", nchar("ghost ~ dbern(0.5) is code"))))
    expect_identical(.bt_jags_syntax_graph(paste0(model, trailing))$stochastic, "z")
  }
})

test_that("JAGS special operators are opaque to structural syntax scans", {

  operators <- c("%#%", "%/*%", "%}%", "%model{%", "%PET_1%", "%]~(y)%", "%*%",
    "%%#%", "%%%#%", "%%%%", paste0("%x\fy%"), paste0("%x", intToUtf8(0xA0), "y%"))
  for(line_end in c("\n", "\r\n", "\r")){
    for(operator in operators){
      hash <- paste0("# ", operator, " hidden }")
      block <- paste0("/* ", operator, " hidden } */")
      source <- paste0(hash, line_end, "model { x <- a ", operator, " b ", block,
        line_end, "y ~ dnorm(x, 1) }")
      opaque <- strrep("%", nchar(operator))
      expect_identical(.JAGS_mask_comments(source), paste0(strrep(" ", nchar(hash)),
        line_end, "model { x <- a ", opaque, " b ", strrep(" ", nchar(block)),
        line_end, "y ~ dnorm(x, 1) }"))
      graph <- .bt_jags_syntax_graph(source)
      expect_identical(graph$deterministic$x, c("a", "b"))
      expect_identical(graph$stochastic, "y")
      expect_identical(.bt_sampler_tokens(source), c("model", "x", "a", "b", "y", "dnorm", "x"))
    }
  }
  # When no final percent follows the body, the longest valid token ends
  # within the leading run. The following hash begins an ordinary comment.
  expect_identical(.JAGS_mask_comments("a %%%# gone\r b"), "a %%%      \r b")

  source <- paste(
    "data { precision <- a %model{% b %}% c }",
    "model {",
    "  for(i in 1:2) { g[i] ~ dbern(p[i %]~(% 1]) }",
    "}", sep = "\n"
  )
  opening <- regexpr("model {", source, fixed = TRUE)[1L] + nchar("model {") - 1L
  expect_identical(.JAGS_model_block_opening(source), as.integer(opening))
  graph <- .bt_jags_syntax_graph(source)
  expect_identical(graph$data, "precision")
  expect_identical(graph$stochastic, "g")
  marked <- .bt_sampler_insert_markers(source, "mark <- 0")
  expect_identical(marked, sub("\\}$", "  mark <- 0\n}", source))
})

test_that("special operators keep real gate dependencies and tag calls visible", {

  prior_list <- list(bias = prior_mixture(list(prior_none(),
    prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE)))
  gate_plan <- function(node){
    syntax <- paste0("model { precision <- 96 %PET_1% 4; y ~ dnorm(", node, ", precision) }")
    .bt_sampler_declarations(JAGS_add_priors(syntax, prior_list), list(y = 10),
      prior_list, NULL, sampler_gates = TRUE)
  }
  expect_identical(gate_plan("PET")$declarations$gates$node, "PET_1")
  expect_equal(nrow(gate_plan("PET_1")$declarations$gates), 0L)

  declaration <- function(prior_list, formula_design, data){
    list(markers = list(list(name = "private_terminal", tag = "package_zero", arguments = list("x"))))
  }
  plan <- function(operator, relation){
    legacy <- list(model = paste0("model { x ~ dnorm(0, 1); precision <- 96 ", operator, " 4; ", relation, " }"))
    .bt_sampler_apply_plan(legacy, NULL, NULL, NULL, NULL, declaration, NULL, "Package")
  }
  for(operator in c("%#%", "%/*%", "%%#%", "%}%")){
    expect_error(plan(operator, "scientific <- package_zero(x)"),
      "tag function 'package_zero' is already called by the model", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
    expect_error(plan(operator, "y ~ dnorm(x + private_terminal, precision)"),
      "name 'private_terminal' is already used or reserved", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }
  expect_match(plan("%private_terminal%", "y ~ dnorm(x, precision)")$model,
    "private_terminal <- package_zero(x)", fixed = TRUE)
})

test_that("explicit probability requests reserve the private finite functions", {

  plan <- function(syntax, nodes){
    legacy <- list(model = syntax, data = list(y = 0.2))
    .bt_sampler_apply_plan(legacy, NULL, NULL, NULL, NULL, NULL, nodes)
  }
  for(call in c("bt_finite_witness_tag(z, w)", "bt_finite_tag(z2)")){
    syntax <- paste0("model { z ~ dbern(0.5); y ~ dnorm(z, 1); extra <- ", call, " }")
    expect_error(plan(syntax, "z"), "The model calls the private finite readout function",
      fixed = TRUE, class = "BayesTools_sampler_declaration_error")
    expect_identical(plan(syntax, NULL)$model, syntax)
  }
  syntax <- "model { z ~ dbern(0.5)\n# extra <- bt_finite_tag(z2)\n y ~ dnorm(z, 1) }"
  expect_match(plan(syntax, "z")$model, "BayesTools_finite_1 <- bt_finite_tag(z)", fixed = TRUE)
})

test_that("parallel merge checks probability schema by indicator and reports mixed gate ownership", {

  make_worker <- function(ids, start, gate_sampler){
    draws <- function(values) coda::mcmc.list(lapply(ids, function(i){
      coda::mcmc(values(i), start = start, thin = 2L)
    }))
    list(
      mcmc = draws(function(i) cbind(x = rep(i, 3L))),
      end_state = lapply(ids, function(id) list(x = id)),
      samplers = data.frame(Index.No = c(1L, 2L, 3L),
        Sampler = c(gate_sampler, "BayesTools::GatedSlice", "BayesTools::Finite"),
        Node = c("x", "w", "b")),
      sampler_probabilities = list(b = draws(function(i){
        p <- rep(i / 4, 3L)
        cbind("p[0]" = p, "p[1]" = 1 - p, "log_p[0]" = log(p),
          "log_p[1]" = log1p(-p), selected = 0, lower = 0)
      }))
    )
  }
  groups <- list(c(1L, 3L), 2L)
  workers <- list(make_worker(groups[[1L]], 101L, "BayesTools::GatedSlice"),
                  make_worker(groups[[2L]], 51L, "base::RealSlicer"))
  plan <- data.frame(node = "b")
  fit <- .bt_rjags_merge_workers(workers, groups, 2L, plan)
  expect_identical(fit$samplers, workers[[1L]]$samplers)
  expect_identical(lapply(fit$samplers_by_worker, "[[", "chains"), groups)
  expect_identical(lapply(fit$samplers_by_worker, "[[", "samplers"), lapply(workers, "[[", "samplers"))
  expect_identical(vapply(fit$mcmc, function(chain) chain[1L, "x"], integer(1)), 1:3)
  expect_identical(vapply(fit$sampler_probabilities$b, function(chain) 4 * chain[1L, "p[0]"], numeric(1)), c(1, 2, 3))
  expect_identical(lapply(fit$mcmc, stats::time), lapply(fit$sampler_probabilities$b, stats::time))

  gates <- data.frame(marker = c("first", "second"), node = c("x", "w"),
    branch = NA_integer_, prior = c("x", "w"))
  claims <- .bt_sampler_fit_gate_claims(gates, NULL, fit)
  expect_identical(claims$status, c("mixed", "claimed"))
  expect_identical(claims$sampler, c(NA_character_, "BayesTools::GatedSlice"))
  fit_without_readouts <- .bt_rjags_merge_workers(workers, groups, 2L, NULL)
  expect_identical(.bt_sampler_fit_gate_claims(gates, NULL, fit_without_readouts), claims)

  incompatible <- workers
  colnames(incompatible[[2L]]$sampler_probabilities$b[[1L]])[1L] <- "p[3]"
  expect_error(.bt_rjags_merge_workers(incompatible, groups, 2L, plan),
    "incompatible probability traces for 'b'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
  incompatible[[2L]]$sampler_probabilities["b"] <- list(NULL)
  expect_error(.bt_rjags_merge_workers(incompatible, groups, 2L, plan),
    "incompatible probability traces for 'b'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
})

test_that("readout decoding distinguishes an unowned node from an incomplete trace", {

  readout <- c(alias = "BayesTools_finite_readout_1", declined = "BayesTools_finite_readout_2")
  row <- c(0.25, 0.75, log(0.25), log(0.75), 1, 0, -1)
  values <- matrix(rep(row, 3L), nrow = 3L, byrow = TRUE,
    dimnames = list(NULL, c(paste0(readout[[1L]], "[", 1:6, "]"), readout[[2L]])))
  draws <- coda::mcmc.list(coda::mcmc(values, start = 20L, thin = 2L))
  decoded <- .bt_rjags_probability_draws(draws, readout)
  expect_identical(names(decoded), names(readout))
  expect_s3_class(decoded$alias, "mcmc.list")
  expect_null(decoded$declined)
  expect_identical(stats::time(decoded$alias[[1L]]), stats::time(draws[[1L]]))
  values[2L, 7L] <- 1
  invalid <- coda::mcmc.list(coda::mcmc(values))
  expect_error(.bt_rjags_probability_draws(invalid, readout),
    "incomplete probability trace for 'declined'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
})
