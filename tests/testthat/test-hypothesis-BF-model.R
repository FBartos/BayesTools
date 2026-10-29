skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Fit-level hypothesis Bayes factors (the engine, the provider, and
# the density backend registry)
# ============================================================================ #
#
# The fits are synthetic: a compiled formula with a scaled predictor and a
# treatment factor, with posterior draws sampled from known normals. Cached
# MCMC fits are used in test-hypothesis-BF-model-fixture.R.
#
# ============================================================================ #

.hypothesis_model_test_fit <- function(seed = 3L, n = 400L){

  set.seed(seed)
  data <- data.frame(x = c(-2, -1, 0, 3, 5, 7), f = factor(rep(c("A", "B", "C"), 2)))
  priors <- list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1)),
    f         = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  compiled <- JAGS_formula(~ 1 + x + f, "mu", data, priors, formula_scale = TRUE)
  columns  <- c("mu_intercept", "mu_x", "mu_f[1]", "mu_f[2]")
  draws <- cbind(
    matrix(stats::rnorm(n * length(columns), 0.3, 0.5), n, length(columns),
           dimnames = list(NULL, columns)),
    sigma = abs(stats::rnorm(n, 1, 0.2))
  )
  all_priors <- c(compiled$prior_list,
                  list(sigma = prior("normal", list(0, 1), list(0, Inf))))

  .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(draws)), all_priors,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale)
  )
}

# The extraction a user would write by hand for a fixed coefficient.
.hypothesis_model_manual_coefficient <- function(fit, parameter, hypothesis,
                                                 seed = NULL, n = 10000L, ...){

  samples <- as_mixed_posteriors(fit, parameter, n_prior_samples = n)
  posterior <- marginal_posterior(
    samples, parameter, prior_samples = TRUE, use_formula = FALSE, n_samples = n
  )
  hypothesis_BF.default(posterior, hypothesis = hypothesis, parameter = parameter,
                        seed = seed, ...)
}

# A toy ordinate backend: the exact normal ordinate of the draws.
.hypothesis_model_toy_backend <- function(name = "toyNormal", state = new.env(),
                                          ...){

  state$attached <- 0L
  state$requests <- list()
  hypothesis_density_backend(
    name = name, type = "ordinate",
    attach = function(posterior, request, control, cache){
      state$attached <- state$attached + 1L
      state$requests[[length(state$requests) + 1L]] <- request
      sample <- if(identical(request$kind, "level")) posterior[[request$level]] else posterior
      draws  <- as.numeric(sample)
      for(value in request$values){
        if(isTRUE(value == 99)){
          stop(structure(class = c("toy_failure", "error", "condition"),
                         list(message = "toy failure", call = NULL)))
        }
        entry <- posterior_ordinate_attribute(
          value = value,
          ordinate = stats::dnorm(value, mean(draws), stats::sd(draws)),
          method = "toy", density_method = name,
          diagnostics = list(BF_error_percent = 1,
                             warning = paste0("toy warning at ", value))
        )
        posterior_metadata(sample, "posterior_ordinate") <- posterior_ordinate_append(
          posterior_metadata(sample, "posterior_ordinate"), entry
        )
      }
      if(identical(request$kind, "level")){
        posterior[[request$level]] <- sample
        posterior
      }else{
        sample
      }
    },
    ...
  )
}

with_toy_backend <- function(backend, code, env = parent.frame()){

  register_hypothesis_density_backend(backend, overwrite = TRUE)
  withr::defer(unregister_hypothesis_density_backend(backend$name), envir = env)
  force(code)
}


test_that("hypothesis_BF is a generic and the draw-based method is unchanged", {

  set.seed(1)
  draws <- data.frame(theta = stats::rnorm(500, 0.4, 0.3))
  prior_draws <- data.frame(theta = stats::rnorm(500, 0, 1))
  direct  <- hypothesis_BF(draws, prior_draws, "theta > 0", seed = 1)
  by_name <- hypothesis_BF.default(posterior = draws, prior = prior_draws,
                                   hypothesis = "theta > 0", seed = 1)
  expect_identical(attr(direct, "raw_log_BF"), attr(by_name, "raw_log_BF"))
  expect_error(
    hypothesis_BF(draws, prior_draws, "theta > 0", densty_method = "KDE"),
    "unused argument"
  )
  expect_true(is.function(getS3method("hypothesis_BF", "BayesTools_fit")))
  expect_true(is.function(getS3method("hypothesis_BF", "default")))

})


test_that("the fit method equals the three-step extraction of a coefficient", {

  fit <- .hypothesis_model_test_fit()
  engine <- hypothesis_BF(fit, "x = 0.1", seed = 1)
  manual <- .hypothesis_model_manual_coefficient(fit, "mu_x", "mu_x = 0.1", seed = 1)
  expect_s3_class(engine, "BayesTools_hypothesis_BF")
  expect_identical(attr(engine, "raw_log_BF"), attr(manual, "raw_log_BF"))
  expect_identical(attr(engine, "raw_BF"), attr(manual, "raw_BF"))
  expect_identical(engine$Null, "x = 0.1")
  expect_identical(rownames(engine), "(mu) x")

  # public names need no backticks, with the same result under every alias
  for(name in c("mu_x", "(mu) x")){
    again <- hypothesis_BF(fit, paste0("`", name, "` = 0.1"), seed = 1)
    expect_identical(attr(again, "raw_log_BF"), attr(engine, "raw_log_BF"), info = name)
  }
  scalar <- hypothesis_BF(fit, "sigma > 1", seed = 1)
  expect_identical(rownames(scalar), "sigma")

})


test_that("transform_scaled tests coefficients on the original scale", {

  fit <- .hypothesis_model_test_fit()
  slope <- JAGS_formula_coefficient_transform(fit, "mu")$matrix["mu_x", "mu_x"]
  original <- hypothesis_BF(fit, "x = 0.1", transform_scaled = TRUE, seed = 1)
  fitted   <- hypothesis_BF(fit, paste0("x = ", format(0.1 / slope, digits = 17)), seed = 1)
  expect_equal(attr(original, "raw_log_BF"), attr(fitted, "raw_log_BF"),
               tolerance = 1e-6)
  default <- hypothesis_BF(fit, "x = 0.1", seed = 1)
  expect_false(isTRUE(all.equal(attr(default, "raw_log_BF"), attr(original, "raw_log_BF"))))

})


test_that("factor levels, vs statements, and several statements bind in order", {

  fit <- .hypothesis_model_test_fit()
  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 10000)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 10000)
  engine <- hypothesis_BF(fit, "f[B] = 0.1 vs f[B] > 0.1", seed = 1)
  manual <- hypothesis_BF.default(levels, hypothesis = "mu_f[B] = 0.1 vs mu_f[B] > 0.1",
                                  parameter = "mu_f", seed = 1)
  expect_identical(attr(engine, "raw_log_BF"), attr(manual, "raw_log_BF"))
  expect_identical(engine$Alternative, "f[B] = 0.1")
  expect_identical(engine$Null, "f[B] > 0.1")

  statements <- c("x > 0", "f[B] > f[C]", "x < 1", "sigma > 1", "f[B] > 0")
  table <- hypothesis_BF(fit, statements, seed = 1)
  expect_identical(nrow(table), 5L)
  expect_identical(rownames(table), c("(mu) x (1)", "f (2)", "(mu) x (3)", "sigma", "f (5)"))
  single <- vapply(statements, function(statement){
    attr(hypothesis_BF(fit, statement, seed = 1), "raw_log_BF")
  }, numeric(1))
  expect_equal(unname(attr(table, "raw_log_BF")), unname(single))
  expect_identical(attr(table, "hypothesis_ast")$statements |> length(), 5L)

  logged <- hypothesis_BF(fit, statements, logBF = TRUE, BF01 = TRUE, seed = 1)
  expect_equal(as.numeric(logged$BF), -unname(single), tolerance = 1e-12)

})


test_that("refused statements stop with the classes of the engine", {

  fit <- .hypothesis_model_test_fit()
  condition <- function(code) tryCatch(code, condition = identity)

  unknown <- condition(hypothesis_BF(fit, "nonsense = 0"))
  expect_s3_class(unknown, "BayesTools_hypothesis_statement")
  expect_s3_class(unknown, "BayesTools_parameter_not_found")
  expect_true(all(c("alias", "available") %in% names(unclass(unknown))))

  multiple <- condition(hypothesis_BF(fit, "x > f[B]"))
  expect_s3_class(multiple, "BayesTools_hypothesis_ambiguous")
  expect_s3_class(multiple, "BayesTools_hypothesis_statement")

  fixed <- condition(hypothesis_BF(fit, "f[A] = 0"))
  expect_s3_class(fixed, "BayesTools_hypothesis_fixed")
  expect_s3_class(fixed, "BayesTools_hypothesis_unavailable")

  direct <- condition(hypothesis_BF(fit, "2 * x = 0"))
  expect_s3_class(direct, "BayesTools_hypothesis_statement")

  contrast <- condition(hypothesis_BF(fit, "mu_f{1} = 0"))
  expect_s3_class(contrast, "BayesTools_hypothesis_statement")

  mismatch <- condition(hypothesis_BF(fit, "x > 0", namespace = "model"))
  expect_s3_class(mismatch, "BayesTools_hypothesis_component_mismatch")

  expect_error(hypothesis_BF(fit, "x = 0", density_method = "precomputed"),
               "must be one of 'KDE', 'normal'")
  expect_error(hypothesis_BF(fit, "x = 0", density_method = "nothing"),
               "must be one of")
  expect_error(hypothesis_BF(fit, "x = 0", density_control = list(a = 1)),
               "not used by density_method = 'KDE'")
  expect_error(hypothesis_BF(fit, "x = 0", namespace = "mu", component = "A"),
               "not both")
  expect_warning(
    unused <- hypothesis_BF(fit, "x = 0", seed = 1, nonsense = 1),
    "Unused argument in hypothesis_BF\\(\\): 'nonsense'"
  )
  expect_identical(
    attr(unused, "raw_log_BF"),
    attr(hypothesis_BF(fit, "x = 0", seed = 1), "raw_log_BF")
  )

})


test_that("hypothesis_stop raises the classes of the refusals", {

  for(type in c("target", "fixed", "method")){
    error <- tryCatch(hypothesis_stop("reason", type), error = identity)
    expect_s3_class(error, paste0("BayesTools_hypothesis_", type))
    expect_s3_class(error, "BayesTools_hypothesis_unavailable")
    expect_identical(conditionMessage(error), "reason")
  }
  statement <- tryCatch(hypothesis_stop("restate", "statement", fields = list(a = 1)),
                        error = identity)
  expect_identical(class(statement)[1:2], c("BayesTools_hypothesis_statement", "error"))
  expect_identical(statement$a, 1)
  ambiguous <- tryCatch(hypothesis_stop("which", "ambiguous", class = "extra"),
                        error = identity)
  expect_identical(class(ambiguous)[1:4], c("BayesTools_hypothesis_ambiguous",
    "BayesTools_hypothesis_statement", "extra", "error"))
  expect_error(hypothesis_stop("reason", "nonsense"))

})


test_that("the registry keeps built-ins and validates registrations", {

  backends <- hypothesis_density_backends()
  expect_identical(backends$name[1:3], c("KDE", "normal", "precomputed"))
  expect_identical(backends$type[1:3], c("draws", "draws", "stored"))
  expect_identical(backends$advertised[1:3], c(TRUE, FALSE, FALSE))

  toy <- .hypothesis_model_toy_backend()
  expect_s3_class(toy, "BayesTools_hypothesis_density_backend")
  expect_error(register_hypothesis_density_backend("toy"), "hypothesis_density_backend")
  expect_error(hypothesis_density_backend("1bad", "ordinate", attach = function(...) NULL),
               "start with a letter")
  expect_error(hypothesis_density_backend("noattach", "ordinate"), "needs an 'attach'")
  expect_error(hypothesis_density_backend("draw", "draws", attach = function(...) NULL),
               "Only an ordinate backend")
  expect_error(hypothesis_density_backend("draw", "draws", method = "other"))
  expect_error(hypothesis_density_backend("x", "ordinate", attach = function(...) NULL,
                                          reuse = 1), "must be a function")

  register_hypothesis_density_backend(toy)
  withr::defer(unregister_hypothesis_density_backend("toyNormal"))
  expect_true("toyNormal" %in% hypothesis_density_backends()$name)
  expect_error(register_hypothesis_density_backend(toy), "registered already")
  expect_s3_class(register_hypothesis_density_backend(toy, overwrite = TRUE),
                  "BayesTools_hypothesis_density_backend")
  expect_error(
    register_hypothesis_density_backend(.hypothesis_model_toy_backend("TOYNORMAL")),
    "differs from the registered density backend 'toyNormal' by case only"
  )
  expect_error(
    register_hypothesis_density_backend(
      hypothesis_density_backend("KDE", "draws", method = "KDE"), overwrite = TRUE),
    "built-in"
  )
  expect_error(unregister_hypothesis_density_backend("KDE"), "built-in")
  expect_true(unregister_hypothesis_density_backend("toyNormal"))
  expect_false(unregister_hypothesis_density_backend("toyNormal"))
  expect_false("toyNormal" %in% hypothesis_density_backends()$name)

})


test_that("density_method matching ignores case and accepts partial names", {

  fit <- .hypothesis_model_test_fit()
  toy <- .hypothesis_model_toy_backend()
  with_toy_backend(toy, {
    exact   <- hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal", seed = 1)
    for(name in c("TOYNORMAL", "toy", "T")){
      again <- hypothesis_BF(fit, "x = 0.1", density_method = name, seed = 1)
      expect_identical(attr(again, "raw_log_BF"), attr(exact, "raw_log_BF"), info = name)
    }
    expect_identical(
      attr(hypothesis_BF(fit, "x = 0.1", density_method = "kde", seed = 1), "raw_log_BF"),
      attr(hypothesis_BF(fit, "x = 0.1", seed = 1), "raw_log_BF")
    )
    # 'n' selects 'normal' only while it is the only registered name that starts with n
    expect_error(hypothesis_BF(fit, "x = 0.1", density_method = "zz"), "must be one of")
  })

})


test_that("an ordinate backend attaches ordinates that the draws evaluate", {

  fit <- .hypothesis_model_test_fit()
  state <- new.env()
  toy <- .hypothesis_model_toy_backend(state = state, policy = list(error = c(0.025, 0.05)),
                                       diagnostics = function(posterior, request){
                                         data.frame(method = "toy", n = 1L)
                                       })
  with_toy_backend(toy, {
    expect_identical(hypothesis_density_backends()$policy[[4]]$error, c(0.025, 0.05))
    toy_table <- hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal", seed = 1)
    normal    <- hypothesis_BF(fit, "x = 0.1", density_method = "normal", seed = 1)
    # the exact normal ordinate is the normal approximation of the draws
    expect_equal(attr(toy_table, "raw_log_BF"), attr(normal, "raw_log_BF"), tolerance = 1e-8)
    expect_identical(state$attached, 1L)
    request <- state$requests[[1L]]
    expect_identical(request$kind, "scalar")
    expect_identical(request$route, "scalar")
    expect_identical(request$parameter, "mu_x")
    expect_identical(request$label, "(mu) x")
    expect_identical(request$values, 0.1)
    expect_identical(request$density_method, "toyNormal")
    expect_null(request$conditional)
    expect_identical(request$spec$type, "primitive")

    # warnings and diagnostics are lifted onto the rows of the statements
    expect_identical(unname(attr(toy_table, "warnings")), "toy warning at 0.1")
    expect_identical(names(attr(toy_table, "warnings")), "(mu) x")
    expect_identical(attr(toy_table, "density_diagnostics")$method, "toy")
    expect_identical(as.numeric(toy_table$BF_error), 1)

    # region statements are method-free: nothing is attached
    state$attached <- 0L
    region <- hypothesis_BF(fit, "x > 0", density_method = "toyNormal", seed = 1)
    expect_identical(state$attached, 0L)
    expect_identical(attr(region, "raw_log_BF"),
                     attr(hypothesis_BF(fit, "x > 0", seed = 1), "raw_log_BF"))
    expect_null(attr(region, "density_diagnostics"))

    # levels and linear combinations of levels
    levels <- hypothesis_BF(fit, c("f[B] = 0", "f[B] = 0.2", "f[C] = 0"),
                            density_method = "toyNormal", seed = 1)
    expect_setequal(
      paste(names(attr(levels, "warnings")), attr(levels, "warnings")),
      c("f (1) toy warning at 0", "f (2) toy warning at 0.2", "f (3) toy warning at 0")
    )
    expect_equal(
      attr(levels, "raw_log_BF"),
      attr(hypothesis_BF(fit, c("f[B] = 0", "f[B] = 0.2", "f[C] = 0"),
                         density_method = "normal", seed = 1), "raw_log_BF"),
      tolerance = 1e-8
    )
    combination <- hypothesis_BF(fit, "f[B] - 2 * f[C] = 0.1",
                                 density_method = "toyNormal", seed = 1)
    expect_equal(
      attr(combination, "raw_log_BF"),
      attr(hypothesis_BF(fit, "f[B] - 2 * f[C] = 0.1", density_method = "normal",
                         seed = 1), "raw_log_BF"),
      tolerance = 1e-8
    )
    last <- state$requests[[length(state$requests)]]
    expect_identical(last$route, "linear")
    expect_identical(last$kind, "scalar")
    expect_identical(last$spec$type, "linear")

    # one statement per group: the groups bind into one table in order
    mixed <- hypothesis_BF(fit, c("x = 0.1", "f[B] = 0.1", "x > 0"),
                           density_method = "toyNormal", seed = 1)
    expect_identical(rownames(mixed), c("(mu) x (1)", "f", "(mu) x (3)"))
    expect_identical(names(attr(mixed, "warnings")), c("(mu) x (1)", "f"))
    expect_identical(nrow(attr(mixed, "density_diagnostics")), 2L)
  })

})


test_that("a backend failure propagates and its settings are validated by it", {

  fit <- .hypothesis_model_test_fit()
  toy <- .hypothesis_model_toy_backend(control = function(control){
    if(is.null(control)) return(list(scale = 1))
    if(!is.list(control) || !identical(names(control), "scale")){
      stop("'density_control' needs 'scale'.", call. = FALSE)
    }
    control
  })
  with_toy_backend(toy, {
    failure <- tryCatch(hypothesis_BF(fit, "x = 99", density_method = "toyNormal"),
                        error = identity)
    expect_s3_class(failure, "toy_failure")
    expect_error(
      hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal",
                    density_control = list(other = 1)),
      "needs 'scale'"
    )
    expect_s3_class(
      hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal",
                    density_control = list(scale = 2)),
      "BayesTools_hypothesis_BF"
    )
  })

})


test_that("backend availability is evaluated on every call", {

  fit <- .hypothesis_model_test_fit()
  state <- new.env()
  state$refuse <- FALSE
  state$seen <- list()
  toy <- .hypothesis_model_toy_backend(
    state = state,
    available = function(object, plan, control){
      state$seen[[length(state$seen) + 1L]] <- list(object = object, plan = plan, control = control)
      if(state$refuse){
        list(reason = "toy backend refuses this fit.", class = "toy_refusal")
      }
    },
    control = function(control) if(is.null(control)) list(scale = 1) else control
  )
  with_toy_backend(toy, {
    ok <- hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal")
    expect_s3_class(ok, "BayesTools_hypothesis_BF")
    expect_s3_class(state$seen[[1L]]$object, "BayesTools_fit")
    expect_identical(state$seen[[1L]]$plan$kind, "scalar")
    expect_identical(state$seen[[1L]]$control, list(scale = 1))

    state$refuse <- TRUE
    refused <- tryCatch(hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal"),
                        error = identity)
    expect_s3_class(refused, "toy_refusal")
    expect_s3_class(refused, "BayesTools_hypothesis_method")
    expect_s3_class(refused, "BayesTools_hypothesis_unavailable")
    expect_identical(conditionMessage(refused), "toy backend refuses this fit.")
    # a region statement needs no density method
    expect_s3_class(hypothesis_BF(fit, "x > 0", density_method = "toyNormal"),
                    "BayesTools_hypothesis_BF")
    # the refusal is not cached with the plan
    state$refuse <- FALSE
    expect_s3_class(hypothesis_BF(fit, "x = 0.1", density_method = "toyNormal"),
                    "BayesTools_hypothesis_BF")
    # the quantities table reports the refusal per backend
    state$refuse <- TRUE
    quantities <- hypothesis_quantities(fit)
    expect_false(any(grepl("toyNormal", quantities$point_test_methods)))
    expect_true(any(grepl("toy backend refuses", quantities$reason)))
    state$refuse <- FALSE
    quantities <- hypothesis_quantities(fit)
    expect_true(all(grepl("toyNormal", quantities$point_test_methods[quantities$point_test])))
  })

})


test_that("stored ordinates are replaced unless the backend reuses them", {

  fit <- .hypothesis_model_test_fit()
  run <- function(reuse){
    state <- new.env()
    toy <- .hypothesis_model_toy_backend(state = state, reuse = reuse)
    provider <- hypothesis_provider(fit)
    samples  <- provider$mixed_posteriors("mu_x", NULL, "AND", FALSE, 1000)
    posterior <- marginal_posterior(samples, "mu_x", prior_samples = TRUE,
                                    use_formula = FALSE, n_samples = 1000)
    stale <- posterior_ordinate_attribute(
      value = 0.1, ordinate = 1000, method = "stale", density_method = "toyNormal",
      diagnostics = list(BF_error_percent = 1), parameter = "mu_x")
    other <- posterior_ordinate_attribute(
      value = 0.3, ordinate = 5, method = "other", density_method = "toyNormal",
      diagnostics = list(BF_error_percent = 1), parameter = "mu_x")
    posterior_metadata(posterior, "posterior_ordinate") <- posterior_ordinate_append(stale, other)
    request <- BayesTools:::.hypothesis_model_request(
      provider, toy, list(kind = "scalar"), kind = "scalar", parameter = "mu_x",
      label = "mu_x", values = 0.1)
    attached <- BayesTools:::.hypothesis_model_attach_request(
      BayesTools:::.hypothesis_model_keep_requested_ordinates(
        posterior, data.frame(level = NA_character_, value = 0.1)),
      toy, request, NULL, new.env())
    list(
      attached = state$attached,
      entries  = BayesTools:::.hypothesis_model_ordinate_entries(
        posterior_metadata(attached, "posterior_ordinate"))
    )
  }
  replaced <- run(NULL)
  expect_identical(replaced$attached, 1L)
  expect_length(replaced$entries, 1L)
  expect_identical(replaced$entries[[1L]]$method, "toy")
  reused <- run(function(entry, request, control) identical(entry$method, "stale"))
  expect_identical(reused$attached, 0L)
  expect_length(reused$entries, 1L)
  expect_identical(reused$entries[[1L]]$method, "stale")

})


test_that("a provider adds a quantity through the catalog extension", {

  fit <- .hypothesis_model_test_fit()
  model <- structure(list(fit = fit), class = "bt_hypothesis_test_model")
  calls <- new.env(); calls$sample <- character(); calls$mixed <- 0L
  provider_of <- function(object, ...){
    provider <- hypothesis_provider(object$fit)
    base_metadata <- provider$metadata
    quantity <- parameter_catalog(object$fit)$quantities
    quantity <- quantity[quantity$canonical_name == "mu_x", , drop = FALSE]
    quantity$quantity_id    <- "TestModel::my_x"
    quantity$canonical_name <- "my_x"
    quantity$provider       <- "TestModel"
    quantity$status         <- "derived"
    quantity$term           <- "my_x"
    quantity$component      <- "mu"
    quantity$display_label  <- "my_x"
    quantity$label_parts    <- I(list(NULL))
    quantity$extraction_key <- I(list(list(type = "test_alias", dependencies = "mu_x")))
    modifyList(provider, list(
      object = object,
      metadata = function(){
        metadata <- base_metadata()
        metadata$catalog <- parameter_catalog_extend(
          metadata$catalog, quantity,
          data.frame(alias = "my_x", quantity_id = "TestModel::my_x", namespace = "mu",
                     component = "mu", simplified = FALSE, stringsAsFactors = FALSE),
          provider = "TestModel")
        entry <- metadata$entries[metadata$entries$parameter == "mu_x", , drop = FALSE]
        entry$quantity_id <- "TestModel::my_x"
        entry$parameter   <- "my_x"
        entry$term        <- "my_x"
        entry$aliases     <- I(list("my_x"))
        metadata$entries  <- rbind(metadata$entries, entry)
        metadata
      },
      sample_parameters = function(parameter){
        calls$sample <- c(calls$sample, parameter)
        if(identical(parameter, "my_x")) "mu_x" else parameter
      },
      mixed_posteriors = function(parameters, conditional, conditional_rule,
                                  transform_scaled, n_prior_samples){
        calls$mixed <- calls$mixed + 1L
        provider$mixed_posteriors(parameters, conditional, conditional_rule,
                                  transform_scaled, n_prior_samples)
      }
    ))
  }
  registerS3method("hypothesis_provider", "bt_hypothesis_test_model", provider_of,
                   envir = asNamespace("BayesTools"))

  result <- hypothesis_BF_model(model, "my_x = 0.1", seed = 1)
  expect_identical(calls$sample, "my_x")
  expect_identical(calls$mixed, 1L)
  fitted <- hypothesis_BF(fit, "x = 0.1", seed = 1)
  expect_equal(attr(result, "raw_log_BF"), attr(fitted, "raw_log_BF"), tolerance = 1e-12)
  expect_identical(rownames(result), "my_x")
  expect_s3_class(hypothesis_BF_model(model, "my_x > 0 vs my_x < 0", seed = 1),
                  "BayesTools_hypothesis_BF")
  quantities <- hypothesis_quantities(model)
  expect_true("my_x" %in% quantities$alias)

  broken <- structure(list(fit = fit), class = "bt_hypothesis_test_broken")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_broken",
                   function(object, ...) list(kind = "model"),
                   envir = asNamespace("BayesTools"))
  expect_error(hypothesis_BF_model(broken, "x = 0"), "Invalid hypothesis provider")
  expect_error(hypothesis_BF_model(1, "x = 0"), "No hypothesis provider")

})


test_that("a marginal-means provider is tested on its stored marginal inference", {

  fit <- .hypothesis_model_test_fit()
  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 2000)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 2000)
  model <- structure(list(), class = "bt_hypothesis_test_marginal")
  state <- new.env(); state$requests <- 0L
  registerS3method("hypothesis_provider", "bt_hypothesis_test_marginal",
    function(object, ...){
      provider <- list(
        kind = "marginal_means", object = object, build = "BayesTools",
        refit_class = NULL, caller = "toy()", ordinate_label = "toy ordinates",
        density_methods = c("KDE", "toyNormal"),
        inference = list(averaged = list(f = levels), conditional = list(f = levels)),
        model_averaged = FALSE, capability_object = "source",
        select = function(hypothesis, parameter){
          list(parameter = "f", aliases = list(f = "f"))
        },
        density_method_default = function() "KDE",
        density_control_default = function(density_method) NULL,
        metadata = function(){
          list(term_map = data.frame(parameter = "f", term = "f", stringsAsFactors = FALSE))
        }
      )
      class(provider) <- c("BayesTools_hypothesis_provider", "list")
      provider
    }, envir = asNamespace("BayesTools"))

  region <- hypothesis_BF_model(model, "f[B] > 0", seed = 1)
  expect_s3_class(region, "BayesTools_hypothesis_BF")
  direct <- hypothesis_BF.default(levels, hypothesis = "f[B] > 0", parameter = "f", seed = 1)
  expect_identical(attr(region, "raw_log_BF"), attr(direct, "raw_log_BF"))
  expect_identical(rownames(hypothesis_BF_model(model, c("f[B] > 0", "f[C] < 1"), seed = 1)),
                   c("f (1)", "f (2)"))
  expect_error(hypothesis_BF_model(model, "f[A] = 0"), class = "BayesTools_hypothesis_fixed")
  expect_error(hypothesis_BF_model(model, "f[D] = 0"), class = "BayesTools_parameter_not_found")

  toy <- .hypothesis_model_toy_backend(state = state, reuse = function(entry, request, control) TRUE)
  with_toy_backend(toy, {
    point <- hypothesis_BF_model(model, "f[B] = 0.2", density_method = "toy")
    expect_identical(state$attached, 1L)
    expect_identical(state$requests[[1L]]$kind, "level")
    expect_identical(state$requests[[1L]]$route, "marginal_mean")
    expect_identical(state$requests[[1L]]$object, "source")
    expect_identical(state$requests[[1L]]$level, "B")
    expect_equal(
      attr(point, "raw_log_BF"),
      attr(hypothesis_BF.default(levels, hypothesis = "f[B] = 0.2", parameter = "f",
                                 density_method = "normal"), "raw_log_BF"),
      tolerance = 1e-8
    )
  })
  expect_identical(hypothesis_quantities(model)$bracket, c("f[A]", "f[B]", "f[C]"))

})


test_that("the exported helpers of the engine describe formula coefficients", {

  fit <- .hypothesis_model_test_fit()
  route <- JAGS_formula_coefficient_route(fit, "mu", "mu_x")
  transform <- JAGS_formula_coefficient_transform(fit, "mu")
  expect_identical(route$type, "affine")
  expect_equal(unname(route$weights[["mu_x"]]), unname(transform$matrix["mu_x", "mu_x"]))
  expect_identical(route$formula_parameter, "mu")
  expect_error(JAGS_formula_coefficient_route(fit, "mu", "mu_absent"),
               class = "BayesTools_refit_required")

  sample <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 500)
  levels <- marginal_posterior(sample, "mu_f", prior_samples = TRUE, use_formula = FALSE,
                               n_samples = 500)
  expect_identical(posterior_linear_weight_space(levels$B), "coefficient")
  expect_true(posterior_formula_measure_check(levels$B, "atoms"))
  weights <- posterior_metadata(levels$B, "linear_weights")
  expect_equal(posterior_fitted_linear_weights(levels$B, weights), weights)
  expect_null(posterior_linear_weight_space(as_mixed_posteriors(fit, "mu_x")$mu_x))

  expect_null(hypothesis_cache_clear())
  expect_s3_class(hypothesis_BF(fit, "x > 0", seed = 1), "BayesTools_hypothesis_BF")

})


test_that("the eligibility table lists what the engine tests", {

  fit <- .hypothesis_model_test_fit()
  quantities <- hypothesis_quantities(fit)
  expect_identical(names(quantities), c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ))
  factor_rows <- quantities[quantities$parameter == "mu_f", ]
  expect_true(all(factor_rows$bracket == "mu_f[level]"))
  expect_true(all(factor_rows$contrast_test))
  expect_true(all(quantities$point_test_methods == "KDE"))
  expect_identical(sort(unique(quantities$component)), c("model", "mu"))
  expect_true(all(c("x", "(mu) x", "mu_x") %in% quantities$alias))
  expect_warning(hypothesis_quantities(fit, nonsense = 1), "Unused argument")
  scaled <- hypothesis_quantities(fit, transform_scaled = TRUE)
  expect_identical(scaled$alias, quantities$alias)

})
