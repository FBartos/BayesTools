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

.hypothesis_model_test_fit <- function(seed = 3L, n = 400L, formula_scale = TRUE){

  set.seed(seed)
  data <- data.frame(x = c(-2, -1, 0, 3, 5, 7), f = factor(rep(c("A", "B", "C"), 2)))
  priors <- list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1)),
    f         = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  compiled <- JAGS_formula(~ 1 + x + f, "mu", data, priors, formula_scale = formula_scale)
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
  state$controls <- list()
  hypothesis_density_backend(
    name = name, type = "ordinate",
    attach = function(posterior, request, control, cache){
      state$attached <- state$attached + 1L
      state$requests[[length(state$requests) + 1L]] <- request
      state$controls[[length(state$controls) + 1L]] <- control
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
  hypothesis_cache_clear()
  set.seed(11)
  engine <- hypothesis_BF(fit, "x = 0.1", seed = 1)
  engine_seed <- .Random.seed
  set.seed(11)
  manual <- .hypothesis_model_manual_coefficient(fit, "mu_x", "mu_x = 0.1", seed = 1)
  manual_seed <- .Random.seed
  expect_s3_class(engine, "BayesTools_hypothesis_BF")
  expect_identical(attr(engine, "raw_log_BF"), attr(manual, "raw_log_BF"))
  expect_identical(attr(engine, "raw_BF"), attr(manual, "raw_BF"))
  # the one-call form leaves the random-number generator where the three steps do
  expect_identical(engine_seed, manual_seed)
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
    registered <- hypothesis_density_backends()
    expect_identical(registered$policy[[match("toyNormal", registered$name)]]$error,
                     c(0.025, 0.05))
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
    kept <- BayesTools:::.hypothesis_model_keep_requested_ordinates(
      posterior, data.frame(level = NA_character_, value = 0.1))
    request <- BayesTools:::.hypothesis_model_request(
      provider, toy, list(kind = "scalar"), kind = "scalar", parameter = "mu_x",
      label = "mu_x", values = 0.1, sample = kept)
    attached <- BayesTools:::.hypothesis_model_attach_request(
      kept, toy, request, NULL, new.env())
    list(
      attached = state$attached,
      entries  = BayesTools:::.hypothesis_model_ordinate_entries(
        posterior_metadata(attached, "posterior_ordinate")),
      kept     = kept,
      request  = request
    )
  }
  replaced <- run(NULL)
  expect_identical(replaced$attached, 1L)
  expect_length(replaced$entries, 1L)
  expect_identical(replaced$entries[[1L]]$method, "toy")

  # reuse() decides for the stored ordinate at the requested value with the
  # draws that carry it and the request, whose condition is theirs
  seen <- NULL
  reused <- run(function(entry, request, control, posterior){
    seen <<- list(entry = entry, request = request, posterior = posterior)
    identical(entry$method, "stale")
  })
  expect_identical(reused$attached, 0L)
  expect_length(reused$entries, 1L)
  expect_identical(reused$entries[[1L]]$method, "stale")
  expect_identical(seen$entry$method, "stale")
  expect_identical(seen$posterior, reused$kept)
  expect_identical(
    BayesTools:::.hypothesis_model_ordinate_entries(
      posterior_metadata(seen$posterior, "posterior_ordinate"))[[1L]]$method,
    "stale"
  )
  expect_identical(seen$request$condition,
                   BayesTools:::.hypothesis_model_condition_metadata(seen$posterior))
  expect_identical(seen$request$condition$condition_key, "<averaged>")

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

  toy <- .hypothesis_model_toy_backend(state = state,
                                       reuse = function(entry, request, control, posterior) TRUE)
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


# A fit whose two formula parameters share the predictors and the treatment
# factor 'f': the references to 'f' are ambiguous without a 'namespace'.
.hypothesis_model_shared_factor_fit <- function(seed = 5L, n = 400L){

  set.seed(seed)
  data <- data.frame(x = c(-2, -1, 0, 3, 5, 7), f = factor(rep(c("A", "B", "C"), 2)))
  compile <- function(parameter){
    JAGS_formula(~ 1 + x + f, parameter, data, list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1)),
      f         = prior_factor("normal", list(0, 1), contrast = "treatment")
    ), formula_scale = TRUE)
  }
  mu        <- compile("mu")
  log_sigma <- compile("log_sigma")
  columns   <- c(paste0("mu_", c("intercept", "x", "f[1]", "f[2]")),
                 paste0("log_sigma_", c("intercept", "x", "f[1]", "f[2]")))
  draws <- matrix(stats::rnorm(n * length(columns), 0.3, 0.5), n, length(columns),
                  dimnames = list(NULL, columns))

  .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(draws)), c(mu$prior_list, log_sigma$prior_list),
    list(mu = mu$formula_design, log_sigma = log_sigma$formula_design),
    list(mu = mu$formula_scale, log_sigma = log_sigma$formula_scale)
  )
}

# The draws 'sample' of a level of 'mu_f' as the contribution of a formula term
# multiplied by 'multiplier' (a number, or the name of a parameter); 'scaled'
# declares a scaling of the coefficients that the weights cannot be remapped
# through.
.hypothesis_model_formula_contribution <- function(sample, multiplier, scaled = FALSE){

  posterior_metadata(sample, "linear_weight_space") <- "formula_contribution"
  context <- posterior_metadata(sample, "prior_context")
  attr(context$prior_list$mu_f, "multiply_by") <- multiplier
  if(scaled){
    context$formula_scale <- list(mu = list(x = 1))
  }
  posterior_metadata(sample, "prior_context") <- context

  sample
}

# A model with a marginal-means provider that stores the marginal posterior of
# the levels of 'f' ('levels_hook' changes them); '...' replace elements of the
# provider.
.hypothesis_model_toy_marginal <- function(fit, class, levels_hook = NULL, ...){

  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 2000)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 2000)
  if(is.function(levels_hook)){
    levels <- levels_hook(levels)
  }
  overrides <- list(...)
  registerS3method("hypothesis_provider", class, function(object, ...){
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
    provider <- utils::modifyList(provider, overrides)
    class(provider) <- c("BayesTools_hypothesis_provider", "list")
    provider
  }, envir = asNamespace("BayesTools"))

  structure(list(), class = class)
}


test_that("a level token selects the level of a factor term", {

  fit <- .hypothesis_model_test_fit()
  # the references to the term stand for its level, as if it were written in
  # the statement
  cases <- list(
    c("f = 0.1", "f[B] = 0.1"),
    c("f = 0.1 vs f > 0.1", "f[B] = 0.1 vs f[B] > 0.1"),
    c("f > 0.1", "f[B] > 0.1"),
    c("f > 0.1 & f < 0.4", "f[B] > 0.1 & f[B] < 0.4")
  )
  for(case in cases){
    selected <- hypothesis_BF(fit, case[[1L]], component = "B", seed = 1)
    explicit <- hypothesis_BF(fit, case[[2L]], seed = 1)
    expect_identical(nrow(selected), 1L, info = case[[1L]])
    expect_identical(attr(selected, "raw_log_BF"), attr(explicit, "raw_log_BF"),
                     info = case[[1L]])
    expect_identical(rownames(selected), rownames(explicit), info = case[[1L]])
  }
  # without the token the whole term is tested, including its fixed level
  expect_error(hypothesis_BF(fit, "f > 0.1", seed = 1), class = "BayesTools_hypothesis_fixed")
  # several statements on the level
  several <- hypothesis_BF(fit, c("f = 0.1", "f > 0"), component = "B", seed = 1)
  expect_identical(
    attr(several, "raw_log_BF"),
    attr(hypothesis_BF(fit, c("f[B] = 0.1", "f[B] > 0"), seed = 1), "raw_log_BF")
  )
  expect_identical(nrow(several), 2L)

  # the level is a fixed reference level, or no level of the term
  expect_error(hypothesis_BF(fit, "f = 0.1", component = "A"),
               class = "BayesTools_hypothesis_fixed")
  expect_error(hypothesis_BF(fit, "f = 0.1", component = "Z"),
               class = "BayesTools_parameter_not_found")
  # a level in the statement must be the token's
  expect_error(hypothesis_BF(fit, "f[C] = 0.1", component = "B"),
               class = "BayesTools_hypothesis_component_mismatch")

  # the level reaches the density backends
  state <- new.env()
  toy <- .hypothesis_model_toy_backend(state = state)
  with_toy_backend(toy, {
    point <- hypothesis_BF(fit, "f = 0.1", component = "B",
                           density_method = "toyNormal", seed = 1)
    expect_identical(state$requests[[1L]]$kind, "level")
    expect_identical(state$requests[[1L]]$level, "B")
    expect_identical(state$attached, 1L)
    expect_equal(
      attr(point, "raw_log_BF"),
      attr(hypothesis_BF(fit, "f[B] = 0.1", density_method = "toyNormal", seed = 1),
           "raw_log_BF"),
      tolerance = 1e-12
    )
  })

})


test_that("the namespace and the level token combine", {

  fit <- .hypothesis_model_shared_factor_fit()
  ambiguous <- tryCatch(hypothesis_BF(fit, "f = 0.1", component = "B"), error = identity)
  expect_s3_class(ambiguous, "BayesTools_hypothesis_ambiguous")
  expect_match(conditionMessage(ambiguous), "Set 'namespace' to 'mu', 'log_sigma'", fixed = TRUE)

  # 'namespace' selects the formula parameter that shares the factor, 'component'
  # its level
  for(namespace in c("mu", "log_sigma")){
    selected <- hypothesis_BF(fit, "f = 0.1", namespace = namespace, component = "B", seed = 1)
    explicit <- hypothesis_BF(fit, paste0(namespace, "_f[B] = 0.1"), seed = 1)
    expect_identical(nrow(selected), 1L, info = namespace)
    expect_identical(attr(selected, "raw_log_BF"), attr(explicit, "raw_log_BF"),
                     info = namespace)
    expect_identical(rownames(selected), rownames(explicit), info = namespace)

    # also for a level written in the statement, and for the whole term
    leveled <- hypothesis_BF(fit, "f[B] = 0.1", namespace = namespace, component = "B", seed = 1)
    expect_identical(attr(leveled, "raw_log_BF"), attr(explicit, "raw_log_BF"),
                     info = namespace)
    region <- hypothesis_BF(fit, "f > 0.1", namespace = namespace, component = "C", seed = 1)
    expect_identical(
      attr(region, "raw_log_BF"),
      attr(hypothesis_BF(fit, paste0(namespace, "_f[C] > 0.1"), seed = 1), "raw_log_BF"),
      info = namespace
    )
    expect_error(hypothesis_BF(fit, "f = 0.1", namespace = namespace, component = "A"),
                 class = "BayesTools_hypothesis_fixed", info = namespace)
  }
  expect_false(isTRUE(all.equal(
    attr(hypothesis_BF(fit, "f = 0.1", namespace = "mu", component = "B", seed = 1), "raw_log_BF"),
    attr(hypothesis_BF(fit, "f = 0.1", namespace = "log_sigma", component = "B", seed = 1), "raw_log_BF")
  )))

  # a statement on another formula parameter than the namespace's, and a
  # reference that has no such level, are refused
  expect_error(hypothesis_BF(fit, "log_sigma_f[B] = 0.1", namespace = "mu", component = "B"),
               class = "BayesTools_hypothesis_component_mismatch")
  expect_error(hypothesis_BF(fit, "x = 0.1", namespace = "mu", component = "B"),
               class = "BayesTools_parameter_not_found")

})


test_that("the eligibility table has its columns without tested quantities", {

  fit <- .hypothesis_model_test_fit()
  columns <- names(hypothesis_quantities(fit))
  types   <- vapply(hypothesis_quantities(fit), function(column) class(column)[1L], character(1))
  expect_empty <- function(table, info){
    expect_s3_class(table, "data.frame")
    expect_identical(nrow(table), 0L, info = info)
    expect_identical(names(table), columns, info = info)
    expect_identical(vapply(table, function(column) class(column)[1L], character(1)),
                     types, info = info)
  }

  # a model whose provider lists no entries, and one whose entries are all refused
  model <- structure(list(fit = fit), class = "bt_hypothesis_test_empty")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_empty", function(object, ...){
    provider <- hypothesis_provider(object$fit)
    base     <- provider$metadata
    utils::modifyList(provider, list(
      object   = object,
      metadata = function(){
        metadata <- base()
        metadata$entries <- metadata$entries[0L, , drop = FALSE]
        metadata
      }
    ))
  }, envir = asNamespace("BayesTools"))
  expect_empty(hypothesis_quantities(model), "no entries")
  expect_empty(hypothesis_quantities(model, transform_scaled = TRUE), "no entries, scaled")

  refused <- structure(list(fit = fit), class = "bt_hypothesis_test_refused")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_refused", function(object, ...){
    provider <- hypothesis_provider(object$fit)
    provider$object <- object
    provider$components$refused <- c(mu = "Refused.", model = "Refused.")
    provider
  }, envir = asNamespace("BayesTools"))
  expect_empty(hypothesis_quantities(refused), "all refused")

  # stored marginal means without terms
  empty <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_empty",
    metadata = function(){
      list(term_map = data.frame(parameter = character(), term = character(),
                                 stringsAsFactors = FALSE))
    }
  )
  expect_empty(hypothesis_quantities(empty), "no terms")

})


test_that("marginal-means availability sees the control the backend normalized", {

  fit   <- .hypothesis_model_test_fit()
  state <- new.env()
  state$seen <- list()
  toy <- .hypothesis_model_toy_backend(
    state = state,
    control = function(control){
      if(isTRUE(control$bad)){
        stop(structure(class = c("toy_bad_control", "error", "condition"),
                       list(message = "bad toy control", call = NULL)))
      }
      utils::modifyList(list(ready = TRUE), if(is.null(control)) list() else control)
    },
    available = function(object, plan, control){
      state$seen[[length(state$seen) + 1L]] <- control
      if(isTRUE(control$refuse)){
        return(list(reason = "toy capability refusal", class = "toy_capability"))
      }
      if(!isTRUE(control$ready) && !isTRUE(control$bad)){
        list(reason = "the control was not normalized", class = "toy_unnormalized")
      }
    }
  )
  reset <- function(){
    state$seen <- list(); state$controls <- list(); state$attached <- 0L
  }
  # stored settings that the backend normalizes: an empty list
  model <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_control",
    density_methods = c("KDE", "toyNormal"),
    density_control_default = function(density_method) list()
  )
  with_toy_backend(toy, {

    # the inherited stored settings, an empty control, and settings of the call
    # reach available() and attach() normalized
    for(control in list(NULL, list(), list(extra = 1))){
      reset()
      result <- hypothesis_BF_model(model, "f[B] = 0.2", density_method = "toyNormal",
                                    density_control = control)
      expected <- utils::modifyList(list(ready = TRUE), if(is.null(control)) list() else control)
      expect_s3_class(result, "BayesTools_hypothesis_BF")
      expect_identical(state$seen[[1L]], expected, info = deparse(control))
      expect_identical(state$controls[[1L]], expected, info = deparse(control))
    }

    # a refusal of the backend precedes the failure of the settings, which stops
    # before anything is attached; available() then sees the settings as given
    reset()
    refused <- tryCatch(
      hypothesis_BF_model(model, "f[B] = 0.2", density_method = "toyNormal",
                          density_control = list(bad = TRUE, refuse = TRUE)),
      error = identity
    )
    expect_s3_class(refused, "toy_capability")
    expect_identical(state$seen[[1L]], list(bad = TRUE, refuse = TRUE))
    reset()
    bad <- tryCatch(
      hypothesis_BF_model(model, "f[B] = 0.2", density_method = "toyNormal",
                          density_control = list(bad = TRUE)),
      error = identity
    )
    expect_s3_class(bad, "toy_bad_control")
    expect_identical(state$attached, 0L)
    expect_identical(state$seen[[1L]], list(bad = TRUE))
    # a statement the plans refuse is reported before the settings
    expect_error(
      hypothesis_BF_model(model, "f[A] = 0.2", density_method = "toyNormal",
                          density_control = list(bad = TRUE)),
      class = "BayesTools_hypothesis_fixed"
    )
    # a region statement needs no density method: the settings of the call are
    # still validated, the stored ones are not normalized
    expect_error(
      hypothesis_BF_model(model, "f[B] > 0.2", density_method = "toyNormal",
                          density_control = list(bad = TRUE)),
      class = "toy_bad_control"
    )
    expect_s3_class(
      hypothesis_BF_model(model, "f[B] > 0.2", density_method = "toyNormal"),
      "BayesTools_hypothesis_BF"
    )
  })

  # stored settings that fail to normalize stop the point statements only
  failing <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_failing",
    density_methods = c("KDE", "toyNormal"),
    density_control_default = function(density_method) list(bad = TRUE)
  )
  with_toy_backend(toy, {
    expect_s3_class(
      hypothesis_BF_model(failing, "f[B] > 0.2", density_method = "toyNormal"),
      "BayesTools_hypothesis_BF"
    )
    expect_error(
      hypothesis_BF_model(failing, "f[B] = 0.2", density_method = "toyNormal"),
      class = "toy_bad_control"
    )
  })

})


test_that("the provider checks the settings of a density method first", {

  fit   <- .hypothesis_model_test_fit()
  calls <- new.env()
  calls$checks <- list()
  calls$order  <- character()
  check <- function(density_method, density_control){
    calls$checks[[length(calls$checks) + 1L]] <- list(density_method, density_control)
    calls$order <- c(calls$order, "provider")
    if(density_method %in% c("KDE", "normal") && !is.null(density_control)){
      stop(structure(class = c("toy_draws_control", "error", "condition"),
                     list(message = "the draws methods take no settings", call = NULL)))
    }
  }
  model <- structure(list(fit = fit), class = "bt_hypothesis_test_check")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_check", function(object, ...){
    provider <- hypothesis_provider(object$fit)
    provider$object <- object
    provider$check_density_control <- check
    provider
  }, envir = asNamespace("BayesTools"))
  toy <- .hypothesis_model_toy_backend(control = function(control){
    calls$order <- c(calls$order, "backend")
    if(is.null(control)) list(scale = 1) else control
  })

  # without the hook the built-in backends accept an empty list
  expect_s3_class(hypothesis_BF(fit, "x = 0.1", density_control = list(), seed = 1),
                  "BayesTools_hypothesis_BF")
  expect_error(hypothesis_BF(fit, "x = 0.1", density_control = list(a = 1)),
               "not used by density_method = 'KDE'")

  # with it, the provider's own conditions: for any settings, an empty list too
  for(control in list(list(), list(a = 1))){
    error <- tryCatch(
      hypothesis_BF_model(model, "x = 0.1", density_method = "KDE", density_control = control),
      error = identity
    )
    expect_s3_class(error, "toy_draws_control")
  }
  expect_identical(calls$checks[[1L]], list("KDE", list()))
  # no settings, no check
  calls$checks <- list()
  expect_s3_class(hypothesis_BF_model(model, "x = 0.1", seed = 1), "BayesTools_hypothesis_BF")
  expect_length(calls$checks, 0L)

  # an ordinate backend: the provider's check precedes its control() on every call
  with_toy_backend(toy, {
    calls$order <- character()
    hypothesis_BF_model(model, "x = 0.1", density_method = "toyNormal", seed = 1)
    expect_identical(calls$order, c("provider", "backend"))
    expect_identical(calls$checks[[1L]], list("toyNormal", NULL))
  })

  # stored marginal means
  calls$checks <- list()
  marginal <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_check",
    check_density_control = check
  )
  expect_error(
    hypothesis_BF_model(marginal, "f[B] > 0", density_method = "KDE", density_control = list()),
    class = "toy_draws_control"
  )
  expect_length(calls$checks, 1L)
  # after the checks of the plans
  expect_error(
    hypothesis_BF_model(marginal, "f[A] = 0", density_method = "KDE", density_control = list()),
    class = "BayesTools_hypothesis_fixed"
  )
  calls$checks <- list()
  expect_s3_class(hypothesis_BF_model(marginal, "f[B] > 0", seed = 1), "BayesTools_hypothesis_BF")
  expect_length(calls$checks, 0L)

  expect_error(
    {
      broken <- structure(list(fit = fit), class = "bt_hypothesis_test_check_broken")
      registerS3method("hypothesis_provider", "bt_hypothesis_test_check_broken",
                       function(object, ...){
        provider <- hypothesis_provider(object$fit)
        provider$object <- object
        provider$check_density_control <- "not a function"
        provider
      }, envir = asNamespace("BayesTools"))
      hypothesis_BF_model(broken, "x = 0.1")
    },
    "'check_density_control' must be a function or NULL"
  )

})


test_that("refusals of the density target of a random-effect quantity keep their type and cause", {

  refusal_of <- function(target, parameter = "sd1", evaluation = list(parameter = NULL)){
    provider <- list(
      random = list(density_target = function(parameter, transform_scaled, operation){
        if(is.function(target)) target(parameter) else target
      })
    )
    plan <- list(conditional = FALSE, parameter = parameter, transform_scaled = FALSE)
    BayesTools:::.hypothesis_model_plan_random_precomputed(
      plan, provider, evaluation, label = parameter
    )$refusals$ordinate
  }
  classes_of <- function(refusal){
    backend <- hypothesis_density_backend(
      "toyRefusal", attach = function(...) NULL,
      refusal_classes = function(cause) c(paste0("toy_density_", cause), "toy_density_unavailable")
    )
    plan <- list(refusal = NULL, point = TRUE, targets = list(),
                 method_refusals = list(ordinate = refusal))
    BayesTools:::.hypothesis_model_plan_status(plan, list(), backend)$class
  }

  # the provider's type and cause
  original <- refusal_of(list(reason = "allocated SD", type = "target", cause = "original_scale"))
  expect_identical(original$reason, "allocated SD")
  expect_identical(original$cause, "original_scale")
  expect_identical(classes_of(original), c(
    "BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable",
    "toy_density_original_scale", "toy_density_unavailable"
  ))
  # the defaults are a target without a supported coordinate
  default <- refusal_of(list(reason = "no coordinate"))
  expect_identical(default$cause, "random_target")
  expect_identical(classes_of(default), c(
    "BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable",
    "toy_density_random_target", "toy_density_unavailable"
  ))
  as_method <- refusal_of(list(reason = "r", type = "method", cause = "original_scale"))
  expect_identical(classes_of(as_method)[1:2], c(
    "BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable"
  ))
  expect_error(refusal_of(list(reason = "r", cause = "other")), "Invalid hypothesis provider")
  expect_error(refusal_of(list(reason = "r", type = "statement")), "Invalid hypothesis provider")

  # a variance is evaluated through its SD: its own reason, the SD's type and cause
  variance <- refusal_of(
    function(parameter){
      if(identical(parameter, "sd1")){
        list(reason = "SD refused", type = "target", cause = "original_scale")
      }else{
        list(reason = "variance refused", cause = "random_target")
      }
    },
    parameter = "var1", evaluation = list(parameter = "sd1")
  )
  expect_identical(variance$reason, "variance refused")
  expect_identical(variance$cause, "original_scale")

  # a conditional statement keeps its own cause and type
  conditional <- BayesTools:::.hypothesis_model_plan_random_precomputed(
    list(conditional = TRUE, parameter = "sd1", transform_scaled = FALSE),
    list(random = list()), list(parameter = NULL), label = "sd1"
  )$refusals$ordinate
  expect_identical(conditional$cause, "conditional_random")
  expect_identical(classes_of(conditional)[1:2], c(
    "BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable"
  ))

})


test_that("every callback of a backend receives the request of its group", {

  fit   <- .hypothesis_model_test_fit()
  state <- new.env()
  state$warnings <- list()
  state$diagnostics <- list()
  state$reused <- list()
  toy <- .hypothesis_model_toy_backend(
    state = state,
    warnings = function(posterior, request){
      state$warnings[[length(state$warnings) + 1L]] <- request
      BayesTools:::.hypothesis_model_ordinate_warning_records(posterior)
    },
    diagnostics = function(posterior, request){
      state$diagnostics[[length(state$diagnostics) + 1L]] <- request
      NULL
    }
  )
  reset <- function(){
    state$requests <- list(); state$warnings <- list(); state$diagnostics <- list()
  }
  # the request of the first attach call of the group, as it was built
  expect_first <- function(info){
    expect_length(state$warnings, 1L)
    expect_length(state$diagnostics, 1L)
    expect_false(is.null(state$warnings[[1L]]), info = info)
    expect_identical(state$warnings[[1L]], state$requests[[1L]], info = info)
    expect_identical(state$diagnostics[[1L]], state$requests[[1L]], info = info)
    expect_identical(state$warnings[[1L]]$density_method, "toyNormal", info = info)
    expect_true(is.list(state$warnings[[1L]]$condition), info = info)
    expect_identical(state$warnings[[1L]]$condition$condition_key, "<averaged>", info = info)
  }

  with_toy_backend(toy, {
    for(case in list(
      list("x = 0.1", "scalar", "scalar", NA_character_),
      list("f[B] = 0.1", "level", "linear", "B"),
      list("f[B] - 2 * f[C] = 0.1", "scalar", "linear", NA_character_)
    )){
      reset()
      hypothesis_BF(fit, case[[1L]], density_method = "toyNormal", seed = 1)
      expect_first(case[[1L]])
      expect_identical(state$warnings[[1L]]$kind, case[[2L]], info = case[[1L]])
      expect_identical(state$warnings[[1L]]$route, case[[3L]], info = case[[1L]])
      expect_identical(state$warnings[[1L]]$level, case[[4L]], info = case[[1L]])
    }

    # a group attached level by level: the first request, which differs from the
    # others in its level and values only
    reset()
    hypothesis_BF(fit, c("f[B] = 0.1", "f[C] = 0"), density_method = "toyNormal", seed = 1)
    expect_length(state$requests, 2L)
    expect_first("levels")
    expect_identical(state$warnings[[1L]]$level, "B")
    expect_identical(state$requests[[2L]]$level, "C")

    # every group of a call
    reset()
    hypothesis_BF(fit, c("x = 0.1", "f[B] = 0.1"), density_method = "toyNormal", seed = 1)
    expect_length(state$warnings, 2L)
    expect_length(state$diagnostics, 2L)
    expect_identical(vapply(state$warnings, `[[`, character(1), "kind"), c("scalar", "level"))
    expect_false(any(vapply(c(state$warnings, state$diagnostics), is.null, logical(1))))
  })

  # stored marginal means
  marginal <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_request",
    density_control_default = function(density_method) list()
  )
  with_toy_backend(toy, {
    for(statement in c("f[B] = 0.2", "f[B] - f[C] = 0.1")){
      reset()
      hypothesis_BF_model(marginal, statement, density_method = "toyNormal", seed = 1)
      expect_first(statement)
    }
    expect_identical(state$warnings[[1L]]$object, "source")
    expect_identical(state$warnings[[1L]]$route, "marginal_mean")
    reset()
    hypothesis_BF_model(marginal, c("f[B] = 0.2", "f[C] = 0.1"),
                        density_method = "toyNormal", seed = 1)
    expect_length(state$requests, 2L)
    expect_first("marginal levels")
    expect_identical(state$warnings[[1L]]$kind, "level")
    expect_identical(state$warnings[[1L]]$level, "B")
    expect_identical(state$requests[[2L]]$level, "C")
  })

})


test_that("the fitted linear weights of stochastic products raise their own condition", {

  fit     <- .hypothesis_model_test_fit(formula_scale = FALSE)
  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 500)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 500)
  weights <- posterior_metadata(levels$B, "linear_weights")
  multiplied <- function(multiplier){
    .hypothesis_model_formula_contribution(levels$B, multiplier)
  }

  # a constant multiplier is folded in
  expect_equal(posterior_fitted_linear_weights(multiplied(3), weights), c("mu_f[1]" = 3))

  # a stochastic one has no fitted linear projection: its own condition, with
  # the cause for which the density backends refuse the target, and no class
  # of a refusal of a hypothesis
  condition <- tryCatch(posterior_fitted_linear_weights(multiplied("sigma"), weights),
                        error = identity)
  expect_identical(class(condition), c(
    "BayesTools_fitted_linear_weights_unavailable",
    "BayesTools_formula_transform_unavailable", "error", "condition"
  ))
  expect_identical(
    conditionMessage(condition),
    paste0("Ordinate density methods are unavailable for stochastic formula ",
           "contribution products. Use density_method = 'KDE'.")
  )
  expect_identical(condition$reason, "stochastic_multiplier")
  expect_identical(condition$cause, "original_scale")
  expect_identical(names(unclass(condition)), c("message", "call", "reason", "cause"))
  expect_false(any(grepl("^BayesTools_hypothesis_", class(condition))))
  labelled <- tryCatch(
    posterior_fitted_linear_weights(multiplied("sigma"), weights, ordinate_label = "Toy ordinates"),
    error = identity
  )
  expect_match(conditionMessage(labelled), "^Toy ordinates are unavailable")

})


test_that("the weight helpers name the build and the refit class of the caller", {

  fit     <- .hypothesis_model_test_fit(formula_scale = FALSE)
  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 500)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 500)
  weights <- posterior_metadata(levels$B, "linear_weights")
  # contribution weights of a scaled formula cannot be remapped: refit
  stale <- .hypothesis_model_formula_contribution(levels$B, 3, scaled = TRUE)

  default <- tryCatch(posterior_linear_weight_space(stale), error = identity)
  expect_s3_class(default, "BayesTools_refit_required")
  expect_match(conditionMessage(default), "with the current BayesTools build.", fixed = TRUE)
  expect_false(inherits(default, "toy_refit_required"))

  for(helper in list(
    function(...) posterior_linear_weight_space(stale, ...),
    function(...) posterior_fitted_linear_weights(stale, weights, ...)
  )){
    named <- tryCatch(helper(build = "Toy/BayesTools", refit_class = "toy_refit_required"),
                      error = identity)
    expect_identical(class(named)[1:2], c("toy_refit_required", "BayesTools_refit_required"))
    expect_match(conditionMessage(named), "with the current Toy/BayesTools build.", fixed = TRUE)
    expect_s3_class(tryCatch(helper(), error = identity), "BayesTools_refit_required")
    expect_error(helper(build = 1), "build")
    expect_error(helper(refit_class = 1), "refit_class")
  }
  # without a refit request the arguments change nothing
  expect_identical(posterior_linear_weight_space(levels$B, build = "Toy", refit_class = "toy"),
                   "coefficient")

})


test_that("marginal means with stochastic products refuse the ordinate backends by their cause", {

  fit <- .hypothesis_model_test_fit(formula_scale = FALSE)
  toy <- .hypothesis_model_toy_backend(
    refusal_classes = function(cause) c(paste0("toy_density_", cause), "toy_density_unavailable")
  )
  model <- .hypothesis_model_toy_marginal(
    fit, "bt_hypothesis_test_marginal_stochastic",
    levels_hook = function(levels){
      for(level in c("B", "C")){
        levels[[level]] <- .hypothesis_model_formula_contribution(levels[[level]], "sigma")
      }
      levels
    },
    density_control_default = function(density_method) list()
  )
  expected <- c(
    "BayesTools_hypothesis_method", "BayesTools_hypothesis_unavailable",
    "toy_density_original_scale", "toy_density_unavailable", "error", "condition"
  )
  with_toy_backend(toy, {
    refused <- tryCatch(
      hypothesis_BF_model(model, "f[B] = 0.2", density_method = "toyNormal"),
      error = identity
    )
    # the classes of the refusal and of the backend for the cause, and no fields
    expect_identical(class(refused), expected)
    expect_identical(names(unclass(refused)), c("message", "call"))
    expect_identical(
      conditionMessage(refused),
      paste0("toy ordinates are unavailable for stochastic formula contribution ",
             "products. Use density_method = 'KDE'.")
    )

    # a linear combination of the levels, whose weights refuse the same way
    local_mocked_bindings(
      .hypothesis_model_fitted_linear_weights = function(provider, sample, weights){
        BayesTools:::.hypothesis_model_stop_stochastic_weights(
          "toy ordinates are unavailable for stochastic formula contribution products."
        )
      },
      .package = "BayesTools"
    )
    combination <- tryCatch(
      hypothesis_BF_model(model, "f[B] - f[C] = 0.1", density_method = "toyNormal"),
      error = identity
    )
    expect_identical(class(combination), expected)
    expect_identical(names(unclass(combination)), c("message", "call"))
    # the draws and the region statements do not need the fitted weights
    expect_s3_class(hypothesis_BF_model(model, "f[B] = 0.2", density_method = "KDE"),
                    "BayesTools_hypothesis_BF")
    expect_s3_class(hypothesis_BF_model(model, "f[B] > 0.2", density_method = "toyNormal"),
                    "BayesTools_hypothesis_BF")
  })

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
  # the methods of the backends of BayesTools (backends of other packages that
  # are registered in the session are listed next to them)
  registered <- hypothesis_density_backends()
  own        <- registered$name[!is.na(registered$package) & registered$package == "BayesTools"]
  listed     <- lapply(strsplit(quantities$point_test_methods, ", ", fixed = TRUE), intersect, own)
  expect_true(all(vapply(listed, identical, logical(1), y = "KDE")))
  expect_identical(sort(unique(quantities$component)), c("model", "mu"))
  expect_true(all(c("x", "(mu) x", "mu_x") %in% quantities$alias))
  expect_warning(hypothesis_quantities(fit, nonsense = 1), "Unused argument")
  scaled <- hypothesis_quantities(fit, transform_scaled = TRUE)
  expect_identical(scaled$alias, quantities$alias)

})


# ---------------------------------------------------------------------------- #
# The session cache of plans
# ---------------------------------------------------------------------------- #

# The cache of an object is found by the hash of the object's whole content. The
# objects here are plain lists: the cache does not read what it is a cache of.
# The provider of an object names the object the cache belongs to.
.hypothesis_model_test_cache <- function(object){

  .hypothesis_model_plan_cache(list(object = object))
}


test_that("an object has one cache, found by its content", {

  hypothesis_cache_clear()
  withr::defer(hypothesis_cache_clear())
  object  <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "p")
  copy    <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "p")
  other   <- list(fit = list(draws = c(1, 2, 4)), data = "a", priors = "p")
  changed <- list(fit = list(draws = c(1, 2, 3)), data = "a", priors = "q")

  cache <- .hypothesis_model_test_cache(object)
  expect_true(is.environment(cache))
  expect_identical(.hypothesis_model_test_cache(object), cache)
  expect_identical(.hypothesis_model_test_cache(copy), cache)
  # other draws, or other priors, are other objects
  expect_false(identical(.hypothesis_model_test_cache(other), cache))
  expect_false(identical(.hypothesis_model_test_cache(changed), cache))
  # without an object the cache is that of the call
  expect_false(identical(.hypothesis_model_plan_cache(), .hypothesis_model_plan_cache()))
  expect_false(isTRUE(.hypothesis_model_plan_cache()[[".bounded"]]))
  expect_true(isTRUE(cache[[".bounded"]]))

  # an object that cannot be hashed has the cache of its call only
  local_mocked_bindings(
    hash = function(x, ...) stop("cannot hash"),
    .package = "rlang"
  )
  expect_false(identical(.hypothesis_model_test_cache(object),
                         .hypothesis_model_test_cache(object)))

})


test_that("a plan is computed once per object, statement and arguments", {

  hypothesis_cache_clear()
  withr::defer(hypothesis_cache_clear())
  builds <- 0L
  build  <- function(){
    builds <<- builds + 1L
    list(value = builds)
  }
  cache <- .hypothesis_model_test_cache(list(a = 1))
  expect_identical(.hypothesis_model_plan_cached(cache, "one", build), list(value = 1L))
  expect_identical(.hypothesis_model_plan_cached(cache, "one", build), list(value = 1L))
  expect_identical(.hypothesis_model_plan_cached(cache, "two", build), list(value = 2L))
  expect_identical(builds, 2L)
  # the cache of another object has its own entries, and a failing computation
  # keeps nothing
  other <- .hypothesis_model_test_cache(list(a = 2))
  expect_identical(.hypothesis_model_plan_cached(other, "one", build), list(value = 3L))
  expect_error(.hypothesis_model_plan_cached(cache, "failing", function(){
    stop("no plan")
  }), "no plan")
  expect_false(exists("failing", envir = cache, inherits = FALSE))
  # without a cache the plan is computed each time
  expect_identical(.hypothesis_model_plan_cached(NULL, "one", build), list(value = 4L))
  expect_identical(.hypothesis_model_plan_cached(NULL, "one", build), list(value = 5L))

})


test_that("the cache of an object keeps a bounded number of entries and no entry above its size", {

  hypothesis_cache_clear()
  withr::defer(hypothesis_cache_clear())
  cache <- .hypothesis_model_test_cache(list(a = 1))
  limit <- .hypothesis_model_plan_cache_limit()
  for(i in seq_len(limit + 5L)){
    .hypothesis_model_plan_cached(cache, paste0("entry-", i), function() list(i))
  }
  kept <- setdiff(ls(cache, all.names = TRUE), c(".bounded", ".order", ".sizes"))
  expect_length(kept, limit)
  # the oldest entries are the ones dropped
  expect_false(exists("entry-1", envir = cache, inherits = FALSE))
  expect_true(exists(paste0("entry-", limit + 5L), envir = cache, inherits = FALSE))
  expect_identical(cache[[".order"]], paste0("entry-", 6:(limit + 5L)))

  # entries above the size bound are computed again, not kept
  local_mocked_bindings(
    .hypothesis_model_plan_cache_max_bytes = function() 4000,
    .package = "BayesTools"
  )
  cache  <- .hypothesis_model_test_cache(list(a = 2))
  builds <- 0L
  large  <- function(){
    builds <<- builds + 1L
    numeric(1000L)
  }
  .hypothesis_model_plan_cached(cache, "large", large)
  .hypothesis_model_plan_cached(cache, "large", large)
  expect_identical(builds, 2L)
  expect_false(exists("large", envir = cache, inherits = FALSE))
  # the entries together stay under the size bound, the oldest out
  small <- function() numeric(200L)
  for(i in 1:6){
    .hypothesis_model_plan_cached(cache, paste0("small-", i), small)
  }
  expect_lte(sum(cache[[".sizes"]]), 4000)
  expect_true(exists("small-6", envir = cache, inherits = FALSE))
  expect_false(exists("small-1", envir = cache, inherits = FALSE))

})


test_that("the package keeps the caches of a few objects, the oldest out", {

  hypothesis_cache_clear()
  withr::defer(hypothesis_cache_clear())
  limit    <- .hypothesis_model_plan_registry_limit()
  objects  <- lapply(seq_len(limit + 2L), function(i) list(draws = i))
  caches   <- lapply(objects, .hypothesis_model_test_cache)
  registry <- setdiff(ls(.hypothesis_model_plan_registry, all.names = TRUE), ".order")
  expect_length(registry, limit)
  # the last objects' caches are kept and found again; the first is new again
  expect_identical(.hypothesis_model_test_cache(objects[[limit + 2L]]), caches[[limit + 2L]])
  expect_false(identical(.hypothesis_model_test_cache(objects[[1L]]), caches[[1L]]))

  hypothesis_cache_clear()
  expect_length(ls(.hypothesis_model_plan_registry, all.names = TRUE), 0L)

})


# The number of plans `expr` builds.
.hypothesis_model_test_plan_builds <- function(expr){

  builds   <- 0L
  original <- .hypothesis_model_plan
  value    <- with_mocked_bindings(
    expr,
    .hypothesis_model_plan = function(...){
      builds <<- builds + 1L
      original(...)
    },
    .package = "BayesTools"
  )

  list(value = value, builds = builds)
}


test_that("hypothesis_BF() and hypothesis_quantities() build the plan of a statement once per fit", {

  fit      <- .hypothesis_model_test_fit()
  provider <- hypothesis_provider(fit)
  metadata <- provider$metadata()
  entries  <- metadata$entries
  entry    <- as.list(entries[entries$parameter == "mu_x", setdiff(names(entries), "aliases"),
                              drop = FALSE])
  statement <- paste0("`", .hypothesis_model_quantities_reference_root(metadata, entry), "` = 0")
  plan_of   <- function(...){
    .hypothesis_model_plans(provider, statement, component = entry$component,
                            metadata = metadata, ...)[[1L]]
  }
  hypothesis_cache_clear()
  withr::defer(hypothesis_cache_clear())

  first <- .hypothesis_model_test_plan_builds(plan_of())
  expect_identical(first$builds, 1L)
  again <- .hypothesis_model_test_plan_builds(plan_of())
  expect_identical(again$builds, 0L)
  expect_identical(again$value, first$value)
  # another argument of the plan is another plan of the statement
  expect_identical(.hypothesis_model_test_plan_builds(plan_of(n_samples = 500L))$builds, 1L)
  expect_identical(.hypothesis_model_test_plan_builds(plan_of(n_samples = 500L))$builds, 0L)
  expect_identical(.hypothesis_model_test_plan_builds(plan_of(transform_scaled = TRUE))$builds, 1L)

  # hypothesis_quantities() plans every statement of the fit once, the statement
  # above among them, and once again for the next call
  hypothesis_cache_clear()
  quantities <- .hypothesis_model_test_plan_builds(hypothesis_quantities(fit))
  expect_gt(quantities$builds, 1L)
  expect_identical(.hypothesis_model_test_plan_builds(plan_of())$builds, 0L)
  repeated <- .hypothesis_model_test_plan_builds(hypothesis_quantities(fit))
  expect_identical(repeated$builds, 0L)
  expect_identical(repeated$value, quantities$value)

  # and so does hypothesis_BF(), which reads the plan the method above kept
  call <- .hypothesis_model_test_plan_builds(suppressWarnings(tryCatch(
    hypothesis_BF(fit, statement, namespace = entry$component, seed = 1),
    error = function(e) e
  )))
  expect_identical(call$builds, 0L)
  expect_s3_class(call$value, "BayesTools_hypothesis_BF")

})


test_that("the results of a cache are those of a call without one", {

  fit      <- .hypothesis_model_test_fit()
  metadata <- hypothesis_provider(fit)$metadata()
  entries  <- metadata$entries
  run <- function(cached){
    hypothesis_cache_clear()
    with_mocked_bindings(
      {
        # each statement twice: the second time from the cache
        statements <- lapply(seq_len(nrow(entries)), function(i){
          entry     <- as.list(entries[i, setdiff(names(entries), "aliases"), drop = FALSE])
          statement <- paste0(
            "`", .hypothesis_model_quantities_reference_root(metadata, entry), "` > 0"
          )
          call <- function(){
            suppressWarnings(tryCatch(
              hypothesis_BF(fit, statement, namespace = entry$component, seed = 1),
              error = function(e) e
            ))
          }
          list(call(), call())
        })
        list(statements = statements, quantities = hypothesis_quantities(fit))
      },
      .hypothesis_model_plan_cache = if(cached){
        .hypothesis_model_plan_cache
      }else{
        function(provider = NULL) new.env(parent = emptyenv())
      },
      .package = "BayesTools"
    )
  }
  withr::defer(hypothesis_cache_clear())
  cached   <- run(TRUE)
  uncached <- run(FALSE)
  expect_identical(cached, uncached)
  expect_gt(length(cached$statements), 1L)
  for(pair in cached$statements){
    expect_identical(pair[[1L]], pair[[2L]])
  }

})


# ---------------------------------------------------------------------------- #
# Plans of formula coefficients on the original scale
# ---------------------------------------------------------------------------- #

# Draws with the declared metadata of a posterior ('...' are the fields of
# posterior_metadata()).
.hypothesis_model_test_draws <- function(x, ..., class = NULL){

  if(!is.null(class)){
    class(x) <- class
  }
  fields <- list(...)
  for(field in names(fields)){
    posterior_metadata(x, field) <- fields[[field]]
  }

  x
}


# The averaged, atom-free mixed posterior of a scalar target; 'mixture' marks a
# mixture of components, 'source' the structure that declares the atoms.
.hypothesis_model_test_exp_affine_sample <- function(mixture = FALSE,
                                                     source = "single_model_structure"){

  condition_event <- structure(list(
    conditional      = character(),
    conditional_rule = "AND",
    families         = list(),
    condition_key    = "<averaged>"
  ), class = "BayesTools_condition_event")
  sample_class <- c("mixed_posteriors.simple", "mixed_posteriors")
  if(isTRUE(mixture)){
    sample_class <- c("mixed_posteriors.mixture", sample_class)
  }

  .hypothesis_model_test_draws(
    c(.15, .25, .35),
    class     = sample_class,
    condition = list(
      conditional              = character(),
      condition_key            = "<averaged>",
      resolved_condition_event = condition_event,
      averaged                 = TRUE
    ),
    atoms     = posterior_atom_attribute(source = source)
  )
}


.hypothesis_model_test_exp_affine_prior_density <- function(){

  structure(list(
    points = data.frame(x = numeric(), p = numeric())
  ), class = c("prior_linear_density", "prior_density"))
}


test_that("exp-affine targets are certified only as atom-free, averaged scalar posteriors", {

  prior_density <- .hypothesis_model_test_exp_affine_prior_density()
  # atom-free averaged posteriors are certified, whether or not they mix components
  mixture <- .hypothesis_model_test_exp_affine_sample(
    mixture = TRUE, source = "product_space_structure"
  )
  posterior <- .hypothesis_model_test_exp_affine_sample()
  expect_null(.hypothesis_model_plan_exp_affine_certify(mixture, prior_density))
  expect_null(.hypothesis_model_plan_exp_affine_certify(posterior, prior_density))

  refusal <- function(sample, density = prior_density){
    out <- .hypothesis_model_plan_exp_affine_certify(sample, density)
    expect_identical(out[["class"]],
                     c("BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable"))
    out[["reason"]]
  }
  # declared point masses of the posterior
  atomic <- mixture
  posterior_metadata(atomic, "atoms") <- posterior_atom_attribute(
    point_masses = data.frame(x = 0, mass = .1)
  )
  expect_match(refusal(atomic), "posterior is atom-free", fixed = TRUE)
  atomic <- posterior
  posterior_metadata(atomic, "atoms") <- posterior_atom_attribute(
    point_masses = data.frame(x = 0, mass = 0.25)
  )
  expect_match(refusal(atomic), "posterior is atom-free", fixed = TRUE)
  # unconditioned draws are declared by 'averaged', not by the condition key
  keyed <- posterior
  posterior_metadata(keyed, "condition") <- list(
    conditional   = character(),
    condition_key = "<averaged>",
    averaged      = FALSE
  )
  expect_match(refusal(keyed), "structural evidence for an unconditional posterior", fixed = TRUE)
  unproven <- posterior
  posterior_metadata(unproven, "condition") <- NULL
  expect_match(refusal(unproven), "structural evidence for an unconditional posterior", fixed = TRUE)
  # point masses of the prior, and draws that are no certified scalar posterior
  atomic_prior <- prior_density
  atomic_prior[["points"]] <- data.frame(x = 0, p = 0.25)
  expect_match(refusal(posterior, atomic_prior), "prior is atom-free", fixed = TRUE)
  expect_match(
    refusal(as.numeric(posterior)),
    "require a certified scalar mixed posterior",
    fixed = TRUE
  )

})


test_that("exp-affine targets are tested on their own scale with the certified prior density", {

  density <- .hypothesis_model_test_exp_affine_prior_density()
  # the target's marginal posterior (marginal_posterior() in the plan): its
  # draws with the prior density and the exact support of the map
  target  <- .hypothesis_model_test_draws(
    c(.15, .25, .35),
    class         = c("marginal_posterior.simple", "marginal_posterior"),
    prior_density = density,
    support       = posterior_support_attribute(c(0, Inf))
  )
  captured <- NULL
  # the engine evaluates a group with the draw-based method directly
  local_mocked_bindings(
    hypothesis_BF.default = function(posterior, hypothesis, ...){
      captured <<- list(posterior = posterior, hypothesis = hypothesis, ...)
      structure(
        data.frame(BF = 2, check.names = FALSE),
        hypothesis_ast = hypothesis, raw_BF = 2,
        class = c("BayesTools_hypothesis_BF", "data.frame")
      )
    },
    .package = "BayesTools"
  )
  plan <- list(
    parameter   = "log_tau_intercept",
    point       = TRUE,
    conditional = FALSE,
    draws       = list(posterior = target),
    targets     = list()
  )
  hypothesis <- hypothesis_parse("log_tau_intercept = 0.2")
  .hypothesis_model_plan_execute_scalar(
    plans           = list(plan),
    hypothesis      = hypothesis,
    provider        = list(),
    backend         = .hypothesis_model_backend("KDE"),
    logBF           = FALSE,
    BF01            = FALSE,
    seed            = 1,
    density_control = NULL,
    columns         = "all"
  )

  # no prior draws: the posterior carries the prior density, and the
  # statement keeps its original-scale values
  expect_null(captured[["prior"]])
  expect_identical(captured[["density_method"]], "KDE")
  expect_identical(posterior_metadata(captured[["posterior"]], "prior_density"), density)
  expect_equal(as.numeric(captured[["posterior"]]), c(.15, .25, .35))
  expect_identical(captured[["hypothesis"]], hypothesis)

})


test_that("scalar and exp-affine plans catch only declared measure refusal conditions", {

  condition <- tryCatch(
    .bt_formula_density_stop("Declared law unavailable.", reason = "numerical_scale_unavailable"),
    error = identity
  )
  plan <- list(parameter = "mu_x", selected = list(entry = list()), transform_scaled = TRUE)
  plan_of <- function(exp){
    if(exp){
      .hypothesis_model_plan_exp_affine(plan, list(), list(), list())
    }else{
      .hypothesis_model_plan_scalar(plan, list(), list())
    }
  }
  for(exp in c(FALSE, TRUE)){
    with_mocked_bindings({
      result <- plan_of(exp)
      expect_identical(result$refusal, .hypothesis_model_refusal_from_condition(condition))
      expect_identical(result$kind, if(exp) "exp_affine" else "scalar")
    }, .hypothesis_model_plan_scalar_available = function(...) stop(condition),
       .hypothesis_model_plan_exp_affine_available = function(...) stop(condition),
       .package = "BayesTools")
    for(class in c("BayesTools_refit_required", "BayesTools_formula_transform_unavailable",
                   "simpleError")){
      error <- structure(list(message = "Strict owner or programming failure.", call = NULL),
                         class = c(class, "error", "condition"))
      with_mocked_bindings({
        expect_error(plan_of(exp), class = class)
      }, .hypothesis_model_plan_scalar_available = function(...) stop(error),
         .hypothesis_model_plan_exp_affine_available = function(...) stop(error),
         .package = "BayesTools")
    }
  }

})


test_that("point values outside the open support of a coefficient are refused", {

  description <- "transformed coefficient 'log_tau_intercept'"
  refusal <- .hypothesis_model_plan_support_refusal(0, c(0, Inf), description)
  expect_identical(
    refusal[["reason"]],
    paste0("Point-null value 0 is outside or on the boundary of the open support ",
           "for transformed coefficient 'log_tau_intercept'.")
  )
  expect_identical(refusal[["class"]],
                   c("BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable"))
  for(value in c(0, -1, Inf, NA)){
    expect_match(
      .hypothesis_model_plan_support_refusal(value, c(0, Inf), description)[["reason"]],
      "outside or on the boundary",
      fixed = TRUE,
      info = value
    )
  }
  expect_match(.hypothesis_model_plan_support_refusal(1, c(0, 1), description)[["reason"]],
               "outside or on the boundary", fixed = TRUE)
  expect_null(.hypothesis_model_plan_support_refusal(0.5, c(0, Inf), description))
  expect_null(.hypothesis_model_plan_support_refusal(1, c(0, Inf), description))
  expect_null(.hypothesis_model_plan_support_refusal(-1, c(-Inf, Inf), description))
  expect_null(.hypothesis_model_plan_support_refusal(0, NULL, description))

})


# A fit of a formula parameter with the draws of its coefficients.
.hypothesis_model_test_route_fit <- function(parameter, coefficients, formula, prior_list, data){

  compiled <- JAGS_formula(
    formula = formula, parameter = parameter, data = data,
    prior_list = prior_list, formula_scale = TRUE
  )
  draws <- do.call(cbind, coefficients)
  colnames(draws) <- paste0(parameter, "_", names(coefficients))

  .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(draws)), compiled$prior_list,
    stats::setNames(list(compiled$formula_design), parameter),
    stats::setNames(list(compiled$formula_scale), parameter)
  )
}


test_that("the route of a formula coefficient follows its fitted transform", {

  # synthetic compiler rows, not posterior-inference evidence
  location <- .hypothesis_model_test_route_fit(
    "mu", list(intercept = c(-.2, 0, .2), x = c(-.1, 0, .1)), ~ 1 + x,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    data.frame(x = c(2, 4, 6))
  )
  intercept <- JAGS_formula_coefficient_route(location, "mu", "mu_intercept")
  expect_identical(intercept[["type"]], "affine")
  expect_identical(intercept[["formula_parameter"]], "mu")
  expect_identical(intercept[["target"]], "mu_intercept")
  expect_identical(intercept[["weights"]], c(mu_intercept = 1, mu_x = -2))
  expect_identical(intercept[["support"]], c(-Inf, Inf))
  slope <- JAGS_formula_coefficient_route(location, "mu", "mu_x")
  expect_identical(slope[["type"]], "affine")
  expect_identical(slope[["weights"]], c(mu_x = 0.5))

  # a unit log-intercept is an identity map on its positive support, which
  # keeps the primitive target of the coefficient itself
  formula <- ~ 1
  attr(formula, "log(intercept)") <- TRUE
  unit <- .hypothesis_model_test_route_fit(
    "log_tau", list(intercept = c(.2, .3, .4)), formula,
    list(intercept = prior("gamma", list(2, 2))), data.frame(row = 1:3)
  )
  route <- JAGS_formula_coefficient_route(unit, "log_tau", "log_tau_intercept")
  expect_identical(route[["type"]], "identity")
  expect_identical(route[["weights"]], c(log_tau_intercept = 1))
  expect_identical(route[["support"]], c(0, Inf))
  description <- "transformed coefficient 'log_tau_intercept'"
  for(value in c(0, -1)){
    expect_match(
      .hypothesis_model_plan_support_refusal(value, route[["support"]], description)[["reason"]],
      "outside or on the boundary",
      fixed = TRUE,
      info = value
    )
  }
  expect_null(.hypothesis_model_plan_support_refusal(1, route[["support"]], description))

  # a log-intercept with a slope is the exponential of a linear combination
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  scale <- .hypothesis_model_test_route_fit(
    "log_tau", list(intercept = c(.2, .3, .4), x = c(-.1, 0, .1)), formula,
    list(intercept = prior("gamma", list(2, 2)), x = prior("normal", list(0, 1))),
    data.frame(x = c(2, 4, 6))
  )
  exp_affine <- JAGS_formula_coefficient_route(scale, "log_tau", "log_tau_intercept")
  expect_identical(exp_affine[["type"]], "exp_affine")
  expect_identical(exp_affine[["support"]], c(0, Inf))

  # a coefficient fixed by the model has no route
  compiled <- JAGS_formula(
    formula = ~ 1 + x, parameter = "mu", data = data.frame(x = c(2, 4, 6)),
    prior_list = list(intercept = prior("point", list(3)), x = prior("point", list(0))),
    formula_scale = TRUE
  )
  fixed <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(matrix(numeric(), 3L, 0L))), compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale)
  )
  expect_identical(
    JAGS_formula_coefficient_route(fixed, "mu", "mu_x")[["reason"]],
    "The fitted coefficient 'mu_x' is structurally fixed and has no posterior hypothesis route."
  )

})


test_that("the route of a formula coefficient follows the declared map type of the transform", {

  fit <- .hypothesis_model_test_route_fit(
    "mu", list(intercept = c(-.2, 0, .2), x = c(-.1, 0, .1)), ~ 1 + x,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    data.frame(x = c(2, 4, 6))
  )
  declared <- JAGS_formula_coefficient_transform(fit, "mu")
  # artificial consumer-branch test: copies of the genuine descriptor deliberately
  # contradict or omit producer metadata; the consumer must follow the declared
  # map type and the raw weights
  declared[["matrix"]]["mu_intercept", ] <- c(1, 1)
  declared[["source_transforms"]][["mu_x"]] <- "log"
  route <- function(targets, caller = "hypothesis_BF()"){
    transform <- declared
    transform[["targets"]] <- targets
    local_mocked_bindings(
      JAGS_formula_coefficient_transform = function(...) transform,
      .package = "BayesTools"
    )
    JAGS_formula_coefficient_route(fit, "mu", "mu_intercept", caller = caller)
  }
  targets <- declared[["targets"]]
  targets[["map_type"]][targets[["target"]] == "mu_intercept"] <- "unsupported"
  unsupported <- route(targets, caller = "toy()")
  expect_identical(unsupported[["type"]], "unsupported")
  expect_match(unsupported[["reason"]], "not supported by toy()", fixed = TRUE)

  targets[["map_type"]][targets[["target"]] == "mu_intercept"] <- "affine"
  affine <- route(targets)
  expect_identical(affine[["type"]], "affine")
  expect_identical(affine[["weights"]], c(mu_intercept = 1, mu_x = 1))
  expect_identical(affine[["support"]], c(-Inf, Inf))

  # transforms without the declared map types are not certified
  for(targets in list(NULL, data.frame(target = "mu_intercept"))){
    uncertified <- route(targets)
    expect_identical(uncertified[["type"]], "unsupported")
    expect_match(uncertified[["reason"]], "lacks the certified structural metadata", fixed = TRUE)
  }

})


# ---------------------------------------------------------------------------- #
# Selection, rewriting, and the references of point statements
# ---------------------------------------------------------------------------- #

test_that("fitted hypotheses resolve and rewrite the structured AST", {

  fit      <- .hypothesis_model_test_fit()
  provider <- hypothesis_provider(fit)
  ast      <- hypothesis_parse("x > abs(x)")
  # the AST is parsed once: the selection and the rewrite take it as it is
  local_mocked_bindings(
    hypothesis_parse = function(...) stop("hypothesis was parsed twice"),
    .package = "BayesTools"
  )

  selected <- .hypothesis_model_select_parameter(
    provider   = provider,
    hypothesis = ast,
    component  = "mu"
  )
  rewritten <- .hypothesis_model_rewrite(
    hypothesis = ast,
    aliases    = selected[["aliases"]],
    parameter  = selected[["parameter"]]
  )

  expect_identical(selected[["parameter"]], "mu_x")
  expect_s3_class(selected[["resolution"]], "BayesTools_hypothesis_resolution")
  expect_s3_class(rewritten, "BayesTools_hypothesis_ast")
  expect_identical(hypothesis_render(rewritten), "mu_x > abs(mu_x)")

})


test_that("hypothesis plans rewrite aliases independently by statement", {

  fit      <- .hypothesis_model_test_fit()
  provider <- hypothesis_provider(fit)
  plans <- .hypothesis_model_plans(provider, c("x > 0", "mu_x < 0"), component = "mu")

  expect_identical(
    vapply(plans, function(plan) hypothesis_render(plan[["hypothesis"]]), character(1L)),
    c("mu_x > 0", "mu_x < 0")
  )

})


test_that("hypothesis alias rewriting preserves trailing level-reference backticks", {

  rewritten <- .hypothesis_model_rewrite(
    hypothesis = "alloc > alloc[random]",
    aliases    = list(alloc = "mu_alloc", mu_alloc = "mu_alloc"),
    parameter  = "mu_alloc"
  )

  expect_s3_class(rewritten, "BayesTools_hypothesis_ast")
  expect_equal(hypothesis_render(rewritten), "mu_alloc > `mu_alloc[random]`")

})


test_that("point-null references must be direct", {

  expect_error(
    .hypothesis_model_point_refs(hypothesis_parse("2 * mu = 0"), parameter = "mu"),
    "direct parameter or level reference"
  )
  expect_equal(
    .hypothesis_model_point_refs(hypothesis_parse("mu = 0"), parameter = "mu")[["value"]],
    0
  )
  precise <- hypothesis_parse("mu = 0.701406683025")
  expect_identical(
    .hypothesis_model_point_refs(precise, parameter = "mu")[["value"]],
    precise[["statements"]][[1L]][["left"]][["value"]]
  )
  expect_equal(
    nrow(.hypothesis_model_point_refs(hypothesis_parse("2 * mu > 0"), parameter = "mu")),
    0L
  )

})


test_that("only linear combinations of factor levels bypass the direct guard", {

  plan_of <- function(statement){
    ast <- hypothesis_parse(statement)
    .hypothesis_model_plan_new(
      statement = ast, hypothesis = ast, parameter = "mu", label = "mu",
      component = "mu"
    )
  }
  # scalar parameters need a direct point reference; factor levels are linear
  # targets whose combinations are compiled (.hypothesis_model_plan_linear())
  # the statement has to be restated: a statement error, not an unavailable test
  refusal <- .hypothesis_model_plan_direct_refusal(plan_of("2 * mu = 0"))
  expect_identical(refusal[["class"]], "BayesTools_hypothesis_statement")
  expect_match(refusal[["reason"]], "direct parameter or level reference", fixed = TRUE)
  expect_null(.hypothesis_model_plan_direct_refusal(plan_of("mu = 0")))
  expect_null(.hypothesis_model_plan_direct_refusal(plan_of("mu > 0")))

})


test_that("resolution refusals of the references are statement errors by their class", {

  provider <- hypothesis_provider(.hypothesis_model_test_fit())
  metadata <- list(catalog = NULL, entries = data.frame())
  select   <- function(hypothesis){
    tryCatch(
      .hypothesis_model_select_parameter(provider, hypothesis, "auto", metadata),
      error = identity
    )
  }
  resolution <- function(class, ...){
    structure(
      class = c(class, "BayesTools_parameter_resolution_error", "error", "condition"),
      list(message = "The statement was refused.", call = NULL, ...)
    )
  }
  # a statement without parameter symbols and a level of another component are
  # statement errors by their class, with the message and fields of the
  # resolution, whatever the statement's form
  refusals <- list(
    no_parameters      = resolution("BayesTools_hypothesis_no_parameters"),
    component_mismatch = resolution(
      "BayesTools_hypothesis_component_mismatch",
      symbol = "mu_g[a]", component = "mods"
    )
  )
  for(name in names(refusals)){
    local_mocked_bindings(
      hypothesis_resolve = function(...) stop(refusals[[name]]),
      .package = "BayesTools"
    )
    expected <- refusals[[name]]
    class(expected) <- c("BayesTools_hypothesis_statement", class(expected))
    for(statement in c("1 > 0", "mu = 0")){
      expect_identical(select(statement), expected, info = paste(name, statement))
    }
  }
  # other errors of the resolution are no statement refusals
  unclassed <- simpleError("The resolution failed.")
  local_mocked_bindings(
    hypothesis_resolve = function(...) stop(unclassed),
    .package = "BayesTools"
  )
  expect_identical(select("1 > 0"), unclassed)

})


test_that("metadata of an older build are a request to refit with the build of the provider", {

  fit      <- .hypothesis_model_test_fit()
  provider <- utils::modifyList(hypothesis_provider(fit), list(
    build       = "Toy/BayesTools",
    refit_class = "toy_refit_required"
  ))
  refit <- c("toy_refit_required", "BayesTools_refit_required", "error", "condition")

  # a catalog without entries
  local_mocked_bindings(
    hypothesis_resolve = function(...) list(occurrences = data.frame(quantity_id = character())),
    .hypothesis_model_entries_for_quantities = function(...) data.frame(),
    .package = "BayesTools"
  )
  expect_error(
    .hypothesis_model_select_parameter(
      provider, "mu = 0", "auto", list(catalog = NULL, entries = data.frame())
    ),
    paste0("Resolved hypothesis metadata are unavailable. Refit the model with ",
           "the current Toy/BayesTools build."),
    fixed = TRUE,
    class = "toy_refit_required"
  )

  # every stop on missing or unsupported fitted metadata has the class of the
  # provider with its parent "BayesTools_refit_required", and no other class
  selected <- list(
    parameter = "mu_x",
    component = "mu",
    entry     = list(formula_parameter = "mu", role = "fixed_coefficient")
  )
  transform <- structure(
    list(schema_version = 1L, target_names = "mu_x"),
    class = "BayesTools_formula_coefficient_transform"
  )
  local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    parameter_coordinates = function(...) data.frame(
      coordinate_name = "mu_z", formula_parameter = "mu", term = "z"
    ),
    .package = "BayesTools"
  )
  # a coefficient transform of an unsupported schema
  error <- tryCatch(
    .hypothesis_model_formula_coefficient_target(provider, selected),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
  expect_match(conditionMessage(error), "are unsupported", fixed = TRUE)
  # a resolved coefficient absent from the fitted transform
  transform[["schema_version"]] <- 4L
  transform[["target_names"]]   <- "mu_z"
  error <- tryCatch(
    .hypothesis_model_formula_coefficient_target(provider, selected),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_match(conditionMessage(error), "is absent from the fitted coefficient transformation",
               fixed = TRUE)
  # weights on coordinates absent from the fitted coordinate table
  error <- tryCatch(
    .hypothesis_model_target_prior_parameters(provider, c(mu_x = 1)),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
  # a factor level without fitted linear weights is stale fitted metadata, not
  # a target refusal of a current fit
  error <- tryCatch(
    .hypothesis_model_plan_level_target(
      plan     = list(draws = list(posterior = list(a = NULL)), label = "g"),
      provider = provider,
      level    = "a",
      value    = 0,
      label    = "g[a]"
    ),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_identical(
    conditionMessage(error),
    paste0("The linear weights of factor level 'g[a]' on the fitted coefficients ",
           "are unavailable. Refit the model with the current Toy/BayesTools build.")
  )

})


# ---------------------------------------------------------------------------- #
# Variances, labels, and tables of results
# ---------------------------------------------------------------------------- #

test_that("variance point statements keep their Bayes factor orientation through the SD", {

  parameter <- "(mu) study: tau2(intercept)"
  sd        <- "(mu) study: tau(intercept)"
  probabilities <- (seq_len(2001L) - 0.5) / 2001
  posterior <- stats::setNames(
    data.frame(stats::qgamma(probabilities, shape = 6, rate = 20)),
    sd
  )
  prior <- stats::setNames(
    data.frame(stats::qgamma(probabilities, shape = 2, rate = 4)),
    sd
  )

  evaluate <- function(statement){

    original <- hypothesis_rewrite(hypothesis_parse(statement), c(theta = parameter))
    transformed <- .hypothesis_model_plan_sd_hypothesis(original, sd)
    expect_equal(
      transformed[["statements"]][[1L]][["left"]][["value"]],
      sqrt(original[["statements"]][[1L]][["left"]][["value"]])
    )
    # draw-only priors have no structural prior density: the point ordinate is
    # estimated from the prior draws, which is warned about
    warned <- FALSE
    out <- withCallingHandlers(
      hypothesis_BF(
        posterior      = posterior,
        prior          = prior,
        hypothesis     = transformed,
        parameter      = sd,
        seed           = 1,
        density_method = "KDE",
        columns        = "all"
      ),
      BayesTools_inexact_ordinate = function(warning){
        warned <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    expect_true(warned, info = statement)
    restored <- .hypothesis_model_plan_sd_restore(out, original)
    # the densities are those of the variance: the SD densities over 2 * sd
    expect_equal(
      as.numeric(restored[["prior"]]),
      as.numeric(out[["prior"]]) / (2 * 0.3),
      tolerance = 1e-12
    )
    .hypothesis_model_restore_hypothesis_labels(out = restored, hypothesis = original)
  }

  implicit_equal     <- evaluate("theta = 0.09")
  explicit_not_equal <- evaluate("theta != 0.09 vs theta = 0.09")
  explicit_equal     <- evaluate("theta = 0.09 vs theta != 0.09")

  expect_identical(implicit_equal["Alternative"], explicit_not_equal["Alternative"])
  expect_identical(implicit_equal["Null"], explicit_not_equal["Null"])
  expect_identical(implicit_equal[["Alternative"]], "(mu) study: tau2(intercept) != 0.09")
  expect_identical(implicit_equal[["Null"]], "(mu) study: tau2(intercept) = 0.09")
  expect_identical(explicit_equal[["Alternative"]], "(mu) study: tau2(intercept) = 0.09")
  expect_equal(
    attr(implicit_equal, "raw_BF"),
    attr(explicit_not_equal, "raw_BF"),
    tolerance = 1e-12
  )
  expect_equal(
    attr(implicit_equal, "raw_BF") * attr(explicit_equal, "raw_BF"),
    1,
    tolerance = 1e-12
  )

})


test_that("hypothesis results restore public interaction labels", {

  probabilities <- (seq_len(2001L) - 0.5) / 2001
  posterior <- stats::setNames(
    data.frame(stats::qnorm(probabilities, mean = -0.2, sd = 0.8)),
    "mu_interaction"
  )
  prior <- stats::setNames(
    data.frame(stats::qnorm(probabilities)),
    "mu_interaction"
  )
  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = c("mu_interaction < 0", "mu_interaction < 1"),
    parameter  = "mu_interaction"
  )
  raw_BF <- attr(out, "raw_BF", exact = TRUE)

  display_hypothesis <- hypothesis_parse(c(
    "alloc:ablat[random] < 0",
    "alloc:ablat[random] < alloc:ablat[systematic]"
  ))
  out <- .hypothesis_model_restore_hypothesis_labels(
    out                = out,
    hypothesis         = display_hypothesis,
    parameter_label    = "alloc:ablat",
    producer_parameter = "mu_interaction",
    route              = "scalar"
  )

  expect_identical(out[["Alternative"]], c(
    "alloc:ablat[random] < 0",
    "alloc:ablat[random] < alloc:ablat[systematic]"
  ))
  expect_identical(out[["Null"]], c(
    "alloc:ablat[random] >= 0",
    "alloc:ablat[random] >= alloc:ablat[systematic]"
  ))
  expect_identical(rownames(out), c("alloc:ablat (1)", "alloc:ablat (2)"))
  expect_false(attr(out, "rownames", exact = TRUE))
  expect_identical(attr(out, "raw_BF", exact = TRUE), raw_BF)
  expect_identical(attr(out, "hypothesis_ast", exact = TRUE), display_hypothesis)
  printed <- capture.output(print(out))
  expect_true(any(grepl("alloc:ablat[random] < 0", printed, fixed = TRUE)))
  expect_false(any(grepl("alloc:ablat1", printed, fixed = TRUE)))
  expect_false(any(grepl("alloc:ablat (", printed, fixed = TRUE)))
  expect_false(any(grepl("mu_interaction", printed, fixed = TRUE)))

})


test_that("hypothesis labels retain expanded factor levels", {

  display     <- hypothesis_parse("Preregistered > 0")
  transformed <- hypothesis_parse("mu_Preregistered > 0")
  out <- data.frame(
    Alternative = rep("mu_Preregistered > 0", 2L),
    Null        = rep("mu_Preregistered <= 0", 2L),
    BF          = c(2, 3),
    row.names   = c(
      "mu_Preregistered[Not Pre-Registered]",
      "mu_Preregistered[Pre-Registered]"
    ),
    check.names = FALSE
  )
  attr(out, "hypothesis_ast") <- transformed

  restored <- .hypothesis_model_restore_hypothesis_labels(
    out                = out,
    hypothesis         = display,
    parameter_label    = "Preregistered",
    producer_parameter = "mu_Preregistered",
    route              = "levels"
  )

  expect_identical(restored[["Alternative"]], rep("Preregistered > 0", 2L))
  expect_identical(restored[["Null"]], rep("Preregistered <= 0", 2L))
  expect_identical(
    rownames(restored),
    c("Preregistered[Not Pre-Registered]", "Preregistered[Pre-Registered]")
  )
  expect_identical(attr(restored, "hypothesis_ast", exact = TRUE), display)

})


test_that("parameter groups bind producer logs in statement order and export visible values", {

  statements <- c("a = 356 vs a != 356", "b > 0", "a = 0")
  make_table <- function(parameter, hypothesis){
    hypothesis_BF(
      c(1, 2, 3), prior("normal", list(0, 1)), hypothesis,
      parameter = parameter, density_method = "normal", columns = "all"
    )
  }
  first  <- make_table("a", statements[c(1L, 3L)])
  second <- make_table("b", statements[2L])
  out <- .hypothesis_model_bind_parameter_results(
    list(first, second), list(c(1L, 3L), 2L), hypothesis_parse(statements)
  )
  expected <- c(attr(first, "raw_log_BF")[1L], attr(second, "raw_log_BF"),
                attr(first, "raw_log_BF")[2L])
  expect_equal(expected[1L], 710, tolerance = 1e-12)
  expect_identical(attr(out, "raw_log_BF"), expected)
  expect_length(attr(out, "numerical_diagnostics"), 3L)
  expect_identical(attr(out$BF_error, "name"), "error%(BF)")
  expect_true(any(grepl("error%(BF)", capture.output(print(out)), fixed = TRUE)))
  expect_identical(
    as.numeric(out$BF_error),
    c(as.numeric(first$BF_error)[1L], as.numeric(second$BF_error), as.numeric(first$BF_error)[2L])
  )
  expect_length(attr(out, "prior_numerical_diagnostics"), 3L)
  inverse <- update(out, BF01 = TRUE)
  expect_equal(as.numeric(inverse$BF)[1L], exp(-710), tolerance = 1e-14)

})


test_that("numerical model probability diagnostics survive the transport of a refusal", {

  draws <- structure(seq_len(5) / 5, class = c("mixed_posteriors.simple", "numeric"))
  attr(draws, "parameter") <- "theta"
  diagnostics <- list(
    model_indices = 1:2, log_prior_probabilities = log(c(.5, .5)),
    log_posterior_probabilities = c(0, -1000), stage = "posterior"
  )
  unavailable <- data.frame(
    column = "theta", measure = "atoms", reason = "Positive continuous mass is unavailable.",
    cause = "numerical_model_probability_unavailable"
  )
  unavailable$diagnostics <- list(diagnostics)
  posterior_metadata(draws, "measure_unavailable") <- unavailable

  condition <- tryCatch(posterior_atoms_free(draws), error = identity)
  refusal   <- .hypothesis_model_refusal_from_condition(condition)
  expect_identical(refusal$fields$diagnostics, diagnostics)
  expect_identical(refusal$reason, condition$message)
  returned <- tryCatch(.hypothesis_model_stop(refusal), error = identity)
  expect_identical(returned$diagnostics, diagnostics)
  malformed <- unavailable
  malformed$diagnostics[[1L]]$model_indices <- c(1L, 1L)
  expect_error(posterior_metadata(draws, "measure_unavailable") <- malformed, "invalid")

})


# ---------------------------------------------------------------------------- #
# Ordinates, warnings, and fixed draws
# ---------------------------------------------------------------------------- #

test_that("hypothesis ordinate metadata is limited to requested values", {

  ordinate <- function(value){
    posterior_ordinate_attribute(
      value          = value,
      ordinate       = value + 1,
      method         = "toy",
      density_method = "toy",
      diagnostics    = list(marker = value)
    )
  }
  posterior <- stats::rnorm(10)
  posterior_metadata(posterior, "posterior_ordinate") <- posterior_ordinate_append(
    ordinate(0), ordinate(1)
  )
  refs <- data.frame(level = NA_character_, value = 1)

  out     <- .hypothesis_model_keep_requested_ordinates(posterior, refs)
  entries <- .hypothesis_model_ordinate_entries(posterior_metadata(out, "posterior_ordinate"))

  expect_length(entries, 1L)
  expect_equal(entries[[1L]][["value"]], 1)
  expect_identical(entries[[1L]][["diagnostics"]][["marker"]], 1)

})


test_that("warning ownership distinguishes exact references and expressions", {

  statements <- c("tau = 0", "tau_total = 0", "2 * tau = 0", "g[A] = 0", "g[AA] = 0")
  ast   <- hypothesis_parse(statements)
  table <- data.frame(
    Alternative = statements,
    Null = c("tau != 0", "tau_total != 0", "2 * tau != 0", "g[A] != 0", "g[AA] != 0"),
    row.names = c("tau", "tau_total", "scaled_tau", "g[A]", "g[AA]"),
    check.names = FALSE
  )
  attr(table, "hypothesis_ast") <- ast
  record <- function(parameter, label = parameter, level = NA_character_, value = 0){
    data.frame(parameter = parameter, label = label, level = level, value = value,
               warning = "Declared ordinate warning.", stringsAsFactors = FALSE)
  }
  warning_rows <- .hypothesis_model_warning_rows
  expect_identical(warning_rows(table, record("tau")), "tau")
  expect_identical(warning_rows(table, record("tau_total")), "tau_total")
  expect_identical(warning_rows(table, record("2 * tau")), "scaled_tau")
  expect_identical(warning_rows(table, record("g", "g[A]", "A")), "g[A]")
  expect_identical(warning_rows(table, record("tau", value = .1)), character())
  expanded <- table[rep(4L, 4L), , drop = FALSE]
  rownames(expanded) <- c("g[A] (1)", "g[B] (1)", "g[A] (2)", "g[B] (2)")
  expanded$Alternative <- "g = 0"
  expanded$Null <- "g != 0"
  attr(expanded, "hypothesis_ast") <- hypothesis_parse(c("g = 0", "g = 0"))
  expect_identical(warning_rows(expanded, record("g", "g[A]", "A")), c("g[A] (1)", "g[A] (2)"))

})


test_that("primitive AST symbols match the warning records of ordinates", {

  record <- data.frame(
    label     = "mu_intercept",
    parameter = "mu_intercept",
    value     = 0,
    stringsAsFactors = FALSE
  )
  side <- hypothesis_parse("mu_intercept = 0")$statements[[1L]][["left"]]
  expect_true(.hypothesis_model_warning_side_matches(side, record))
  side <- hypothesis_parse("tau = 0")$statements[[1L]][["left"]]
  expect_false(.hypothesis_model_warning_side_matches(side, record))

})


test_that("multiple declared atomic locations are not one fixed scalar", {

  cell <- .hypothesis_model_test_draws(
    rep(.2, 3),
    atoms = posterior_atom_attribute(data.frame(x = c(.2, .4), mass = c(.4, .6)))
  )
  expect_false(.hypothesis_model_plan_draws_fixed(cell))
  duplicate <- posterior_atom_attribute(data.frame(x = c(.2, .2), mass = c(.4, .6)))
  fixed <- .hypothesis_model_test_draws(rep(.2, 3), atoms = duplicate)
  expect_true(.hypothesis_model_plan_draws_fixed(fixed))
  expect_error(posterior_atom_attribute(data.frame(x = c(-.2, .2), mass = c(0, 1))))
  expect_false(.hypothesis_model_plan_draws_fixed(numeric()))

})


test_that("region prior probabilities are evaluated for densities with provenance", {

  # a combination without a structural ordinate route (three gamma terms) has
  # provenance: its region probabilities are evaluated, although its ordinate
  # method is "unsupported_provenance" like that of a grid without provenance;
  # the method field does not route regions
  gammas <- .prior_linear_combination_density(
    list(
      x = prior("gamma", list(2, 1)),
      y = prior("gamma", list(2, 1)),
      z = prior("gamma", list(2, 1))
    ),
    c(x = 1, y = 1, z = 1)
  )
  expect_identical(prior_density_ordinate(gammas, 6)[["method"]], "unsupported_provenance")
  expect_true(.hypothesis_model_plan_density_has_provenance(gammas))
  expect_false(.hypothesis_model_plan_density_has_provenance(NULL))

})


# ---------------------------------------------------------------------------- #
# Groups of statements
# ---------------------------------------------------------------------------- #

# The Bayes factors of the rows of one call are those of the statements run one
# at a time.
.hypothesis_model_test_expect_combined_equals_single <- function(call, statements, info){

  combined <- call(statements)
  expect_s3_class(combined, "BayesTools_hypothesis_BF")
  expect_identical(nrow(combined), length(statements), info = info)
  for(i in seq_along(statements)){
    expect_equal(
      as.numeric(combined[["BF"]][[i]]), as.numeric(call(statements[[i]])[["BF"]]),
      tolerance = 1e-8, info = paste(info, statements[[i]])
    )
  }

  invisible(combined)
}


test_that("different level expressions with equal fitted weights are separate groups", {

  fit      <- .hypothesis_model_test_fit()
  provider <- hypothesis_provider(fit)
  groups <- function(statements){
    vapply(.hypothesis_model_plans(provider, statements), `[[`, character(1), "group")
  }
  # the reference level A is fixed: f[B] - f[A] and f[B] - 2 * f[A] have the
  # fitted weights of f[B], but they are different linear combinations of levels
  equal_weights <- c("f[B] - f[A] = 0.1", "f[B] - 2 * f[A] = 0.1")
  plans <- .hypothesis_model_plans(provider, equal_weights)
  expect_identical(
    lapply(plans, function(plan) plan[["targets"]][[1L]][["spec"]][["weights"]]),
    list(c("mu_f[1]" = 1), c("mu_f[1]" = 1))
  )
  expect_length(unique(groups(equal_weights)), 2L)
  # one level expression, two thresholds: one group
  expect_length(unique(groups(c("f[B] - f[A] = 0.1", "f[B] - f[A] = 0.2"))), 1L)
  # a multiple of an expression is another expression, with other weights
  expect_length(unique(groups(c("f[B] - f[C] = 0.1", "2 * f[B] - 2 * f[C] = 0.1"))), 2L)
  # the groups are those of a request in any order
  expect_identical(groups(rev(equal_weights)), rev(groups(equal_weights)))

  call <- function(statements) hypothesis_BF(fit, statements, seed = 1, columns = "all")
  .hypothesis_model_test_expect_combined_equals_single(call, equal_weights, "fit")
  same <- .hypothesis_model_test_expect_combined_equals_single(
    call, c("f[B] - f[A] = 0.1", "f[B] - f[A] = 0.2"), "same expression"
  )
  expect_identical(nrow(same), 2L)

  # the stored marginal means of the same model
  marginal <- .hypothesis_model_toy_marginal(fit, "bt_hypothesis_test_marginal_groups")
  marginal_plans <- function(provider, statements){
    lapply(statements, function(statement){
      .hypothesis_model_plan_marginal_means(provider, statement, NULL, NULL)
    })
  }
  marginal_provider <- hypothesis_provider(marginal)
  expect_length(unique(vapply(marginal_plans(marginal_provider, equal_weights),
                              `[[`, character(1), "group")), 2L)
  expect_length(unique(vapply(marginal_plans(marginal_provider,
                                             c("f[B] - f[A] = 0.1", "f[B] - f[A] = 0.2")),
                              `[[`, character(1), "group")), 1L)
  marginal_call <- function(statements){
    hypothesis_BF_model(marginal, statements, seed = 1, columns = "all")
  }
  .hypothesis_model_test_expect_combined_equals_single(marginal_call, equal_weights,
                                                       "marginal means")

})


# ---------------------------------------------------------------------------- #
# Scale of random-effect quantities
# ---------------------------------------------------------------------------- #

# The values of the predictor of the fits below: 'x' (standardized by the fit)
# has the sample SD of these values, so that the original-scale SD of a random
# slope is its fitted-scale SD times 1 / sd(x_values).
.hypothesis_model_test_x_values <- function() rep(c(2, 4, 6, 9), each = 3L)

# A fit whose random block has standardized predictors; its draws are samples
# of the priors (the posterior is the prior, the Bayes factors are near 1).
.hypothesis_model_test_scaled_fit <- function(prior_random, formula, z = FALSE,
                                              x_prior = prior("normal", list(0, 1)),
                                              formula_scale = TRUE, n = 2000L,
                                              seed = 42L, priors_hook = NULL,
                                              intercept_prior = prior("normal", list(0, 1))){

  data <- data.frame(x = .hypothesis_model_test_x_values(),
                     z = rep(c(1, 5, 2), 4L), g = factor(rep(1:4, each = 3L)))
  priors <- list(intercept = intercept_prior, x = x_prior)
  if(z) priors$z <- prior("normal", list(0, 1))
  compiled <- JAGS_formula(formula, "mu", data, priors, prior_random = prior_random,
                           formula_scale = formula_scale)
  # 'priors_hook' replaces priors of the compiled prior list
  if(is.function(priors_hook)) compiled$prior_list <- priors_hook(compiled)
  draws <- .generate_prior_sample_matrix(compiled$prior_list, n, seed = seed)

  .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(draws)), compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale)
  )
}

# A runjags-class fit of a standardized LKJ(2) block of two columns, with the
# correlation of the two columns drawn from its exact Beta marginal.
.hypothesis_model_test_scaled_lkj_fit <- function(formula, z, n = 2000L, seed = 5L){

  set.seed(seed)
  data <- data.frame(x = .hypothesis_model_test_x_values(),
                     z = rep(c(1, 5, 2), 4L), g = factor(rep(1:4, each = 3L)))
  priors <- list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))
  if(z) priors$z <- prior("normal", list(0, 1))
  compiled <- JAGS_formula(
    formula, "mu", data, priors, formula_scale = TRUE,
    prior_random = prior_random(g = random_block(
      sd = prior("uniform", list(1, 2)), cor = prior_lkj(eta = 2),
      monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
    ))
  )
  term     <- compiled$formula_design$random_effects[[1L]]
  r        <- 2 * stats::rbeta(n, 2, 2) - 1
  sd_names <- unique(term$sd_parameter_names)
  cor_names <- as.vector(outer(1:2, 1:2, function(row, column){
    paste0(term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
  }))
  chol_names <- as.vector(.bt_random_effect_cholesky_names(term, 2L))
  fixed <- cbind(mu_intercept = stats::rnorm(n), mu_x = stats::rnorm(n))
  if(z) fixed <- cbind(fixed, mu_z = stats::rnorm(n))
  sds <- matrix(stats::runif(n * length(sd_names), 1, 2), n, length(sd_names),
                dimnames = list(NULL, sd_names))
  correlation <- cbind(1, r, r, 1)
  colnames(correlation) <- cor_names
  cholesky <- cbind(1, r, 0, sqrt(1 - r^2))
  colnames(cholesky) <- chol_names
  draws <- cbind(fixed, sds, correlation, cholesky)

  chains  <- coda::mcmc.list(coda::mcmc(draws))
  carrier <- .parameter_catalog_test_fit(
    chains, compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale)
  )
  fit <- structure(
    list(mcmc = chains, sample = n, summary.pars = list(mutate = NULL),
         monitor = colnames(draws)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  for(name in c("prior_list", "formula_design", "formula_scale", "parameter_map")){
    attr(fit, name) <- attr(carrier, name)
  }
  fit <- .bt_attach_draw_geometry(fit)

  .bt_attach_fit_contract(fit)
}

# The log Bayes factor of the point statement 'q = value' for the draws of a
# quantity with the prior object that declares its law and its bounds: the
# computation of the draws, which reads no support of the catalog.
.hypothesis_model_test_prior_reference <- function(draws, prior, value){

  attr(hypothesis_BF(
    as.numeric(draws), prior, paste0("q = ", format(value, digits = 17)),
    parameter = "q", seed = 1
  ), "raw_log_BF")
}

.hypothesis_model_test_select <- function(fit, name){

  parameter_catalog_resolve(parameter_catalog(fit), name)
}

# The fit as a model whose density methods are those of the backends of
# BayesTools: the refusals of backends that other packages register in the
# session (the fixture tests load RoBMA) are not listed in its eligibility
# table.
.hypothesis_model_test_own_methods <- function(fit){

  registered <- hypothesis_density_backends()
  own <- registered$name[!is.na(registered$package) & registered$package == "BayesTools" &
                           registered$type != "stored"]
  model <- structure(list(fit = fit), class = "bt_hypothesis_test_own_methods")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_own_methods", function(object, ...){
    utils::modifyList(hypothesis_provider(object$fit),
                      list(object = object, density_methods = own))
  }, envir = asNamespace("BayesTools"))

  model
}


test_that("the support and the Bayes factors of random-effect SDs and variances are those of the requested scale", {

  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g)
  )
  original_to_fitted <- 1 / stats::sd(.hypothesis_model_test_x_values())

  for(transform_scaled in c(FALSE, TRUE)){
    k <- if(transform_scaled) original_to_fitted else 1
    info <- paste("transform_scaled =", transform_scaled)
    # the bounds of the Uniform(1, 2) prior of the SD on the scale of the draws
    sd_draws  <- parameter_mixed_posterior(
      fit, .hypothesis_model_test_select(fit, "(mu) sd(x)"), transform_scaled = transform_scaled
    )
    var_draws <- parameter_mixed_posterior(
      fit, .hypothesis_model_test_select(fit, "(mu) var(x)"), transform_scaled = transform_scaled
    )
    expect_equal(posterior_metadata(sd_draws, "support")$bounds, k * c(1, 2), info = info)
    expect_equal(posterior_metadata(var_draws, "support")$bounds, k^2 * c(1, 4), info = info)
    expect_true(all(sd_draws >= k & sd_draws <= 2 * k), info = info)

    # the Savage-Dickey ratio, with the reflection of the posterior density at
    # the bounds of the prior (the bounds themselves only where they are exact)
    values <- if(transform_scaled) c(1.05, 1.1, 1.5, 1.95) else c(1, 1.05, 1.1, 1.5, 1.95, 2)
    for(value in values){
      reference <- .hypothesis_model_test_prior_reference(
        sd_draws, prior("uniform", list(k, 2 * k)), k * value
      )
      expect_equal(
        attr(hypothesis_BF(fit, paste0("(mu) sd(x) = ", format(k * value, digits = 17)),
                           transform_scaled = transform_scaled, seed = 1), "raw_log_BF"),
        reference, tolerance = 1e-10, info = paste(info, "SD", value)
      )
      # a variance is tested through its SD
      expect_equal(
        attr(hypothesis_BF(fit, paste0("(mu) var(x) = ", format((k * value)^2, digits = 17)),
                           transform_scaled = transform_scaled, seed = 1), "raw_log_BF"),
        reference, tolerance = 1e-8, info = paste(info, "variance", value)
      )
    }
  }

  # the discovery table tests both quantities on either scale
  for(transform_scaled in list(NULL, FALSE, TRUE)){
    quantities <- hypothesis_quantities(.hypothesis_model_test_own_methods(fit),
                                        transform_scaled = transform_scaled)
    random <- quantities[quantities$parameter %in% c("(mu) sd(x)", "(mu) var(x)"), ]
    expect_identical(random$point_test, c(TRUE, TRUE), info = format(transform_scaled))
    expect_identical(random$reason, c("", ""), info = format(transform_scaled))
  }

  # a tampered selection is refused on either scale
  selection <- .hypothesis_model_test_select(fit, "(mu) sd(x)")
  selection$quantities$support[[1L]]$bounds <- c(0, 1)
  for(transform_scaled in c(FALSE, TRUE)){
    expect_error(parameter_mixed_posterior(fit, selection, transform_scaled = transform_scaled),
                 "stale or does not belong", info = transform_scaled)
  }

  # a fit without scaled predictors is on one scale
  plain <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g), formula_scale = FALSE
  )
  plain_selection <- .hypothesis_model_test_select(plain, "(mu) sd(x)")
  expect_identical(
    parameter_mixed_posterior(plain, plain_selection, transform_scaled = FALSE),
    parameter_mixed_posterior(plain, plain_selection, transform_scaled = TRUE)
  )

})


test_that("the SD of a standardized intercept column has its prior density on the fitted scale", {

  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (1 + x || g)
  )
  selection <- .hypothesis_model_test_select(fit, "(mu) sd(intercept)")
  draws <- parameter_mixed_posterior(fit, selection, transform_scaled = FALSE)
  expect_equal(posterior_metadata(draws, "support")$bounds, c(1, 2))
  expect_false(is.null(posterior_metadata(draws, "prior_density")))
  for(value in c(1.1, 1.5)){
    expect_equal(
      attr(hypothesis_BF(fit, paste0("(mu) sd(intercept) = ", value), seed = 1), "raw_log_BF"),
      .hypothesis_model_test_prior_reference(draws, prior("uniform", list(1, 2)), value),
      tolerance = 1e-10, info = value
    )
  }
  quantities <- hypothesis_quantities(fit)
  expect_true(quantities$point_test[quantities$parameter == "(mu) sd(intercept)"])

})


test_that("the support of conditioned and gated random-effect SDs is that of the requested scale", {

  slab <- prior_spike_and_slab(
    prior("uniform", list(1, 2)), prior_inclusion = prior("spike", list(0.5))
  )
  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = slab)), ~ 1 + x + (0 + x || g), n = 3000L
  )
  original_to_fitted <- 1 / stats::sd(.hypothesis_model_test_x_values())
  selection <- .hypothesis_model_test_select(fit, "(mu) sd(x)")

  for(transform_scaled in c(FALSE, TRUE)){
    k <- if(transform_scaled) original_to_fitted else 1
    info <- paste("transform_scaled =", transform_scaled)
    # without conditioning: the point at 0 and the slab
    averaged <- parameter_mixed_posterior(fit, selection, transform_scaled = transform_scaled)
    expect_equal(posterior_metadata(averaged, "support")$bounds, c(0, 2 * k), info = info)
    expect_equal(as.numeric(posterior_metadata(averaged, "atoms")$locations), 0, info = info)
    # within the inclusion event: the slab
    conditioned <- parameter_mixed_posterior(
      fit, selection, conditional = TRUE, transform_scaled = transform_scaled
    )
    expect_equal(posterior_metadata(conditioned, "support")$bounds, k * c(1, 2), info = info)
    expect_true(all(conditioned >= k & conditioned <= 2 * k), info = info)
    for(value in c(1.1, 1.5, 1.95)){
      expect_equal(
        attr(hypothesis_BF(fit, paste0("(mu) sd(x) = ", format(k * value, digits = 17)),
                           conditional = TRUE, transform_scaled = transform_scaled, seed = 1),
             "raw_log_BF"),
        .hypothesis_model_test_prior_reference(
          conditioned, prior("uniform", list(k, 2 * k)), k * value
        ),
        tolerance = 1e-10, info = paste(info, value)
      )
    }
  }

})


test_that("the correlation of a standardized block is tested on either scale", {

  fit <- .hypothesis_model_test_scaled_lkj_fit(~ 1 + x + z + (0 + x + z | g), z = TRUE)
  selection <- .hypothesis_model_test_select(fit, "(mu) cor(x,z)")
  for(transform_scaled in c(FALSE, TRUE)){
    draws <- parameter_mixed_posterior(fit, selection, transform_scaled = transform_scaled)
    expect_equal(posterior_metadata(draws, "support")$bounds, c(-1, 1),
                 info = transform_scaled)
    # (cor + 1) / 2 has the Beta(2, 2) law of an LKJ(2) correlation of two columns
    for(value in c(-0.8, 0, 0.3, 0.9)){
      expect_equal(
        attr(hypothesis_BF(fit, paste0("(mu) cor(x,z) = ", value),
                           transform_scaled = transform_scaled, seed = 1), "raw_log_BF"),
        .hypothesis_model_test_prior_reference(
          (as.numeric(draws) + 1) / 2, prior("beta", list(2, 2)), (value + 1) / 2
        ),
        tolerance = 1e-8, info = paste(transform_scaled, value)
      )
    }
  }
  quantities <- hypothesis_quantities(fit)
  expect_true(quantities$point_test[quantities$parameter == "(mu) cor(x,z)"])

  # the correlation of the intercept and a standardized slope on the fitted scale
  fit <- .hypothesis_model_test_scaled_lkj_fit(~ 1 + x + (1 + x | g), z = FALSE)
  draws <- parameter_mixed_posterior(
    fit, .hypothesis_model_test_select(fit, "(mu) cor(intercept,x)"), transform_scaled = FALSE
  )
  expect_equal(posterior_metadata(draws, "support")$bounds, c(-1, 1))
  expect_equal(
    attr(hypothesis_BF(fit, "(mu) cor(intercept,x) = 0.3", seed = 1), "raw_log_BF"),
    .hypothesis_model_test_prior_reference(
      (as.numeric(draws) + 1) / 2, prior("beta", list(2, 2)), 0.65
    ),
    tolerance = 1e-8
  )

})


test_that("the support and the Bayes factors of a standardized coefficient are those of the requested scale", {

  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g), x_prior = prior("uniform", list(-1, 2))
  )
  original_to_fitted <- 1 / stats::sd(.hypothesis_model_test_x_values())

  for(transform_scaled in c(FALSE, TRUE)){
    k <- if(transform_scaled) original_to_fitted else 1
    info <- paste("transform_scaled =", transform_scaled)
    draws <- as_mixed_posteriors(fit, "mu_x", transform_scaled = transform_scaled)$mu_x
    expect_equal(posterior_metadata(draws, "support")$bounds, k * c(-1, 2), info = info)
    for(value in c(-0.9, 0.2, 1.8)){
      expect_equal(
        attr(hypothesis_BF(fit, paste0("x = ", format(k * value, digits = 17)),
                           transform_scaled = transform_scaled, seed = 1), "raw_log_BF"),
        .hypothesis_model_test_prior_reference(draws, prior("uniform", list(-k, 2 * k)), k * value),
        tolerance = 1e-8, info = paste(info, value)
      )
    }
  }

})


test_that("method refusals follow the evaluation method of a registered draws backend", {

  # the exponential of a linear combination of the fitted coefficients is
  # tested with KDE only
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  set.seed(9)
  fit <- .hypothesis_model_test_route_fit(
    "log_tau", list(intercept = stats::rgamma(400L, 2, 2), x = stats::rnorm(400L, 0, 0.3)),
    formula,
    list(intercept = prior("gamma", list(2, 2)), x = prior("normal", list(0, 1))),
    data.frame(x = c(2, 4, 6))
  )
  alias <- hypothesis_density_backend("normalAlias", type = "draws", method = "normal")
  with_toy_backend(alias, {
    refusal <- function(method){
      tryCatch(
        hypothesis_BF(fit, "intercept = 0.5", transform_scaled = TRUE,
                      density_method = method, seed = 1),
        error = identity
      )
    }
    builtin <- refusal("normal")
    expect_s3_class(builtin, "BayesTools_hypothesis_method")
    registered <- refusal("normalAlias")
    expect_s3_class(registered, "BayesTools_hypothesis_method")
    expect_identical(conditionMessage(registered), conditionMessage(builtin))
    expect_s3_class(
      hypothesis_BF(fit, "intercept = 0.5", transform_scaled = TRUE, seed = 1),
      "BayesTools_hypothesis_BF"
    )

    # the table lists the methods that evaluate a quantity
    quantities <- hypothesis_quantities(fit, transform_scaled = TRUE)
    methods    <- quantities$point_test_methods[quantities$alias == "intercept"]
    expect_match(methods, "KDE", fixed = TRUE)
    expect_false(grepl("normalAlias", methods, fixed = TRUE))

    # a statement that the normal approximation evaluates is evaluated by the alias
    expect_identical(
      attr(hypothesis_BF(fit, "x = 0.1", density_method = "normalAlias", seed = 1), "raw_log_BF"),
      attr(hypothesis_BF(fit, "x = 0.1", density_method = "normal", seed = 1), "raw_log_BF")
    )

    # and the stored marginal means refuse the normal approximation by any name
    marginal <- .hypothesis_model_toy_marginal(
      .hypothesis_model_test_fit(), "bt_hypothesis_test_marginal_alias",
      density_methods = c("KDE", "normal", "normalAlias")
    )
    for(method in c("normal", "normalAlias")){
      expect_error(
        hypothesis_BF_model(marginal, "f[B] = 0.2", density_method = method),
        "'normal' density_method is not supported for marginal-means hypotheses.",
        fixed = TRUE, class = "BayesTools_hypothesis_method", info = method
      )
    }
  })

})


test_that("the fit method takes the statements as character strings", {

  fit <- .hypothesis_model_test_fit()
  parsed <- hypothesis_parse("x = 0.1")
  for(call in list(
    function() hypothesis_BF(fit, parsed, seed = 1),
    function() hypothesis_BF_model(fit, parsed, seed = 1)
  )){
    condition <- tryCatch(call(), error = identity)
    expect_s3_class(condition, "BayesTools_hypothesis_statement")
    expect_match(conditionMessage(condition), "must be character strings", fixed = TRUE)
    expect_match(conditionMessage(condition), "Pass the statements as character strings",
                 fixed = TRUE)
  }
  # the rendered statement is accepted, and the draw-based method takes the parsed one
  rendered <- hypothesis_render(parsed)
  expect_type(rendered, "character")
  expect_identical(
    attr(hypothesis_BF(fit, rendered, seed = 1), "raw_log_BF"),
    attr(hypothesis_BF(fit, "x = 0.1", seed = 1), "raw_log_BF")
  )
  set.seed(1)
  draws <- stats::rnorm(500, 0.4, 0.3)
  normal <- prior("normal", list(0, 1))
  expect_identical(
    attr(hypothesis_BF(draws, normal, parsed, parameter = "x", seed = 1), "raw_log_BF"),
    attr(hypothesis_BF(draws, normal, "x = 0.1", parameter = "x", seed = 1), "raw_log_BF")
  )

})


test_that("the messages name the argument that the provider selects the component with", {

  # three formula parameters: a statement on two of them has no component of the third
  set.seed(5)
  data <- data.frame(x = c(-2, -1, 0, 3, 5, 7))
  compile <- function(parameter){
    JAGS_formula(~ 1 + x, parameter, data, list(
      intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))
    ), formula_scale = TRUE)
  }
  parameters <- c("mu", "log_sigma", "log_tau")
  compiled   <- stats::setNames(lapply(parameters, compile), parameters)
  columns    <- paste0(rep(parameters, each = 2L), "_", c("intercept", "x"))
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(matrix(
      stats::rnorm(400L * length(columns), 0.3, 0.5), 400L, length(columns),
      dimnames = list(NULL, columns)
    ))),
    do.call(c, unname(lapply(compiled, `[[`, "prior_list"))),
    lapply(compiled, `[[`, "formula_design"), lapply(compiled, `[[`, "formula_scale")
  )
  model <- structure(list(fit = fit), class = "bt_hypothesis_test_vocabulary")
  registerS3method("hypothesis_provider", "bt_hypothesis_test_vocabulary", function(object, ...){
    modifyList(hypothesis_provider(object$fit),
               list(object = object, component_argument = "component"))
  }, envir = asNamespace("BayesTools"))
  expect_identical(hypothesis_provider(fit)$component_argument, "namespace")
  expect_identical(hypothesis_provider(model)$component_argument, "component")

  mismatch <- function(object, statement){
    condition <- tryCatch(
      hypothesis_BF_model(object, statement, namespace = "log_tau"), error = identity
    )
    expect_s3_class(condition, "BayesTools_hypothesis_component_mismatch")
    expect_s3_class(condition, "BayesTools_hypothesis_statement")
    conditionMessage(condition)
  }
  # a statement of one component
  expect_identical(
    mismatch(fit, "mu_x = 0.1"),
    "The 'hypothesis' argument selects namespace = 'mu' but 'namespace' was set to 'log_tau'."
  )
  expect_identical(
    mismatch(model, "mu_x = 0.1"),
    "The 'hypothesis' argument selects component = 'mu' but 'component' was set to 'log_tau'."
  )
  # a statement of several components
  expect_identical(mismatch(fit, "mu_x = log_sigma_x"),
                   "The hypothesis does not resolve to namespace = 'log_tau'.")
  expect_identical(mismatch(model, "mu_x = log_sigma_x"),
                   "The hypothesis does not resolve to component = 'log_tau'.")
  # the advice of an ambiguous statement is the provider's own
  ambiguous <- tryCatch(hypothesis_BF(fit, "x = 0.1", seed = 1), error = identity)
  expect_s3_class(ambiguous, "BayesTools_hypothesis_ambiguous")
  expect_match(conditionMessage(ambiguous), "Set 'namespace' to", fixed = TRUE)

  # the field is validated
  for(value in list(1, NA_character_, "", c("a", "b"))){
    broken <- structure(list(fit = fit), class = "bt_hypothesis_test_vocabulary_broken")
    registerS3method("hypothesis_provider", "bt_hypothesis_test_vocabulary_broken",
      function(object, ...){
        modifyList(hypothesis_provider(object$fit),
                   list(object = object, component_argument = value))
      }, envir = asNamespace("BayesTools"))
    expect_error(hypothesis_BF_model(broken, "mu_x = 0.1"),
                 "'component_argument' must be a string or NULL", fixed = TRUE)
  }

})


# A runjags-class fit of a standardized block of the nested variance
# allocations of the study and drug terms; it monitors the Cholesky factor of
# the correlation of the study block and not the primitive LKJ variables.
.hypothesis_model_test_scaled_allocation_fit <- function(n = 1000L, seed = 42L){

  compiled <- JAGS_formula(
    ~ 1 + random(1 + x | study, name = "study", covariance = "us") +
      random(1 | drug, name = "drug", covariance = "diag"),
    "mu",
    data.frame(study = factor(rep(1:4, each = 3L)), drug = factor(rep(1:3, 4L)),
               x = rep(c(2, 4, 6), 4L)),
    list(intercept = prior("normal", list(0, 1))), formula_scale = TRUE,
    prior_random = prior_random(
      root  = random_variance_allocation(
        name = "root", terms = c(study = "study", drug = "drug"),
        sd = prior("gamma", list(2, 2))
      ),
      split = random_variance_allocation(
        name = "split", terms = "study", target = "sd_component",
        parent = allocation_ref("root", "study")
      ),
      study = random_block(cor = prior_lkj(2))
    )
  )
  draws <- .generate_prior_sample_matrix(compiled$prior_list, n, seed = seed)
  term  <- compiled$formula_design$random_effects[[1L]]
  r     <- 2 * stats::rbeta(n, 2, 2) - 1
  cholesky <- cbind(1, r, 0, sqrt(1 - r^2))
  colnames(cholesky) <- as.vector(.bt_random_effect_cholesky_names(term, 2L))
  draws <- cbind(
    draws[, !colnames(draws) %in% term$correlation$primitive_names, drop = FALSE],
    cholesky[, !colnames(cholesky) %in% colnames(draws), drop = FALSE]
  )

  chains  <- coda::mcmc.list(coda::mcmc(draws))
  carrier <- .parameter_catalog_test_fit(
    chains, compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale)
  )
  fit <- structure(
    list(mcmc = chains, sample = n, summary.pars = list(mutate = NULL),
         monitor = colnames(draws)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  for(name in c("prior_list", "formula_design", "formula_scale", "parameter_map")){
    attr(fit, name) <- attr(carrier, name)
  }
  fit <- .bt_attach_draw_geometry(fit)

  .bt_attach_fit_contract(fit)
}


test_that("the prior draws of a standardized block are extracted with the catalog of the fitted scale", {

  fit <- .hypothesis_model_test_scaled_allocation_fit()
  selection <- .hypothesis_model_test_select(fit, "(mu) study: sd(x)")
  scaled <- .bt_parameter_fitted_scale(fit, selection)
  # the original-scale SD of the standardized slope needs the Cholesky factor
  # (not monitored); the fitted-scale SD does not
  expect_true(any(grepl("_CORx_L", selection$quantities$extraction_key[[1L]]$dependencies)))
  expect_false(any(grepl("_CORx_L", scaled$selection$quantities$extraction_key[[1L]]$dependencies)))

  # region hypotheses on a quantity whose prior has no density with provenance
  # take their prior probabilities from prior draws
  result <- hypothesis_BF(fit, "(mu) study: sd(x) > 0.5", seed = 1, n_samples = 1000L)
  expect_s3_class(result, "BayesTools_hypothesis_BF")
  raw <- transform_prior_samples(fit, n_samples = 1000L, seed = 1, formula_scale = list())
  prior_values <- as.numeric(as.matrix(parameter_draws(
    JAGS_with_draws(scaled$fit, coda::mcmc.list(coda::mcmc(raw))), scaled$selection
  )))
  posterior <- parameter_mixed_posterior(fit, selection, transform_scaled = FALSE)
  expect_identical(
    attr(result, "raw_log_BF"),
    attr(hypothesis_BF.default(
      posterior = as.numeric(posterior), prior = prior_values,
      hypothesis = "q > 0.5", parameter = "q", seed = 1
    ), "raw_log_BF")
  )
  # the provider callback is the same extraction
  provider <- hypothesis_provider(fit)
  selected <- provider$random$select("(mu) study: sd(x)", FALSE)
  expect_identical(
    as.numeric(provider$random$prior_draws("(mu) study: sd(x)", selected, FALSE, 1000L, 1)$samples),
    prior_values
  )

})


test_that("a quantity that is fixed on the fitted scale is refused as fixed", {

  # the fitted SD of the intercept is fixed and the SD of the slope is sampled:
  # the original-scale SD of the intercept combines both and is sampled
  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("spike", list(1.5)))),
    ~ 1 + x + (1 + x || g), n = 1000L, seed = 2L,
    priors_hook = function(compiled){
      sd_names <- unique(compiled$formula_design$random_effects[[1L]]$sd_parameter_names)
      compiled$prior_list[[sd_names[[2L]]]] <- prior("uniform", list(1, 2))
      compiled$prior_list
    }
  )
  selection <- .hypothesis_model_test_select(fit, "(mu) sd(intercept)")
  expect_false(identical(selection$quantities$status, "structural"))
  expect_identical(
    .bt_parameter_fitted_scale(fit, selection)$selection$quantities$status, "structural"
  )
  for(statement in c("(mu) sd(intercept) > 1", "(mu) var(intercept) > 1")){
    expect_error(hypothesis_BF(fit, statement, seed = 1),
                 class = "BayesTools_hypothesis_fixed", info = statement)
  }
  quantities <- hypothesis_quantities(fit)
  fixed <- quantities[quantities$parameter %in% c("(mu) sd(intercept)", "(mu) var(intercept)"), ]
  expect_identical(fixed$point_test, c(FALSE, FALSE))
  expect_identical(fixed$direction_test, c(FALSE, FALSE))
  expect_match(fixed$reason, "is fixed by the fitted model", fixed = TRUE)
  # the sampled SD of the slope is tested
  expect_s3_class(hypothesis_BF(fit, "(mu) sd(x) > 1.2", seed = 1), "BayesTools_hypothesis_BF")
  expect_true(quantities$point_test[quantities$parameter == "(mu) sd(x)"])

})


test_that("the catalog of the fitted scale describes the fitted coordinates of a standardized block", {

  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g)
  )
  k <- 1 / stats::sd(.hypothesis_model_test_x_values())
  original <- parameter_catalog(fit)
  fitted   <- parameter_catalog(fit, transform_scaled = FALSE)
  expect_identical(parameter_catalog(fit, transform_scaled = TRUE), original)
  expect_s3_class(fitted, "BayesTools_parameter_catalog")
  # the same quantities with the same identifiers, names, and aliases
  expect_identical(fitted$quantities$quantity_id, original$quantities$quantity_id)
  expect_identical(fitted$quantities$canonical_name, original$quantities$canonical_name)
  expect_identical(fitted$aliases, original$aliases)
  expect_error(parameter_catalog(fit, transform_scaled = NA), "transform_scaled")
  expect_error(parameter_catalog(fit, transform_scaled = "no"), "transform_scaled")

  row <- function(catalog, name) catalog$quantities[catalog$quantities$canonical_name == name, ]
  for(name in c("(mu) sd(x)", "(mu) var(x)")){
    power <- if(identical(name, "(mu) sd(x)")) 1 else 2
    bounds <- c(1, 2)^power
    expect_equal(row(original, name)$support[[1L]]$bounds, k^power * bounds, info = name)
    expect_equal(row(fitted, name)$support[[1L]]$bounds, bounds, info = name)
    expect_true(row(fitted, name)$support[[1L]]$exact, info = name)
  }
  expect_identical(row(original, "(mu) sd(x)")$source_type, "one_to_one_transform")
  expect_identical(row(fitted, "(mu) sd(x)")$source_type, "identity")
  expect_equal(row(original, "(mu) sd(x)")$extraction_key[[1L]]$source_scale, k)
  expect_equal(row(fitted, "(mu) sd(x)")$extraction_key[[1L]]$source_scale, 1)
  expect_identical(row(fitted, "(mu) sd(x)")$status, "sampled")
  # quantities that the scales do not touch are shared
  expect_identical(row(fitted, "mu_x"), row(original, "mu_x"))

  # a status that differs: the intercept SD is fixed on the fitted scale and
  # combined with the sampled slope SD on the original scale
  fixed <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("spike", list(1.5)))),
    ~ 1 + x + (1 + x || g), n = 1000L, seed = 2L,
    priors_hook = function(compiled){
      sd_names <- unique(compiled$formula_design$random_effects[[1L]]$sd_parameter_names)
      compiled$prior_list[[sd_names[[2L]]]] <- prior("uniform", list(1, 2))
      compiled$prior_list
    }
  )
  expect_identical(row(parameter_catalog(fixed), "(mu) sd(intercept)")$status, "sampled")
  fixed_row <- row(parameter_catalog(fixed, transform_scaled = FALSE), "(mu) sd(intercept)")
  expect_identical(fixed_row$status, "structural")
  expect_equal(fixed_row$fixed_value, 1.5)

  # a fit without scaled predictors has one catalog, also with the scales
  # of its formula blanked
  plain <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g), formula_scale = FALSE
  )
  expect_identical(parameter_catalog(plain, transform_scaled = FALSE), parameter_catalog(plain))
  blanked <- fit
  attr(blanked, "formula_scale") <- list()
  expect_identical(parameter_catalog(blanked), original)
  expect_identical(parameter_catalog(blanked, transform_scaled = FALSE), fitted)

})


test_that("selections of either catalog are extracted on the scale of their catalog", {

  fit <- .hypothesis_model_test_scaled_fit(
    prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
    ~ 1 + x + (0 + x || g)
  )
  k <- 1 / stats::sd(.hypothesis_model_test_x_values())
  original <- .hypothesis_model_test_select(fit, "(mu) sd(x)")
  fitted   <- parameter_catalog_resolve(
    parameter_catalog(fit, transform_scaled = FALSE), "(mu) sd(x)"
  )
  expect_false(identical(original$quantities, fitted$quantities))
  expect_identical(original$quantity_id, fitted$quantity_id)

  # draws, prior densities, and transforms
  original_draws <- parameter_draws(fit, original)
  fitted_draws   <- parameter_draws(fit, fitted)
  expect_true(all(as.matrix(original_draws) >= k & as.matrix(original_draws) <= 2 * k))
  expect_true(all(as.matrix(fitted_draws) >= 1 & as.matrix(fitted_draws) <= 2))
  expect_equal(as.numeric(as.matrix(original_draws)) / as.numeric(as.matrix(fitted_draws)),
               rep(k, nrow(as.matrix(fitted_draws))), tolerance = 1e-10)
  expect_equal(posterior_metadata(fitted_draws, "support")[[1L]]$bounds, c(1, 2))
  for(case in list(list(selection = original, k = k), list(selection = fitted, k = 1))){
    density <- parameter_prior_density(fit, case$selection)
    expect_equal(range(density$density$x), case$k * c(1, 2), tolerance = 1e-6)
    transform <- parameter_transform(fit, case$selection)
    expect_equal(parameter_transform_forward(c(1, 2), transform), case$k * c(1, 2))
  }
  # the same quantity as the mixed posterior of the fitted scale, from either catalog
  mixed <- parameter_mixed_posterior(fit, original, transform_scaled = FALSE)
  expect_identical(parameter_mixed_posterior(fit, fitted, transform_scaled = FALSE), mixed)
  expect_equal(as.numeric(mixed), as.numeric(as.matrix(fitted_draws)))
  expect_equal(posterior_metadata(mixed, "support")$bounds, c(1, 2))
  # the fit with the scales blanked (the extraction of a downstream package)
  blanked <- fit
  attr(blanked, "formula_scale") <- list()
  expect_equal(as.numeric(as.matrix(parameter_draws(blanked, original))),
               as.numeric(as.matrix(fitted_draws)))
  expect_equal(as.numeric(as.matrix(parameter_draws(blanked, fitted))),
               as.numeric(as.matrix(fitted_draws)))
  expect_no_error(random_effects_marginal_update_plan(fit, fitted))
  expect_no_error(parameter_gate_states(fit, fitted))

  # the selection of the fitted-scale catalog is not a selection of the original scale
  expect_error(parameter_mixed_posterior(fit, fitted, transform_scaled = TRUE),
               "belongs to the catalog of the fitted scale", fixed = TRUE,
               class = "BayesTools_parameter_resolution_error")
  expect_error(parameter_mixed_posterior(fit, fitted),
               "use 'transform_scaled = FALSE'", fixed = TRUE)
  expect_equal(posterior_metadata(parameter_mixed_posterior(fit, original), "support")$bounds,
               k * c(1, 2))

  # selections of another catalog are refused by every function
  tampered <- fitted
  tampered$quantities$support[[1L]]$bounds <- c(0, 1)
  foreign <- parameter_catalog_resolve(
    parameter_catalog(.hypothesis_model_test_fit(), transform_scaled = FALSE), "sigma"
  )
  for(selection in list(tampered, foreign)){
    for(call in list(
      function() parameter_draws(fit, selection),
      function() parameter_prior_density(fit, selection),
      function() parameter_transform(fit, selection),
      function() parameter_gate_states(fit, selection),
      function() parameter_mixed_posterior(fit, selection, transform_scaled = FALSE),
      function() parameter_mixed_posterior(fit, selection, transform_scaled = TRUE)
    )){
      expect_error(call(), "does not belong to this fitted catalog", class = "BayesTools_parameter_resolution_error")
    }
  }
  # a selection of several quantities is refused where one is needed
  both <- parameter_catalog_resolve(
    parameter_catalog(fit, transform_scaled = FALSE), "mu_x"
  )
  both$quantities <- rbind(both$quantities, fitted$quantities)
  both$quantity_id <- both$quantities$quantity_id
  expect_error(parameter_mixed_posterior(fit, both, transform_scaled = FALSE),
               "exactly one parameter quantity", fixed = TRUE)

})


test_that("prior draws of a standardized block are extracted with a selection of the fitted-scale catalog", {

  fit <- .hypothesis_model_test_scaled_allocation_fit()
  raw <- transform_prior_samples(fit, n_samples = 200L, seed = 1, formula_scale = list())
  view <- JAGS_with_draws(fit, coda::mcmc.list(coda::mcmc(raw)))
  original <- .hypothesis_model_test_select(fit, "(mu) study: sd(x)")
  fitted   <- parameter_catalog_resolve(
    parameter_catalog(fit, transform_scaled = FALSE), "(mu) study: sd(x)"
  )
  # the original-scale SD needs the Cholesky coordinates, which are not monitored
  expect_error(parameter_draws(view, original),
               class = "BayesTools_draws_view_coordinate_unavailable")
  draws <- parameter_draws(view, fitted)
  expect_equal(dim(as.matrix(draws)), c(200L, 1L))
  expect_true(all(is.finite(as.matrix(draws))))
  expect_equal(posterior_metadata(draws, "support")[[1L]]$bounds, c(0, Inf))

})


test_that("the fitted-scale catalog keeps the names of the original scale", {

  # a log-intercept of a scaled formula is named by its original-scale label,
  # 'exp(intercept)', on both scales: with a standardized random slope the
  # quantities of the two catalogs differ, without it they are the same
  for(random in c(TRUE, FALSE)){
    formula <- if(random) ~ 1 + x + (0 + x || g) else ~ 1 + x
    attr(formula, "log(intercept)") <- TRUE
    fit <- .hypothesis_model_test_scaled_fit(
      if(random) prior_random(g = random_block(sd = prior("uniform", list(1, 2)))),
      formula, intercept_prior = prior("lognormal", list(0, 0.5))
    )
    original <- parameter_catalog(fit)
    fitted   <- parameter_catalog(fit, transform_scaled = FALSE)
    expect_identical(identical(fitted$quantities, original$quantities), !random)
    expect_true(all(c("exp(intercept)", "(mu) exp(intercept)") %in% original$aliases$alias))
    expect_identical(fitted$aliases, original$aliases)
    intercept <- original$quantities$quantity_id[original$quantities$canonical_name == "mu_intercept"]
    for(catalog in list(original, fitted)){
      for(alias in c("exp(intercept)", "(mu) exp(intercept)")){
        expect_identical(
          parameter_catalog_resolve(catalog, alias, namespace = "mu")$quantity_id, intercept,
          info = paste(random, alias)
        )
      }
    }
    if(random){
      expect_identical(
        parameter_catalog_resolve(fitted, "(mu) sd(x)")$quantity_id,
        parameter_catalog_resolve(original, "(mu) sd(x)")$quantity_id
      )
    }
  }

})
