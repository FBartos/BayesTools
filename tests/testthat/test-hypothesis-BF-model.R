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

# A model with a marginal-means provider that stores the marginal posterior of
# the levels of 'f'; 'overrides' replace elements of the provider.
.hypothesis_model_toy_marginal <- function(fit, class, ...){

  samples <- as_mixed_posteriors(fit, "mu_f", n_prior_samples = 2000)
  levels  <- marginal_posterior(samples, "mu_f", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 2000)
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
