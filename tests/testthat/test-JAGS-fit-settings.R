skip_if_not_test_profile("unit")

test_that("parallel JAGS checks the actual workers' package builds", {

  native <- getLoadedDLLs()[["BayesTools"]][["path"]]
  expect_identical(
    unname(.JAGS_package_builds("BayesTools")$BayesTools$dll),
    unname(tools::md5sum(native))
  )
  expected <- .JAGS_package_builds("stats")
  actual <- expected
  testthat::local_mocked_bindings(
    clusterCall = function(cl, fun, packages) list(actual, actual),
    .package = "parallel"
  )
  expect_identical(unname(.JAGS_require_packages("stats", cl = list())), TRUE)
  message <- paste0(
    "Parallel JAGS fitting is unavailable with mismatched package versions, R code, or native builds: 'stats'. ",
    "Install the parent-session builds into a library and set 'R_LIBS_USER' ",
    "to that library before starting R and its workers."
  )
  actual$stats$version <- "0.0.0"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  actual <- expected
  actual$stats$dll[] <- "different-build"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  actual <- expected
  actual$stats$r_code <- "different-code"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  expect_error(.JAGS_require_packages("stats", cl = list(),
    operation = "Parallel zplot density computation"),
    sub("Parallel JAGS fitting", "Parallel zplot density computation", message, fixed = TRUE),
    fixed = TRUE)
  actual$stats <- NULL
  expect_error(
    .JAGS_require_packages("stats", cl = list()),
    "Required packages are not available: 'stats'.", fixed = TRUE
  )
})

test_that("JAGS fit settings reject missing and non-finite controls", {
  valid <- list(
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    thin = 1,
    autofit = FALSE,
    parallel = FALSE,
    cores = 1,
    silent = TRUE,
    seed = 1
  )
  invalid <- list(
    chains = NA_real_,
    adapt = Inf,
    burnin = NA_real_,
    sample = Inf,
    thin = NA_real_,
    autofit = NA,
    parallel = NA,
    cores = Inf,
    silent = NA,
    seed = NaN
  )

  expect_silent(do.call(JAGS_check_and_list_fit_settings, valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_check_and_list_fit_settings, settings),
      paste0("'", name, "'")
    )
  }
})

test_that("JAGS autofit settings reject missing and non-finite controls", {
  valid <- list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = FALSE,
    monitor = c("mu", "sigma"),
    allow_not_assessable = FALSE
  )
  invalid <- list(
    max_Rhat = Inf,
    min_ESS = Inf,
    max_error = NaN,
    max_SD_error = NA_real_,
    max_time = list(time = Inf, unit = "mins"),
    sample_extend = Inf,
    restarts = NA_real_,
    max_extend = Inf,
    check_indicators = NA,
    monitor = NA_character_,
    allow_not_assessable = NA
  )

  expect_silent(JAGS_check_and_list_autofit_settings(valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      JAGS_check_and_list_autofit_settings(settings),
      paste0("'", if(name == "max_time") "max_time:time" else name, "'")
    )
  }

  invalid_unit <- valid
  invalid_unit$max_time$unit <- NA_character_
  expect_error(
    JAGS_check_and_list_autofit_settings(invalid_unit),
    "'max_time:unit'"
  )

  empty_monitor <- valid
  empty_monitor$monitor <- character()
  expect_error(
    JAGS_check_and_list_autofit_settings(empty_monitor),
    "'monitor' argument must select at least one parameter"
  )
})

test_that("JAGS chain seeds do not overlap across adjacent seeds and keep the initial values", {

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )
  chain_seeds <- function(seed, chains){
    vapply(
      JAGS_get_inits(prior_list, chains = chains, seed = seed),
      function(inits) inits[[".RNG.seed"]],
      numeric(1)
    )
  }

  # The documented derivation: R's generator after set.seed(seed).
  expected_seeds <- local({
    set.seed(11)
    sample.int(.Machine$integer.max, 3)
  })
  expect_identical(chain_seeds(11, 3), as.numeric(expected_seeds))

  # 'seed + chain' gave chain k + 1 of seed s the seed of chain k of seed s + 1.
  for(seed in c(1:50, 665, 666666)){
    seeds      <- chain_seeds(seed, 8)
    next_seeds <- chain_seeds(seed + 1, 8)
    expect_false(any(seeds[-1] == next_seeds[-8]), info = paste("seed", seed))
    expect_length(intersect(seeds, next_seeds), 0)
    expect_length(unique(seeds), 8)
  }

  # A chain's seed does not depend on the number of chains.
  for(seed in c(1, 42, 666)){
    all_seeds <- chain_seeds(seed, 8)
    for(chains in 1:7){
      expect_identical(chain_seeds(seed, chains), all_seeds[seq_len(chains)])
    }
  }

  # Initial values are those of 40b07c5 (drawn right after set.seed(seed));
  # only '.RNG.seed' changed. The caller's generator state is restored.
  set.seed(5)
  state_before_inits <- .Random.seed
  inits <- JAGS_get_inits(prior_list, chains = 3, seed = 11)
  expect_identical(.Random.seed, state_before_inits)
  set.seed(11)
  expected_inits <- lapply(1:3, function(chain){
    list(
      mu    = rng(prior_list[["mu"]], 1),
      sigma = rng(prior_list[["sigma"]], 1)
    )
  })
  expect_identical(
    lapply(inits, function(chain_inits) chain_inits[c("mu", "sigma")]),
    expected_inits
  )
  expect_identical(
    vapply(inits, function(chain_inits) chain_inits[[".RNG.name"]], character(1)),
    rep("base::Super-Duper", 3)
  )
  # Values printed from 40b07c5 with 17 significant digits. The generator's
  # normal deviates (qnorm()) differ in the last bit where the compiler fuses
  # multiply-adds (Apple silicon), as can R's reading of 17-digit literals
  # where its long double is double, so the stream is pinned to 1e-13.
  expect_equal(
    vapply(inits, function(chain_inits) chain_inits[["mu"]], numeric(1)),
    c(-0.59103110258436842, -1.516553097081865, -1.1590583625324684),
    tolerance = 1e-13
  )
  expect_equal(
    vapply(inits, function(chain_inits) chain_inits[["sigma"]], numeric(1)),
    c(1.5327481321762879, 0.29530333655215796, 1.3229858936207395),
    tolerance = 1e-13
  )
})

test_that("JAGS_get_inits leaves the caller's random-number state unchanged", {

  # The JAGS-free part of the scoped-RNG checks of the seeded public functions;
  # those that fit or read a fit are live tests in test-00-model-fits.R.
  withr::local_preserve_seed()
  withr::defer(RNGkind("default", "default", "default"))
  priors <- list(
    mu = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )

  # A Mersenne-Twister or an L'Ecuyer-CMRG caller state, or no '.Random.seed'.
  .expect_scoped_rng("JAGS_get_inits", function(){
    JAGS_get_inits(priors, chains = 2, seed = 11)
  })

  # The seeded values are those of 'set.seed(seed)', as before the scoping.
  set.seed(11)
  expect_identical(
    vapply(JAGS_get_inits(priors, chains = 2, seed = 11), `[[`, numeric(1), ".RNG.seed"),
    as.numeric(local({ set.seed(11); sample.int(.Machine$integer.max, 2) }))
  )

  # Unseeded calls take their seed from the caller's stream: one draw.
  set.seed(9)
  JAGS_get_inits(priors, chains = 2, seed = NULL)
  after_inits <- .Random.seed
  set.seed(9)
  sample(666666, 1)
  expect_identical(after_inits, .Random.seed)
})

test_that("JAGS restart seeds come from their own stream of the seed", {

  restart_seeds <- function(seed, restarts){
    as.numeric(.JAGS_restart_seeds(seed, restarts))
  }

  # The documented derivation: the first L'Ecuyer-CMRG substream of the seed.
  expected_seeds <- withr::with_preserve_seed({
    set.seed(7, kind = "L'Ecuyer-CMRG")
    assign(
      ".Random.seed",
      parallel::nextRNGStream(get(".Random.seed", envir = globalenv())),
      envir = globalenv()
    )
    sample.int(.Machine$integer.max, 4)
  })
  expect_identical(restart_seeds(7, 4), as.numeric(expected_seeds))

  # 'seed + i' made restart i of seed s the first attempt of seed s + i.
  for(seed in c(1:50, 665, 666666)){
    seeds <- restart_seeds(seed, 10)
    expect_false(any(seeds == seed + seq_len(10)), info = paste("seed", seed))
    expect_length(intersect(seeds, restart_seeds(seed + 1, 10)), 0)
    expect_length(intersect(seeds, .JAGS_chain_seeds(seed, 10)), 0)
    expect_length(unique(seeds), 10)
  }

  # A restart's seed does not depend on the number of restarts.
  for(seed in c(1, 42, 666)){
    all_seeds <- restart_seeds(seed, 10)
    for(restarts in 1:9){
      expect_identical(restart_seeds(seed, restarts), all_seeds[seq_len(restarts)])
    }
  }

  # The caller's RNG kind and state are restored.
  set.seed(3)
  state <- get(".Random.seed", envir = globalenv())
  kind  <- RNGkind()
  restart_seeds(3, 5)
  expect_identical(get(".Random.seed", envir = globalenv()), state)
  expect_identical(RNGkind(), kind)

  # Without a caller '.Random.seed', none is left behind and the kind is kept
  # (removing '.Random.seed' alone would keep the L'Ecuyer-CMRG kind).
  withr::with_preserve_seed({
    RNGkind("Mersenne-Twister")
    rm(".Random.seed", envir = globalenv())
    restart_seeds(3, 5)
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    expect_identical(RNGkind()[[1]], "Mersenne-Twister")
  })
})

test_that("JAGS_fit restarts change the chain seeds and retain generated observations", {

  recorded_inits <- list()
  recorded_data <- list()
  recovered_fit <- structure(
    list(mcmc = list(matrix(0, nrow = 2, ncol = 1,
                            dimnames = list(NULL, "mu")))),
    class = "BayesTools_fit"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      recorded_inits[[length(recorded_inits) + 1L]] <<- list(...)[["inits"]]
      recorded_data[[length(recorded_data) + 1L]] <<- list(...)[["data"]]$y
      if(length(recorded_inits) < 3L){
        stop("backend exploded")
      }
      recovered_fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )
  prior_list <- list(mu = prior("normal", list(0, 1)))

  JAGS_fit(
    model_syntax = "data { y ~ dnorm(0, 1) } model { observed <- y }",
    prior_list = prior_list,
    chains = 2,
    adapt = 50,
    burnin = 50,
    sample = 100,
    autofit_control = list(
      max_Rhat = NULL,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL,
      max_time = list(time = 60, unit = "secs"),
      sample_extend = 1,
      restarts = 3,
      max_extend = 1,
      check_indicators = FALSE
    ),
    seed = 5
  )

  restart_seeds <- .JAGS_restart_seeds(5, 2)
  expect_length(recorded_inits, 3)
  expect_length(recorded_data, 3)
  expect_true(is.finite(recorded_data[[1L]]))
  expect_identical(recorded_data, rep(recorded_data[1L], 3L))
  expect_identical(recorded_inits[[1]], JAGS_get_inits(prior_list, chains = 2, seed = 5))
  for(i in 1:2){
    expect_identical(
      recorded_inits[[i + 1]],
      JAGS_get_inits(prior_list, chains = 2, seed = restart_seeds[[i]])
    )
    expect_false(identical(
      recorded_inits[[i + 1]],
      JAGS_get_inits(prior_list, chains = 2, seed = 5 + i)
    ))
  }
})

test_that("JAGS_extend validates runtime controls before extension", {
  fit <- contract_test_backend_fit()
  invalid <- list(
    parallel = NA,
    cores = Inf,
    silent = NA,
    adapt = 49
  )

  for(name in names(invalid)){
    arguments <- list(fit = fit)
    arguments[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_extend, arguments),
      paste0("'", name, "'")
    )
  }
})

.jags_extend_test_fit <- function(){

  fit <- contract_test_backend_fit()
  attr(fit, "prior_list") <- list()
  attr(fit, "model_syntax") <- "model{}"
  attr(fit, "required_packages") <- character()
  attr(fit, "jags_modules") <- character()
  attr(fit, "add_parameters") <- character()
  attr(fit, "adaptation") <- .bt_adaptation_record(500L)
  fit
}

.jags_extend_test_control <- function(max_time = list(time = 60, unit = "secs")){
  list(
    max_Rhat = NULL,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    max_time = max_time,
    sample_extend = 1,
    restarts = 1,
    max_extend = 1,
    check_indicators = FALSE
  )
}

test_that("JAGS_extend rejects stale fitted metadata before backend work", {

  stale_contract <- .bt_attach_fit_contract(.jags_extend_test_fit())
  attr(stale_contract, "fit_contract")$formula_design_version <- 2L
  package_calls <- 0L
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...){
      package_calls <<- package_calls + 1L
    },
    .package = "BayesTools"
  )

  expect_error(
    JAGS_extend(
      stale_contract,
      autofit_control = .jags_extend_test_control()
    ),
    "missing or unsupported 'formula_design' metadata",
    fixed = TRUE
  )
  expect_identical(package_calls, 0L)

  formula_result <- JAGS_formula(
    ~ 1,
    "mu",
    data.frame(row = 1:2),
    list(intercept = prior("normal", list(0, 1)))
  )
  stale_design <- formula_result$formula_design
  stale_design$schema_version <- 2L
  actual_mismatch <- .jags_extend_test_fit()
  attr(actual_mismatch, "formula_design") <- list(mu = stale_design)
  actual_mismatch <- .bt_attach_fit_contract(actual_mismatch)

  expect_error(
    JAGS_extend(
      actual_mismatch,
      autofit_control = .jags_extend_test_control()
    ),
    "cannot replay this fitted formula",
    fixed = TRUE
  )
  expect_identical(package_calls, 0L)
})

test_that("arguments added since 0.3.0 follow the released arguments", {

  released <- list(
    JAGS_fit = c(
      "model_syntax", "data", "prior_list", "formula_list", "formula_data_list",
      "formula_prior_list", "formula_scale_list", "chains", "adapt", "burnin",
      "sample", "thin", "autofit", "autofit_control", "parallel", "cores",
      "silent", "seed", "add_parameters", "required_packages"
    ),
    JAGS_extend = c("fit", "autofit_control", "parallel", "cores", "silent"),
    JAGS_bridgesampling = c(
      "fit", "log_posterior", "data", "prior_list", "formula_list",
      "formula_data_list", "formula_prior_list", "formula_scale_list",
      "add_parameters", "add_bounds", "maxiter", "silent"
    ),
    # 'transform_orthonormal' (after 'transform_factors') was removed
    runjags_estimates_table = c(
      "fit", "transformations", "title", "footnotes", "warnings", "conditional",
      "probs", "remove_spike_0", "transform_factors", "formula_prefix",
      "remove_inclusion", "remove_parameters", "remove_formulas",
      "keep_parameters", "keep_formulas", "return_samples", "transform_scaled",
      "remove_diagnostics", "diagnostic_columns"
    )
  )
  released$JAGS_estimates_table <- released$runjags_estimates_table

  for(name in names(released)){
    formal_names <- names(formals(get(name, envir = asNamespace("BayesTools"))))
    expect_identical(formal_names[seq_along(released[[name]])], released[[name]], info = name)
  }
  for(name in c("JAGS_fit", "JAGS_bridgesampling")){
    expect_identical(utils::tail(names(formals(get(name, envir = asNamespace("BayesTools")))), 1L), "...", info = name)
  }
})

test_that("JAGS_extend has no seed argument; the chains continue their stored states", {

  # 0.3.0 only set R's seed here and never changed the draws, so the argument
  # was removed (the extension continues the stored generator states)
  expect_false("seed" %in% names(formals(JAGS_extend)))
  expect_identical(
    names(formals(JAGS_extend)),
    c("fit", "autofit_control", "parallel", "cores", "silent", "runtime_setup", "runtime_cache", "worker_output", "adapt", "runtime_report")
  )
  expect_error(
    JAGS_extend(.jags_extend_test_fit(), autofit_control = .jags_extend_test_control(), seed = 1),
    "unused argument (seed = 1)",
    fixed = TRUE
  )
})

test_that("JAGS_extend preserves the last valid fit after a backend error", {

  fit <- .jags_extend_test_fit()
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(...){
      stop("backend exploded")
    },
    .package = "BayesTools"
  )

  expect_false("seed" %in% names(formals(JAGS_extend)))
  expect_warning(
    result <- JAGS_extend(
      fit,
      autofit_control = .jags_extend_test_control()
    ),
    "returning the last valid fit.*backend exploded"
  )
  expected_fit <- fit
  attr(expected_fit, "warnings") <-
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  expect_identical(result, expected_fit)
  expect_identical(
    attr(result, "warnings"),
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  )
})

test_that("JAGS_fit autofit preserves the last valid fit after a backend error", {

  initial_fit <- structure(
    list(mcmc = list(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "mu")))),
    class = "BayesTools_fit"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...) initial_fit,
    .bt_rjags_extend = function(...) stop("backend exploded"),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...) FALSE,
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )

  expect_warning(
    result <- JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit = TRUE,
      autofit_control = list(
        max_Rhat = NULL,
        min_ESS = NULL,
        max_error = NULL,
        max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"),
        sample_extend = 1,
        restarts = 1,
        max_extend = 1,
        check_indicators = FALSE
      ),
      silent = TRUE,
      seed = 1
    ),
    "returning the last valid fit.*backend exploded"
  )
  expect_s3_class(result, "BayesTools_fit")
  expect_false(inherits(result, "error"))
  expect_identical(result$mcmc, initial_fit$mcmc)
  expect_identical(
    attr(result, "warnings"),
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  )
})

test_that("fits leave JAGS_fit and JAGS_extend without the compiled rjags model", {

  # The autofit extensions of one call continue it; the returned fit drops it.
  compiled_model <- function(){
    model_data <- stats::runif(1e4)
    structure(list(ptr = function() model_data), class = "jags")
  }
  backend_fit <- function(){
    structure(
      list(
        mcmc = list(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "mu"))),
        live_model = compiled_model()
      ),
      class = "BayesTools_fit"
    )
  }
  # a fitted object (built before the fitted-metadata helpers are mocked)
  # that still carries a compiled model
  extended_input <- .jags_extend_test_fit()
  extended_input$live_model <- backend_fit()$live_model
  received_models <- list()
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...) backend_fit(),
    .bt_rjags_extend = function(fit, ...){
      received_models[[length(received_models) + 1L]] <<- fit["live_model"]
      fit$live_model <- compiled_model()
      fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...) FALSE,
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )

  fit <- JAGS_fit(
    model_syntax = "model{ mu ~ dnorm(0, 1) }",
    prior_list = list(mu = prior("normal", list(0, 1))),
    chains = 1, adapt = 50, burnin = 50, sample = 100,
    autofit = TRUE,
    autofit_control = .jags_extend_test_control(),
    silent = TRUE,
    seed = 1
  )
  # the extension continued the compiled model of the initial run
  expect_length(received_models, 1L)
  expect_s3_class(received_models[[1L]]$live_model, "jags")
  expect_null(fit$live_model)

  received_models <- list()
  control <- .jags_extend_test_control()
  control$max_extend <- 2
  extended <- JAGS_extend(extended_input, autofit_control = control)
  # the first extension recompiles from the stored states; the second
  # continues the model the first one compiled
  expect_length(received_models, 2L)
  expect_null(received_models[[1L]]$live_model)
  expect_s3_class(received_models[[2L]]$live_model, "jags")
  expect_null(extended$live_model)
})

test_that("JAGS build checks fingerprint loaded R definitions", {

  original <- .JAGS_package_builds("stats")$stats
  sd <- stats::sd
  body(sd) <- quote(stop("modified implementation"))
  testthat::local_mocked_bindings(sd = sd, .package = "stats")
  changed <- .JAGS_package_builds("stats")$stats
  expect_identical(changed$version, original$version)
  expect_identical(changed$dll, original$dll)
  expect_false(identical(changed$r_code, original$r_code))
})

test_that("JAGS build checks distinguish adjacent numeric literals", {

  probe <- function(value){
    fn <- function() NULL
    body(fn) <- value
    environment(fn) <- asNamespace("BayesTools")
    fn
  }
  testthat::local_mocked_bindings(
    .prior_linear_density_default_grid = probe(1), .package = "BayesTools")
  original <- .JAGS_package_builds("BayesTools")$BayesTools
  testthat::local_mocked_bindings(
    .prior_linear_density_default_grid = probe(1 + .Machine$double.eps),
    .package = "BayesTools")
  changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_false(identical(changed$r_code, original$r_code))
  expect_identical(changed[c("version", "dll")], original[c("version", "dll")])

  testthat::local_mocked_bindings(
    .prior_linear_density_default_grid = compiler::cmpfun(probe(1)),
    .package = "BayesTools")
  compiled <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_identical(compiled, original)
})

test_that("JAGS_fit reports and records backend errors recovered by a restart", {

  backend_calls <- 0L
  recovered_fit <- structure(
    list(mcmc = list(matrix(0, nrow = 2, ncol = 1,
                            dimnames = list(NULL, "mu")))),
    class = "BayesTools_fit"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      backend_calls <<- backend_calls + 1L
      if(backend_calls == 1L){
        stop("backend exploded")
      }
      recovered_fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )

  fit_once <- function(silent){

    JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit_control = list(
        max_Rhat = NULL,
        min_ESS = NULL,
        max_error = NULL,
        max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"),
        sample_extend = 1,
        restarts = 2,
        max_extend = 1,
        check_indicators = FALSE
      ),
      silent = silent,
      seed = 1
    )
  }
  expected_message <-
    "JAGS fitting attempt 1 failed and was restarted: backend exploded."

  for(silent in c(FALSE, TRUE)){
    backend_calls <- 0L
    emitted <- character()
    result <- withCallingHandlers(fit_once(silent), warning = function(w){
      expect_identical(backend_calls, 1L)
      expect_null(conditionCall(w))
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

    expect_identical(emitted, if(silent) character() else expected_message)
    expect_identical(backend_calls, 2L)
    expect_s3_class(result, "BayesTools_fit")
    expect_identical(attr(result, "warnings"), expected_message)
  }
})

test_that("the final failed fitting attempt is not announced as a restart", {

  backend_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      backend_calls <<- backend_calls + 1L
      stop("backend exploded")
    },
    .package = "BayesTools"
  )
  inits_drawn <- 0L
  drawn_inits <- BayesTools:::.JAGS_fit_inits
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .JAGS_fit_inits = function(...){
      inits_drawn <<- inits_drawn + 1L
      drawn_inits(...)
    },
    .package = "BayesTools"
  )
  emitted <- character()
  result <- withCallingHandlers(
    JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit_control = list(
        max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"), sample_extend = 1,
        restarts = 3, max_extend = 1, check_indicators = FALSE
      ),
      silent = FALSE,
      seed = 1
    ),
    warning = function(w){
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "error")
  expect_identical(backend_calls, 3L)
  expect_identical(emitted, c(
    "JAGS fitting attempt 1 failed and was restarted: backend exploded.",
    "JAGS fitting attempt 2 failed and was restarted: backend exploded.",
    "The model estimation failed with the following error: backend exploded"
  ))
  # the initial values of the first attempt and of the two restarts only
  expect_identical(inits_drawn, 3L)
})

test_that("JAGS_extend resets its time budget for every call", {

  fit <- .jags_extend_test_fit()
  extension_calls <- 0L
  clock_calls <- 0L
  clock_values <- as.POSIXct(
    c(
      "2026-01-01 00:00:00",
      "2026-01-01 00:00:01",
      "2026-01-01 01:00:00",
      "2026-01-01 01:00:01"
    ),
    tz = "UTC"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      extension_calls <<- extension_calls + 1L
      fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_jags_extend_time = function(){
      clock_calls <<- clock_calls + 1L
      clock_values[[clock_calls]]
    },
    .bt_check_convergence = function(...) TRUE,
    .package = "BayesTools"
  )

  control <- .jags_extend_test_control(
    max_time = list(time = 2, unit = "secs")
  )
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_equal(extension_calls, 2L)
  expect_equal(clock_calls, 4L)
})

.jags_extend_test_fit_with_draws <- function(){

  fit <- .jags_extend_test_fit()
  set.seed(71)
  draws <- cbind(mu = stats::rnorm(50), theta = stats::rnorm(50))
  fit$mcmc <- coda::mcmc.list(coda::mcmc(draws), coda::mcmc(draws + 0.1))
  fit <- add_test_chain_states(fit)
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)))
  attr(fit, "add_parameters") <- "theta"
  attach_test_parameter_map(fit)
}

test_that("JAGS_extend forwards explicit convergence monitor policy", {

  fit <- .jags_extend_test_fit_with_draws()
  convergence_arguments <- NULL
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...){
      convergence_arguments <<- list(...)
      TRUE
    },
    .package = "BayesTools"
  )

  control <- .jags_extend_test_control()
  control$monitor <- "mu"
  control$allow_not_assessable <- TRUE
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_identical(convergence_arguments$monitor, "mu")
  expect_true(convergence_arguments$allow_not_assessable)
})

test_that("JAGS_extend keeps the fit's warnings after successful extensions", {

  fit <- .jags_extend_test_fit()
  previous <- "JAGS fitting attempt 1 failed and was restarted: backend exploded."
  attr(fit, "warnings") <- previous
  # The backend returns new draws without BayesTools attributes.
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      structure(
        list(mcmc = fit$mcmc),
        class = "BayesTools_fit"
      )
    },
    .package = "BayesTools"
  )
  converged <- TRUE
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...) converged,
    .package = "BayesTools"
  )

  extended <- JAGS_extend(fit, autofit_control = .jags_extend_test_control())
  expect_identical(attr(extended, "warnings"), previous)

  converged <- FALSE
  control <- .jags_extend_test_control()
  control$max_extend <- 2
  stopped <- JAGS_extend(fit, autofit_control = control)
  expect_identical(
    attr(stopped, "warnings"),
    c(
      previous,
      "The automatic model fitting was terminated due to the 'max_extend' constraint."
    )
  )
})

test_that("JAGS_extend resolves explicit convergence monitors before extending", {

  fit <- .jags_extend_test_fit_with_draws()
  extension_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      extension_calls <<- extension_calls + 1L
      fit
    },
    .package = "BayesTools"
  )

  control <- .jags_extend_test_control()
  control$monitor <- "theta[2]"
  expect_error(
    JAGS_extend(fit, autofit_control = control),
    "The requested convergence monitor 'theta[2]' is not available in the fitted model.",
    fixed = TRUE
  )
  expect_identical(extension_calls, 0L)

  # An additional monitor named explicitly is checked after the extension.
  control$monitor <- "theta"
  control$min_ESS <- 1
  extended <- JAGS_extend(fit, autofit_control = control)
  expect_identical(extension_calls, 1L)
  expect_s3_class(extended, "BayesTools_fit")
})

test_that("a random-block SD source that prior_list samples is monitored once", {

  # F318 (4): the block-level random_sd_source("tau") made the formula monitor
  # 'tau' through add_parameters, which JAGS_fit() refused as a repetition of
  # the 'tau' prior of prior_list
  captured <- NULL
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      captured <<- list(...)
      stop("mocked backend", call. = FALSE)
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .package = "BayesTools"
  )
  data <- data.frame(g = factor(rep(letters[1:4], each = 3L)))
  fit_with <- function(syntax, prior_list, add_parameters = NULL){
    captured <<- NULL
    JAGS_fit(
      model_syntax = syntax,
      data = list(y = rep(0, 12L)),
      prior_list = prior_list,
      formula_list = list(mu = ~ 1 + random(1 | g, name = "g", covariance = "diag")),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
      formula_random_prior_list = list(mu = prior_random(
        g = random_block(sd_source = random_sd_source("tau"))
      )),
      add_parameters = add_parameters,
      chains = 1, adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE
    )
  }
  syntax <- "model{
  for(i in 1:N_mu){
    y[i] ~ dnorm(mu[i], 1)
  }
}"
  tau_prior <- list(tau = prior("normal", list(0, 1), list(0, Inf)))

  # the source sampled through prior_list is one parameter, monitored through its prior
  expect_no_error(fit_with(syntax, tau_prior))
  expect_identical(sum(captured$monitor == "tau"), 1L)
  expect_true(all(c("mu_intercept", "tau", "mu__xREx__g_xRE_Zx") %in% captured$monitor))

  # a source that the model syntax defines is monitored by the formula
  expect_no_error(fit_with(
    "model{
  tau ~ dnorm(0, 1) T(0, )
  for(i in 1:N_mu){
    y[i] ~ dnorm(mu[i], 1)
  }
}",
    NULL
  ))
  expect_identical(sum(captured$monitor == "tau"), 1L)

  # the user's own repetition of a prior is still refused
  expect_error(
    fit_with(syntax, tau_prior, add_parameters = "tau"),
    "The 'add_parameters' argument must not include parameters already monitored through 'prior_list': tau.",
    fixed = TRUE
  )
})

test_that("JAGS_fit rejects unknown convergence monitors before sampling", {

  backend_calls <- 0L
  set.seed(72)
  draws <- cbind(mu = stats::rnorm(50), theta = stats::rnorm(50))
  sampled_fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(draws), coda::mcmc(draws + 0.1))
    ),
    class = "BayesTools_fit"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      backend_calls <<- backend_calls + 1L
      sampled_fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )
  fit_with_monitor <- function(monitor){
    JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1)\n theta ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      add_parameters = "theta",
      chains = 2,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit = TRUE,
      autofit_control = list(
        max_Rhat = NULL,
        min_ESS = 1,
        max_error = NULL,
        max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"),
        sample_extend = 1,
        restarts = 1,
        max_extend = 1,
        check_indicators = FALSE,
        monitor = monitor
      ),
      silent = TRUE,
      seed = 1
    )
  }

  expect_error(
    fit_with_monitor("thetaa"),
    "The requested convergence monitor 'thetaa' is not monitored by the model.",
    fixed = TRUE
  )
  expect_identical(backend_calls, 0L)

  fit <- fit_with_monitor("theta")
  expect_identical(backend_calls, 1L)
  expect_s3_class(fit, "BayesTools_fit")
  expect_null(attr(fit, "warnings"))
})

test_that("autofit treats deterministic user monitors of constants as structural", {

  # 'theta' is a deterministic model constant; 'eta' is sampled. Both are
  # user monitors, as in BayesTools 0.3.0. Their convergence roles come from
  # the model syntax and are stored in the parameter map.
  sampled_fit <- function(eta_constant){
    set.seed(76)
    chain <- function(){
      draws <- cbind(mu = stats::rnorm(200), theta = 2, eta = stats::rnorm(200))
      if(eta_constant){
        draws[, "eta"] <- .5
      }
      coda::mcmc(draws)
    }
    structure(
      list(mcmc = coda::mcmc.list(chain(), chain())),
      class = "BayesTools_fit"
    )
  }
  control <- list(
    max_Rhat = 1.2, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    max_time = list(time = 60, unit = "secs"), sample_extend = 1,
    restarts = 1, max_extend = 1, check_indicators = FALSE
  )
  prior_list <- list(mu = prior("normal", list(0, 1)))
  max_extend_warning <-
    "The automatic model fitting was terminated due to the 'max_extend' constraint."
  extension_calls <- 0L

  fit_autofit <- function(eta_constant){
    sampled <- add_test_chain_states(sampled_fit(eta_constant))
    testthat::local_mocked_bindings(
      .bt_rjags_fit = function(model, data, ...){
        sampled$model <- model
        sampled["data"] <- list(data)
        sampled
      },
      .bt_rjags_extend = function(...){
        extension_calls <<- extension_calls + 1L
        sampled
      },
      .package = "BayesTools"
    )
    testthat::local_mocked_bindings(
      .JAGS_require_packages = function(...) invisible(NULL),
      .JAGS_load_modules = function(...) invisible(NULL),
      .package = "BayesTools"
    )
    extension_calls <<- 0L
    JAGS_fit(
      model_syntax = "model{\n theta <- 2\n eta ~ dnorm(0, 1)\n}",
      prior_list = prior_list,
      add_parameters = c("theta", "eta"),
      chains = 2, adapt = 50, burnin = 50, sample = 100,
      autofit = TRUE, autofit_control = control, silent = TRUE, seed = 1
    )
  }
  extend_autofit <- function(fit, eta_constant){
    sampled <- add_test_chain_states(sampled_fit(eta_constant))
    testthat::local_mocked_bindings(
      .bt_rjags_extend = function(...) sampled,
      .package = "BayesTools"
    )
    testthat::local_mocked_bindings(
      .JAGS_require_packages = function(...) invisible(NULL),
      .JAGS_load_modules = function(...) invisible(NULL),
      .package = "BayesTools"
    )
    attr(JAGS_extend(fit, autofit_control = control), "warnings")
  }

  # The constant deterministic 'theta' does not block convergence and is
  # reported as structural by the post-fit check.
  fit <- fit_autofit(eta_constant = FALSE)
  expect_identical(extension_calls, 0L)
  expect_null(attr(fit, "warnings"))
  expect_identical(
    parameter_coordinates(fit)$convergence_role,
    c("sampled", "structural", "sampled")
  )
  check <- JAGS_check_convergence(
    fit, prior_list, max_Rhat = 1.2, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL
  )
  expect_true(check)
  expect_equal(
    stats::setNames(attr(check, "diagnostics")$state, attr(check, "diagnostics")$parameter),
    c(mu = "assessable", theta = "structural_constant", eta = "assessable")
  )
  expect_null(extend_autofit(fit, eta_constant = FALSE))

  # A constant sampled 'eta' is a stuck sampler and keeps autofit running.
  fit <- fit_autofit(eta_constant = TRUE)
  expect_identical(extension_calls, 1L)
  expect_identical(attr(fit, "warnings"), max_extend_warning)
  expect_identical(extend_autofit(fit, eta_constant = TRUE), max_extend_warning)
})

.convergence_state_of <- function(result, parameter){

  diagnostics <- attr(result, "diagnostics")
  diagnostics$state[diagnostics$parameter == parameter]
}

test_that("autofit checks user and sampled generated monitors but not derived ones", {

  stem <- "mu__xREx__g_xRE_"
  columns <- c(
    "mu_intercept", "theta", "mu__xREx__g_intercept", "mu__xREx__g_x",
    paste0(stem, "CORx_L[", c("1,1", "1,2", "2,1", "2,2"), "]"),
    paste0(stem, "CORx_R[", c("1,1", "1,2", "2,1", "2,2"), "]"),
    paste0(stem, "CORx_lkj_u[1]"),
    paste0(stem, "Zx[1,", 1:2, "]")
  )
  # Two chains; 'stuck' columns differ by 10 SD between chains (R-hat >> 1.2).
  # Correlation diagonals and the upper Cholesky entry are exact constants.
  sampled_fit <- function(stuck){
    set.seed(74)
    chain <- function(shift){
      draws <- matrix(
        stats::rnorm(200 * length(columns)), ncol = length(columns),
        dimnames = list(NULL, columns)
      )
      draws[, grepl("_L\\[1,1\\]|_R\\[(1,1|2,2)\\]", columns)] <- 1
      draws[, grepl("_L\\[1,2\\]", columns)] <- 0
      sd_columns <- c("mu__xREx__g_intercept", "mu__xREx__g_x")
      draws[, sd_columns] <- abs(draws[, sd_columns])
      draws[, paste0(stem, "CORx_lkj_u[1]")] <- stats::pnorm(
        draws[, paste0(stem, "CORx_lkj_u[1]")]
      )
      draws[, stuck] <- draws[, stuck] + shift
      coda::mcmc(draws)
    }
    structure(
      list(
        mcmc = coda::mcmc.list(chain(0), chain(10))
      ),
      class = "BayesTools_fit"
    )
  }
  control <- list(
    max_Rhat = 1.2, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    max_time = list(time = 60, unit = "secs"), sample_extend = 1,
    restarts = 1, max_extend = 1, check_indicators = FALSE
  )
  set.seed(75)
  formula_data <- data.frame(
    g = factor(rep(paste0("g", 1:3), each = 4)),
    x = stats::rnorm(12)
  )
  max_extend_warning <-
    "The automatic model fitting was terminated due to the 'max_extend' constraint."
  extension_calls <- 0L

  fit_autofit <- function(stuck){
    sampled <- add_test_chain_states(sampled_fit(stuck))
    testthat::local_mocked_bindings(
      .bt_rjags_fit = function(model, data, ...){
        sampled$model <- model
        sampled["data"] <- list(data)
        sampled
      },
      .bt_rjags_extend = function(...){
        extension_calls <<- extension_calls + 1L
        sampled
      },
      .package = "BayesTools"
    )
    testthat::local_mocked_bindings(
      .JAGS_require_packages = function(...) invisible(NULL),
      .JAGS_load_modules = function(...) invisible(NULL),
      .package = "BayesTools"
    )
    extension_calls <<- 0L
    JAGS_fit(
      model_syntax = "model{\n theta <- 2 * mu_intercept\n}",
      formula_list = list(
        mu = ~ 1 + random(1 + x | g, name = "g", covariance = "us")
      ),
      formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(g = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf))
      ))),
      add_parameters = "theta",
      chains = 2, adapt = 50, burnin = 50, sample = 100,
      autofit = TRUE, autofit_control = control, silent = TRUE, seed = 1
    )
  }

  # A stuck derived Cholesky entry does not block convergence, and the
  # default post-fit check (as run by RoBMA) uses the same selection.
  fit <- fit_autofit(paste0(stem, "CORx_L[2,1]"))
  expect_identical(extension_calls, 0L)
  expect_null(attr(fit, "warnings"))
  default_check <- JAGS_check_convergence(
    fit, max_Rhat = 1.2, min_ESS = NULL, max_error = NULL, max_SD_error = NULL
  )
  expect_true(default_check)
  expect_identical(
    unname(.convergence_state_of(default_check, paste0(stem, "CORx_L[2,1]"))),
    "not_requested"
  )
  requested <- JAGS_check_convergence(
    fit, max_Rhat = 1.2, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    monitor = paste0(stem, "CORx_L")
  )
  expect_false(requested)
  expect_match(attr(requested, "errors"), "CORx_L[2,1]", fixed = TRUE)

  # A stuck user node, LKJ primitive, or latent effect is checked, as the
  # user's add_parameters were on BayesTools 0.3.0.
  for(stuck in c("theta", paste0(stem, "CORx_lkj_u[1]"), paste0(stem, "Zx[1,2]"))){
    fit <- fit_autofit(stuck)
    expect_identical(extension_calls, 1L, info = stuck)
    expect_identical(attr(fit, "warnings"), max_extend_warning, info = stuck)
  }

  # JAGS_extend applies the same selection to the stored fit.
  extend_autofit <- function(fit, stuck){
    sampled <- add_test_chain_states(sampled_fit(stuck))
    testthat::local_mocked_bindings(
      .bt_rjags_extend = function(...) sampled,
      .package = "BayesTools"
    )
    testthat::local_mocked_bindings(
      .JAGS_require_packages = function(...) invisible(NULL),
      .JAGS_load_modules = function(...) invisible(NULL),
      .package = "BayesTools"
    )
    attr(JAGS_extend(fit, autofit_control = control), "warnings")
  }
  fit <- fit_autofit(paste0(stem, "CORx_L[2,1]"))
  expect_null(extend_autofit(fit, paste0(stem, "CORx_L[2,1]")))
  for(stuck in c("theta", paste0(stem, "CORx_lkj_u[1]"), paste0(stem, "Zx[1,2]"))){
    expect_identical(
      extend_autofit(fit, stuck), max_extend_warning, info = stuck
    )
  }
})


test_that("loaded numerical constants and immutable attributes enter build parity", {

  original <- .JAGS_package_builds("BayesTools")$BayesTools
  nested <- structure(list(orders = list(c(3L, 7L), c(15L, 31L))),
    source = list(kind = "quadrature", version = 1L))
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_identical(changed$version, original$version)
  expect_identical(changed$dll, original$dll)
  expect_false(identical(changed$r_code, original$r_code))
  nested$orders[[2L]][2L] <- 63L
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  value_changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_false(identical(value_changed$r_code, changed$r_code))
  attr(nested, "source")$version <- 2L
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  attribute_changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_false(identical(attribute_changed$r_code, value_changed$r_code))

  cache <- new.env(parent = emptyenv())
  cache$value <- 1
  mutable <- structure(list(orders = c(3L, 7L)), metadata = list(cache = cache))
  testthat::local_mocked_bindings(.bt_parameter_map_version = mutable, .package = "BayesTools")
  first <- .JAGS_package_builds("BayesTools")$BayesTools
  cache$value <- 2
  second <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_identical(second, first)
})


test_that("post-fit cleanup errors identify the actual parallel operation", {

  testthat::local_mocked_bindings(
    stopCluster = function(cl) stop("socket is unavailable"),
    .package = "parallel"
  )
  expect_warning(.JAGS_finish_runtime_setup(
    function(context) stop("finish must not run"), chains = 2L, cl = list(),
    operation = "Parallel zplot density computation"),
    paste0("Parallel zplot density computation worker cleanup failed: socket is unavailable. ",
      "The runtime finish callback was not run."), fixed = TRUE)
})

test_that("seeds outside the integer range or missing are refused instead of coerced to NA", {

  limit <- .Machine$integer.max
  settings <- function(seed) JAGS_check_and_list_fit_settings(
    chains = 4, adapt = 50, burnin = 50, sample = 100, thin = 1, autofit = FALSE,
    parallel = FALSE, cores = 1, silent = TRUE, seed = seed
  )
  expect_silent(settings(limit))
  expect_silent(settings(-limit))
  for(seed in list(limit + 1, -limit - 1, 2^40)){
    expect_error(settings(seed), "must be equal or lower than|must be equal or higher than", info = seed)
    expect_error(
      JAGS_get_inits(list(mu = prior("normal", list(0, 1))), chains = 1, seed = seed),
      "must be equal or lower than|must be equal or higher than", info = seed
    )
  }
  for(seed in list(NA_real_, Inf, -Inf)){
    expect_error(
      JAGS_get_inits(list(mu = prior("normal", list(0, 1))), chains = 1, seed = seed),
      "The 'seed' ", fixed = TRUE, info = seed
    )
  }
  expect_length(JAGS_get_inits(list(mu = prior("normal", list(0, 1))), chains = 2, seed = limit), 2L)
})


# The extension arguments that reach BayesTools:::.bt_rjags_extend(), one list per call.
.jags_extension_recorder <- function(live_model = FALSE, warning_message = NULL){

  calls <- list()
  extend <- function(fit, ...){
    calls[[length(calls) + 1L]] <<- c(list(fit = fit), list(...))
    if(!is.null(warning_message)){
      warning(warningCondition(warning_message, class = "BayesTools_adaptation_incomplete"))
    }
    fit$live_model <- if(live_model) structure(list(), class = "jags")
    fit
  }
  list(extend = extend, calls = function() calls)
}

.jags_adaptation_test_record <- function(initial, adapt = integer(), rebuilt = logical()){
  list(
    initial = as.integer(initial),
    extensions = data.frame(adapt = as.integer(adapt), rebuilt = rebuilt)
  )
}

test_that("JAGS_extend inherits the initial adaptation and applies an override to its call only", {

  fit <- .jags_extend_test_fit()
  attr(fit, "adaptation") <- .jags_adaptation_test_record(5000)
  recorder <- .jags_extension_recorder()
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  testthat::local_mocked_bindings(.bt_check_convergence = function(...) TRUE, .package = "BayesTools")
  control <- .jags_extend_test_control()

  omitted  <- JAGS_extend(fit, autofit_control = control)
  explicit <- JAGS_extend(fit, autofit_control = control, adapt = 5000)
  override <- JAGS_extend(fit, autofit_control = control, adapt = 20000)
  after    <- JAGS_extend(override, autofit_control = control)

  # what reached the driver: the adaptation, no burn-in, the stored chain states
  calls <- recorder$calls()
  expect_identical(lapply(calls, `[[`, "adapt"), list(5000L, 5000L, 20000L, 5000L))
  expect_true(all(vapply(calls, function(call) !"burnin" %in% names(call), logical(1))))
  expect_identical(lapply(calls, `[[`, "sample"), rep(list(1), 4))
  expect_true(all(vapply(calls, function(call) !"inits" %in% names(call), logical(1))))
  expect_identical(calls[[1L]]$fit$end_state, fit$end_state)

  # an explicit value equal to the initial one is the omitted value
  expect_identical(omitted, explicit)
  expect_identical(attr(omitted, "adaptation"), .jags_adaptation_test_record(5000, 5000, TRUE))
  # the override is recorded for its extension and does not replace the
  # initial setting, which the next extension without 'adapt' uses again
  expect_identical(attr(override, "adaptation"), .jags_adaptation_test_record(5000, 20000, TRUE))
  expect_identical(
    attr(after, "adaptation"),
    .jags_adaptation_test_record(5000, c(20000, 5000), c(TRUE, TRUE))
  )
  expect_identical(attr(fit, "adaptation"), .jags_adaptation_test_record(5000))
})

test_that("JAGS_extend adapts only when the backend rebuilds the sampler", {

  fit <- .jags_extend_test_fit()
  attr(fit, "adaptation") <- .jags_adaptation_test_record(800)
  control <- .jags_extend_test_control()
  control$max_extend <- 2
  testthat::local_mocked_bindings(.bt_check_convergence = function(...) FALSE, .package = "BayesTools")

  # serially, the second automatic extension continues the compiled model of
  # the first; the returned fit keeps no compiled model
  recorder <- .jags_extension_recorder(live_model = TRUE)
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  serial <- suppressWarnings(JAGS_extend(fit, autofit_control = control, adapt = 900))
  calls <- recorder$calls()
  expect_identical(lapply(calls, `[[`, "adapt"), list(900L, 900L))
  expect_null(calls[[1L]]$fit$live_model)
  expect_s3_class(calls[[2L]]$fit$live_model, "jags")
  expect_identical(
    attr(serial, "adaptation"),
    .jags_adaptation_test_record(800, c(900, 900), c(TRUE, FALSE))
  )
  expect_null(serial$live_model)

  # parallel workers rebuild the sampler for every extension
  recorder <- .jags_extension_recorder(live_model = TRUE)
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  testthat::local_mocked_bindings(
    .JAGS_runtime_cluster_start = function(...) list("worker-1", "worker-2"),
    .JAGS_runtime_cluster_stop = function(cl) invisible(NULL),
    .package = "BayesTools"
  )
  parallel <- suppressWarnings(JAGS_extend(fit, autofit_control = control, parallel = TRUE))
  calls <- recorder$calls()
  expect_identical(vapply(calls, function(call) length(call$cl), integer(1)), c(2L, 2L))
  expect_identical(lapply(calls, `[[`, "adapt"), list(800L, 800L))
  expect_true(all(vapply(calls, function(call) !"burnin" %in% names(call), logical(1))))
  expect_identical(
    attr(parallel, "adaptation"),
    .jags_adaptation_test_record(800, c(800, 800), c(TRUE, TRUE))
  )
})

test_that("JAGS_extend validates 'adapt' as JAGS_fit does", {

  fit <- .jags_extend_test_fit()
  extension_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      extension_calls <<- extension_calls + 1L
      fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(.bt_check_convergence = function(...) TRUE, .package = "BayesTools")
  fit_settings <- function(adapt) JAGS_check_and_list_fit_settings(
    chains = 1, adapt = adapt, burnin = 50, sample = 100, thin = 1,
    autofit = FALSE, parallel = FALSE, cores = 1, silent = TRUE, seed = 1
  )

  invalid <- list(49, 0, -1, 50.5, NA, NA_real_, Inf, "100", c(100, 200), TRUE,
    numeric(0), integer(0), character(0), .Machine$integer.max + 1)
  for(adapt in invalid){
    initial_error <- tryCatch(fit_settings(adapt), error = conditionMessage)
    extend_error  <- tryCatch(
      JAGS_extend(fit, autofit_control = .jags_extend_test_control(), adapt = adapt),
      error = conditionMessage
    )
    expect_true(is.character(initial_error), info = deparse(adapt))
    expect_identical(extend_error, initial_error, info = deparse(adapt))
  }
  expect_identical(extension_calls, 0L)

  expect_identical(fit_settings(50)$adapt, 50)
  extended <- JAGS_extend(fit, autofit_control = .jags_extend_test_control(), adapt = 50)
  expect_identical(attr(extended, "adaptation"), .jags_adaptation_test_record(500, 50, TRUE))
})

test_that("JAGS_extend requires a current fit, chain states and adaptation record", {

  fit <- .jags_extend_test_fit()
  legacy <- missing_adaptation <- malformed <- missing_state <- fit
  class(legacy) <- c("runjags", "BayesTools_fit")
  attr(missing_adaptation, "adaptation") <- NULL
  attr(malformed, "adaptation") <- list(initial = "500")
  missing_state$end_state <- NULL
  for (object in list(legacy, missing_adaptation, malformed, missing_state)) {
    expect_error(JAGS_extend(object, autofit_control = .jags_extend_test_control()),
                 class = "BayesTools_refit_required")
  }
})

test_that("JAGS_extend issues and keeps incomplete-adaptation warnings", {

  fit <- .jags_extend_test_fit()
  message <- "The JAGS adaptation phase was incomplete after 500 iterations."
  recorder <- .jags_extension_recorder(warning_message = message)
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  testthat::local_mocked_bindings(.bt_check_convergence = function(...) TRUE, .package = "BayesTools")

  expect_warning(
    extended <- JAGS_extend(fit, autofit_control = .jags_extend_test_control()),
    class = "BayesTools_adaptation_incomplete"
  )
  expect_identical(attr(extended, "warnings"), message)
})

test_that("JAGS_fit forwards its adaptation to automatic extensions", {

  backend_fit <- function(...){
    structure(
      list(
        mcmc = list(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "mu"))),
        live_model = structure(list(), class = "jags")
      ),
      class = "BayesTools_fit"
    )
  }
  run_calls <- list()
  recorder <- .jags_extension_recorder(live_model = TRUE)
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      run_calls[[length(run_calls) + 1L]] <<- list(...)
      backend_fit()
    },
    .bt_rjags_extend = recorder$extend,
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...) FALSE,
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .JAGS_runtime_cluster_start = function(...) list("worker-1", "worker-2"),
    .JAGS_runtime_cluster_stop = function(cl) invisible(NULL),
    .JASP_progress_bar_start = function(...) invisible(NULL),
    .JASP_progress_bar_tick = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )
  control <- .jags_extend_test_control()
  control$max_extend <- 2
  fit_with <- function(...) suppressWarnings(JAGS_fit(
    model_syntax = "model{ mu ~ dnorm(0, 1) }",
    prior_list = list(mu = prior("normal", list(0, 1))),
    chains = 2, adapt = 750, burnin = 50, sample = 100,
    autofit_control = control, silent = TRUE, seed = 1, ...
  ))

  # serial autofit continues the compiled model of the initial run
  serial <- fit_with(autofit = TRUE)
  expect_identical(run_calls[[1L]]$adapt, 750)
  expect_identical(lapply(recorder$calls(), `[[`, "adapt"), list(750L, 750L))
  expect_true(all(vapply(recorder$calls(), function(call) !"burnin" %in% names(call), logical(1))))
  expect_identical(
    attr(serial, "adaptation"),
    .jags_adaptation_test_record(750, c(750, 750), c(FALSE, FALSE))
  )

  # parallel autofit rebuilds the sampler for every extension
  recorder <- .jags_extension_recorder()
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  parallel <- fit_with(autofit = TRUE, parallel = TRUE, cores = 2)
  expect_identical(vapply(recorder$calls(), function(call) length(call$cl), integer(1)), c(2L, 2L))
  expect_identical(lapply(recorder$calls(), `[[`, "adapt"), list(750L, 750L))
  expect_identical(
    attr(parallel, "adaptation"),
    .jags_adaptation_test_record(750, c(750, 750), c(TRUE, TRUE))
  )

  # Progress is delivered by the initial backend run.
  recorder <- .jags_extension_recorder()
  testthat::local_mocked_bindings(.bt_rjags_extend = recorder$extend, .package = "BayesTools")
  jasp <- fit_with(autofit = FALSE, parallel = TRUE, cores = 2, is_JASP = TRUE)
  expect_length(recorder$calls(), 0L)
  expect_identical(run_calls[[3L]]$sample, 100)
  expect_identical(attr(jasp, "adaptation"), .jags_adaptation_test_record(750))

  # a fit without extensions records only the initial setting
  plain <- fit_with(autofit = FALSE)
  expect_identical(attr(plain, "adaptation"), .jags_adaptation_test_record(750))
})
