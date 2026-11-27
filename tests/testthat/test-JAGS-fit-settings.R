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
    monitor = c("mu", "sigma")
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
    monitor = NA_character_
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

.jags_extend_test_fit <- function(chains = 1L){

  fit <- contract_test_backend_fit(chains = chains)
  fit$thin <- 1L
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
    .bt_autofit_clock = function(){
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

.jags_extend_test_fit_with_draws <- function(chains = 2L){

  fit <- .jags_extend_test_fit()
  set.seed(71)
  draws <- cbind(mu = stats::rnorm(50), theta = stats::rnorm(50))
  fit$mcmc <- coda::mcmc.list(lapply(seq_len(chains), function(i) coda::mcmc(draws + 0.1 * (i - 1L))))
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

  # the check of automatic fitting ignores diagnostics that are not
  # assessable: they do not extend a fit
  control <- .jags_extend_test_control()
  control$monitor <- "mu"
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_identical(convergence_arguments$monitor, "mu")
  expect_true(convergence_arguments$allow_not_assessable)
})

test_that("autofit_control no longer accepts allow_not_assessable", {

  # removed after 0.4.0, without a deprecation: automatic fitting never
  # extends a fit for diagnostics that are not assessable, and
  # JAGS_check_convergence() keeps the argument
  expect_false("allow_not_assessable" %in% names(formals(JAGS_fit)$autofit_control))
  expect_false("allow_not_assessable" %in% names(formals(JAGS_extend)$autofit_control))
  expect_identical(formals(JAGS_check_convergence)$allow_not_assessable, FALSE)
  expect_false("allow_not_assessable" %in% names(JAGS_check_and_list_autofit_settings(list(sample_extend = 1000))))

  message <- "The 'allow_not_assessable' objects are not recognized by the 'autofit_control' argument."
  for(allow in c(FALSE, TRUE)){
    expect_error(
      JAGS_check_and_list_autofit_settings(list(allow_not_assessable = allow)),
      message, fixed = TRUE
    )
    expect_error(
      JAGS_fit(
        model_syntax = "model{ mu ~ dnorm(0, 1) }", data = list(),
        prior_list = list(mu = prior("normal", list(0, 1))),
        autofit_control = list(allow_not_assessable = allow)
      ),
      message, fixed = TRUE
    )
    control <- .jags_extend_test_control()
    control$allow_not_assessable <- allow
    expect_error(
      JAGS_extend(.jags_extend_test_fit(), autofit_control = control),
      message, fixed = TRUE
    )
  }
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

  # A constant sampled 'eta' is a stuck sampler: not assessable, which the
  # post-fit check reports, but more draws do not change that, so autofit
  # stops once the assessable parameters converge.
  fit <- fit_autofit(eta_constant = TRUE)
  expect_identical(extension_calls, 0L)
  expect_null(attr(fit, "warnings"))
  check <- JAGS_check_convergence(
    fit, prior_list, max_Rhat = 1.2, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(check)
  expect_identical(attr(check, "errors"), "R-hat diagnostic is not assessable for 'eta'.")
  expect_null(extend_autofit(fit, eta_constant = TRUE))
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

test_that("the incomplete-adaptation warning belongs to the fit diagnostics", {

  testthat::local_mocked_bindings(adapt = function(...) FALSE, .package = "rjags")
  condition <- tryCatch(.bt_rjags_group(
    inits = list(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1L)),
    model = "model{ mu ~ dnorm(0, 1)\n y ~ dnorm(mu, 1) }", data = list(y = 0.5), monitor = "mu",
    adapt = 10L, burnin = 0L, sample = 10L, thin = 1L, silent = TRUE, keep_live = FALSE
  ), warning = identity)
  expect_identical(class(condition)[1:2], c("BayesTools_adaptation_incomplete", "BayesTools_fit_diagnostic"))
  expect_identical(conditionMessage(condition), "The JAGS adaptation phase was incomplete after 10 iterations.")

  # the iterations are written out also where R would print them in
  # scientific notation: 500 draws at thin = 200 are 1e+05 iterations
  condition <- tryCatch(.bt_rjags_group(
    inits = list(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1L)),
    model = "model{ mu ~ dnorm(0, 1)\n y ~ dnorm(mu, 1) }", data = list(y = 0.5), monitor = "mu",
    adapt = 500, burnin = 0, sample = 100, thin = 200, silent = TRUE, keep_live = FALSE
  ), warning = identity)
  expect_identical(conditionMessage(condition), "The JAGS adaptation phase was incomplete after 100000 iterations.")
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
  jasp <- fit_with(autofit = FALSE, parallel = TRUE, cores = 2, is_JASP = TRUE,
    is_JASP_prefix = "Subgroup A")
  expect_length(recorder$calls(), 0L)
  expect_identical(run_calls[[3L]]$sample, 100)
  expect_identical(run_calls[[3L]]$jasp_prefix, "Subgroup A: ")
  expect_false(run_calls[[3L]]$retain_workers)
  expect_identical(attr(jasp, "adaptation"), .jags_adaptation_test_record(750))

  plain <- fit_with(autofit = FALSE)
  expect_identical(attr(plain, "adaptation"), .jags_adaptation_test_record(750))

  # During JASP autofit the coordinator identifies the live worker models
  # with TRUE. All extensions reuse the same worker group.
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(retain_workers, ...){
      expect_true(retain_workers)
      fit <- backend_fit()
      fit$live_model <- TRUE
      fit
    },
    .bt_rjags_extend = function(fit, ...) fit,
    .package = "BayesTools"
  )
  jasp_auto <- fit_with(autofit = TRUE, parallel = TRUE, cores = 2, is_JASP = TRUE)
  expect_identical(attr(jasp_auto, "adaptation"),
    .jags_adaptation_test_record(750, c(750, 750), c(FALSE, FALSE)))

})

# ---------------------------------------------------------------------------- #
# sample_extend = "auto"
# ---------------------------------------------------------------------------- #
test_that("sample_extend accepts a positive whole number or 'auto'", {

  control <- function(sample_extend){
    list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05,
      max_time = list(time = 60, unit = "mins"), sample_extend = sample_extend,
      restarts = 10, max_extend = 10)
  }
  expect_identical(JAGS_check_and_list_autofit_settings(control("auto"))$sample_extend, "auto")
  expect_identical(JAGS_check_and_list_autofit_settings(control(1000))$sample_extend, 1000)
  # a named "auto" is stored as the plain string
  expect_identical(JAGS_check_and_list_autofit_settings(control(c(a = "auto")))$sample_extend, "auto")
  invalid <- list(
    list(value = "Auto",            message = "The 'Auto' values are not recognized by the 'sample_extend' argument."),
    list(value = c("auto", "auto"), message = "The 'sample_extend' argument must have length '1'."),
    list(value = NA_character_,     message = "The 'sample_extend' argument cannot contain NA/NaN values."),
    list(value = 0,                 message = "The 'sample_extend' must be equal or higher than 1."),
    list(value = 1.5,               message = "The 'sample_extend' argument must be an integer vector."),
    list(value = TRUE,              message = "The 'sample_extend' argument must be a numeric vector."),
    list(value = NULL,              message = "The 'sample_extend' argument cannot be NULL.")
  )
  fit <- .jags_extend_test_fit()
  extension_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, ...){
      extension_calls <<- extension_calls + 1L
      fit
    },
    .package = "BayesTools"
  )
  for(case in invalid){
    expect_error(JAGS_check_and_list_autofit_settings(control(case$value)), case$message,
      fixed = TRUE, info = deparse(case$value))
    expect_error(JAGS_extend(fit, autofit_control = control(case$value)), case$message,
      fixed = TRUE, info = deparse(case$value))
  }
  expect_identical(extension_calls, 0L)
  expect_null(JAGS_check_and_list_autofit_settings(control(NULL), skip_sample_extend = TRUE)$sample_extend)
})

test_that("automatic extensions are sized from the failed convergence check", {

  # one row per checked quantity, as .bt_check_convergence() reports them
  diagnostics <- function(state = "assessable", Rhat = 1.01, ESS = 1000,
                          MCMC_error = 0.005, MCMC_SD_error = 0.02){
    data.frame(parameter = paste0("p", seq_along(state)), state = state, Rhat = Rhat,
      ESS = ESS, MCMC_error = MCMC_error, MCMC_SD_error = MCMC_SD_error,
      stringsAsFactors = FALSE)
  }
  size <- function(rows, n = 4000, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01,
                   max_SD_error = 0.05){
    .bt_autofit_extension_size(rows, n, max_Rhat, min_ESS, max_error, max_SD_error)
  }
  # n = 4000 draws per chain: n * (f - 1) * 1.1 rounded up to a multiple of
  # 100, then at most 8000, then at least 1000
  cases <- list(
    # each criterion alone: f = 2 (ESS 250 of 500), 2.25 ((0.015 / 0.01)^2),
    # 1.44 ((0.06 / 0.05)^2: 1936 draws, rounded up), 1.5625 (ESS 320: 2475)
    list(rows = diagnostics(ESS = 250),                               size = 4400),
    list(rows = diagnostics(MCMC_error = 0.015),                      size = 5500),
    list(rows = diagnostics(MCMC_SD_error = 0.06),                    size = 2000),
    list(rows = diagnostics(ESS = 320),                               size = 2500),
    # several at once, in one quantity and across quantities: the largest
    list(rows = diagnostics(ESS = 250, MCMC_error = 0.015, MCMC_SD_error = 0.06), size = 5500),
    list(rows = diagnostics(state = c("assessable", "assessable"), ESS = c(250, 1000),
      MCMC_error = c(0.005, 0.015)),                                   size = 5500),
    # R-hat: at least 2, also when the quantity misses another criterion by
    # less (ESS 495 or 400 of 500), and a larger growth of another criterion
    # (ESS 200: f = 2.5) wins
    list(rows = diagnostics(Rhat = 1.2),                              size = 4400),
    list(rows = diagnostics(Rhat = 1.1, ESS = 495),                   size = 4400),
    list(rows = diagnostics(Rhat = 1.2, ESS = 400),                   size = 4400),
    list(rows = diagnostics(Rhat = 1.2, ESS = 200),                   size = 6600),
    list(rows = diagnostics(state = c("assessable", "assessable"), Rhat = c(1.2, 1.01),
      ESS = c(1000, 400)),                                             size = 4400),
    # caps: a small shortfall takes the minimum, a large one (or ESS 0) 2 n
    list(rows = diagnostics(ESS = 450),                               size = 1000),
    list(rows = diagnostics(ESS = 10),                                size = 8000),
    list(rows = diagnostics(ESS = 0),                                 size = 8000),
    list(rows = diagnostics(MCMC_error = Inf, MCMC_SD_error = Inf),   size = 8000),
    # not assessable: what its computed diagnostics require, the undefined
    # ones add nothing; alone or next to a failing quantity, only the
    # computed diagnostics count (ESS 400 of 500: f = 1.25, 1100)
    list(rows = diagnostics(state = "not_assessable", Rhat = NA, ESS = NA,
      MCMC_error = NA, MCMC_SD_error = NA),                            size = 1000),
    list(rows = diagnostics(state = "not_assessable", Rhat = NA, ESS = 100),
                                                                       size = 8000),
    list(rows = diagnostics(state = c("assessable", "not_assessable"), Rhat = c(1.01, NA),
      ESS = c(400, NA), MCMC_error = c(0.005, NA), MCMC_SD_error = c(0.02, NA)),
                                                                       size = 1100),
    # quantities that were not checked and constants add nothing
    list(rows = diagnostics(state = c("not_checked", "structural_constant", "not_requested"),
      ESS = c(10, NA, 10)),                                             size = 1000),
    list(rows = diagnostics()[0L, ],                                   size = 1000)
  )
  for(i in seq_along(cases)){
    expect_identical(size(cases[[i]]$rows), cases[[i]]$size, info = i)
  }
  # criteria that are not set add nothing; a zero target is met
  expect_identical(size(diagnostics(ESS = 10, Rhat = 2), min_ESS = NULL, max_Rhat = NULL), 1000)
  expect_identical(size(diagnostics(MCMC_error = 1, MCMC_SD_error = 1), max_error = NULL,
    max_SD_error = NULL), 1000)
  expect_identical(size(diagnostics(ESS = 0), min_ESS = 0), 1000)
  # 2 n between multiples of 100 is rounded down (3050 to 3000; 1002 to 1000);
  # below the minimum, the minimum holds
  expect_identical(size(diagnostics(ESS = 10), n = 1525), 3000)
  expect_identical(size(diagnostics(ESS = 10), n = 501), 1000)
  expect_identical(size(diagnostics(ESS = 10), n = 300), 1000)
  expect_identical(size(diagnostics(ESS = 0), n = 300), 1000)
  # below the cap, the size is rounded up before the cap applies: n = 1525
  # with f = 2 needs 1677.5 draws, 1700 (cap 3000)
  expect_identical(size(diagnostics(ESS = 250), n = 1525), 1700)
  # whole growth stays exact: 6000 * 1 * 1.1 = 6600, not 6700
  expect_identical(size(diagnostics(ESS = 500 / 2), n = 6000), 6600)
  # also when f - 1 is inexact: 500 / 450 - 1 exceeds 1 / 9, but 9000 and
  # 45000 draws need exactly 1100 and 5500
  expect_identical(size(diagnostics(ESS = 450), n = 9000), 1100)
  expect_identical(size(diagnostics(ESS = 450), n = 45000), 5500)
  # sizes are whole doubles, also above the integer range (a saturated
  # deficit at 1.1e9 draws per chain), until the callers cap them at the
  # largest multiple of 100 whose draws * thin sampling iterations the backend
  # takes as a signed integer
  expect_identical(size(diagnostics(ESS = 10), n = 1.1e9), 2.2e9)
  expect_identical(.bt_autofit_iteration_cap(1), 2147483600)
  expect_identical(min(size(diagnostics(ESS = 10), n = 1.1e9), .bt_autofit_iteration_cap(1)),
    2147483600)
  # n = 5000 with thin = 250000: a saturated 10000 is capped at 8500
  # (2147483647 / 250000 = 8589.9)
  expect_identical(.bt_autofit_iteration_cap(250000), 8500)
  expect_identical(min(size(diagnostics(ESS = 10), n = 5000), .bt_autofit_iteration_cap(250000)),
    8500)
})

# Mocked extensions append 'sample' rows; the mocked check reports 'mu', whose
# ESS is 4% of the draws per chain and passes from 500, and with
# 'not_assessable' also 'theta', whose diagnostics are not assessable (as for
# a component that no chain visits); with 'assessable = FALSE' only 'theta'.
# With 'theta_ESS', the R-hat of 'theta' is undefined but its ESS is that
# share of the draws per chain. Its result follows the failures of the real
# check under the allow_not_assessable that it is given.
.autofit_size_test_backend <- function(env = parent.frame(), not_assessable = FALSE,
                                       assessable = TRUE, theta_ESS = NA_real_){

  state <- new.env(parent = emptyenv())
  state$samples <- integer()
  state$checks <- list()
  grow <- function(fit, sample){
    fit$mcmc <- coda::mcmc.list(lapply(fit$mcmc, function(chain){
      coda::mcmc(rbind(as.matrix(chain),
        matrix(0, nrow = sample, ncol = ncol(chain), dimnames = list(NULL, colnames(chain)))))
    }))
    fit
  }
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, sample, ...){
      state$samples <- c(state$samples, as.integer(sample))
      grow(fit, sample)
    },
    .bt_check_convergence = function(fit, ..., max_Rhat, fail_fast, allow_not_assessable){
      n <- nrow(fit$mcmc[[1L]])
      state$checks[[length(state$checks) + 1L]] <- list(n = n, fail_fast = fail_fast,
        Rhat = !is.null(max_Rhat), allow_not_assessable = allow_not_assessable)
      rows <- data.frame(parameter = c("mu", "theta"),
        state = c("assessable", "not_assessable"), Rhat = NA_real_,
        ESS = c(0.04 * n, theta_ESS * n), MCMC_error = NA_real_, MCMC_SD_error = NA_real_)
      rows <- rows[c(assessable, not_assessable), , drop = FALSE]
      fails <- .bt_convergence_failures(rows, max_Rhat = NULL, min_ESS = 500,
        max_error = NULL, max_SD_error = NULL, allow_not_assessable = allow_not_assessable)
      structure(length(fails) == 0L, diagnostics = rows, errors = fails)
    },
    .package = "BayesTools",
    .env = env
  )
  state
}

.autofit_size_test_control <- function(sample_extend = "auto", max_extend = 10,
                                       max_time = NULL){
  list(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL,
    max_time = max_time, sample_extend = sample_extend, restarts = 1, max_extend = max_extend)
}

test_that("JAGS_extend sizes each extension from the last check and stops at max_time", {

  # (12500 - n) * 1.1 draws for n draws per chain: from the supplied fit's 2
  # draws the minimum, then 2 n rounded down (2004 to 2000, 6004 to 6000), and
  # 3847.8 rounded up
  backend <- .autofit_size_test_backend()
  fit <- .jags_extend_test_fit()
  extended <- JAGS_extend(fit, autofit_control = .autofit_size_test_control())
  expect_identical(backend$samples, c(1000L, 2000L, 6000L, 3900L))
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(2L, 1002L, 3002L, 9002L, 12902L))
  expect_false(any(vapply(backend$checks, `[[`, logical(1), "fail_fast")))
  # the one-chain fit: the check of the supplied fit, which only sizes the
  # first extension, skips R-hat without the one-chain warning of the others
  expect_identical(vapply(backend$checks, `[[`, logical(1), "Rhat"), c(FALSE, rep(TRUE, 4L)))
  expect_identical(nrow(extended$mcmc[[1L]]), 12902L)
  expect_null(attr(extended, "warnings"))
  # a named "auto" runs as "auto"
  backend <- .autofit_size_test_backend()
  named <- JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    sample_extend = c(step = "auto")))
  expect_identical(backend$samples, c(1000L, 2000L, 6000L, 3900L))
  expect_identical(nrow(named$mcmc[[1L]]), 12902L)
  # at most the draws whose sampling iterations the backend can run, to which
  # the minimum yields: with thin 2.5e6, 800 (2147483647 / 2.5e6 = 859.0);
  # with an adaptation of 50 draws, whose iterations the rebuilt sampler
  # counts with those of the extension ((50 + 800) * 2.5e6 is in range)
  backend <- .autofit_size_test_backend()
  thinned <- fit
  thinned$thin <- 2500000L
  suppressWarnings(JAGS_extend(thinned, autofit_control = .autofit_size_test_control(
    max_extend = 1), adapt = 50))
  expect_identical(backend$samples, 800L)

  # an extension is sized only when it runs: max_extend = 2 sizes two
  size_from <- .bt_autofit_extension_size_from
  sized <- 0L
  testthat::local_mocked_bindings(
    .bt_autofit_extension_size_from = function(...){
      sized <<- sized + 1L
      size_from(...)
    },
    .package = "BayesTools"
  )
  backend <- .autofit_size_test_backend()
  limited <- suppressWarnings(JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    max_extend = 2)))
  expect_identical(backend$samples, c(1000L, 2000L))
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(2L, 1002L, 3002L))
  expect_identical(sized, 2L)
  expect_identical(attr(limited, "warnings"),
    "The automatic model fitting was terminated due to the 'max_extend' constraint.")

  # max_time: an extension is at most the draws per chain that the time left
  # allows at the rate of this call's draws so far, rounded down to a multiple
  # of 100 (also below the minimum of 1000), and a cap below 100 stops the
  # fitting with the max_time warning. The first extension has no rate yet.
  # Clock from 0 with a 25.55 s limit: 1000 draws in 10 s leave 15.55 s for
  # 1555 draws, so 1500 instead of 2000; 2500 draws in 24 s leave 1.55 s for
  # 161 draws, so 100 instead of 5000; 2600 draws in 25.5 s leave 5: stop
  clock <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC") + c(0, 1, 10, 24, 25.5)
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_autofit_clock = function(){
      calls <<- calls + 1L
      clock[[calls]]
    },
    .package = "BayesTools"
  )
  backend <- .autofit_size_test_backend()
  stopped <- JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    max_time = list(time = 25.55, unit = "secs")))
  expect_identical(backend$samples, c(1000L, 1500L, 100L))
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(2L, 1002L, 2502L, 2602L))
  expect_identical(calls, 5L)
  expect_identical(attr(stopped, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")
  # the same in minutes
  calls <- 0L
  backend <- .autofit_size_test_backend()
  stopped <- JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    max_time = list(time = 25.55 / 60, unit = "mins")))
  expect_identical(backend$samples, c(1000L, 1500L, 100L))

  # past the limit at the first extension, the supplied fit is neither
  # checked nor sized
  calls <- 0L
  clock <- clock[[1L]] + c(0, 2.5)
  sized <- 0L
  backend <- .autofit_size_test_backend()
  stopped <- JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    max_time = list(time = 2, unit = "secs")))
  expect_length(backend$samples, 0L)
  expect_length(backend$checks, 0L)
  expect_identical(sized, 0L)
  expect_identical(attr(stopped, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")

  # no time left (max_time = 0) on a clock that has not advanced: no rate is
  # observed, but the cap is 0 and the fitting stops before any extension
  calls <- 0L
  clock <- rep(clock[[1L]], 20L)
  backend <- .autofit_size_test_backend()
  stopped <- JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    max_time = list(time = 0, unit = "secs")))
  expect_length(backend$samples, 0L)
  expect_length(backend$checks, 0L)
  expect_identical(sized, 0L)
  expect_identical(calls, 2L)
  expect_identical(attr(stopped, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")

  # the fixed step is unchanged: no check of the supplied fit, fail_fast checks
  calls <- 0L
  clock <- rep(clock[[1L]], 20L)
  backend <- .autofit_size_test_backend()
  fixed <- suppressWarnings(JAGS_extend(fit, autofit_control = .autofit_size_test_control(
    sample_extend = 1000, max_extend = 3)))
  expect_identical(backend$samples, c(1000L, 1000L, 1000L))
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(1002L, 2002L, 3002L))
  expect_true(all(vapply(backend$checks, `[[`, logical(1), "fail_fast")))
  expect_identical(attr(fixed, "warnings"),
    "The automatic model fitting was terminated due to the 'max_extend' constraint.")
})

test_that("JAGS_extend sizes the first automatic extension of a multi-chain fit with R-hat", {

  # two chains of 1525 draws, an R-hat above max_Rhat and the other criteria
  # met: the check of the supplied fit requests R-hat, so the first extension
  # is 1.1 n = 1677.5 rounded up to 1700 (without R-hat, the minimum of 1000)
  fit <- .jags_extend_test_fit(chains = 2L)
  fit$mcmc <- coda::mcmc.list(lapply(1:2, function(i){
    coda::mcmc(matrix(0, nrow = 1525L, ncol = 1L, dimnames = list(NULL, "x")))
  }))
  fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(fit))
  samples <- integer()
  checks <- list()
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, sample, ...){
      samples <<- c(samples, as.integer(sample))
      fit$mcmc <- coda::mcmc.list(lapply(fit$mcmc, function(chain){
        coda::mcmc(rbind(as.matrix(chain),
          matrix(0, nrow = sample, ncol = ncol(chain), dimnames = list(NULL, colnames(chain)))))
      }))
      fit
    },
    .bt_check_convergence = function(fit, ..., max_Rhat, fail_fast){
      n <- nrow(fit$mcmc[[1L]])
      checks[[length(checks) + 1L]] <<- list(n = n, Rhat = !is.null(max_Rhat))
      rows <- data.frame(parameter = "x", state = "assessable",
        Rhat = if(is.null(max_Rhat)) NA_real_ else 1.2, ESS = 1e6,
        MCMC_error = NA_real_, MCMC_SD_error = NA_real_)
      structure(n > 1525L, diagnostics = rows)
    },
    .package = "BayesTools"
  )
  extended <- JAGS_extend(fit, autofit_control = .autofit_size_test_control())
  expect_identical(samples, 1700L)
  expect_identical(vapply(checks, `[[`, integer(1), "n"), c(1525L, 3225L))
  expect_identical(vapply(checks, `[[`, logical(1), "Rhat"), c(TRUE, TRUE))
  expect_identical(nrow(extended$mcmc[[2L]]), 3225L)
})

test_that("an automatic JAGS_extend gives the one-chain warning as often as the fixed step", {

  # the real convergence check of a one-chain fit warns that R-hat cannot be
  # computed; the check of the supplied fit, which only sizes the first
  # "auto" extension, does not repeat it, so two extensions give two
  # warnings with either step
  fit <- .jags_extend_test_fit_with_draws(chains = 1L)
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(fit, sample, ...){
      fit$mcmc <- coda::mcmc.list(lapply(fit$mcmc, function(chain){
        draws <- matrix(stats::rnorm(sample * ncol(chain)), ncol = ncol(chain),
          dimnames = list(NULL, colnames(chain)))
        coda::mcmc(rbind(as.matrix(chain), draws))
      }))
      fit
    },
    .package = "BayesTools"
  )
  one_chain_warnings <- function(sample_extend){
    count <- 0L
    withCallingHandlers(
      JAGS_extend(fit, autofit_control = list(max_Rhat = 1.05, min_ESS = 1e6,
        max_error = NULL, max_SD_error = NULL, max_time = NULL,
        sample_extend = sample_extend, restarts = 1, max_extend = 2), silent = FALSE),
      warning = function(w){
        if(grepl("Only one chain was run", conditionMessage(w), fixed = TRUE)){
          count <<- count + 1L
        }
        invokeRestart("muffleWarning")
      }
    )
    count
  }
  set.seed(1)
  expect_identical(one_chain_warnings(1000), 2L)
  expect_identical(one_chain_warnings("auto"), 2L)
})

# A JAGS_fit() whose initial run returns 'sample' rows of each of 'chains'
# chains (one by default), for the mocked check of
# .autofit_size_test_backend().
.autofit_size_test_fitter <- function(env = parent.frame()){

  n_chains <- 1L
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(sample, ...){
      structure(list(
        mcmc = coda::mcmc.list(lapply(seq_len(n_chains), function(chain){
          coda::mcmc(matrix(0, nrow = sample, ncol = 1, dimnames = list(NULL, "mu")))
        })),
        live_model = structure(list(), class = "jags")
      ), class = "BayesTools_fit")
    },
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_build_parameter_coordinates = function(...) NULL,
    .package = "BayesTools",
    .env = env
  )
  function(control, sample = 1000, thin = 1, chains = 1L, silent = TRUE){
    n_chains <<- chains
    JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = chains, adapt = 100, burnin = 50, sample = sample, thin = thin,
      autofit = TRUE, autofit_control = control, silent = silent, seed = 1
    )
  }
}

test_that("JAGS_fit sizes automatic extensions from each failed check", {

  backend <- .autofit_size_test_backend()
  fit_with <- .autofit_size_test_fitter()
  # (12500 - n) * 1.1 from 1000 draws: 2 n, 2 n, then 3850 rounded up
  fit <- fit_with(.autofit_size_test_control())
  expect_identical(backend$samples, c(2000L, 6000L, 3900L))
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(1000L, 3000L, 9000L, 12900L))
  expect_false(any(vapply(backend$checks, `[[`, logical(1), "fail_fast")))
  expect_identical(nrow(fit$mcmc[[1L]]), 12900L)
  expect_identical(attr(fit, "adaptation")$extensions$adapt, rep(100L, 3L))

  # at most the draws whose sampling iterations the backend can run: 4500
  # draws with ESS 180 need 8800, and thin 250000 allows 8500
  backend <- .autofit_size_test_backend()
  thinned <- suppressWarnings(fit_with(.autofit_size_test_control(max_extend = 1),
    sample = 4500, thin = 250000))
  expect_identical(backend$samples, 8500L)

  # max_time = 0 stops before the first extension, also on a clock that has
  # not advanced since the start: the time left is 0
  testthat::local_mocked_bindings(
    .bt_autofit_clock = function() as.POSIXct("2026-01-01 00:00:00", tz = "UTC"),
    .package = "BayesTools"
  )
  backend <- .autofit_size_test_backend()
  stopped <- fit_with(.autofit_size_test_control(max_time = list(time = 0, unit = "secs")))
  expect_length(backend$samples, 0L)
  expect_length(backend$checks, 1L)
  expect_identical(attr(stopped, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")

  # the fixed step is unchanged
  backend <- .autofit_size_test_backend()
  fixed <- fit_with(.autofit_size_test_control(sample_extend = 1000, max_extend = 2))
  expect_identical(backend$samples, c(1000L, 1000L))
  expect_true(all(vapply(backend$checks, `[[`, logical(1), "fail_fast")))

  # max_time caps each extension at the draws per chain that the time left
  # allows at the rate of all draws so far, on the clock that starts before
  # the initial run. With a 40 s limit: 1000 draws in 10 s leave 30 s for
  # 3000, so 2000 runs; 3000 draws in 30 s leave 10 s for 1000 (instead of
  # 6000); 4000 draws in 39.5 s leave 0.5 s for 50: stop
  clock <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC") + c(0, 10, 30, 39.5, 41)
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_autofit_clock = function(){
      calls <<- calls + 1L
      clock[[calls]]
    },
    .package = "BayesTools"
  )
  backend <- .autofit_size_test_backend()
  capped <- fit_with(.autofit_size_test_control(max_time = list(time = 40, unit = "secs")))
  expect_identical(backend$samples, c(2000L, 1000L))
  expect_identical(calls, 4L)
  expect_identical(attr(capped, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")
  # the fixed step under the same clock: 1000 each, until the limit passes
  calls <- 0L
  backend <- .autofit_size_test_backend()
  fixed <- fit_with(.autofit_size_test_control(sample_extend = 1000,
    max_time = list(time = 40, unit = "secs")))
  expect_identical(backend$samples, c(1000L, 1000L, 1000L))
  expect_identical(calls, 5L)
  expect_identical(attr(fixed, "warnings"),
    "The automatic model fitting was terminated due to the 'max_time' constraint.")

  # equivalent limits give the same cap: 3500 draws in 3500 s leave 100 s for
  # exactly 100 draws under 3600 secs, 60 mins and 1 hours; those run, and
  # the limit has passed at the next gate
  clock <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC") + c(0, 3500, 3601)
  for(limit_unit in list(list(time = 3600, unit = "secs"), list(time = 60, unit = "mins"),
                         list(time = 1, unit = "hours"))){
    calls <- 0L
    backend <- .autofit_size_test_backend()
    boundary <- fit_with(.autofit_size_test_control(max_time = limit_unit), sample = 3500)
    expect_identical(backend$samples, 100L, info = limit_unit$unit)
    expect_identical(calls, 3L, info = limit_unit$unit)
    expect_identical(attr(boundary, "warnings"),
      "The automatic model fitting was terminated due to the 'max_time' constraint.",
      info = limit_unit$unit)
  }
})

# The value of 'expr' and the warnings that it signals.
.autofit_collect_warnings <- function(expr){

  warnings <- list()
  value <- withCallingHandlers(expr, warning = function(w){
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

# The warning of an automatic fit that stopped with diagnostics that are not
# assessable: one warning, of the documented classes, with the 'message';
# with silent = TRUE, no warning at all.
.expect_autofit_not_assessable_warning <- function(run, message, silent, info = NULL){

  if(silent){
    expect_length(run$warnings, 0L)
    return(invisible())
  }
  expect_length(run$warnings, 1L)
  if(length(run$warnings) != 1L){
    return(invisible())
  }
  expect_identical(class(run$warnings[[1L]]), c("BayesTools_autofit_not_assessable",
    "BayesTools_fit_diagnostic", "warning", "condition"), info = info)
  expect_identical(conditionMessage(run$warnings[[1L]]), message, info = info)
}

.autofit_not_assessable_message <- function(not_assessable){
  paste0(not_assessable, " Automatic fitting did not extend the fit for these diagnostics;",
    " check the fit with 'JAGS_check_convergence()'.")
}

test_that("automatic fitting extends only for quantities that can be assessed", {

  fit_with <- .autofit_size_test_fitter()
  extensions <- function(fit) nrow(attr(fit, "adaptation")$extensions)
  # one chain: the check names the ESS of 'theta' only, as the public check
  # skips R-hat
  theta_ESS <- .autofit_not_assessable_message("ESS diagnostic is not assessable for 'theta'.")

  # JAGS_fit(): 'mu' meets min_ESS with the initial 12500 draws (ESS 500) and
  # 'theta' is not assessable: the fitting stops after the initial sampling
  # with either step, after one check that ignores 'theta'. Unless silent, it
  # warns about 'theta' once; the warning is not stored with the fit
  for(control in list(.autofit_size_test_control(), .autofit_size_test_control(sample_extend = 1000))){
    for(silent in c(TRUE, FALSE)){
      info <- paste(control$sample_extend, silent)
      backend <- .autofit_size_test_backend(not_assessable = TRUE)
      run <- .autofit_collect_warnings(fit_with(control, sample = 12500, silent = silent))
      fit <- run$value
      .expect_autofit_not_assessable_warning(run, theta_ESS, silent, info = info)
      expect_length(backend$samples, 0L)
      expect_identical(vapply(backend$checks, `[[`, logical(1), "allow_not_assessable"), TRUE,
        info = info)
      expect_identical(nrow(fit$mcmc[[1L]]), 12500L, info = info)
      expect_identical(extensions(fit), 0L, info = info)
      expect_null(attr(fit, "warnings"))
    }
  }

  # a failing quantity next to one that is not assessable: the extension is
  # what the failing quantity alone needs. From 10000 draws, ESS 400 of 500
  # (f = 1.25) needs 2750 draws, 2800; 12800 draws (ESS 512) then stop. A
  # factor of 2 for 'theta' would have run 11000 and gone on to max_extend.
  backend <- .autofit_size_test_backend(not_assessable = TRUE)
  run <- .autofit_collect_warnings(fit_with(.autofit_size_test_control(), sample = 10000,
    silent = FALSE))
  expect_identical(backend$samples, 2800L)
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(10000L, 12800L))
  .expect_autofit_not_assessable_warning(run, theta_ESS, silent = FALSE)
  expect_null(attr(run$value, "warnings"))

  # an undefined R-hat next to a computed ESS that fails: the ESS of 'theta'
  # still extends the fit and sizes the extension. Two chains of 10000 draws,
  # 'mu' with ESS 400 (f = 1.25) and 'theta' with ESS 200 (f = 2.5): 'theta'
  # needs 16500 draws; 26500 draws (ESS 1060 and 530) then stop, with the
  # warning about the R-hat of 'theta'. Without the ESS of 'theta', 2800
  # draws would have run (as above).
  for(sample_extend in list("auto", 1000)){
    backend <- .autofit_size_test_backend(not_assessable = TRUE, theta_ESS = 0.02)
    run <- .autofit_collect_warnings(fit_with(.autofit_size_test_control(
      sample_extend = sample_extend, max_extend = 20), sample = 10000, chains = 2L,
      silent = FALSE))
    if(identical(sample_extend, "auto")){
      expect_identical(backend$samples, 16500L)
      expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(10000L, 26500L))
    }else{
      # the fixed step runs until 'theta' reaches ESS 500 at 25000 draws
      expect_identical(backend$samples, rep(1000L, 15L))
      expect_identical(backend$checks[[length(backend$checks)]]$n, 25000L)
    }
    expect_identical(nrow(run$value$mcmc[[2L]]), backend$checks[[length(backend$checks)]]$n)
    .expect_autofit_not_assessable_warning(run, .autofit_not_assessable_message(
      "R-hat diagnostic is not assessable for 'theta'."), silent = FALSE,
      info = paste(sample_extend))
    expect_null(attr(run$value, "warnings"))
  }

  # nothing assessable: the fitting stops after the initial sampling
  backend <- .autofit_size_test_backend(not_assessable = TRUE, assessable = FALSE)
  run <- .autofit_collect_warnings(fit_with(.autofit_size_test_control(), silent = FALSE))
  expect_length(backend$samples, 0L)
  expect_identical(extensions(run$value), 0L)
  .expect_autofit_not_assessable_warning(run, theta_ESS, silent = FALSE)
  expect_null(attr(run$value, "warnings"))

  # everything assessable: the stop signals nothing, also when not silent
  backend <- .autofit_size_test_backend()
  run <- .autofit_collect_warnings(fit_with(.autofit_size_test_control(), sample = 12500,
    silent = FALSE))
  expect_length(backend$samples, 0L)
  expect_length(run$warnings, 0L)
  expect_null(attr(run$value, "warnings"))

  # a fit that stops at max_extend keeps its termination warning, stored and
  # signalled, and gives no warning about 'theta'
  backend <- .autofit_size_test_backend(not_assessable = TRUE)
  run <- .autofit_collect_warnings(fit_with(.autofit_size_test_control(max_extend = 1),
    silent = FALSE))
  termination <- "The automatic model fitting was terminated due to the 'max_extend' constraint."
  expect_identical(vapply(run$warnings, conditionMessage, character(1)), termination)
  expect_false(any(vapply(run$warnings, inherits, logical(1), "BayesTools_autofit_not_assessable")))
  expect_identical(attr(run$value, "warnings"), termination)
})

test_that("automatic JAGS_extend extends only for quantities that can be assessed", {

  theta_ESS <- .autofit_not_assessable_message("ESS diagnostic is not assessable for 'theta'.")
  # with 'theta' next to 'mu', the extensions are those without 'theta'
  # (1000, 2000, 6000, 3900 to 12902 draws, as in the test above), and the
  # fixed step stops at the first check that 'mu' passes; unless silent, the
  # stop warns about 'theta' once, without storing the warning
  extend_fit <- .jags_extend_test_fit()
  for(silent in c(TRUE, FALSE)){
    backend <- .autofit_size_test_backend(not_assessable = TRUE)
    run <- .autofit_collect_warnings(JAGS_extend(extend_fit,
      autofit_control = .autofit_size_test_control(), silent = silent))
    expect_identical(backend$samples, c(1000L, 2000L, 6000L, 3900L))
    expect_true(all(vapply(backend$checks, `[[`, logical(1), "allow_not_assessable")))
    expect_identical(nrow(run$value$mcmc[[1L]]), 12902L)
    .expect_autofit_not_assessable_warning(run, theta_ESS, silent, info = paste(silent))
    expect_null(attr(run$value, "warnings"))
    backend <- .autofit_size_test_backend(not_assessable = TRUE)
    run <- .autofit_collect_warnings(JAGS_extend(extend_fit,
      autofit_control = .autofit_size_test_control(sample_extend = 5000, max_extend = 10),
      silent = silent))
    # 15002 draws give an ESS of 600
    expect_identical(backend$samples, c(5000L, 5000L, 5000L))
    expect_true(all(vapply(backend$checks, `[[`, logical(1), "fail_fast")))
    .expect_autofit_not_assessable_warning(run, theta_ESS, silent, info = paste(silent))
    expect_null(attr(run$value, "warnings"))
  }
  # nothing assessable: one extension, sized from no failure (the minimum)
  backend <- .autofit_size_test_backend(not_assessable = TRUE, assessable = FALSE)
  run <- .autofit_collect_warnings(JAGS_extend(extend_fit,
    autofit_control = .autofit_size_test_control(), silent = FALSE))
  expect_identical(backend$samples, 1000L)
  expect_identical(vapply(backend$checks, `[[`, integer(1), "n"), c(2L, 1002L))
  .expect_autofit_not_assessable_warning(run, theta_ESS, silent = FALSE)
  expect_null(attr(run$value, "warnings"))
  # everything assessable: the stop signals nothing, also when not silent
  backend <- .autofit_size_test_backend()
  run <- .autofit_collect_warnings(JAGS_extend(extend_fit,
    autofit_control = .autofit_size_test_control(), silent = FALSE))
  expect_identical(backend$samples, c(1000L, 2000L, 6000L, 3900L))
  expect_length(run$warnings, 0L)
  expect_null(attr(run$value, "warnings"))
  # a fit that stops at max_extend keeps its termination warning only
  backend <- .autofit_size_test_backend(not_assessable = TRUE)
  run <- .autofit_collect_warnings(JAGS_extend(extend_fit,
    autofit_control = .autofit_size_test_control(max_extend = 1), silent = FALSE))
  termination <- "The automatic model fitting was terminated due to the 'max_extend' constraint."
  expect_identical(vapply(run$warnings, conditionMessage, character(1)), termination)
  expect_identical(attr(run$value, "warnings"), termination)
})

test_that("the time cap of automatic extensions follows the observed rate", {

  limit <- list(time = 60, unit = "secs")
  # 1000 draws in 10 s leave 50 s for 5000; 999 draws leave 4995, rounded down
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 10, draws = 1000), 5000)
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 10, draws = 999), 4900)
  # elapsed is in seconds, and the limit is converted to seconds: equivalent
  # limits give the same cap (3500 draws in 3500 s leave 100 s for 100)
  expect_identical(.bt_autofit_time_cap(list(time = 2, unit = "mins"), elapsed = 30, draws = 1000), 3000)
  for(limit_unit in list(list(time = 3600, unit = "secs"), list(time = 60, unit = "mins"),
                         list(time = 1, unit = "hours"))){
    expect_identical(.bt_autofit_time_cap(limit_unit, elapsed = 3500, draws = 3500), 100,
      info = limit_unit$unit)
  }
  # no limit, or time left but no draws or no time on the clock yet: no cap
  expect_identical(.bt_autofit_time_cap(NULL, elapsed = 10, draws = 1000), Inf)
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 10, draws = 0), Inf)
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 0, draws = 1000), Inf)
  # no time left: 0, also without an observed rate
  zero <- list(time = 0, unit = "secs")
  expect_identical(.bt_autofit_time_cap(zero, elapsed = 0, draws = 1000), 0)
  expect_identical(.bt_autofit_time_cap(zero, elapsed = 0, draws = 0), 0)
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 60, draws = 0), 0)
  expect_identical(.bt_autofit_time_cap(limit, elapsed = 61, draws = 1000), 0)
  # less than 100 draws left: below 100, so the caller stops
  expect_lt(.bt_autofit_time_cap(limit, elapsed = 59.5, draws = 1000), 100)
})

test_that("an automatic fit of a slowly mixing model reaches its ESS target in fewer extensions with 'auto'", {

  skip_if_not_installed("rjags")
  # y[i] ~ N(a + b, 1) with 100 observations: the posterior correlation of a
  # and b is -100 / 101, and their Gibbs updates mix slowly (an ESS of about
  # 1% of the draws per chain)
  set.seed(2)
  data <- list(y = stats::rnorm(100, 0.5), N = 100L)
  fit_with <- function(sample_extend, max_time = NULL){
    JAGS_fit("model{ for(i in 1:N){ y[i] ~ dnorm(a + b, 1) } }", data = data,
      prior_list = list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1))),
      chains = 2, adapt = 100, burnin = 100, sample = 1000, seed = 3, silent = TRUE,
      autofit = TRUE, autofit_control = list(max_Rhat = NULL, min_ESS = 500, max_error = NULL,
        max_SD_error = NULL, max_time = max_time, sample_extend = sample_extend, max_extend = 50))
  }
  extensions <- function(fit) nrow(attr(fit, "adaptation")$extensions)
  reached <- function(fit) isTRUE(JAGS_check_convergence(fit, max_Rhat = NULL, min_ESS = 500,
    max_error = NULL, max_SD_error = NULL))

  fixed <- fit_with(1000)
  auto <- fit_with("auto")
  # without time, the readiness check also stops extending burn-in, and the
  # chains (which disagree after 100 sweeps) start sampling with its warning
  expect_warning(
    stopped_fixed <- fit_with(1000, max_time = list(time = 0, unit = "secs")),
    class = "BayesTools_burnin_not_converged"
  )
  expect_false(reached(stopped_fixed))
  expect_true(reached(fixed))
  expect_true(reached(auto))
  expect_null(attr(fixed, "warnings"))
  expect_null(attr(auto, "warnings"))
  # the fixed step needs about 24 extensions, "auto" about four
  expect_gt(extensions(fixed), 10L)
  expect_lte(extensions(auto), 5L)
  expect_lt(extensions(auto), extensions(fixed))

  # max_time stops "auto" as it stops the fixed step: here before any extension
  expect_warning(
    stopped <- fit_with("auto", max_time = list(time = 0, unit = "secs")),
    "Increase 'adapt', 'burnin' or 'thin'.", fixed = TRUE, class = "BayesTools_burnin_not_converged"
  )
  expect_identical(extensions(stopped), 0L)
  expect_identical(nrow(stopped$mcmc[[1L]]), 1000L)
  expect_identical(attr(stopped, "burnin_check")$stop_reason, "max_time")
  expect_identical(attr(stopped, "warnings"), c(
    "The automatic model fitting was terminated due to the 'max_time' constraint.",
    BayesTools:::.bt_burnin_warning(attr(stopped, "burnin_check"))
  ))
})

test_that("an automatic fit stops at once when only a component that no chain visits is not assessable", {

  skip_if_not_installed("rjags")
  # 'theta' mixes a spike at 0 with a normal component of prior weight 1e-12,
  # which no chain visits: 'theta' stays 0 and its diagnostics are not
  # assessable, while 'mu' converges in the initial sampling. The fitting
  # does not wait for a visit: it stops after the initial sampling, well
  # before max_time (before, it extended until max_time), and warns about
  # 'theta' unless silent
  set.seed(1)
  data <- list(y = stats::rnorm(50), N = 50L)
  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    theta = prior_mixture(list(
      prior("spike", list(0), prior_weights = 1),
      prior("normal", list(0, 1), prior_weights = 1e-12)
    ))
  )
  max_time <- list(time = 30, unit = "secs")
  fit_with <- function(sample_extend, silent){
    JAGS_fit("model{ for(i in 1:N){ y[i] ~ dnorm(mu + theta, 1) } }", data = data,
      prior_list = prior_list, chains = 2, adapt = 100, burnin = 100, sample = 1000,
      seed = 1, silent = silent, autofit = TRUE, autofit_control = list(max_Rhat = 1.05,
        min_ESS = 500, max_error = NULL, max_SD_error = NULL, max_time = max_time,
        sample_extend = sample_extend, max_extend = NULL))
  }
  check <- function(fit, allow_not_assessable){
    JAGS_check_convergence(fit, max_Rhat = 1.05, min_ESS = 500, max_error = NULL,
      max_SD_error = NULL, allow_not_assessable = allow_not_assessable)
  }
  # "auto" runs with silent = FALSE and warns; the fixed step runs silently
  for(sample_extend in list("auto", 1000)){
    info <- paste(sample_extend)
    silent <- !identical(sample_extend, "auto")
    start <- Sys.time()
    run <- .autofit_collect_warnings(fit_with(sample_extend, silent = silent))
    fit <- run$value
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    expect_lt(elapsed, max_time$time / 3)
    expect_identical(nrow(attr(fit, "adaptation")$extensions), 0L, info = info)
    expect_identical(nrow(fit$mcmc[[1L]]), 1000L, info = info)
    .expect_autofit_not_assessable_warning(run, .autofit_not_assessable_message(
      "R-hat, ESS diagnostics are not assessable for 'theta'."), silent, info = info)
    expect_null(attr(fit, "warnings"))
    expect_identical(unique(unlist(lapply(fit$mcmc, function(chain) chain[, "theta"]))), 0,
      info = info)
    # the check of the returned fit reports 'theta' unless allow_not_assessable
    for(allow in c(FALSE, TRUE)){
      checked <- check(fit, allow_not_assessable = allow)
      diagnostics <- attr(checked, "diagnostics")
      expect_identical(diagnostics$state[match(c("mu", "theta"), diagnostics$parameter)],
        c("assessable", "not_assessable"), info = info)
      if(allow){
        expect_true(checked)
        expect_null(attr(checked, "errors"))
      }else{
        expect_false(checked)
        expect_identical(attr(checked, "errors"),
          "R-hat, ESS diagnostics are not assessable for 'theta'.", info = info)
      }
    }
  }
})

# ============================================================================ #
# 'adapt', 'burnin' and 'sample' count draws at the fit's thinning (the live
# fits are in test-certification-fit-assertions.R)
# ============================================================================ #
test_that("the settings are refused when the last draw would be numbered beyond the integer range", {

  settings <- function(adapt, burnin, sample, thin){
    JAGS_check_and_list_fit_settings(chains = 1, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
      autofit = FALSE, parallel = FALSE, cores = 1, silent = TRUE, seed = 1)
  }
  limit <- .Machine$integer.max
  message <- function(last){
    paste0("The settings 'adapt', 'burnin', 'sample' and 'thin' would number the last retained draw ", last,
      ", beyond the largest iteration number (2147483647): a fit runs ('adapt' + 'burnin' + 'sample') times ",
      "'thin' iterations. Reduce 'adapt', 'burnin', 'sample' or 'thin'.")
  }
  # the last draw is iteration (adapt + burnin + sample - 1) * thin + 1
  expect_silent(settings(50, 50, limit - 100, 1))
  expect_error(settings(50, 50, limit - 99, 1), message("2147483648"), fixed = TRUE)
  expect_silent(settings(100, 100, 1073741624, 2))
  expect_error(settings(100, 100, 1073741625, 2), message("2147483649"), fixed = TRUE)
  expect_error(settings(1073741824, 50, 100, 2), message("2147483947"), fixed = TRUE)
  expect_error(settings(50, 1073741824, 100, 2L), message("2147483947"), fixed = TRUE)
  # a single setting beyond the integer range is named by its own check
  expect_error(settings(limit + 1, 50, 100, 1), "The 'adapt' must be equal or lower than 2147483647.", fixed = TRUE)
})
