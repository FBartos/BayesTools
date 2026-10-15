skip_if_not_test_profile("unit")

test_that("manual native load ownership survives unavailable symbol inventory", {
  original_private <- .BayesTools_private
  private_fields <- as.list(original_private, all.names = TRUE)
  loader_original <- .BayesTools_load_native_routines
  unload_original <- getFromNamespace(".onUnload", "BayesTools")
  for(scenario in c("standard", "manual", "missing_symbols", "failed_manual")){
    isolated <- new.env(parent = asNamespace("BayesTools"))
    state <- new.env(parent = emptyenv())
    directory <- withr::local_tempdir()
    module <- file.path(directory, paste0("BayesTools", .Platform$dynlib.ext))
    writeBin(charToRaw("inert fixture; loading is mocked"), module)
    state$module_location <- directory; state$lib_name <- tempdir()
    state$native_routines_loaded <- FALSE; state$native_dll_path <- NULL; state$native_dll_loaded_manually <- FALSE
    standard_calls <- manual_calls <- standard_unload <- manual_unload <- 0L
    isolated$.BayesTools_private <- state
    isolated$.BayesTools_native_routines_loaded <- function(pkgname){
      (scenario == "standard" && standard_calls > 0L) || (scenario == "manual" && manual_calls > 0L)
    }
    isolated$library.dynam <- function(...){standard_calls <<- standard_calls + 1L; if(scenario != "standard") stop("isolated standard-loader failure")}
    isolated$dyn.load <- function(...){manual_calls <<- manual_calls + 1L; if(scenario == "failed_manual") stop("isolated manual-loader failure")}
    isolated$requireNamespace <- function(...) FALSE
    isolated$library.dynam.unload <- function(...){standard_unload <<- standard_unload + 1L}
    isolated$dyn.unload <- function(...){manual_unload <<- manual_unload + 1L}
    loader <- loader_original; environment(loader) <- isolated
    unload <- unload_original; environment(unload) <- isolated
    expect_identical(body(loader), body(loader_original))
    expect_identical(formals(loader), formals(loader_original))
    condition <- NULL
    result <- withCallingHandlers(loader(pkgname = "BayesTools", libname = tempdir(), warn = TRUE),
      warning = function(w){condition <<- w; invokeRestart("muffleWarning")})
    if(scenario %in% c("missing_symbols", "failed_manual")){
      expect_s3_class(condition, "warning")
      expect_match(conditionMessage(condition), "BayesTools native routines failed to load", fixed = TRUE)
    }else expect_null(condition)
    expect_identical(result, scenario %in% c("standard", "manual"))
    expect_identical(state$native_dll_loaded_manually, scenario %in% c("manual", "missing_symbols"))
    if(scenario %in% c("manual", "missing_symbols")) expect_identical(state$native_dll_path, normalizePath(module, winslash = "/"))
    else expect_null(state$native_dll_path)
    expect_identical(standard_calls, 1L)
    expect_identical(manual_calls, if(scenario == "standard") 0L else 1L)
    unload(tempdir())
    expect_identical(standard_unload, 1L)
    expect_identical(manual_unload, if(scenario %in% c("manual", "missing_symbols")) 1L else 0L)
    expect_false(state$native_dll_loaded_manually)
    expect_null(state$native_dll_path)
  }
  expect_identical(as.list(original_private, all.names = TRUE), private_fields)
})

test_that("native loader inventory matches registered call routines", {

  dll <- getLoadedDLLs()[["BayesTools"]]
  expect_false(is.null(dll))

  registered <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expect_setequal(.BayesTools_native_symbols(), names(registered))
  expect_true(all(vapply(
    .BayesTools_native_symbols(),
    is.loaded,
    logical(1),
    PACKAGE = "BayesTools"
  )))
})

test_that("native routine guards name their routine family", {

  testthat::local_mocked_bindings(
    .BayesTools_native_routines_loaded = function(pkgname = "BayesTools") FALSE,
    .BayesTools_load_native_routines = function(...) invisible(FALSE),
    .package = "BayesTools"
  )
  expect_error(.BayesTools_require_native_lkj(),
    "BayesTools native LKJ routines are not loaded.", fixed = TRUE)
  expect_error(.BayesTools_require_native_nonlocal(),
    "BayesTools native nonlocal-prior routines are not loaded.", fixed = TRUE)
  expect_error(.BayesTools_require_native_invgamma(),
    "BayesTools native inverse-gamma routines are not loaded.", fixed = TRUE)
  expect_error(
    .bt_random_effect_native_structured_cholesky("cs", rho = 0.5, coordinates = c(1, 1)),
    "BayesTools native structured-Cholesky routines are not loaded.", fixed = TRUE)
})

test_that("unloading the package unloads its JAGS module without a message", {

  seen <- new.env(parent = emptyenv())
  seen$quiet <- NULL
  testthat::local_mocked_bindings(
    list.modules = function() "BayesTools",
    unload.module = function(name, quiet = FALSE){
      seen$name <- name
      seen$quiet <- quiet
      if(!quiet) cat("Module", name, "unloaded\n")
      invisible()
    },
    .package = "rjags"
  )
  isolated <- new.env(parent = asNamespace("BayesTools"))
  isolated$.BayesTools_private <- list2env(list(
    native_routines_loaded = TRUE, native_dll_path = NULL, native_dll_loaded_manually = FALSE
  ), parent = emptyenv())
  isolated$library.dynam.unload <- function(...) invisible(NULL)
  unload <- getFromNamespace(".onUnload", "BayesTools")
  environment(unload) <- isolated
  output <- capture.output(unload(tempdir()))
  expect_identical(seen$name, "BayesTools")
  expect_true(isTRUE(seen$quiet))
  expect_identical(output, character())
})

test_that("the JAGS module loads and registers its six distribution and function names (F065)", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  expect_true(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))
  expect_true("BayesTools" %in% rjags::list.modules())

  # the generated syntax names the four distributions and the two functions;
  # JAGS compiles (and initializes, without sampling) the model only if the
  # module registered every one of them
  registered <- c("dbt_lkj_cpc", "dbt_invgamma", "dbt_moment", "dbt_invmoment",
                  "bt_lkj_cholesky", "bt_lkj_corr")
  syntax <- JAGS_add_priors(
    paste0("model{
", JAGS_lkj_corr_cholesky("Sigma", K = 3, eta = 1)$syntax, "
}"),
    list(
      a = prior("moment", list(location = 0, tau = 1, order = 2)),
      b = prior("invmoment", list(location = 0, tau = 1, order = 2, df = 3)),
      c = prior("invgamma", list(2, 1))
    )
  )
  for(name in registered){
    expect_true(grepl(name, syntax, fixed = TRUE), info = name)
  }
  model <- rjags::jags.model(textConnection(syntax), n.chains = 1, n.adapt = 0, quiet = TRUE)
  expect_s3_class(model, "jags")
})
