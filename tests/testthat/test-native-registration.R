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
