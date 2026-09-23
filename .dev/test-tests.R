.bayestools_test_project_root <- local({

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    ofile <- frame[["ofile"]]
    if (is.null(ofile)) "" else ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  candidates <- c(
    dirname(normalizePath(ofiles, winslash = "/", mustWork = FALSE)),
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  )

  for (candidate in unique(candidates)) {
    repeat {
      description <- file.path(candidate, "DESCRIPTION")
      if (file.exists(description)) {
        package <- tryCatch(
          read.dcf(description, fields = "Package")[[1L]],
          error = function(error) ""
        )
        if (identical(package, "BayesTools")) {
          return(candidate)
        }
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate)) {
        break
      }
      candidate <- parent
    }
  }

  stop(
    "Cannot locate the BayesTools package root from the sourced file or working directory.",
    call. = FALSE
  )
})

.bayestools_test_files_dir <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
if (!nzchar(.bayestools_test_files_dir)) {
  .bayestools_test_files_dir <- file.path(
    tempdir(),
    "BayesTools_test_files"
  )
}
dir.create(
  .bayestools_test_files_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
.bayestools_test_files_dir <- normalizePath(
  .bayestools_test_files_dir,
  winslash = "/",
  mustWork = TRUE
)
Sys.setenv(BAYESTOOLS_TEST_FILES_DIR = .bayestools_test_files_dir)

if (!exists("reference_table_candidates", mode = "function")) {
  source(
    file.path(
      .bayestools_test_project_root,
      "tests",
      "testthat",
      "helper-00-reference-table-review.R"
    ),
    local = TRUE
  )
}
if (!exists("bayestools_quiet_llm_reporter", mode = "function")) {
  source(
    file.path(
      .bayestools_test_project_root,
      "tests",
      "testthat",
      "helper-test-profiles.R"
    ),
    local = TRUE
  )
}
.bayestools_clean_test_cache <- local({

  cache_root <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
  if (!nzchar(cache_root)) {
    cache_root <- file.path(tempdir(), "BayesTools_test_files")
  }
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_root <- normalizePath(cache_root, winslash = "/", mustWork = TRUE)
  unsafe_roots <- normalizePath(
    c(path.expand("~"), .bayestools_test_project_root),
    winslash = "/",
    mustWork = TRUE
  )

  function() {

    if (identical(dirname(cache_root), cache_root) ||
        cache_root %in% unsafe_roots) {
      stop("Refusing to use a broad directory as the BayesTools test cache.",
           call. = FALSE)
    }

    cache_names <- c("fits", "margliks", "temp")
    targets <- file.path(cache_root, cache_names)
    parents <- normalizePath(
      dirname(targets),
      winslash = "/",
      mustWork = TRUE
    )
    if (length(parents) != length(targets) || any(parents != cache_root)) {
      stop("Refusing to clean paths outside the BayesTools test cache.", call. = FALSE)
    }

    for (target in targets) {
      unlink(target, recursive = TRUE, force = TRUE)
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
    }
    message("Cleaned BayesTools test cache: ", cache_root)

    invisible(TRUE)
  }
})


.bayestools_test_validate_flag <- function(value, name) {

  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }

  value
}


.bayestools_test_rscript <- function() {

  command <- Sys.which("Rscript")
  if (!nzchar(command)) {
    command <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
  }

  command
}


.run_bayestools_test_profile <- function(profile, filter = NULL) {

  profile_script <- file.path(
    .bayestools_test_project_root,
    "tools",
    "test-profile.R"
  )
  arguments <- c(shQuote(profile_script), profile)
  if (!is.null(filter)) {
    arguments <- c(arguments, shQuote(filter))
  }

  message("Running BayesTools test profile: ", profile)
  status <- system2(.bayestools_test_rscript(), args = arguments)
  if (status != 0L) {
    stop("BayesTools test profile failed: ", profile, call. = FALSE)
  }

  invisible(TRUE)
}


review_test_snapshots <- function(root = file.path(
                                    .bayestools_test_project_root,
                                    "tests",
                                    "testthat"
                                  ),
                                  force = FALSE) {

  force <- .bayestools_test_validate_flag(force, "force")
  root  <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!interactive()) {
    message("Cached test snapshots require an interactive session for review.")
    return(invisible(FALSE))
  }

  snapshot_dir <- file.path(root, "_snaps")
  candidates <- if (dir.exists(snapshot_dir)) {
    list.files(
      snapshot_dir,
      pattern    = "[.]new[.]",
      recursive  = TRUE,
      full.names = TRUE
    )
  } else {
    character()
  }
  if (force || length(candidates) > 0L) {
    required  <- c("shiny", "diffviewer")
    available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
    if (all(available)) {
      testthat::snapshot_review(path = root)
    } else {
      message(
        "Interactive visual-snapshot review is unavailable; install: ",
        paste(required[!available], collapse = ", "), "."
      )
    }
  } else {
    message("No visual snapshots to update.")
  }

  review_reference_tables(root = file.path(dirname(root), "results"))
  invisible(TRUE)
}


# Run BayesTools' profile lanes from an interactive development session.
# A filtered refit runs the centralized fit lane before the selected tests.
test_tests <- function(filter = NULL, reporter = "progress", refit = FALSE,
                       update = FALSE, update_timings = FALSE,
                       regenerate = FALSE,
                       load_package = TRUE, stop_on_failure = FALSE,
                       root = file.path(
                         .bayestools_test_project_root, "tests", "testthat"
                       )) {

  refit           <- .bayestools_test_validate_flag(refit, "refit")
  update          <- .bayestools_test_validate_flag(update, "update")
  update_timings  <- .bayestools_test_validate_flag(
    update_timings,
    "update_timings"
  )
  regenerate      <- .bayestools_test_validate_flag(regenerate, "regenerate")
  load_package    <- .bayestools_test_validate_flag(load_package, "load_package")
  stop_on_failure <- .bayestools_test_validate_flag(
    stop_on_failure,
    "stop_on_failure"
  )
  if (!is.null(filter) &&
      (!is.character(filter) || length(filter) != 1L || is.na(filter))) {
    stop("'filter' must be NULL or one regular expression.", call. = FALSE)
  }
  if (!is.character(reporter) || length(reporter) != 1L ||
      is.na(reporter) || !nzchar(reporter)) {
    stop("'reporter' must be one reporter name.", call. = FALSE)
  }
  if (update_timings) {
    stop(
      "BayesTools tests have no timing baselines; routine non-refit tests ",
      "share a 15-minute suite target.",
      call. = FALSE
    )
  }

  refit  <- refit || regenerate
  update <- update || regenerate
  root         <- normalizePath(root, winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(
    file.path(root, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(project_root, .bayestools_test_project_root)) {
    stop("'root' must be BayesTools' 'tests/testthat' directory.", call. = FALSE)
  }

  if (load_package) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("The devtools package is required to load BayesTools.", call. = FALSE)
    }
    devtools::load_all(project_root, quiet = TRUE)
  }

  environment_names <- c(
    "AGENT",
    "BAYESTOOLS_TEST_PROFILE",
    "BAYESTOOLS_TEST_SKIP_REFIT",
    "BAYESTOOLS_TEST_REPORTER",
    "BAYESTOOLS_TEST_STOP_ON_FAILURE",
    "BAYESTOOLS_TEST_QUIET_SKIPS",
    "NOT_CRAN",
    "VDIFFR_RUN_TESTS"
  )
  old_environment <- Sys.getenv(environment_names, unset = NA_character_)
  on.exit({
    for (i in seq_along(environment_names)) {
      name  <- environment_names[[i]]
      value <- old_environment[[i]]
      if (is.na(value)) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(value), name))
      }
    }
  }, add = TRUE)
  on.exit({
    if (interactive()) {
      try(review_test_snapshots(root = root, force = update))
    }
  }, add = TRUE)

  Sys.setenv(
    BAYESTOOLS_TEST_PROFILE           = "all",
    BAYESTOOLS_TEST_SKIP_REFIT        = if (refit) "FALSE" else "TRUE",
    BAYESTOOLS_TEST_REPORTER          = reporter,
    BAYESTOOLS_TEST_STOP_ON_FAILURE   = if (stop_on_failure) "TRUE" else "FALSE",
    BAYESTOOLS_TEST_QUIET_SKIPS       = if (identical(reporter, "llm")) {
      "TRUE"
    } else {
      "FALSE"
    },
    NOT_CRAN                          = "true",
    VDIFFR_RUN_TESTS                  = "true"
  )
  if (identical(reporter, "llm")) {
    Sys.setenv(AGENT = "1")
  } else {
    Sys.unsetenv("AGENT")
  }

  results <- list()
  if (is.null(filter)) {
    if (refit) {
      .bayestools_clean_test_cache()
    }
    results[["all"]] <- .run_bayestools_test_profile("all")
  } else {
    if (refit) {
      .bayestools_clean_test_cache()
      results[["fit"]] <- .run_bayestools_test_profile("fit")
      Sys.setenv(BAYESTOOLS_TEST_SKIP_REFIT = "TRUE")
    }
    results[["tests"]] <- .run_bayestools_test_profile(
      "filter",
      filter = filter
    )
  }

  invisible(results)
}
