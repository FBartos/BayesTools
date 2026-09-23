args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) > 0L && nzchar(args[[1L]])) args[[1L]] else Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit")
filter <- if (identical(tolower(profile), "filter")) {
  if (length(args) < 2L || !any(nzchar(args[-1L]))) {
    stop("Profile 'filter' requires one testthat filter.", call. = FALSE)
  }
  paste(args[-1L], collapse = "|")
} else {
  if (length(args) > 1L) {
    stop("Unexpected test-profile arguments: ", paste(args[-1L], collapse = ", "), call. = FALSE)
  }
  NULL
}

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-profile.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

if (Sys.getenv("AGENT") == "" &&
    Sys.getenv("BAYESTOOLS_TEST_REPORTER") == "") {
  Sys.setenv(AGENT = "1")
}

source(file.path("tests", "testthat", "helper-test-profiles.R"))

profile_selection <- if (is.null(filter)) profile else "all"
Sys.setenv(BAYESTOOLS_TEST_PROFILE = profile_selection)

requested_profile_parts <- unique(unlist(strsplit(tolower(profile_selection), "[,;[:space:]]+")))
selected_profiles <- bayestools_normalize_test_profiles(profile_selection)
run_all_profiles <- "all" %in% requested_profile_parts

resolve_test_files_dir <- function() {
  test_files_dir <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
  if (test_files_dir == "") {
    test_files_dir <- file.path(tempdir(), "BayesTools_test_files")
  }
  if (!dir.exists(test_files_dir)) {
    dir.create(test_files_dir, showWarnings = FALSE, recursive = TRUE)
  }
  normalizePath(test_files_dir, winslash = "/", mustWork = TRUE)
}

Sys.setenv(BAYESTOOLS_TEST_FILES_DIR = resolve_test_files_dir())

if (any(c("fit", "visual-fixture") %in% selected_profiles)) {
  Sys.setenv(NOT_CRAN = "true")
}
if (any(c("visual", "visual-fixture") %in% selected_profiles)) {
  Sys.setenv(VDIFFR_RUN_TESTS = "true")
}
if ("fit" %in% selected_profiles && Sys.getenv("BAYESTOOLS_TEST_SKIP_REFIT") == "") {
  Sys.setenv(BAYESTOOLS_TEST_SKIP_REFIT = "false")
}

message(
  "Running BayesTools test profile: ",
  if (is.null(filter)) Sys.getenv("BAYESTOOLS_TEST_PROFILE") else "filter"
)
message("Selected profile lane(s): ", paste(selected_profiles, collapse = ", "))
message("BAYESTOOLS_TEST_FILES_DIR: ", Sys.getenv("BAYESTOOLS_TEST_FILES_DIR"))
message("BAYESTOOLS_TEST_SKIP_REFIT: ", Sys.getenv("BAYESTOOLS_TEST_SKIP_REFIT", unset = "<unset>"))
message("NOT_CRAN: ", Sys.getenv("NOT_CRAN", unset = "<unset>"))
message("VDIFFR_RUN_TESTS: ", Sys.getenv("VDIFFR_RUN_TESTS", unset = "<unset>"))
message("AGENT: ", Sys.getenv("AGENT", unset = "<unset>"))

started <- Sys.time()
on.exit({
  elapsed <- difftime(Sys.time(), started, units = "secs")
  message("Elapsed seconds: ", round(as.numeric(elapsed), 1))
}, add = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The devtools package is required to run tools/test-profile.R.", call. = FALSE)
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The testthat package is required to run tools/test-profile.R.", call. = FALSE)
}
reporter <- Sys.getenv("BAYESTOOLS_TEST_REPORTER", unset = "llm")
if (identical(reporter, "llm") &&
    !exists("LlmReporter", envir = asNamespace("testthat"), inherits = FALSE)) {
  stop("testthat::LlmReporter is required; install testthat >= 3.3.0.", call. = FALSE)
}

test_filter <- if (!is.null(filter)) {
  filter
} else if (run_all_profiles) {
  NULL
} else {
  bayestools_test_profile_filter(selected_profiles)
}

if (!is.null(filter)) {
  message("testthat filter: ", filter)
} else if (!is.null(test_filter)) {
  message("testthat filter selects ", length(bayestools_test_profile_selected_contexts(selected_profiles)), " context(s).")
}

env_flag <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  value <- tolower(trimws(value))
  if (value %in% c("1", "true", "yes", "on")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "no", "off")) {
    return(FALSE)
  }

  stop("Environment variable '", name, "' must be TRUE or FALSE.", call. = FALSE)
}

stop_on_failure <- env_flag("BAYESTOOLS_TEST_STOP_ON_FAILURE", TRUE)
reporter_object <- if (identical(reporter, "llm")) {
  if (env_flag("BAYESTOOLS_TEST_QUIET_SKIPS", FALSE)) {
    bayestools_quiet_llm_reporter()
  } else {
    testthat::LlmReporter$new()
  }
} else {
  reporter
}

message("testthat reporter: ", reporter)
devtools::test(
  filter = test_filter,
  stop_on_failure = stop_on_failure,
  reporter = reporter_object
)
