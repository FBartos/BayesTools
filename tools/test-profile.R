args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) > 0L && nzchar(args[[1L]])) args[[1L]] else Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit")

if (Sys.getenv("AGENT") == "") {
  Sys.setenv(AGENT = "1")
}

source(file.path("tests", "testthat", "helper-test-profiles.R"))

Sys.setenv(BAYESTOOLS_TEST_PROFILE = profile)

requested_profile_parts <- unique(unlist(strsplit(tolower(profile), "[,;[:space:]]+")))
selected_profiles <- bayestools_normalize_test_profiles(profile)
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

message("Running BayesTools test profile: ", Sys.getenv("BAYESTOOLS_TEST_PROFILE"))
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
if (!exists("LlmReporter", envir = asNamespace("testthat"), inherits = FALSE)) {
  stop("testthat::LlmReporter is required; install testthat >= 3.3.0.", call. = FALSE)
}

test_filter <- if (run_all_profiles) NULL else bayestools_test_profile_filter(selected_profiles)

if (!is.null(test_filter)) {
  message("testthat filter selects ", length(bayestools_test_profile_selected_contexts(selected_profiles)), " context(s).")
}

message("testthat reporter: LlmReporter")
devtools::test(
  filter = test_filter,
  stop_on_failure = TRUE,
  reporter = testthat::LlmReporter$new()
)
