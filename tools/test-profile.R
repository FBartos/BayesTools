args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) > 0L && nzchar(args[[1L]])) args[[1L]] else Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit")

known_profiles <- c("unit", "fixture", "visual", "fit", "all")
if (!profile %in% known_profiles && !grepl("[,;[:space:]]", profile)) {
  stop("Unknown profile '", profile, "'. Expected one of: ", paste(known_profiles, collapse = ", "), call. = FALSE)
}

Sys.setenv(BAYESTOOLS_TEST_PROFILE = profile)

profile_parts <- unique(unlist(strsplit(tolower(profile), "[,;[:space:]]+")))
if ("all" %in% profile_parts || "fit" %in% profile_parts || "heavy" %in% profile_parts || "slow" %in% profile_parts) {
  Sys.setenv(NOT_CRAN = "true")
}
if ("all" %in% profile_parts || "visual" %in% profile_parts || "plot" %in% profile_parts ||
    "plots" %in% profile_parts || "snapshot" %in% profile_parts || "snapshots" %in% profile_parts) {
  Sys.setenv(VDIFFR_RUN_TESTS = "true")
}

message("Running BayesTools test profile: ", Sys.getenv("BAYESTOOLS_TEST_PROFILE"))
message("BAYESTOOLS_TEST_FILES_DIR: ", Sys.getenv("BAYESTOOLS_TEST_FILES_DIR", unset = "<test helper default>"))
message("NOT_CRAN: ", Sys.getenv("NOT_CRAN", unset = "<unset>"))
message("VDIFFR_RUN_TESTS: ", Sys.getenv("VDIFFR_RUN_TESTS", unset = "<unset>"))

started <- Sys.time()
on.exit({
  elapsed <- difftime(Sys.time(), started, units = "secs")
  message("Elapsed seconds: ", round(as.numeric(elapsed), 1))
}, add = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The devtools package is required to run tools/test-profile.R.", call. = FALSE)
}

devtools::test(stop_on_failure = TRUE)
