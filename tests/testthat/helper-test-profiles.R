bayestools_known_test_profiles <- c("unit", "fixture", "visual", "fit")

bayestools_normalize_test_profiles <- function(profiles = NULL) {
  if (is.null(profiles) || length(profiles) == 0L) {
    profiles <- Sys.getenv("BAYESTOOLS_TEST_PROFILE", "all")
  }

  profiles <- unlist(strsplit(as.character(profiles), "[,;[:space:]]+"))
  profiles <- tolower(profiles[nzchar(profiles)])

  if (length(profiles) == 0L) {
    profiles <- "all"
  }

  aliases <- c(
    fixtures = "fixture",
    plot = "visual",
    plots = "visual",
    snapshot = "visual",
    snapshots = "visual",
    heavy = "fit",
    slow = "fit"
  )

  alias_matches <- profiles %in% names(aliases)
  profiles[alias_matches] <- unname(aliases[profiles[alias_matches]])

  if ("all" %in% profiles) {
    profiles <- bayestools_known_test_profiles
  }

  invalid <- setdiff(profiles, bayestools_known_test_profiles)
  if (length(invalid) > 0L) {
    stop(
      "Unknown BAYESTOOLS_TEST_PROFILE value(s): ",
      paste(invalid, collapse = ", "),
      ". Expected one or more of: all, ",
      paste(bayestools_known_test_profiles, collapse = ", "),
      call. = FALSE
    )
  }

  unique(profiles)
}

bayestools_test_profile <- function() {
  bayestools_normalize_test_profiles()
}

bayestools_test_profile_label <- function(profiles = bayestools_test_profile()) {
  paste(profiles, collapse = ",")
}

bayestools_test_profile_includes <- function(profiles) {
  any(bayestools_test_profile() %in% bayestools_normalize_test_profiles(profiles))
}

skip_if_not_test_profile <- function(profiles) {
  expected <- bayestools_normalize_test_profiles(profiles)

  if (!any(bayestools_test_profile() %in% expected)) {
    testthat::skip(sprintf(
      "Skipping %s-profile test file because BAYESTOOLS_TEST_PROFILE=%s.",
      paste(expected, collapse = "/"),
      bayestools_test_profile_label()
    ))
  }

  invisible(TRUE)
}

skip_if_not_visual_tests <- function() {
  skip_if_not_test_profile("visual")
  testthat::skip_if_not_installed("vdiffr")
  invisible(TRUE)
}

skip_if_not_heavy_tests <- function() {
  skip_if_not_test_profile("fit")
  invisible(TRUE)
}

skip_if_no_jags_runtime <- function() {
  testthat::skip_if_not_installed("rjags")
  invisible(TRUE)
}
