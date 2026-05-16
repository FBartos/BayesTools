bayestools_known_test_profiles <- c("unit", "fixture", "visual", "visual-fixture", "fit")

bayestools_test_profile_contexts <- list(
  unit = c(
    "distributions-mpoint",
    "distributions-point",
    "distributions-tools",
    "distributions-weightfunctions",
    "factor-interaction-coefficients",
    "fixture-catalog-static",
    "interpret",
    "JAGS-diagnostic-plot-data",
    "JAGS-formula-design-oracles",
    "JAGS-lkj-cholesky",
    "JAGS-marginal-distributions",
    "JAGS-posterior-extraction",
    "marginal-inference-conditioning",
    "marginal-prior-samplers",
    "model-averaging-edge-cases",
    "model-averaging-plots-edge-cases",
    "priors-coverage",
    "priors-density-numeric",
    "priors-informed",
    "priors-linear-density",
    "priors-plot-data",
    "priors-print",
    "priors-tools",
    "selection-kernels",
    "summary-tables-helpers",
    "tools-evaluation",
    "tools-input",
    "weightfunction-plot-analytic",
    "weightfunction-redesign"
  ),
  fixture = c(
    "fixture-integrity",
    "JAGS-ensemble-tables",
    "JAGS-fit",
    "JAGS-formula-scale",
    "JAGS-formula",
    "JAGS-summary-tables",
    "model-averaging",
    "selection-kernels",
    "summary-tables",
    "weightfunction-redesign"
  ),
  visual = c(
    "JAGS-ensemble-plots",
    "marginal-prior-samplers",
    "model-averaging-plots",
    "priors",
    "priors-density",
    "priors-plot",
    "priors-print"
  ),
  `visual-fixture` = c(
    "JAGS-diagnostic-plots",
    "JAGS-ensemble-plots",
    "JAGS-marginal-distributions",
    "model-averaging-plots"
  ),
  fit = c(
    "00-model-fits",
    "fixture-integrity",
    "JAGS-fit-edge-cases",
    "JAGS-fit-lm-oracles",
    "JAGS-lkj-cholesky-fit",
    "JAGS-marglik"
  )
)

bayestools_normalize_test_profiles <- function(profiles = NULL) {
  if (is.null(profiles) || length(profiles) == 0L) {
    profiles <- Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit")
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
    `cached-visual` = "visual-fixture",
    cached_visual = "visual-fixture",
    heavy = "fit",
    `jags-visual` = "visual-fixture",
    jags_visual = "visual-fixture",
    slow = "fit",
    visual_fixture = "visual-fixture",
    `visual-fixtures` = "visual-fixture",
    visual_fixtures = "visual-fixture"
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

bayestools_test_profile_selected_contexts <- function(profiles = NULL) {
  unique(unlist(bayestools_test_profile_contexts[bayestools_normalize_test_profiles(profiles)]))
}

bayestools_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

bayestools_test_profile_filter <- function(profiles = NULL) {
  contexts <- bayestools_test_profile_selected_contexts(profiles)
  paste0("^(", paste(bayestools_escape_regex(contexts), collapse = "|"), ")$")
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

skip_if_not_visual_fixture_tests <- function() {
  skip_if_not_test_profile("visual-fixture")
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
