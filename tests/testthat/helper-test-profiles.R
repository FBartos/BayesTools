bayestools_known_test_profiles <- c("unit", "fixture", "visual", "visual-fixture", "fit")


bayestools_quiet_llm_reporter <- function(...) {

  reporter_class <- R6::R6Class(
    classname = "BayesToolsQuietLlmReporter",
    inherit   = testthat::LlmReporter,
    public    = list(
      add_result = function(context, test, result) {

        if (self$is_full()) {
          return(invisible())
        }
        if (inherits(result, "expectation_skip")) {
          self$n_skip <- self$n_skip + 1L
          return(invisible())
        }

        super$add_result(context, test, result)
      }
    )
  )

  reporter_class$new(...)
}

bayestools_test_profile_contexts <- list(
  unit = c(
    "backend-fingerprint",
    "distributions-mpoint",
    "distributions-point",
    "distributions-tools",
    "distributions-weightfunctions",
    "factor-interaction-coefficients",
    "fixture-catalog-static",
    "hypothesis-ast",
    "hypothesis-BF",
    "hypothesis-BF-parser-adversarial",
    "interpret",
    "interactive-test-runner",
    "JAGS-bridge-compiler",
    "JAGS-bridge-formula-context-validation",
    "JAGS-bridgesampling-wrapper",
    "JAGS-convergence",
    "JAGS-diagnostic-plot-data",
    "JAGS-diagnostics-controls",
    "JAGS-draw-geometry",
    "JAGS-fit-contract",
    "JAGS-fit-settings",
    "JAGS-formula-coefficient-density",
    "JAGS-formula-default-priors",
    "JAGS-formula-design-oracles",
    "JAGS-formula-prediction-targets",
    "JAGS-indexed-parameters",
    "JAGS-lkj-cholesky",
    "JAGS-marginal-distributions",
    "JAGS-parameter-catalog",
    "JAGS-parameter-coordinates",
    "JAGS-posterior-extraction",
    "JAGS-random-effects-compile",
    "JAGS-selection-inits",
    "JAGS-structured-rho-support",
    "JAGS-random-effect-scaling",
    "marginal-inference-conditioning",
    "marginal-prior-rng",
    "marginal-prior-samplers",
    "model-averaging-compatibility-guards",
    "model-averaging-edge-cases",
    "model-averaging-plots-edge-cases",
    "native-registration",
    "prior-density-ordinate",
    "prior-ordered",
    "priors-coverage",
    "priors-density-numeric",
    "priors-informed",
    "priors-linear-density",
    "priors-nonlocal",
    "priors-plot-data",
    "priors-print",
    "priors-tools",
    "precomputed-vignette-cache",
    "random-parameterization",
    "random-effects-correlation-draws",
    "random-effects-independent-scalability",
    "random-effects-marginal-covariance",
    "random-effects-rank-one",
    "random-effects-marginal-update",
    "random-effects-memory",
    "random-effects-structured-local",
    "random-effects-summary-posterior",
    "random-effects-vignette-cache",
    "reference-table-review",
    "selection-kernels",
    "summary-tables-helpers",
    "test-layout-policy",
    "tools-evaluation",
    "tools-input",
    "weightfunction-plot-analytic",
    "weightfunction-redesign"
  ),
  fixture = c(
    "fixture-integrity",
    "hypothesis-BF-bridge",
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
    "fixture-integrity"
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
