skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Summary Tables Helper Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for format_BF, format_estimates, and other summary table
#   formatting utilities.
#
# DEPENDENCIES:
#   - common-functions.R: test_reference_table, REFERENCE_DIR
#
# SKIP CONDITIONS:
#   - None (can run on CRAN - pure R with reference file testing)
#
# TAGS: @evaluation, @summary-tables, @formatting
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "summary-tables-helpers")
source(testthat::test_path("common-functions.R"))


test_that("format_BF works correctly", {

  # Basic usage
  BF <- format_BF(3.5)
  expect_equal(as.numeric(BF), 3.5)
  expect_equal(attr(BF, "name"), "BF")
  expect_false(attr(BF, "logBF"))
  expect_false(attr(BF, "BF01"))

  # With BF01 = TRUE (inverted)
  BF_01 <- format_BF(2, BF01 = TRUE)
  expect_equal(as.numeric(BF_01), 0.5)
  expect_equal(attr(BF_01, "name"), "1/BF")
  expect_true(attr(BF_01, "BF01"))

  # With logBF = TRUE
  BF_log <- format_BF(exp(2), logBF = TRUE)
  expect_equal(as.numeric(BF_log), 2, tolerance = 1e-10)
  expect_match(attr(BF_log, "name"), "log\\(BF\\)")
  expect_true(attr(BF_log, "logBF"))

  # With inclusion = TRUE
  BF_incl <- format_BF(5, inclusion = TRUE)
  expect_equal(attr(BF_incl, "name"), "Inclusion BF")

  # With BF01 = TRUE and inclusion = TRUE
  BF_excl <- format_BF(5, BF01 = TRUE, inclusion = TRUE)
  expect_equal(attr(BF_excl, "name"), "Exclusion BF")

  # Combined logBF and BF01
  BF_both <- format_BF(10, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(BF_both), log(0.1), tolerance = 1e-10)
  expect_match(attr(BF_both, "name"), "log\\(1/BF\\)")

  # Vector input with NA - NA_real_ must be part of numeric vector
  BF_vec_na <- format_BF(c(1, 2, NA_real_))
  expect_equal(as.numeric(BF_vec_na)[1:2], c(1, 2))
  expect_true(is.na(as.numeric(BF_vec_na)[3]))

  # Vector input
  BF_vec <- format_BF(c(1, 2, 3))
  expect_equal(as.numeric(BF_vec), c(1, 2, 3))

})


test_that("format_BF input validation works", {

  expect_error(format_BF(-1), "must be equal or higher than 0")
  expect_error(format_BF("3"), "must be a numeric")
  expect_error(format_BF(3, logBF = "TRUE"), "must be a logical")
  expect_error(format_BF(3, BF01 = "TRUE"), "must be a logical")

})


test_that("BF MC error helper formats relative percentage column", {

  expect_equal(.BF_error_column_name(), "error%(Inclusion BF)")
  expect_equal(.BF_error_column_name(BF01 = TRUE), "error%(Inclusion BF)")

})


test_that("format_BF preserves finite-sample BF bounds across BF scales", {

  BF <- c(9, 1 / 9)
  attr(BF, "bound_operator") <- c(">", "<")

  formatted <- format_BF(BF, inclusion = TRUE)
  expect_equal(as.numeric(formatted), as.numeric(BF))
  expect_equal(attr(formatted, "bound_operator"), c(">", "<"))
  expect_equal(attr(formatted[1], "bound_operator"), ">")

  exclusion <- format_BF(BF, BF01 = TRUE, inclusion = TRUE)
  expect_equal(as.numeric(exclusion), 1 / as.numeric(BF))
  expect_equal(attr(exclusion, "bound_operator"), c("<", ">"))
  expect_equal(attr(exclusion, "name"), "Exclusion BF")

  log_exclusion <- format_BF(BF, logBF = TRUE, BF01 = TRUE, inclusion = TRUE)
  expect_equal(as.numeric(log_exclusion), log(1 / as.numeric(BF)), tolerance = 1e-12)
  expect_equal(attr(log_exclusion, "bound_operator"), c("<", ">"))

})


test_that("JAGS diagnostic column selectors support shortcuts and explicit subsets", {

  expect_equal(
    .normalize_diagnostic_columns(TRUE, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns"),
    c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
  )
  expect_equal(
    .normalize_diagnostic_columns(FALSE, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns"),
    character()
  )
  expect_equal(
    .normalize_diagnostic_columns("all", .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns"),
    c("ESS", "MCMC_error", "BF_error_percent")
  )
  expect_equal(
    .normalize_diagnostic_columns("none", .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns"),
    character()
  )
  expect_equal(
    .normalize_diagnostic_columns(c("R_hat", "ESS", "ESS"), .JAGS_estimates_diagnostic_columns(), "diagnostic_columns"),
    c("R_hat", "ESS")
  )

  expect_error(
    .normalize_diagnostic_columns(c("all", "ESS"), .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns"),
    "only by itself"
  )
  expect_error(
    .normalize_diagnostic_columns("R_hat", .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns"),
    "not recognized"
  )
})


test_that("empty JAGS tables honor diagnostic column options and logical shortcuts", {

  old_options <- options(
    BayesTools.JAGS_estimates_diagnostic_columns = c("ESS", "R_hat"),
    BayesTools.JAGS_BF_diagnostic_columns = "BF_error_percent"
  )
  on.exit(options(old_options), add = TRUE)

  estimates_empty <- runjags_estimates_empty_table()
  expect_equal(colnames(estimates_empty), c("Mean", "SD", "0.025", "0.5", "0.975", "ESS", "R_hat"))
  expect_equal(attr(estimates_empty, "type"), c(rep("estimate", 5), "ESS", "R_hat"))

  estimates_all <- runjags_estimates_empty_table(diagnostic_columns = TRUE)
  expect_equal(colnames(estimates_all), c("Mean", "SD", "0.025", "0.5", "0.975", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))

  estimates_none <- runjags_estimates_empty_table(diagnostic_columns = FALSE)
  expect_equal(colnames(estimates_none), c("Mean", "SD", "0.025", "0.5", "0.975"))

  estimates_removed <- runjags_estimates_empty_table(remove_diagnostics = TRUE, diagnostic_columns = "all")
  expect_equal(colnames(estimates_removed), c("Mean", "SD", "0.025", "0.5", "0.975"))

  inference_empty <- runjags_inference_empty_table()
  expect_equal(colnames(inference_empty), c("prior_prob", "post_prob", "inclusion_BF", "BF_error_percent"))
  expect_equal(attr(inference_empty, "type"), c("prior_prob", "post_prob", "inclusion_BF", "BF_error"))
  expect_equal(attr(inference_empty[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
  expect_null(attr(inference_empty, "footnotes"))

  inference_ess <- runjags_inference_empty_table(BF_diagnostic_columns = "ESS")
  expect_equal(colnames(inference_ess), c("prior_prob", "post_prob", "inclusion_BF", "ESS"))
  expect_equal(attr(inference_ess, "type"), c("prior_prob", "post_prob", "inclusion_BF", "ESS"))
  expect_null(attr(inference_ess, "footnotes"))

  inference_all <- runjags_inference_empty_table(BF_diagnostic_columns = TRUE)
  expect_equal(colnames(inference_all), c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error_percent"))

  inference_none <- runjags_inference_empty_table(BF_diagnostic_columns = FALSE)
  expect_equal(colnames(inference_none), c("prior_prob", "post_prob", "inclusion_BF"))
})


test_that("indicator BF diagnostics use indicator MCSE and handle boundaries", {

  indicator <- list(c(rep(0, 50), rep(1, 50)), c(rep(0, 50), rep(1, 50)))
  diagnostics <- .indicator_BF_diagnostics(indicator, prior_prob = 0.5, BF = 1)

  expect_equal(diagnostics$post_prob, 0.5)
  expect_equal(diagnostics$visits, 100)
  expect_equal(diagnostics$n_samples, 200)
  expect_true(is.finite(diagnostics$MCMC_error))
  expect_true(is.finite(diagnostics$ESS))
  expect_equal(diagnostics$BF_error_percent, 100 * diagnostics$MCMC_error / (diagnostics$post_prob * (1 - diagnostics$post_prob)))
  expect_gt(diagnostics$BF_error_percent, 0)

  boundary <- .indicator_BF_diagnostics(list(rep(0, 20), rep(0, 20)), prior_prob = 0.5, BF = 0)
  expect_true(is.na(boundary$BF_error_percent))
  expect_null(.indicator_BF_warnings("x", boundary))
  boundary_reporting <- .indicator_BF_reporting_value(BF = 0, post_prob = 0, prior_prob = 0.5, n_samples = 40)
  expect_equal(boundary_reporting$value, 1 / 39)
  expect_equal(boundary_reporting$operator, "<")

  boundary_one <- .indicator_BF_diagnostics(list(rep(1, 20), rep(1, 20)), prior_prob = 0.5, BF = Inf)
  expect_true(is.na(boundary_one$BF_error_percent))
  expect_null(.indicator_BF_warnings("x", boundary_one))
  boundary_one_reporting <- .indicator_BF_reporting_value(BF = Inf, post_prob = 1, prior_prob = 0.5, n_samples = 40)
  expect_equal(boundary_one_reporting$value, 39)
  expect_equal(boundary_one_reporting$operator, ">")

  prior_boundary <- .indicator_BF_diagnostics(indicator, prior_prob = 1, BF = 0)
  expect_true(is.na(prior_boundary$BF_error_percent))

})


.runjags_inference_fit_for_test <- function(indicator, prior_prob = 0.5) {

  samples <- coda::mcmc(
    cbind(theta = indicator, theta_indicator = indicator),
    start = 1,
    end   = length(indicator),
    thin  = 1
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(samples), sample = length(indicator)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(prior_prob))
    )
  )

  fit
}


test_that("runjags_inference_table reports BEAST/BSSVS-style bounds at indicator boundaries", {

  skip_if_not_installed("runjags")

  all_included <- runjags_inference_table(.runjags_inference_fit_for_test(rep(1, 10)))
  expect_equal(as.numeric(all_included["theta", "inclusion_BF"]), 9)
  expect_equal(attr(all_included[["inclusion_BF"]], "bound_operator"), ">")
  expect_equal(attr(all_included[["inclusion_BF"]][1], "bound_operator"), ">")
  expect_true(is.na(as.numeric(runjags_inference_table(
    .runjags_inference_fit_for_test(rep(1, 10)),
    BF_diagnostic_columns = "BF_error_percent"
  )["theta", "BF_error_percent"])))
  expect_null(attr(all_included, "footnotes"))
  expect_null(attr(all_included, "warnings"))
  expect_true(any(grepl(">9.000", capture_output_lines(all_included, print = TRUE), fixed = TRUE)))
  expect_match(
    .interpret.BF(all_included[["inclusion_BF"]][1], "theta inclusion", attr(all_included[["inclusion_BF"]], "name")),
    "Inclusion BF > 9.00",
    fixed = TRUE
  )

  all_excluded <- runjags_inference_table(.runjags_inference_fit_for_test(rep(0, 10)))
  expect_equal(as.numeric(all_excluded["theta", "inclusion_BF"]), 1 / 9)
  expect_equal(attr(all_excluded[["inclusion_BF"]], "bound_operator"), "<")
  expect_equal(attr(all_excluded[["inclusion_BF"]][1], "bound_operator"), "<")
  expect_null(attr(all_excluded, "footnotes"))
  expect_null(attr(all_excluded, "warnings"))
  expect_true(any(grepl("<0.111", capture_output_lines(all_excluded, print = TRUE), fixed = TRUE)))

  log_exclusion <- update(all_included, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(log_exclusion["theta", "inclusion_BF"]), log(1 / 9), tolerance = 1e-12)
  expect_equal(attr(log_exclusion[["inclusion_BF"]], "bound_operator"), "<")

  exclusion <- update(all_excluded, BF01 = TRUE)
  expect_equal(as.numeric(exclusion["theta", "inclusion_BF"]), 9)
  expect_equal(attr(exclusion[["inclusion_BF"]], "bound_operator"), ">")

  prior_quarter_included <- runjags_inference_table(
    .runjags_inference_fit_for_test(rep(1, 10), prior_prob = 0.25),
    BF_diagnostic_columns = "BF_error_percent"
  )
  expect_equal(as.numeric(prior_quarter_included["theta", "prior_prob"]), 0.25)
  expect_equal(as.numeric(prior_quarter_included["theta", "post_prob"]), 1)
  expect_equal(as.numeric(prior_quarter_included["theta", "inclusion_BF"]), 27)
  expect_equal(attr(prior_quarter_included[["inclusion_BF"]], "bound_operator"), ">")
  expect_true(is.na(as.numeric(prior_quarter_included["theta", "BF_error_percent"])))

  prior_quarter_excluded <- runjags_inference_table(
    .runjags_inference_fit_for_test(rep(0, 10), prior_prob = 0.25),
    BF_diagnostic_columns = "BF_error_percent"
  )
  expect_equal(as.numeric(prior_quarter_excluded["theta", "prior_prob"]), 0.25)
  expect_equal(as.numeric(prior_quarter_excluded["theta", "post_prob"]), 0)
  expect_equal(as.numeric(prior_quarter_excluded["theta", "inclusion_BF"]), 1 / 3)
  expect_equal(attr(prior_quarter_excluded[["inclusion_BF"]], "bound_operator"), "<")
  expect_true(is.na(as.numeric(prior_quarter_excluded["theta", "BF_error_percent"])))

})


test_that("runjags_inference_table reports undefined BF at prior probability boundaries", {

  skip_if_not_installed("runjags")

  prior_excludes <- runjags_inference_table(
    .runjags_inference_fit_for_test(rep(0, 10), prior_prob = 0),
    BF_diagnostic_columns = "BF_error_percent"
  )
  expect_equal(as.numeric(prior_excludes["theta", "prior_prob"]), 0)
  expect_equal(as.numeric(prior_excludes["theta", "post_prob"]), 0)
  expect_true(is.na(as.numeric(prior_excludes["theta", "inclusion_BF"])))
  expect_true(is.na(attr(prior_excludes[["inclusion_BF"]], "bound_operator")))
  expect_true(is.na(as.numeric(prior_excludes["theta", "BF_error_percent"])))

  prior_includes <- runjags_inference_table(
    .runjags_inference_fit_for_test(rep(1, 10), prior_prob = 1),
    BF_diagnostic_columns = "BF_error_percent"
  )
  expect_equal(as.numeric(prior_includes["theta", "prior_prob"]), 1)
  expect_equal(as.numeric(prior_includes["theta", "post_prob"]), 1)
  expect_true(is.na(as.numeric(prior_includes["theta", "inclusion_BF"])))
  expect_true(is.na(attr(prior_includes[["inclusion_BF"]], "bound_operator")))
  expect_true(is.na(as.numeric(prior_includes["theta", "BF_error_percent"])))

})


test_that("BayesTools table row subsetting keeps printed diagnostics local", {

  table <- data.frame(
    prior_prob       = c(0.5, 0.5),
    post_prob        = c(0.99, 0.4),
    inclusion_BF     = format_BF(c(99, 2/3), inclusion = TRUE),
    BF_error_percent = c(25.253, 1.2)
  )
  rownames(table) <- c("theta", "beta")
  class(table) <- c("BayesTools_table", "BayesTools_runjags_inference", class(table))
  attr(table, "type") <- c("prior_prob", "post_prob", "inclusion_BF", "BF_error")
  attr(table, "rownames") <- TRUE
  attr(table, "parameters") <- c("theta", "beta")
  attr(table, "warnings") <- stats::setNames(
    c(
      "Bayes factor MC error for theta is based on only 10 posterior samples from the less frequent model.",
      "Bayes factor MC error for beta is based on only 15 posterior samples from the less frequent model."
    ),
    c("theta", "beta")
  )
  attr(table[["BF_error_percent"]], "name") <- "error%(Inclusion BF)"

  subset <- table["theta", , drop = FALSE]

  expect_equal(attr(subset[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
  expect_equal(names(attr(subset, "warnings")), "theta")
  output <- capture_output_lines(subset, print = TRUE, width = 150)
  expect_true(any(grepl("error%(Inclusion BF)", output, fixed = TRUE)))
  expect_true(any(grepl("theta", output, fixed = TRUE)))
  expect_false(any(grepl("beta", output, fixed = TRUE)))
})


test_that("update preserves relative BF MC error percentage across BF scales", {

  table <- data.frame(
    prior_prob       = 0.5,
    post_prob        = 0.75,
    inclusion_BF     = format_BF(3, inclusion = TRUE),
    ESS              = 100,
    MCMC_error       = 0.02,
    BF_error_percent = 8
  )
  class(table) <- c("BayesTools_table", class(table))
  attr(table, "type") <- c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error")
  attr(table[["BF_error_percent"]], "name") <- "error%(Inclusion BF)"

  log_table <- update(table, logBF = TRUE)
  expect_equal(as.numeric(log_table[["inclusion_BF"]]), log(3))
  expect_equal(as.numeric(log_table[["BF_error_percent"]]), 8)
  expect_equal(attr(log_table[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  default_again <- update(log_table)
  expect_equal(as.numeric(default_again[["inclusion_BF"]]), 3)
  expect_equal(as.numeric(default_again[["BF_error_percent"]]), 8)
  expect_equal(attr(default_again[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  BF01_table <- update(table, BF01 = TRUE)
  expect_equal(as.numeric(BF01_table[["inclusion_BF"]]), 1/3)
  expect_equal(as.numeric(BF01_table[["BF_error_percent"]]), 8)
  expect_equal(attr(BF01_table[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  BF01_from_log <- update(log_table, BF01 = TRUE)
  expect_equal(as.numeric(BF01_from_log[["inclusion_BF"]]), 1/3)
  expect_equal(as.numeric(BF01_from_log[["BF_error_percent"]]), 8)
  expect_equal(attr(BF01_from_log[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  log_BF01_table <- update(table, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(log_BF01_table[["inclusion_BF"]]), log(1/3))
  expect_equal(as.numeric(log_BF01_table[["BF_error_percent"]]), 8)
  expect_equal(attr(log_BF01_table[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  default_from_log_BF01 <- update(log_BF01_table)
  expect_equal(as.numeric(default_from_log_BF01[["inclusion_BF"]]), 3)
  expect_equal(as.numeric(default_from_log_BF01[["BF_error_percent"]]), 8)
  expect_equal(attr(default_from_log_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

})


.runjags_table_fit_for_test <- function(samples) {

  samples <- coda::mcmc(samples)
  fit <- structure(
    list(
      mcmc = coda::mcmc.list(samples),
      sample = nrow(samples),
      summary.pars = list(mutate = NULL),
      monitor = colnames(samples)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- list(
    theta = prior("normal", list(0, 1)),
    beta  = prior("mnormal", list(mean = 0, sd = 1, K = 1))
  )
  fit
}

.expect_runjags_estimate_values_for_test <- function(table, samples, probs) {

  expected <- cbind(
    Mean = colMeans(samples),
    SD   = apply(samples, 2, stats::sd)
  )
  for(prob in probs){
    expected <- cbind(expected, apply(samples, 2, stats::quantile, probs = prob))
    colnames(expected)[ncol(expected)] <- prob
  }

  expect_equal(rownames(table), colnames(samples))
  expect_equal(colnames(table), colnames(expected))
  expect_equal(unname(as.matrix(table)), unname(expected), tolerance = 1e-12)
  expect_equal(attr(table, "type"), rep("estimate", ncol(table)))
  expect_equal(attr(table, "parameters"), colnames(samples))
  expect_true(attr(table, "rownames"))
}

test_that("runjags_estimates_table reports exact sample summaries without snapshots", {

  skip_if_not_installed("runjags")

  posterior <- matrix(
    c(
      1, 2, 3, 4,
      10, 20, 30, 40
    ),
    nrow = 4,
    dimnames = list(NULL, c("theta", "beta[1]"))
  )
  fit <- .runjags_table_fit_for_test(posterior)

  estimates <- suppressWarnings(runjags_estimates_table(
    fit,
    probs = c(0.25, 0.75),
    remove_diagnostics = TRUE
  ))

  .expect_runjags_estimate_values_for_test(estimates, posterior, c(0.25, 0.75))

  transformed_posterior <- posterior
  transformed_posterior[, "theta"] <- exp(transformed_posterior[, "theta"])

  transformed <- suppressWarnings(runjags_estimates_table(
    fit,
    transformations = list(theta = list(fun = exp)),
    probs = c(0.25, 0.75),
    remove_diagnostics = TRUE
  ))

  .expect_runjags_estimate_values_for_test(transformed, transformed_posterior, c(0.25, 0.75))
})

test_that("raw random-effect correlation aliases include logit-scale rho columns", {

  samples <- matrix(
    seq_len(5L * 5L),
    nrow = 5L,
    dimnames = list(
      NULL,
      c(
        "mu_intercept",
        "mu__xREx__id_sd",
        "mu__xREx__id_rho",
        "mu__xREx__id_rho_z",
        "mu__xREx__id_rho_logit"
      )
    )
  )
  formula_design <- list(
    mu = structure(
      list(
        parameter = "mu",
        random_effects = list(list(parameter_stem = "mu__xREx__id"))
      ),
      class = c("BayesTools_formula_design", "list")
    )
  )

  removed <- BayesTools:::.bt_JAGS_estimates_filter_raw_random_columns(
    model_samples = samples,
    prior_list = list(),
    formula_design = formula_design,
    remove_parameters = "random_correlation"
  )
  expect_true("mu__xREx__id_sd" %in% colnames(removed))
  expect_false("mu__xREx__id_rho" %in% colnames(removed))
  expect_false("mu__xREx__id_rho_z" %in% colnames(removed))
  expect_false("mu__xREx__id_rho_logit" %in% colnames(removed))

  kept <- BayesTools:::.bt_JAGS_estimates_filter_raw_random_columns(
    model_samples = samples,
    prior_list = list(),
    formula_design = formula_design,
    keep_parameters = "random_correlation"
  )
  expect_true("mu_intercept" %in% colnames(kept))
  expect_false("mu__xREx__id_sd" %in% colnames(kept))
  expect_true("mu__xREx__id_rho" %in% colnames(kept))
  expect_true("mu__xREx__id_rho_z" %in% colnames(kept))
  expect_true("mu__xREx__id_rho_logit" %in% colnames(kept))
})


test_that("runjags summary helper keeps only requested diagnostics", {

  posterior <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      10, 20, 30, 40, 50, 60
    ),
    nrow = 6,
    dimnames = list(NULL, c("theta", "beta[1]"))
  )

  summary <- .runjags_summary_fast(
    posterior,
    n_samples = 3,
    n_chains = 2,
    conditional = FALSE,
    probs = 0.5,
    diagnostic_columns = c("ESS", "R_hat")
  )

  expect_equal(colnames(summary), c("Mean", "SD", "0.5", "ESS", "R_hat"))
  expect_equal(attr(summary, "type"), NULL)
  expect_true(all(is.finite(summary[["ESS"]])))
  expect_true(all(is.na(summary[["R_hat"]]) | is.finite(summary[["R_hat"]])))

  no_diagnostics <- .runjags_summary_fast(
    posterior,
    n_samples = 3,
    n_chains = 2,
    conditional = FALSE,
    probs = 0.5,
    diagnostic_columns = FALSE
  )
  expect_equal(colnames(no_diagnostics), c("Mean", "SD", "0.5"))
})

test_that("ensemble_estimates_table preserves parameter and matrix row semantics", {

  samples <- list(
    theta = c(1, 2, 3, 4),
    beta = matrix(
      c(1, 2, 3, 4, 10, 20, 30, 40),
      ncol = 2,
      dimnames = list(NULL, c("beta[1]", "beta[2]"))
    )
  )

  estimates <- ensemble_estimates_table(
    samples = samples,
    parameters = c("theta", "beta"),
    probs = c(0.25, 0.75)
  )

  expect_s3_class(estimates, "BayesTools_table")
  expect_equal(rownames(estimates), c("theta", "beta[1]", "beta[2]"))
  expect_equal(colnames(estimates), c("Mean", "Median", "0.25", "0.75"))
  expect_equal(
    unname(as.matrix(estimates)),
    matrix(
      c(
        2.5,  2.5,  1.75,  3.25,
        2.5,  2.5,  1.75,  3.25,
        25.0, 25.0, 17.50, 32.50
      ),
      nrow = 3,
      byrow = TRUE
    ),
    tolerance = 1e-12
  )
  expect_equal(attr(estimates, "type"), rep("estimate", 4))
  expect_true(attr(estimates, "rownames"))
})

.mock_ensemble_table_model <- function(model_number, prior_list, prior_prob,
                                       post_prob, marglik, inclusion_BF,
                                       fit_summary) {
  fit <- list()
  attr(fit, "prior_list") <- prior_list

  list(
    fit = fit,
    fit_summary = fit_summary,
    inference = list(
      m_number = model_number,
      marglik = marglik,
      prior_prob = prior_prob,
      post_prob = post_prob,
      inclusion_BF = inclusion_BF
    )
  )
}

.mock_ensemble_fit_summary <- function(MCMC_error, MCMC_SD_error, ESS, R_hat) {
  out <- data.frame(
    MCMC_error = MCMC_error,
    MCMC_SD_error = MCMC_SD_error,
    ESS = ESS,
    R_hat = R_hat
  )
  class(out) <- c("BayesTools_table", "BayesTools_runjags_summary", "data.frame")
  attr(out, "type") <- c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
  out
}

test_that("ensemble inference table reports exact inclusion probability algebra", {

  theta <- structure(
    list(
      prior_probs = c(null = 0.25, alt = 0.75),
      post_probs = c(null = 0.10, alt = 0.90),
      BF = 3
    ),
    is_null = c(TRUE, FALSE),
    parameter_name = "theta"
  )
  beta <- structure(
    list(
      prior_probs = c(null = 0.20, alt1 = 0.30, alt2 = 0.50),
      post_probs = c(null = 0.50, alt1 = 0.25, alt2 = 0.25),
      BF = 0.25
    ),
    is_null = c(TRUE, FALSE, FALSE),
    parameter_name = "beta"
  )
  inference <- list(theta = theta, beta = beta)
  attr(inference, "conditional") <- FALSE

  table <- ensemble_inference_table(inference, c("theta", "beta"))

  expect_s3_class(table, "BayesTools_table")
  expect_equal(rownames(table), c("theta", "beta"))
  expect_equal(table$models, c(1, 2))
  expect_equal(table$prior_prob, c(0.75, 0.80), tolerance = 1e-12)
  expect_equal(table$post_prob, c(0.90, 0.50), tolerance = 1e-12)
  expect_equal(as.numeric(table$inclusion_BF), c(3, 0.25), tolerance = 1e-12)
  expect_equal(attr(table, "n_models"), c(2L, 3L))
  expect_equal(attr(table, "type"), c("n_models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_true(attr(table, "rownames"))

  log_BF01 <- ensemble_inference_table(inference, c("theta", "beta"), logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(log_BF01$inclusion_BF), log(c(1 / 3, 4)), tolerance = 1e-12)
  expect_equal(attr(log_BF01$inclusion_BF, "name"), "log(Exclusion BF)")

  attr(inference, "conditional") <- TRUE
  expect_error(
    ensemble_inference_table(inference, c("theta", "beta")),
    "cannot be 'conditional'",
    fixed = TRUE
  )
})

test_that("ensemble summary and diagnostics tables preserve model alignment", {

  models <- list(
    .mock_ensemble_table_model(
      model_number = 3L,
      prior_list = list(theta = prior("point", list(0)), beta = prior("normal", list(0, 1))),
      prior_prob = 0.20,
      post_prob = 0.10,
      marglik = log(2),
      inclusion_BF = 0.5,
      fit_summary = .mock_ensemble_fit_summary(
        MCMC_error = c(0.01, 0.03),
        MCMC_SD_error = c(0.10, 0.20),
        ESS = c(500, 300),
        R_hat = c(1.01, 1.03)
      )
    ),
    .mock_ensemble_table_model(
      model_number = 1L,
      prior_list = list(theta = prior("normal", list(1, 2))),
      prior_prob = 0.30,
      post_prob = 0.70,
      marglik = log(4),
      inclusion_BF = 4,
      fit_summary = .mock_ensemble_fit_summary(
        MCMC_error = c(0.04, NA_real_),
        MCMC_SD_error = c(0.30, NA_real_),
        ESS = c(250, NA_real_),
        R_hat = c(1.05, NA_real_)
      )
    ),
    .mock_ensemble_table_model(
      model_number = 2L,
      prior_list = list(beta = prior("point", list(0))),
      prior_prob = 0.50,
      post_prob = 0.20,
      marglik = log(3),
      inclusion_BF = 1.5,
      fit_summary = .mock_ensemble_fit_summary(
        MCMC_error = c(NA_real_, NA_real_),
        MCMC_SD_error = c(NA_real_, NA_real_),
        ESS = c(NA_real_, NA_real_),
        R_hat = c(NA_real_, NA_real_)
      )
    )
  )

  summary <- ensemble_summary_table(models, c(theta = "theta", beta = "beta"))
  summary_with_spikes <- ensemble_summary_table(
    models,
    c(theta = "theta", beta = "beta"),
    remove_spike_0 = FALSE
  )
  diagnostics <- ensemble_diagnostics_table(models, c(theta = "theta", beta = "beta"))

  expect_equal(summary$Model, c(3, 1, 2))
  expect_equal(summary$prior_prob, c(0.20, 0.30, 0.50), tolerance = 1e-12)
  expect_equal(summary$marglik, log(c(2, 4, 3)), tolerance = 1e-12)
  expect_equal(summary$post_prob, c(0.10, 0.70, 0.20), tolerance = 1e-12)
  expect_equal(as.numeric(summary$inclusion_BF), c(0.5, 4, 1.5), tolerance = 1e-12)
  expect_equal(summary$theta[1], "")
  expect_true(nzchar(summary$theta[2]))
  expect_equal(summary$theta[3], "")
  expect_true(nzchar(summary$beta[1]))
  expect_equal(summary$beta[2], "")
  expect_equal(summary$beta[3], "")
  expect_true(nzchar(summary_with_spikes$theta[1]))
  expect_true(nzchar(summary_with_spikes$beta[3]))
  expect_equal(
    attr(summary, "type"),
    c("integer", "prior", "prior", "prior_prob", "marglik", "post_prob", "inclusion_BF")
  )
  expect_false(attr(summary, "rownames"))

  expect_equal(diagnostics$Model, c(3, 1, 2))
  expect_equal(diagnostics$max_MCMC_error, c(0.03, 0.04, NA_real_), tolerance = 1e-12)
  expect_equal(diagnostics$max_MCMC_SD_error, c(0.20, 0.30, NA_real_), tolerance = 1e-12)
  expect_equal(diagnostics$min_ESS, c(300, 250, NA_real_), tolerance = 1e-12)
  expect_equal(diagnostics$max_R_hat, c(1.03, 1.05, NA_real_), tolerance = 1e-12)
  expect_equal(
    attr(diagnostics, "type"),
    c("integer", "prior", "prior", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat")
  )

  models[[1]]$fit_summary <- models[[1]]$fit_summary[, c("ESS", "R_hat"), drop = FALSE]
  attr(models[[1]]$fit_summary, "type") <- c("ESS", "R_hat")

  missing_diagnostics <- ensemble_diagnostics_table(models, c(theta = "theta", beta = "beta"))
  expect_equal(missing_diagnostics$max_MCMC_error[1], NA_real_)
  expect_equal(missing_diagnostics$max_MCMC_SD_error[1], NA_real_)
  expect_equal(missing_diagnostics$min_ESS[1], 300, tolerance = 1e-12)
  expect_equal(missing_diagnostics$max_R_hat[1], 1.03, tolerance = 1e-12)
})


test_that("add_column works correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Create a proper BayesTools_table with 3 columns (needed for middle position test)
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Add column at end (default)
  result1 <- add_column(test_data, "CI_lower", c(-0.5, 0.8))
  expect_equal(ncol(result1), 4)
  expect_true("CI_lower" %in% names(result1))
  expect_s3_class(result1, "BayesTools_table")
  expect_equal(attr(result1, "type"), c("estimate", "estimate", "estimate", "estimate"))
  test_reference_table(result1, "add_column_end.txt")

  # Add column at specific position
  result2 <- add_column(test_data, "CI_lower", c(-0.5, 0.8), column_position = 2)
  expect_equal(ncol(result2), 4)
  expect_equal(names(result2)[2], "CI_lower")
  test_reference_table(result2, "add_column_position2.txt")

  # Add column at position 1
  result3 <- add_column(test_data, "ID", c(1, 2), column_position = 1)
  expect_equal(ncol(result3), 4)
  expect_equal(names(result3)[1], "ID")
  expect_equal(attr(result3, "type")[1], "integer")
  test_reference_table(result3, "add_column_position1.txt")

  # With specified column type
  result4 <- add_column(test_data, "Prob", c(0.5, 0.8), column_type = "probability")
  expect_equal(ncol(result4), 4)
  expect_equal(attr(result4, "type")[4], "probability")
  test_reference_table(result4, "add_column_probability.txt")

  # Add column with string values (must specify column_type)
  result5 <- add_column(test_data, "Category", c("A", "B"), column_type = "string")
  expect_equal(ncol(result5), 4)
  expect_true("Category" %in% names(result5))
  expect_equal(attr(result5, "type")[4], "string")
  test_reference_table(result5, "add_column_string.txt")

})


test_that("add_column input validation works", {

  # Create a proper BayesTools_table with data (3 columns)
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Wrong table class
  expect_error(add_column(data.frame(a = 1), "b", 2), "must be of class 'BayesTools_table'")

  # Wrong column_values length
  expect_error(add_column(test_data, "new", c(1, 2, 3)), "must be a vector of the same length")

})


test_that("remove_column works correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Create a proper BayesTools_table with data
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1),
    SD = c(0.1, 0.2)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Remove last column (default)
  result0 <- remove_column(test_data)
  expect_equal(ncol(result0), 2)
  expect_false("SD" %in% names(result0))
  expect_s3_class(result0, "BayesTools_table")
  expect_equal(attr(result0, "type"), c("estimate", "estimate"))
  test_reference_table(result0, "remove_column_last.txt")

  # Remove by position
  result2 <- remove_column(test_data, column_position = 2)
  expect_equal(ncol(result2), 2)
  expect_false("Median" %in% names(result2))
  test_reference_table(result2, "remove_column_position2.txt")

})


test_that("remove_column input validation works", {

  # Create a proper BayesTools_table with data
  test_data <- data.frame(
    Mean = c(0.5, 1.2),
    Median = c(0.4, 1.1)
  )
  rownames(test_data) <- c("mu", "sigma")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("estimate", "estimate")
  attr(test_data, "rownames") <- TRUE

  # Wrong table class
  expect_error(remove_column(data.frame(a = 1)), "must be of class 'BayesTools_table'")

  # Invalid column position
  expect_error(remove_column(test_data, column_position = 10), "'column_position'")

})


test_that("column helpers preserve BF columns and table metadata", {

  test_data <- data.frame(
    inclusion_BF = format_BF(c(2, 4), inclusion = TRUE),
    Mean         = c(1.25, 2.50),
    check.names  = FALSE
  )
  rownames(test_data) <- c("theta", "beta")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- c("inclusion_BF", "estimate")
  attr(test_data, "rownames") <- TRUE
  attr(test_data, "title") <- "base title"
  attr(test_data, "footnotes") <- "base footnote"
  attr(test_data, "warnings") <- "base warning"
  attr(test_data, "n_models") <- c(3L, 3L)

  added <- add_column(
    test_data,
    "label",
    c("A", "B"),
    column_position = 2,
    column_type = "string"
  )

  expect_true(is.data.frame(added))
  expect_false("dim" %in% names(attributes(added)))
  expect_false("dimnames" %in% names(attributes(added)))
  expect_s3_class(added, "BayesTools_table")
  expect_equal(colnames(added), c("inclusion_BF", "label", "Mean"))
  expect_equal(attr(added, "type"), c("inclusion_BF", "string", "estimate"))
  expect_equal(attr(added, "title"), "base title")
  expect_equal(attr(added, "footnotes"), "base footnote")
  expect_equal(attr(added, "warnings"), "base warning")
  expect_equal(attr(added, "n_models"), c(3L, 3L))
  expect_equal(as.numeric(added[["inclusion_BF"]]), c(2, 4))
  expect_equal(attr(added[["inclusion_BF"]], "name"), "Inclusion BF")
  expect_false(attr(added[["inclusion_BF"]], "logBF"))
  expect_false(attr(added[["inclusion_BF"]], "BF01"))

  removed <- remove_column(added, column_position = 2)

  expect_true(is.data.frame(removed))
  expect_equal(colnames(removed), colnames(test_data))
  expect_equal(attr(removed, "type"), attr(test_data, "type"))
  expect_equal(attr(removed, "title"), attr(test_data, "title"))
  expect_equal(attr(removed, "footnotes"), attr(test_data, "footnotes"))
  expect_equal(attr(removed, "warnings"), attr(test_data, "warnings"))
  expect_equal(as.numeric(removed[["inclusion_BF"]]), c(2, 4))
  expect_equal(attr(removed[["inclusion_BF"]], "name"), "Inclusion BF")

  updated <- update(
    removed,
    title = "updated title",
    footnotes = "extra footnote",
    warnings = "extra warning",
    remove_parameters = "beta",
    logBF = TRUE
  )

  expect_equal(rownames(updated), "theta")
  expect_equal(as.numeric(updated[["inclusion_BF"]]), log(2), tolerance = 1e-12)
  expect_equal(attr(updated[["inclusion_BF"]], "name"), "log(Inclusion BF)")
  expect_true(attr(updated[["inclusion_BF"]], "logBF"))
  expect_equal(attr(updated, "title"), "updated title")
  expect_equal(attr(updated, "footnotes"), c("base footnote", "extra footnote"))
  expect_equal(attr(updated, "warnings"), c("base warning", "extra warning"))
  expect_equal(attr(updated, "type"), attr(test_data, "type"))
})


test_that("remove_column keeps one-column tables as empty BayesTools tables", {

  test_data <- data.frame(Mean = c(1.25, 2.50))
  rownames(test_data) <- c("theta", "beta")
  class(test_data) <- c("BayesTools_table", "data.frame")
  attr(test_data, "type") <- "estimate"
  attr(test_data, "rownames") <- TRUE
  attr(test_data, "title") <- "single column"

  removed <- remove_column(test_data)

  expect_true(is.data.frame(removed))
  expect_s3_class(removed, "BayesTools_table")
  expect_equal(dim(removed), c(2L, 0L))
  expect_equal(attr(removed, "type"), character())
  expect_equal(attr(removed, "title"), "single column")
  expect_true(attr(removed, "rownames"))
})


test_that("ensemble_estimates_empty_table works correctly", {

  empty_table <- ensemble_estimates_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_estimates_empty.txt")

})


test_that("ensemble_inference_empty_table works correctly", {

  empty_table <- ensemble_inference_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_inference_empty.txt")

})


test_that("ensemble_summary_empty_table works correctly", {

  empty_table <- ensemble_summary_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_summary_empty.txt")

})


test_that("ensemble_diagnostics_empty_table works correctly", {

  empty_table <- ensemble_diagnostics_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "ensemble_diagnostics_empty.txt")

})
