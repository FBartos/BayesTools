# ============================================================================ #
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ============================================================================ #
if (!exists("GENERATE_REFERENCE_FILES")) {
  GENERATE_REFERENCE_FILES <- FALSE
}


test_files_dir <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
if (test_files_dir == "" || !dir.exists(test_files_dir)) {
  test_files_dir <- file.path(tempdir(), "BayesTools_test_files")
}

# Setup directory for saving fitted models
temp_fits_dir    <- file.path(test_files_dir, "fits")
temp_marglik_dir <- file.path(test_files_dir, "margliks")
temp_temp_dir    <- file.path(test_files_dir, "temp")

if (!dir.exists(temp_fits_dir))    dir.create(temp_fits_dir,    showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_marglik_dir)) dir.create(temp_marglik_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_temp_dir))    dir.create(temp_temp_dir,    showWarnings = FALSE, recursive = TRUE)

# Set environment variable so other test files can locate pre-fitted models
Sys.setenv(BAYESTOOLS_TEST_FILES_DIR = test_files_dir)

# Use skip_if_no_fits() for tests that need pre-fitted models.

# ============================================================================ #
# HELPER FUNCTIONS: Reference File Testing
# ============================================================================ #
require("runjags")

# Process reference file: save if GENERATE_REFERENCE_FILES=TRUE, test otherwise
test_reference_table <- function(table, filename, info_msg = NULL,
                                 print_dir = REFERENCE_DIR) {
  if (GENERATE_REFERENCE_FILES) {
    # Save mode
    if (!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    writeLines(capture_output_lines(table, print = TRUE, width = 150),
               file.path(print_dir, filename))
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      actual_output   <- capture_output_lines(table, print = TRUE, width = 150)
      expect_equal(actual_output, expected_output, info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

parse_numeric_table_lines <- function(lines) {

  number_pattern <- "([-+]?Inf|NA|NaN|[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
  row_pattern <- paste0(
    "^\\s*(.*?)\\s+",
    number_pattern, "\\s+",
    number_pattern, "\\s+",
    number_pattern, "\\s+",
    number_pattern, "\\s+",
    number_pattern, "\\s+",
    number_pattern, "\\s*$"
  )

  matches <- regexec(row_pattern, lines, perl = TRUE)
  parsed  <- regmatches(lines, matches)
  is_row  <- lengths(parsed) == 8

  if (!any(is_row)) {
    return(list(
      labels = character(),
      values = matrix(numeric(), nrow = 0, ncol = 6),
      is_row = is_row
    ))
  }

  row_data <- do.call(rbind, lapply(parsed[is_row], function(x) x[-1]))
  values   <- matrix(as.numeric(row_data[, -1]), nrow = nrow(row_data))
  colnames(values) <- c("Mean", "SD", "0.025", "0.5", "0.975", "Inclusion BF")

  list(
    labels = row_data[, 1],
    values = values,
    is_row = is_row
  )
}

test_reference_table_numeric <- function(table, filename, tolerance = 1e-2,
                                         info_msg = NULL, print_dir = REFERENCE_DIR) {
  if (GENERATE_REFERENCE_FILES) {
    test_reference_table(table, filename, info_msg = info_msg, print_dir = print_dir)
  } else {
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      actual_output   <- capture_output_lines(table, print = TRUE, width = 150)

      expected_table <- parse_numeric_table_lines(expected_output)
      actual_table   <- parse_numeric_table_lines(actual_output)

      expect_equal(actual_table$labels, expected_table$labels, info = info_msg)
      expect_equal(actual_table$values, expected_table$values, tolerance = tolerance, info = info_msg)
      expect_equal(actual_output[!actual_table$is_row], expected_output[!expected_table$is_row], info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

test_reference_text <- function(text, filename, info_msg = NULL,
                                print_dir = REFERENCE_DIR) {
  if (GENERATE_REFERENCE_FILES) {
    # Save mode
    if (!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    writeLines(text, file.path(print_dir, filename))
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      expected_output <- paste0(expected_output, collapse = "\n")
      expect_equal(text, expected_output, info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

# Skip if pre-fitted models are not available
skip_if_no_fits <- function() {
  model_registry_file <- file.path(test_files_dir, "model_registry.RDS")
  if (!file.exists(model_registry_file)) {
    skip("Pre-fitted models not found. Run test-00-model-fits.R first.")
  }
}

# ============================================================================ #
# STANDARD TEST FIXTURES: Reusable Prior Definitions
# ============================================================================ #
# These fixtures reduce duplication across test files. Use these instead of
# creating new prior definitions when testing standard functionality.

# Standard log_posterior function for marginal likelihood tests
# Returns 0 (log of 1) for prior-only models
STANDARD_LOG_POSTERIOR <- function(parameters, data) {

  return(0)
}

if (isNamespaceLoaded("BayesTools")) {

  # Standard simple priors (commonly used across tests)
  STANDARD_PRIORS <- list(
    normal          = prior("normal", list(0, 1)),
    normal_trunc    = prior("normal", list(0, 1), list(0, Inf)),
    lognormal       = prior("lognormal", list(0, 0.5)),
    t               = prior("t", list(0, 0.5, 5)),

    cauchy          = prior("Cauchy", list(0, 1)),
    cauchy_trunc    = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    gamma           = prior("gamma", list(2, 1)),
    invgamma        = prior("invgamma", list(3, 2)),
    invgamma_trunc  = prior("invgamma", list(3, 2), list(1, 3)),
    exp             = prior("exp", list(1.5)),
    beta            = prior("beta", list(3, 2)),
    uniform         = prior("uniform", list(0, 1)),
    spike           = prior("spike", list(0)),
    PET             = prior_PET("normal", list(0, 1)),
    PEESE           = prior_PEESE("gamma", list(1, 1))
  )

  # Standard factor priors (for contrast testing)
  STANDARD_FACTOR_PRIORS <- list(
    orthonormal     = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal"),
    meandif         = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif"),
    treatment       = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment"),
    independent     = prior_factor("normal", list(mean = 0, sd = 1), contrast = "independent"),
    orth_cauchy     = prior_factor("mcauchy", list(location = 0, scale = 1), contrast = "orthonormal"),
    orth_spike      = prior_factor("point", list(0), contrast = "orthonormal")
  )

  # Complete prior collections for comprehensive testing
  ALL_SIMPLE_PRIORS <- list(
    p1    = prior("normal", list(0, 1)),
    p2    = prior("normal", list(0, 1), list(1, Inf)),
    p3    = prior("lognormal", list(0, 0.5)),
    p4    = prior("t", list(0, 0.5, 5)),
    p5    = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p6    = prior("gamma", list(2, 1)),
    p7    = prior("invgamma", list(3, 2), list(1, 3)),
    p8    = prior("exp", list(1.5)),
    p9    = prior("beta", list(3, 2)),
    p10   = prior("uniform", list(1, 5)),
    PET   = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
  )

  ALL_VECTOR_PRIORS <- list(
    mnormal = prior("mnormal", list(mean = 0, sd = 1, K = 3)),
    mcauchy = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    mt      = prior("mt", list(location = 2, scale = 0.5, df = 5, K = 2))
  )

  ALL_FACTOR_PRIORS <- list(
    orthonormal = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal"),
    meandif     = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif"),
    treatment   = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment"),
    independent = prior_factor("normal", list(mean = 0, sd = 1), contrast = "independent")
  )
}
# ============================================================================ #
# HELPER FUNCTIONS: Prior Distribution Testing
# ============================================================================ #

#' Test a prior distribution for consistency
#'
#' Validates that a prior distribution's rng, pdf, cdf, quant, mean, and sd

' functions work correctly and are mutually consistent.
#'
#' @param prior A prior object to test
#' @param skip_moments Logical; skip mean/sd validation (for distributions
#'   with undefined moments like Cauchy)
#' @return invisible(); used for visual regression testing
test_prior <- function(prior, skip_moments = FALSE) {
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  if (is.prior.discrete(prior)) {
    barplot(table(samples) / length(samples), main = print(prior, plot = TRUE),
            width = 1 / (max(samples) + 1), space = 0,
            xlim = c(-0.25, max(samples) + 0.25))
  } else if (is.prior.spike_and_slab(prior)) {
    xh <- hist(samples[samples != 0], breaks = 50, plot = FALSE)
    xh$density <- xh$density * mean(samples != 0)
    plot(xh, main = print(prior, plot = TRUE), freq = FALSE)
  } else {
    hist(samples, main = print(prior, plot = TRUE), breaks = 50, freq = FALSE)
  }
  # tests density function
  lines(prior, individual = TRUE)

  # tests quantile function
  if (!is.prior.spike_and_slab(prior) && !is.prior.mixture(prior)) {
    abline(v = quant(prior, 0.5), col = "blue", lwd = 2)
  }
  # tests that cdf(quant(x)) == x

  if (!is.prior.point(prior) && !is.prior.discrete(prior) &&
      !is.prior.spike_and_slab(prior) && !is.prior.mixture(prior)) {
    expect_equal(.25, cdf(prior, quant(prior, 0.25)), tolerance = 1e-4)
    expect_equal(.25, ccdf(prior, quant(prior, 0.75)), tolerance = 1e-4)
  }
  # test mean and sd functions
  if (!skip_moments) {
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples), sd(prior), tolerance = 1e-2)
  }
  return(invisible())
}

#' Test a weight function prior distribution
#'
#' Validates weight function priors with multiple components.
#'
#' @param prior A weight function prior object
#' @param skip_moments Logical; skip moment validation
#' @return invisible()
test_weightfunction <- function(prior, skip_moments = FALSE) {
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 10000)
  set.seed(1)
  densities <- density(prior, individual = TRUE)

  if (!all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))) {
    quantiles <- mquant(prior, 0.5)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
  par(mfcol = c(1, ncol(samples) - 1))

  for (i in 1:(ncol(samples) - 1)) {
    hist(samples[, i], main = print(prior, plot = TRUE), breaks = 50, freq = FALSE)
    lines(densities[[i]])
    if (!all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))) {
      abline(v = quantiles[i], col = "blue", lwd = 2)
    }
    if (!grepl("fixed", prior$distribution) &&
        !all(names(prior$parameters) %in% c("steps", "alpha1", "alpha2"))) {
      expect_equal(.25, mcdf(prior, mquant(prior, 0.25)[, i])[, i], tolerance = 1e-5)
      expect_equal(.25, mccdf(prior, mquant(prior, 0.75)[, i])[, i], tolerance = 1e-5)
    }
    if (!skip_moments) {
      expect_equal(apply(samples, 2, mean), mean(prior), tolerance = 1e-2)
      expect_equal(apply(samples, 2, sd), sd(prior), tolerance = 1e-2)
    }
  }
  return(invisible())
}

#' Test an orthonormal contrast prior
#'
#' Validates orthonormal factor priors.
#'
#' @param prior An orthonormal prior object
#' @param skip_moments Logical; skip moment validation
#' @return invisible()
test_orthonormal <- function(prior, skip_moments = FALSE) {
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  samples <- samples[abs(samples) < 10]
  hist(samples, main = print(prior, plot = TRUE), breaks = 50, freq = FALSE)
  # tests density function
  lines(prior, individual = TRUE)
  # tests quantile function
  abline(v = mquant(prior, 0.5), col = "blue", lwd = 2)
  # tests that mcdf(mquant(x)) == x
  if (!is.prior.point(prior)) {
    expect_equal(.25, mcdf(prior, mquant(prior, 0.25)), tolerance = 1e-5)
    expect_equal(.25, mccdf(prior, mquant(prior, 0.75)), tolerance = 1e-5)
  }
  # test mean and sd functions
  if (!skip_moments) {
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples), sd(prior), tolerance = 1e-2)
  }
  return(invisible())
}

#' Test a mean difference contrast prior
#'
#' Validates meandif factor priors.
#'
#' @param prior A meandif prior object
#' @param skip_moments Logical; skip moment validation
#' @return invisible()
test_meandif <- function(prior, skip_moments = FALSE) {
  set.seed(1)
  # tests rng and print function (for plot)
  samples <- rng(prior, 100000)
  samples <- samples[abs(samples) < 10]
  hist(samples, main = print(prior, plot = TRUE), breaks = 50, freq = FALSE)
  # tests density function
  lines(prior, individual = TRUE)
  # tests quantile function
  abline(v = mquant(prior, 0.5), col = "blue", lwd = 2)
  # tests that mcdf(mquant(x)) == x
  if (!is.prior.point(prior)) {
    expect_equal(.25, mcdf(prior, mquant(prior, 0.25)), tolerance = 1e-5)
    expect_equal(.25, mccdf(prior, mquant(prior, 0.75)), tolerance = 1e-5)
  }
  # test mean and sd functions
  if (!skip_moments) {
    expect_equal(mean(samples), mean(prior), tolerance = 1e-2)
    expect_equal(sd(samples), sd(prior), tolerance = 1e-2)
  }
  return(invisible())
}

# Helper function to save fitted models and register metadata
save_fit <- function(fit, name, marglik = NULL, simple_priors = FALSE, vector_priors = FALSE,
                     factor_priors = FALSE, pub_bias_priors = FALSE,
                     weightfunction_priors = FALSE, spike_and_slab_priors = FALSE,
                     mixture_priors = FALSE, formulas = FALSE,
                     random_effects = FALSE, interactions = FALSE,
                     expression_priors = FALSE, multi_formula = FALSE,
                     autofit = FALSE, parallel = FALSE, thinning = FALSE,
                     add_parameters = FALSE, note = "") {

  saveRDS(fit, file = file.path(temp_fits_dir, paste0(name, ".RDS")))

  # Save marglik if provided
  if (!is.null(marglik)) {
    saveRDS(marglik, file = file.path(temp_marglik_dir, paste0(name, ".RDS")))
  }

  # Return model metadata entry for registry
  list(
    fit = fit,
    marglik = marglik,
    registry_entry = data.frame(
      model_name = name,
      has_marglik = !is.null(marglik),
      simple_priors = simple_priors,
      vector_priors = vector_priors,
      factor_priors = factor_priors,
      pub_bias_priors = pub_bias_priors,
      weightfunction_priors = weightfunction_priors,
      spike_and_slab_priors = spike_and_slab_priors,
      mixture_priors = mixture_priors,
      formulas = formulas,
      random_effects = random_effects,
      interactions = interactions,
      expression_priors = expression_priors,
      multi_formula = multi_formula,
      autofit = autofit,
      parallel = parallel,
      thinning = thinning,
      add_parameters = add_parameters,
      note = note,
      stringsAsFactors = FALSE
    )
  )
}

# Skip model fitting if cached fits exist and ROBMA_TEST_SKIP_REFIT is TRUE
skip_refit_if_cached <- function(name) {
  # refitting settings
  skip_refit <- Sys.getenv("BAYESTOOLS_TEST_SKIP_REFIT")
  skip_refit <- skip_refit != "" && as.logical(skip_refit)

  # fitted indicator
  fitted_indicator <- file.exists(file.path(temp_temp_dir, paste0(name, ".txt")))

  if (skip_refit && fitted_indicator) {
    skip("Skipping model refitting: cached fits exist and BAYESTOOLS_TEST_SKIP_REFIT=TRUE.")
  }

  # tests are not going to be skipped -- add fits done indicator into `temp_temp_dir`
  file.create(file.path(temp_temp_dir, paste0(name, ".txt")))
}

# Clean cached fitted models and margliks
clean_cached_fits <- function(name) {

  if (!missing(name)) {
    # remove only the specific `name`` fitted indicator files side-effects from `temp_temp_dir`
    file.remove(file.path(temp_temp_dir, paste0(name, ".txt")))
  } else {
    # Remove all cached files from test directories
    unlink(temp_fits_dir,    recursive = TRUE)
    unlink(temp_marglik_dir, recursive = TRUE)
    unlink(temp_temp_dir,    recursive = TRUE)

    # Recreate empty directories
    dir.create(temp_fits_dir,    showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_marglik_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_temp_dir,    showWarnings = FALSE, recursive = TRUE)
  }

  return(invisible(TRUE))
}
