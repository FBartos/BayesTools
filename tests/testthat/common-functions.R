# ============================================================================ #
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ============================================================================ #
if (!exists("GENERATE_REFERENCE_FILES")) {
  GENERATE_REFERENCE_FILES <- FALSE
}


test_files_dir <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
if (test_files_dir == "") {
  test_files_dir <- file.path(tempdir(), "BayesTools_test_files")
}
if (!dir.exists(test_files_dir)) {
  dir.create(test_files_dir, showWarnings = FALSE, recursive = TRUE)
}
test_files_dir <- normalizePath(test_files_dir, winslash = "/", mustWork = TRUE)

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

  number_pattern <- "(?:<?[-+]?Inf|NA|NaN|<?[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
  row_pattern <- paste0(
    "^\\s*(.*?)\\s+",
    "(", number_pattern, "(?:\\s+", number_pattern, ")+)",
    "\\s*$"
  )

  matches <- regexec(row_pattern, lines, perl = TRUE)
  parsed  <- regmatches(lines, matches)
  is_row  <- lengths(parsed) == 3

  if (!any(is_row)) {
    return(list(
      labels = character(),
      values = matrix(numeric(), nrow = 0, ncol = 0),
      is_row = is_row
    ))
  }

  row_data <- do.call(rbind, lapply(parsed[is_row], function(x) x[-1]))
  row_values <- lapply(strsplit(trimws(row_data[, 2]), "\\s+"), function(tokens) {
    tokens <- sub("^<", "", tokens)
    tokens[tokens == "NA"] <- NA_character_
    suppressWarnings(as.numeric(tokens))
  })
  value_widths <- lengths(row_values)
  if (length(unique(value_widths)) != 1L) {
    stop("Parsed numeric table rows have inconsistent widths.", call. = FALSE)
  }
  values <- do.call(rbind, row_values)
  colnames(values) <- paste0("V", seq_len(ncol(values)))

  list(
    labels = row_data[, 1],
    values = values,
    is_row = is_row
  )
}

normalize_reference_non_table_lines <- function(lines) {
  gsub("\\s+", " ", trimws(lines))
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
      expect_equal(
        normalize_reference_non_table_lines(actual_output[!actual_table$is_row]),
        normalize_reference_non_table_lines(expected_output[!expected_table$is_row]),
        info = info_msg
      )
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

# Skip if pre-fitted models are not available in lightweight profiles. Explicit
# fixture-bearing profiles fail so missing or stale caches cannot be bypassed.
.test_cache_required_for_active_profile <- function() {
  if (exists("bayestools_test_profile_includes", mode = "function")) {
    return(bayestools_test_profile_includes(c("fixture", "fit", "visual-fixture")))
  }

  profiles <- tolower(Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit"))
  profiles <- unlist(strsplit(profiles, "[,;[:space:]]+"))
  any(profiles %in% c("fixture", "fixtures", "fit", "heavy", "slow", "visual-fixture", "visual_fixture"))
}

.test_cache_unavailable <- function(message) {
  if (.test_cache_required_for_active_profile()) {
    testthat::fail(message)
  } else {
    skip(message)
  }
}

.test_cache_required_catalog <- function() {
  if (exists("bayestools_required_fit_catalog", mode = "function")) {
    return(bayestools_required_fit_catalog())
  }

  NULL
}

skip_if_no_fits <- function() {
  model_registry_file <- file.path(test_files_dir, "model_registry.RDS")
  if (!file.exists(model_registry_file)) {
    .test_cache_unavailable("Pre-fitted models not found. Run test-00-model-fits.R first.")
  }

  catalog <- .test_cache_required_catalog()
  required_fits <- if (is.null(catalog)) NULL else catalog$model_name
  required_margliks <- if (is.null(catalog)) NULL else catalog$model_name[catalog$has_marglik]

  cache_complete <- .test_cache_requirements_complete(
    required_fits = required_fits,
    required_margliks = required_margliks,
    registry_file = model_registry_file
  )
  cache_current <- .test_cache_metadata_current(
    "model-fit",
    required_fits = required_fits,
    required_margliks = required_margliks,
    registry_file = model_registry_file
  )

  if (!cache_complete || !cache_current) {
    .test_cache_unavailable("Pre-fitted model cache is missing required artifacts or has stale metadata. Run test-00-model-fits.R first.")
  }
}

skip_if_missing_fits <- function(names) {
  skip_if_no_fits()

  fit_files <- file.path(temp_fits_dir, paste0(names, ".RDS"))
  missing <- names[!file.exists(fit_files)]

  if (length(missing) > 0L) {
    .test_cache_unavailable(paste0(
      "Required pre-fitted models not found: ",
      paste(missing, collapse = ", "),
      ". Run test-00-model-fits.R first."
    ))
  }
}

skip_if_no_fit_cache <- function() {
  skip_if_no_fits()
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

  quantiles <- mquant(prior, 0.5)

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
  par(mfcol = c(1, ncol(samples) - 1))

  for (i in 2:ncol(samples)) {
    hist(samples[, i], main = print(prior, plot = TRUE), breaks = 50, freq = FALSE)
    lines(densities[[i]])
    abline(v = quantiles[i], col = "blue", lwd = 2)
    if (prior$weights$type != "fixed") {
      expect_equal(.25, mcdf(prior, mquant(prior, 0.25)[, i])[, i], tolerance = 1e-5)
      expect_equal(.25, mccdf(prior, mquant(prior, 0.75)[, i])[, i], tolerance = 1e-5)
    }
    if (!skip_moments) {
      expect_equal(unname(apply(samples, 2, mean)), unname(mean(prior)), tolerance = 1e-2)
      expect_equal(unname(apply(samples, 2, sd)), unname(sd(prior)), tolerance = 1e-2)
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

.fixture_metadata_from_catalog <- function(name) {
  if (!exists("bayestools_expected_fit_catalog", mode = "function")) {
    return(list())
  }

  catalog <- bayestools_expected_fit_catalog()
  row_index <- match(name, catalog$model_name)
  if (is.na(row_index)) {
    return(list())
  }

  row <- catalog[row_index, , drop = FALSE]
  list(
    profile = row$profile,
    model_family = row$model_family,
    formula = row$formula,
    scale_policy = row$scale_policy,
    prior_features = row$prior_features,
    oracle_type = row$oracle_type,
    expected_chains = row$expected_chains,
    expected_iterations = row$expected_iterations,
    expected_parameter_count = row$expected_parameter_count,
    expected_monitor = row$expected_monitor[[1]],
    expected_formula_parameters = row$expected_formula_parameters[[1]],
    expected_formula_scale = row$expected_formula_scale[[1]]
  )
}

.fixture_metadata_from_registry_entry <- function(registry_entry) {
  flag_columns <- setdiff(
    names(registry_entry),
    c("model_name", "has_marglik", "note")
  )
  flags <- as.list(registry_entry[flag_columns])
  flags <- lapply(flags, isTRUE)

  list(
    schema_version = 1L,
    generated_by = "save_fit",
    cache_scope = "model-fit",
    model_name = registry_entry$model_name,
    has_marglik = isTRUE(registry_entry$has_marglik),
    flags = flags,
    note = registry_entry$note,
    catalog = .fixture_metadata_from_catalog(registry_entry$model_name)
  )
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

  registry_entry <- data.frame(
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

  attr(fit, "fixture_metadata") <- .fixture_metadata_from_registry_entry(registry_entry)
  saveRDS(fit, file = file.path(temp_fits_dir, paste0(name, ".RDS")))

  # Save marglik if provided
  if (!is.null(marglik)) {
    saveRDS(marglik, file = file.path(temp_marglik_dir, paste0(name, ".RDS")))
  }

  # Return model metadata entry for registry
  list(
    fit = fit,
    marglik = marglik,
    registry_entry = registry_entry
  )
}

# Fit-cache helpers. The marker is written only after the cache has been
# validated, so a failed fit run cannot make the next run skip refitting.
.test_cache_indicator_file <- function(name) {
  file.path(temp_temp_dir, paste0(name, ".txt"))
}

.test_cache_bool_env <- function(name, default) {
  value <- Sys.getenv(name)
  if (value == "") {
    return(default)
  }

  value <- tolower(value)
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

.test_cache_paths_complete <- function(paths) {
  length(paths) == 0L || all(file.exists(paths))
}

.test_cache_registry_requirements <- function(registry_file = NULL) {
  if (is.null(registry_file)) {
    return(list(fits = character(), margliks = character(), files = character()))
  }

  if (!file.exists(registry_file)) {
    return(list(fits = character(), margliks = character(), files = registry_file))
  }

  registry <- readRDS(registry_file)
  fit_names <- unique(registry$model_name)
  has_marglik <- registry$has_marglik
  has_marglik[is.na(has_marglik)] <- FALSE
  marglik_names <- unique(registry$model_name[has_marglik])

  list(
    fits = fit_names,
    margliks = marglik_names,
    files = registry_file
  )
}

.test_cache_expected_requirements <- function(required_fits = NULL,
                                              required_margliks = NULL,
                                              required_files = NULL,
                                              registry_file = NULL) {
  registry_requirements <- .test_cache_registry_requirements(registry_file)

  list(
    fits = sort(unique(c(required_fits, registry_requirements$fits))),
    margliks = sort(unique(c(required_margliks, registry_requirements$margliks))),
    files = sort(unique(c(required_files, registry_requirements$files)))
  )
}

.test_cache_package_r_files <- function(files = character(), patterns = character()) {
  package_r_dir <- file.path(testthat::test_path("..", ".."), "R")
  package_source_files <- file.path(package_r_dir, files)
  for (pattern in patterns) {
    package_source_files <- c(
      package_source_files,
      list.files(package_r_dir, pattern = pattern, full.names = TRUE)
    )
  }
  package_source_files <- unique(package_source_files)
  package_source_files <- package_source_files[file.exists(package_source_files)]
  names(package_source_files) <- paste0(
    "package_R_",
    tools::file_path_sans_ext(basename(package_source_files))
  )
  package_source_files
}

.test_cache_existing_test_files <- function(files) {
  paths <- testthat::test_path(files)
  paths <- paths[file.exists(paths)]
  names(paths) <- paste0(
    "test_",
    gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(paths)))
  )
  paths
}

.test_cache_source_files <- function(name) {
  generator_sources <- c(
    description = testthat::test_path("..", "..", "DESCRIPTION"),
    .test_cache_package_r_files(
      files = c(
        "JAGS-fit.R",
        "JAGS-formula.R",
        "JAGS-marglik.R",
        "priors.R",
        "priors-tools.R",
        "priors-informed.R",
        "selection-kernels.R",
        "tools.R"
      ),
      patterns = "^distributions-.*\\.R$"
    ),
    test_00_model_fits = testthat::test_path("test-00-model-fits.R")
  )

  fixture_consumer_sources <- c(
    .test_cache_package_r_files(
      files = c(
        "JAGS-fit.R",
        "JAGS-formula.R",
        "model-averaging.R",
        "selection-kernels.R",
        "summary-tables.R",
        "interpret.R"
      )
    ),
    common_functions = testthat::test_path("common-functions.R"),
    expected_fit_catalog = testthat::test_path("helper-expected-fit-catalog.R"),
    semantic_oracles = testthat::test_path("helper-semantic-oracles.R"),
    .test_cache_existing_test_files(c(
      "test-fixture-integrity.R",
      "test-JAGS-ensemble-tables.R",
      "test-JAGS-fit.R",
      "test-JAGS-formula-scale.R",
      "test-JAGS-formula.R",
      "test-JAGS-summary-tables.R",
      "test-model-averaging.R",
      "test-selection-kernels.R",
      "test-summary-tables.R",
      "test-weightfunction-redesign.R"
    ))
  )

  visual_fixture_consumer_sources <- c(
    .test_cache_package_r_files(
      files = c(
        "JAGS-fit.R",
        "model-averaging.R",
        "model-averaging-plots.R",
        "priors-plot.R",
        "summary-tables.R"
      )
    ),
    common_functions = testthat::test_path("common-functions.R"),
    expected_fit_catalog = testthat::test_path("helper-expected-fit-catalog.R"),
    .test_cache_existing_test_files(c(
      "test-JAGS-diagnostic-plots.R",
      "test-JAGS-ensemble-plots.R",
      "test-JAGS-marginal-distributions.R",
      "test-model-averaging-plots.R"
    ))
  )

  sources <- switch(
    name,
    `model-fit` = generator_sources,
    `fixture-consumer` = fixture_consumer_sources,
    `visual-fixture-consumer` = visual_fixture_consumer_sources,
    character()
  )

  sources[file.exists(sources)]
}

.test_cache_file_md5 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(tools::md5sum(path))
}

.test_cache_text_md5 <- function(text) {
  hash_file <- tempfile()
  on.exit(unlink(hash_file), add = TRUE)
  writeLines(text, hash_file, useBytes = TRUE)
  .test_cache_file_md5(hash_file)
}

.test_cache_source_functions <- function(name) {
  switch(
    name,
    `model-fit` = c(
      test_helper_fixture_metadata_from_catalog = ".fixture_metadata_from_catalog",
      test_helper_fixture_metadata_from_registry_entry = ".fixture_metadata_from_registry_entry",
      test_helper_save_fit = "save_fit",
      catalog_registry_flag_columns = "bayestools_registry_flag_columns",
      catalog_registry_schema_columns = "bayestools_registry_schema_columns",
      catalog_optional_fit_requirements = "bayestools_optional_fit_requirements",
      catalog_find_calls = ".bayestools_find_calls",
      catalog_static_string = ".bayestools_static_string",
      catalog_static_logical = ".bayestools_static_logical",
      catalog_save_fit_rows = ".bayestools_save_fit_catalog_rows",
      catalog_semantic_fit_overrides = "bayestools_semantic_fit_catalog_overrides",
      catalog_expected_fit = "bayestools_expected_fit_catalog",
      catalog_required_fit = "bayestools_required_fit_catalog"
    ),
    character()
  )
}

.test_cache_function_hashes <- function(name) {
  source_functions <- .test_cache_source_functions(name)
  if (length(source_functions) == 0L) {
    return(character())
  }

  vapply(source_functions, function(function_name) {
    if (!exists(function_name, mode = "function")) {
      return(NA_character_)
    }

    function_object <- get(function_name, mode = "function")
    .test_cache_text_md5(paste(deparse(function_object), collapse = "\n"))
  }, character(1))
}

.test_cache_source_hashes <- function(name) {
  source_files <- .test_cache_source_files(name)
  file_hashes <- vapply(source_files, .test_cache_file_md5, character(1))
  c(file_hashes, .test_cache_function_hashes(name))
}

.test_cache_marker_values <- function(indicator_file) {
  if (!file.exists(indicator_file)) {
    return(character())
  }

  lines <- readLines(indicator_file, warn = FALSE)
  pieces <- regmatches(lines, regexec("^([^:]+):\\s*(.*)$", lines))
  pieces <- pieces[lengths(pieces) == 3L]
  if (length(pieces) == 0L) {
    return(character())
  }

  values <- vapply(pieces, `[[`, character(1), 3L)
  names(values) <- vapply(pieces, `[[`, character(1), 2L)
  values
}

.test_cache_name_line <- function(names) {
  paste(sort(unique(names)), collapse = ",")
}

.test_cache_metadata_current <- function(name, required_fits = NULL,
                                         required_margliks = NULL,
                                         required_files = NULL,
                                         registry_file = NULL) {
  indicator_file <- .test_cache_indicator_file(name)
  marker <- .test_cache_marker_values(indicator_file)
  if (length(marker) == 0L) {
    return(FALSE)
  }

  marker_value <- function(key) {
    if (!key %in% names(marker)) {
      return(NA_character_)
    }

    unname(marker[[key]])
  }

  marker_test_files_dir <- marker_value("test_files_dir")
  if (is.na(marker_test_files_dir) ||
      !identical(
        normalizePath(marker_test_files_dir, winslash = "/", mustWork = FALSE),
        normalizePath(test_files_dir, winslash = "/", mustWork = TRUE)
      )) {
    return(FALSE)
  }

  requirements <- .test_cache_expected_requirements(
    required_fits = required_fits,
    required_margliks = required_margliks,
    required_files = required_files,
    registry_file = registry_file
  )

  if (!identical(marker_value("required_fits"), .test_cache_name_line(requirements$fits)) ||
      !identical(marker_value("required_margliks"), .test_cache_name_line(requirements$margliks))) {
    return(FALSE)
  }

  if (!is.null(registry_file) && file.exists(registry_file)) {
    registry_md5 <- .test_cache_file_md5(registry_file)
    if (!identical(marker_value("registry_md5"), registry_md5)) {
      return(FALSE)
    }
  }

  source_hashes <- .test_cache_source_hashes(name)
  for (source_name in names(source_hashes)) {
    key <- paste0("source_md5_", source_name)
    if (!identical(marker_value(key), unname(source_hashes[[source_name]]))) {
      return(FALSE)
    }
  }

  TRUE
}

.test_cache_requirements_complete <- function(required_fits = NULL,
                                              required_margliks = NULL,
                                              required_files = NULL,
                                              registry_file = NULL) {
  requirements <- .test_cache_expected_requirements(
    required_fits = required_fits,
    required_margliks = required_margliks,
    required_files = required_files,
    registry_file = registry_file
  )

  fit_files <- file.path(temp_fits_dir, paste0(requirements$fits, ".RDS"))
  marglik_files <- file.path(temp_marglik_dir, paste0(requirements$margliks, ".RDS"))

  .test_cache_paths_complete(fit_files) &&
    .test_cache_paths_complete(marglik_files) &&
    .test_cache_paths_complete(requirements$files)
}

.test_cache_prune_extra_artifacts <- function(requirements) {
  prune_dir <- function(directory, expected_files) {
    actual_files <- list.files(directory, pattern = "\\.RDS$", full.names = TRUE)
    extra_files <- actual_files[!basename(actual_files) %in% expected_files]
    if (length(extra_files) > 0L) {
      unlink(extra_files)
    }
  }

  prune_dir(temp_fits_dir, paste0(requirements$fits, ".RDS"))
  prune_dir(temp_marglik_dir, paste0(requirements$margliks, ".RDS"))
  invisible(TRUE)
}

# Skip model fitting if a validated cache exists and BAYESTOOLS_TEST_SKIP_REFIT is TRUE.
skip_refit_if_cached <- function(name, required_fits = NULL, required_margliks = NULL,
                                 required_files = NULL, registry_file = NULL) {
  skip_refit <- .test_cache_bool_env("BAYESTOOLS_TEST_SKIP_REFIT", TRUE)
  indicator_file <- .test_cache_indicator_file(name)

  cache_complete <- file.exists(indicator_file) &&
    .test_cache_requirements_complete(
      required_fits = required_fits,
      required_margliks = required_margliks,
      required_files = required_files,
      registry_file = registry_file
    ) &&
    .test_cache_metadata_current(
      name = name,
      required_fits = required_fits,
      required_margliks = required_margliks,
      required_files = required_files,
      registry_file = registry_file
    )

  if (skip_refit && cache_complete) {
    skip("Skipping model refitting: validated cache exists and BAYESTOOLS_TEST_SKIP_REFIT=TRUE.")
  }

  if (file.exists(indicator_file)) {
    unlink(indicator_file)
  }

  invisible(FALSE)
}

mark_refit_cache_complete <- function(name, required_fits = NULL, required_margliks = NULL,
                                      required_files = NULL, registry_file = NULL) {
  cache_complete <- .test_cache_requirements_complete(
    required_fits = required_fits,
    required_margliks = required_margliks,
    required_files = required_files,
    registry_file = registry_file
  )

  if (!cache_complete) {
    stop("Cannot mark cache complete because required cache artifacts are missing.", call. = FALSE)
  }

  requirements <- .test_cache_expected_requirements(
    required_fits = required_fits,
    required_margliks = required_margliks,
    required_files = required_files,
    registry_file = registry_file
  )
  source_hashes <- .test_cache_source_hashes(name)
  source_hash_lines <- paste0("source_md5_", names(source_hashes), ": ", unname(source_hashes))
  registry_hash_line <- if (!is.null(registry_file) && file.exists(registry_file)) {
    paste("registry_md5:", .test_cache_file_md5(registry_file))
  } else {
    NULL
  }

  .test_cache_prune_extra_artifacts(requirements)

  writeLines(
    c(
      paste("name:", name),
      paste("completed_at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("test_files_dir:", test_files_dir),
      paste("required_fits:", .test_cache_name_line(requirements$fits)),
      paste("required_margliks:", .test_cache_name_line(requirements$margliks)),
      registry_hash_line,
      source_hash_lines
    ),
    .test_cache_indicator_file(name)
  )

  invisible(TRUE)
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
