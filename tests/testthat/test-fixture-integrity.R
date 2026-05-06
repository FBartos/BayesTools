skip_if_not_test_profile(c("fixture", "fit"))

# ============================================================================ #
# TEST FILE: Fixture Integrity
# ============================================================================ #
#
# PURPOSE:
#   Validates the first expected-fit catalog slice against the generated fit
#   registry, cached fit files, monitored posterior columns, formula metadata,
#   sample dimensions, and marginal-likelihood file policy.
#
# TAGS: @fixture, @catalog, @integrity
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

if (bayestools_test_profile_includes("fit")) {
  skip_if_not_installed("rjags")
}

test_that("expected fit catalog is internally well formed", {
  catalog <- bayestools_expected_fit_catalog()

  expect_expected_fit_catalog_schema(catalog)
  expect_true(any(vapply(catalog$expected_monitor, length, integer(1)) > 0L))
})

test_that("fixture registry matches the expected catalog slice", {
  registry_file <- file.path(test_files_dir, "model_registry.RDS")
  expect_true(file.exists(registry_file), info = paste("Missing fixture registry:", registry_file))

  registry <- readRDS(registry_file)
  expect_registry_matches_expected_catalog(registry)
})

test_that("validated fixture cache marker belongs to the active cache directory", {
  expect_cache_completion_marker("model-fit")
})

test_that("fixture cache metadata rejects stale artifact declarations", {
  registry_file <- file.path(test_files_dir, "model_registry.RDS")
  expect_true(file.exists(registry_file), info = paste("Missing fixture registry:", registry_file))

  marker_name <- paste0("stale-model-fit-", Sys.getpid())
  marker_file <- .test_cache_indicator_file(marker_name)
  withr::defer(unlink(marker_file))

  writeLines(
    c(
      paste("name:", marker_name),
      "completed_at: 2026-05-04 00:00:00 CEST",
      paste("test_files_dir:", test_files_dir),
      "required_fits: fit_that_should_not_exist",
      "required_margliks:",
      paste("registry_md5:", .test_cache_file_md5(registry_file)),
      "source_md5_test_00_model_fits: stale-source-hash"
    ),
    marker_file
  )

  catalog <- bayestools_required_fit_catalog()
  expect_false(
    .test_cache_metadata_current(
      marker_name,
      required_fits = catalog$model_name,
      required_margliks = catalog$model_name[catalog$has_marglik],
      registry_file = registry_file
    )
  )
})

test_that("fixture cache metadata rejects stale source declarations", {
  registry_file <- file.path(test_files_dir, "model_registry.RDS")
  expect_true(file.exists(registry_file), info = paste("Missing fixture registry:", registry_file))

  marker_file <- .test_cache_indicator_file("model-fit")
  original_marker <- if (file.exists(marker_file)) readLines(marker_file, warn = FALSE) else NULL
  withr::defer({
    if (is.null(original_marker)) {
      unlink(marker_file)
    } else {
      writeLines(original_marker, marker_file)
    }
  })

  catalog <- bayestools_required_fit_catalog()
  source_hashes <- .test_cache_source_hashes("model-fit")
  source_hashes[[1L]] <- "stale-source-hash"

  writeLines(
    c(
      "name: model-fit",
      "completed_at: 2026-05-04 00:00:00 CEST",
      paste("test_files_dir:", test_files_dir),
      paste("required_fits:", .test_cache_name_line(catalog$model_name)),
      paste("required_margliks:", .test_cache_name_line(catalog$model_name[catalog$has_marglik])),
      paste("registry_md5:", .test_cache_file_md5(registry_file)),
      paste0("source_md5_", names(source_hashes), ": ", unname(source_hashes))
    ),
    marker_file
  )

  expect_false(
    .test_cache_metadata_current(
      "model-fit",
      required_fits = catalog$model_name,
      required_margliks = catalog$model_name[catalog$has_marglik],
      registry_file = registry_file
    )
  )
})

test_that("fixture artifact validators reject malformed temporary RDS payloads", {
  invalid_rds <- tempfile(fileext = ".RDS")
  malformed_registry_file <- tempfile(fileext = ".RDS")
  malformed_fit_file <- tempfile(fileext = ".RDS")
  malformed_marglik_file <- tempfile(fileext = ".RDS")
  withr::defer(unlink(c(
    invalid_rds,
    malformed_registry_file,
    malformed_fit_file,
    malformed_marglik_file
  )))

  writeLines("not an RDS payload", invalid_rds)
  expect_error(readRDS(invalid_rds))

  saveRDS(data.frame(model_name = "fit_missing_schema"), malformed_registry_file)
  expect_error(
    expect_registry_schema(readRDS(malformed_registry_file)),
    class = "expectation_failure"
  )

  catalog <- bayestools_required_fit_catalog()
  row <- catalog[1L, , drop = FALSE]
  malformed_fit <- structure(
    list(mcmc = NULL),
    class = c("runjags", "BayesTools_fit", "list")
  )
  saveRDS(malformed_fit, malformed_fit_file)
  expect_error(
    expect_fit_declared_metadata(readRDS(malformed_fit_file), row),
    class = "expectation_failure"
  )

  malformed_marglik <- structure(list(logml = NA_real_), class = "bridge")
  saveRDS(malformed_marglik, malformed_marglik_file)
  expect_error(
    expect_marglik_object(readRDS(malformed_marglik_file)),
    class = "expectation_failure"
  )
})

test_that("registry artifact policy matches cached fit and marginal-likelihood files", {
  registry_file <- file.path(test_files_dir, "model_registry.RDS")
  expect_true(file.exists(registry_file), info = paste("Missing fixture registry:", registry_file))

  registry <- readRDS(registry_file)
  expect_registry_artifact_policy(registry)
})

test_that("cataloged fixture files expose expected monitors and metadata", {
  catalog <- bayestools_required_fit_catalog()

  for (model_name in catalog$model_name) {
    row <- catalog[catalog$model_name == model_name, , drop = FALSE]
    fit_file <- expect_fit_file_present(model_name, catalog = catalog)
    marglik_file <- expect_marglik_file_present_or_absent(model_name, catalog = catalog)

    fit <- readRDS(fit_file)
    expect_fit_declared_metadata(fit, row)
    if (length(row$expected_monitor[[1]]) > 0L) {
      expect_fit_monitors(fit, row$expected_monitor[[1]])
      expect_fit_formula_metadata(
        fit,
        expected_formula_parameters = row$expected_formula_parameters[[1]],
        expected_formula_scale = row$expected_formula_scale[[1]]
      )
      expect_fit_sample_dimensions(
        fit,
        expected_chains = row$expected_chains,
        expected_iterations = row$expected_iterations
      )
    }

    if (isTRUE(row$has_marglik)) {
      expect_marglik_object(readRDS(marglik_file))
    }
  }
})
