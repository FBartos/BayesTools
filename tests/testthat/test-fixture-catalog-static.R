skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Static Fixture Catalog
# ============================================================================ #
#
# PURPOSE:
#   Validates source-derived fixture catalog accountability without requiring
#   generated JAGS fixture artifacts.
#
# TAGS: @fixture, @catalog, @unit
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("source-derived fixture catalog covers every generated fit", {
  source_rows <- .bayestools_save_fit_catalog_rows()
  catalog <- bayestools_expected_fit_catalog()
  expect_expected_fit_catalog_schema(catalog)

  expect_equal(nrow(source_rows), 73L)
  expect_equal(sum(source_rows$has_marglik), 22L)
  expect_equal(nrow(catalog), nrow(source_rows))
  expect_equal(catalog$model_name, source_rows$model_name)
  expect_equal(catalog$fit_file, paste0(catalog$model_name, ".RDS"))
  expect_equal(catalog$has_marglik, source_rows$has_marglik)
  expect_equal(catalog$note, source_rows$note)
  expect_false(anyDuplicated(catalog$model_name) > 0L)

  optional_catalog <- bayestools_optional_fit_requirements()
  expect_setequal(
    catalog$model_name[!is.na(catalog$requires_package)],
    optional_catalog$model_name
  )
  required_catalog <- bayestools_required_fit_catalog(catalog)
  unavailable_optional <- optional_catalog$model_name[
    !vapply(optional_catalog$requires_package, requireNamespace, logical(1), quietly = TRUE)
  ]
  expect_setequal(
    setdiff(catalog$model_name, required_catalog$model_name),
    unavailable_optional
  )

  semantic_catalog <- bayestools_semantic_fit_catalog_overrides()
  expect_setequal(
    catalog$model_name[catalog$oracle_type != "registry-metadata"],
    semantic_catalog$model_name
  )
  expect_equal(catalog$model_name[catalog$oracle_type != "registry-metadata"], catalog$model_name)
  expect_true(all(vapply(catalog$expected_monitor, length, integer(1)) > 0L))
})

test_that("source-derived fixture catalog preserves registry schema flags", {
  catalog <- bayestools_expected_fit_catalog()
  flag_cols <- c("has_marglik", bayestools_registry_flag_columns())

  for (flag in flag_cols) {
    expect_type(catalog[[flag]], "logical")
    expect_false(anyNA(catalog[[flag]]), info = paste("Catalog flag:", flag))
  }

  expect_true(all(nzchar(catalog$model_name)))
  expect_true(all(nzchar(catalog$note)))
  expect_type(catalog$requires_package, "character")
  expect_equal(is.na(catalog$marglik_file), !catalog$has_marglik)
  expect_equal(
    catalog$marglik_file[catalog$has_marglik],
    paste0(catalog$model_name[catalog$has_marglik], ".RDS")
  )
})

test_that("legacy cache markers are treated as stale metadata", {
  marker_name <- paste0("legacy-marker-", Sys.getpid())
  marker_file <- .test_cache_indicator_file(marker_name)
  withr::defer(unlink(marker_file), testthat::teardown_env())

  writeLines(
    c(
      paste("name:", marker_name),
      "completed_at: 2026-05-03 23:05:40 CEST",
      paste("test_files_dir:", test_files_dir)
    ),
    marker_file
  )

  expect_false(
    .test_cache_metadata_current(
      marker_name,
      required_fits = c("fit_a", "fit_b"),
      required_margliks = "fit_a"
    )
  )
})

test_that("fixture cache availability policy is profile-aware", {
  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "unit")
  expect_false(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "fixture")
  expect_true(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "fit")
  expect_true(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "visual-fixture")
  expect_true(.test_cache_required_for_active_profile())
})

test_that("model-fit cache marker hashes only fit-generation sources", {
  skip_if_not(.test_cache_package_sources_available(), "Repository R source files are not available in this installed-package test context.")

  source_files <- .test_cache_source_files("model-fit")

  expect_true("description" %in% names(source_files))
  expect_true("package_R_JAGS-fit" %in% names(source_files))
  expect_true("package_R_JAGS-formula" %in% names(source_files))
  expect_true("package_R_JAGS-formula-random" %in% names(source_files))
  expect_true("package_R_JAGS-lkj-cholesky" %in% names(source_files))
  expect_true("package_R_JAGS-marglik" %in% names(source_files))
  expect_true("package_R_random-effects-metadata" %in% names(source_files))
  expect_true("package_R_random-effects-reconstruction" %in% names(source_files))
  expect_true("package_R_random-priors" %in% names(source_files))
  expect_true("package_R_random-effects-summary" %in% names(source_files))
  expect_true("package_src_r_lkj_cc" %in% names(source_files))
  expect_true("package_src_lkj_BTLKJCore_cc" %in% names(source_files))
  expect_true("package_src_functions_BTLKJCholesky_cc" %in% names(source_files))
  expect_true("package_src_distributions_DBTLKJCPC_cc" %in% names(source_files))
  expect_false("package_R_summary-tables" %in% names(source_files))
  expect_false("common_functions" %in% names(source_files))
  expect_false("expected_fit_catalog" %in% names(source_files))
  expect_true("test_00_model_fits" %in% names(source_files))
  expect_true(all(file.exists(source_files)), info = paste(source_files[!file.exists(source_files)], collapse = ", "))

  source_functions <- .test_cache_source_functions("model-fit")
  expect_true("test_helper_save_fit" %in% names(source_functions))
  expect_true("catalog_expected_fit" %in% names(source_functions))
  expect_true("catalog_semantic_fit_overrides" %in% names(source_functions))

  source_hashes <- .test_cache_source_hashes("model-fit")
  expect_setequal(names(source_hashes), c(names(source_files), names(source_functions)))
  expect_false(anyNA(source_hashes))
})

test_that("downstream cache consumer source scopes are separate from model fitting", {
  skip_if_not(.test_cache_package_sources_available(), "Repository R source files are not available in this installed-package test context.")

  model_fit_sources <- .test_cache_source_files("model-fit")
  fixture_sources <- .test_cache_source_files("fixture-consumer")
  visual_fixture_sources <- .test_cache_source_files("visual-fixture-consumer")

  expect_true("common_functions" %in% names(fixture_sources))
  expect_true("expected_fit_catalog" %in% names(fixture_sources))
  expect_true("package_R_JAGS-formula-random" %in% names(fixture_sources))
  expect_true("package_R_random-effects-metadata" %in% names(fixture_sources))
  expect_true("package_R_random-effects-reconstruction" %in% names(fixture_sources))
  expect_true("package_R_random-effects-summary" %in% names(fixture_sources))
  expect_true("package_R_random-priors" %in% names(fixture_sources))
  expect_true("package_R_summary-tables" %in% names(fixture_sources))
  expect_true("package_R_model-averaging" %in% names(fixture_sources))
  expect_false("test_00_model_fits" %in% names(fixture_sources))

  expect_true("common_functions" %in% names(visual_fixture_sources))
  expect_true("expected_fit_catalog" %in% names(visual_fixture_sources))
  expect_true("package_R_JAGS-formula-random" %in% names(visual_fixture_sources))
  expect_true("package_R_random-effects-summary" %in% names(visual_fixture_sources))
  expect_true("package_R_model-averaging-plots" %in% names(visual_fixture_sources))
  expect_true("test_test_JAGS_diagnostic_plots" %in% names(visual_fixture_sources))
  expect_false("test_00_model_fits" %in% names(visual_fixture_sources))

  expect_false("common_functions" %in% names(model_fit_sources))
  expect_false("expected_fit_catalog" %in% names(model_fit_sources))
  expect_false("package_R_summary-tables" %in% names(model_fit_sources))
  expect_false("package_R_model-averaging-plots" %in% names(model_fit_sources))

  expect_false(anyNA(.test_cache_source_hashes("fixture-consumer")))
  expect_false(anyNA(.test_cache_source_hashes("visual-fixture-consumer")))
})
