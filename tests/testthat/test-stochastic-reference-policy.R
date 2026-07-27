skip_if_not_test_profile("unit")


test_that("stochastic reference signatures ignore sampled values only", {

  expected <- c(
    "       Mean   SD 0.025 0.975",
    "theta  1.25 0.50  0.25  2.25",
    "",
    "Conditional summary is based on 79 samples."
  )
  refreshed <- c(
    "       Mean   SD 0.025 0.975",
    "theta  1.30 0.47  0.31  2.18",
    "",
    "Conditional summary is based on 83 samples."
  )

  expect_equal(
    stochastic_reference_signature(refreshed),
    stochastic_reference_signature(expected)
  )

  changed_schema <- refreshed
  changed_schema[[1L]] <- "       Mean   SD 0.050 0.950"
  expect_false(identical(
    stochastic_reference_signature(changed_schema),
    stochastic_reference_signature(expected)
  ))

  changed_row <- refreshed
  changed_row[[2L]] <- sub("theta", "phi", changed_row[[2L]], fixed = TRUE)
  expect_false(identical(
    stochastic_reference_signature(changed_row),
    stochastic_reference_signature(expected)
  ))

  expected_model_summary <- c(
    "log(marglik) -31.94 s ~ Normal(0, 1)[0, Inf]"
  )
  refreshed_model_summary <- c(
    "log(marglik) -29.48 s ~ Normal(0, 1)[0, Inf]"
  )
  expect_equal(
    stochastic_reference_signature(refreshed_model_summary),
    stochastic_reference_signature(expected_model_summary)
  )
})


test_that("stochastic table invariants validate current numeric output", {

  table <- data.frame(
    Mean = c(0.2, 0.8),
    SD = c(0.1, NA),
    "0.025" = c(0.0, 0.7),
    "0.975" = c(0.4, 0.9),
    prior_prob = c(0.5, 0.5),
    post_prob = c(0.4, 0.6),
    ESS = c(100, 200),
    R_hat = c(1.00, 1.01),
    check.names = FALSE
  )
  expect_no_error(expect_stochastic_table_invariants(table))
})


test_that("stochastic references retain presentation structure", {

  reference <- data.frame(
    Mean = 1.25,
    SD = 0.50,
    "0.025" = 0.25,
    "0.975" = 2.25,
    row.names = "theta",
    check.names = FALSE
  )
  refreshed <- reference
  refreshed[1, ] <- c(1.30, 0.47, 0.31, 2.18)

  print_dir <- withr::local_tempdir()
  writeLines(
    capture_output_lines(reference, print = TRUE, width = 150),
    file.path(print_dir, "table.txt")
  )

  expect_no_error(test_reference_table_stochastic(
    refreshed,
    "table.txt",
    print_dir = print_dir
  ))
})


test_that("fitted visual snapshots require the validated fit cache", {

  fitted_visual_files <- c(
    "test-JAGS-diagnostic-plots.R",
    "test-JAGS-ensemble-plots.R",
    "test-JAGS-marginal-distributions.R",
    "test-model-averaging-plots.R"
  )
  for (file in fitted_visual_files) {
    source <- readLines(testthat::test_path(file), warn = FALSE)
    snapshot_lines <- grep("vdiffr::expect_doppelganger", source, fixed = TRUE)
    cache_guard_lines <- grep(
      "^\\s*skip_if_no_fits\\(\\)",
      source,
      perl = TRUE
    )

    expect_true(length(snapshot_lines) > 0L, info = file)
    expect_true(length(cache_guard_lines) > 0L, info = file)
    fitted_snapshot_lines <- snapshot_lines[
      snapshot_lines > min(cache_guard_lines)
    ]
    expect_true(length(fitted_snapshot_lines) > 0L, info = file)
    expect_true(
      all(vapply(
        fitted_snapshot_lines,
        function(line) any(cache_guard_lines < line),
        logical(1)
      )),
      info = paste(file, "must validate the versioned fitted cache before snapshots")
    )
  }
})
