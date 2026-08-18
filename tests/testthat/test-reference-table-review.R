skip_if_not_test_profile("unit")

source(testthat::test_path("common-functions.R"))


test_that("changed reference tables are cached without replacing baselines", {

  print_dir <- withr::local_tempdir()
  expected <- data.frame(value = 1, row.names = "row")
  actual   <- data.frame(value = 2, row.names = "row")
  reference_file <- file.path(print_dir, "table.txt")
  candidate_file <- reference_table_candidate_path(reference_file)
  expected_output <- capture_output_lines(
    expected,
    print = TRUE,
    width = 150
  )
  actual_output <- capture_output_lines(actual, print = TRUE, width = 150)
  writeLines(expected_output, reference_file)

  expect_failure(
    test_reference_table(actual, "table.txt", print_dir = print_dir),
    "Changed reference table cached"
  )
  expect_identical(readLines(reference_file, warn = FALSE), expected_output)
  expect_identical(readLines(candidate_file, warn = FALSE), actual_output)

  expect_no_error(
    test_reference_table(expected, "table.txt", print_dir = print_dir)
  )
  expect_false(file.exists(candidate_file))
})


test_that("stochastic table candidates track structural changes only", {

  print_dir <- withr::local_tempdir()
  expected <- data.frame(
    Mean = "1.25",
    SD = "0.50",
    row.names = "theta"
  )
  refreshed <- expected
  refreshed[1, ] <- c("1.30", "0.47")
  changed_schema <- refreshed
  names(changed_schema)[[2L]] <- "SE"

  reference_file <- file.path(print_dir, "stochastic.txt")
  candidate_file <- reference_table_candidate_path(reference_file)
  writeLines(
    capture_output_lines(expected, print = TRUE, width = 150),
    reference_file
  )
  writeLines("stale", candidate_file)

  expect_no_error(test_reference_table_stochastic(
    refreshed,
    "stochastic.txt",
    print_dir = print_dir
  ))
  expect_false(file.exists(candidate_file))

  expect_failure(
    test_reference_table_stochastic(
      changed_schema,
      "stochastic.txt",
      print_dir = print_dir
    ),
    "Changed reference table cached"
  )
  expect_true(file.exists(candidate_file))
})


test_that("reference-table review applies accept, reject, and skip", {

  root <- withr::local_tempdir()
  names <- c("accept", "reject", "skip")
  for (name in names) {
    directory <- file.path(root, "context")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(name, "old"), file.path(directory, paste0(name, ".txt")))
    writeLines(
      paste(name, "new"),
      file.path(directory, paste0(name, ".new.txt"))
    )
  }
  writeLines("orphan", file.path(root, "orphan.new.txt"))

  changes <- reference_table_candidates(root)
  expect_identical(changes[["name"]], paste0("context/", names))

  review_root <- withr::local_tempdir()
  staged_expected  <- character(nrow(changes))
  staged_candidate <- character(nrow(changes))
  for (i in seq_len(nrow(changes))) {
    staged_expected[[i]] <- file.path(
      review_root,
      paste0(names[[i]], ".txt")
    )
    staged_candidate[[i]] <- reference_table_candidate_path(
      staged_expected[[i]]
    )
    file.copy(changes[["expected"]][[i]], staged_expected[[i]])
    file.copy(changes[["candidate"]][[i]], staged_candidate[[i]])
  }

  file.copy(staged_candidate[[1L]], staged_expected[[1L]], overwrite = TRUE)
  unlink(staged_candidate[[1L]])
  unlink(staged_candidate[[2L]])

  reviewed <- .apply_reference_table_review(
    changes,
    staged_expected,
    staged_candidate
  )
  expect_identical(reviewed[["status"]], c("accepted", "rejected", "skipped"))
  expect_identical(
    readLines(changes[["expected"]][[1L]], warn = FALSE),
    "accept new"
  )
  expect_identical(
    readLines(changes[["expected"]][[2L]], warn = FALSE),
    "reject old"
  )
  expect_false(file.exists(changes[["candidate"]][[1L]]))
  expect_false(file.exists(changes[["candidate"]][[2L]]))
  expect_true(file.exists(changes[["candidate"]][[3L]]))
})
