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

  # 'new' has no baseline: accepting its candidate creates the baseline
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
  writeLines("new baseline", file.path(root, "new.new.txt"))
  names <- c(names, "new")

  changes <- reference_table_candidates(root)
  expect_identical(changes[["name"]], c(paste0("context/", names[1:3]), "new"))

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
    if (file.exists(changes[["expected"]][[i]])) {
      file.copy(changes[["expected"]][[i]], staged_expected[[i]])
    } else {
      writeLines(character(), staged_expected[[i]])
    }
    file.copy(changes[["candidate"]][[i]], staged_candidate[[i]])
  }

  for (accepted in c(1L, 4L)) {
    file.copy(staged_candidate[[accepted]], staged_expected[[accepted]], overwrite = TRUE)
    unlink(staged_candidate[[accepted]])
  }
  unlink(staged_candidate[[2L]])

  reviewed <- .apply_reference_table_review(
    changes,
    staged_expected,
    staged_candidate
  )
  expect_identical(reviewed[["status"]], c("accepted", "rejected", "skipped", "accepted"))
  expect_identical(
    readLines(changes[["expected"]][[1L]], warn = FALSE),
    "accept new"
  )
  expect_identical(
    readLines(changes[["expected"]][[2L]], warn = FALSE),
    "reject old"
  )
  expect_identical(
    readLines(changes[["expected"]][[4L]], warn = FALSE),
    "new baseline"
  )
  expect_false(file.exists(changes[["candidate"]][[1L]]))
  expect_false(file.exists(changes[["candidate"]][[2L]]))
  expect_true(file.exists(changes[["candidate"]][[3L]]))
  expect_false(file.exists(changes[["candidate"]][[4L]]))
})


test_that("text reference candidates are removed on agreement and refreshed on change", {

  print_dir      <- withr::local_tempdir()
  reference_file <- file.path(print_dir, "text.txt")
  candidate_file <- reference_table_candidate_path(reference_file)
  writeLines("stale candidate", candidate_file)
  writeLines(c("first", "second"), reference_file)

  expect_success(test_reference_text("first\nsecond", "text.txt", print_dir = print_dir))
  expect_false(file.exists(candidate_file))

  expect_failure(test_reference_text("latest\noutput", "text.txt", print_dir = print_dir),
                 "Changed reference table cached")
  expect_identical(readLines(candidate_file, warn = FALSE), c("latest", "output"))
  expect_identical(readLines(reference_file, warn = FALSE), c("first", "second"))
})


test_that("snapshot changes from a run become candidates and reviewed files are restored", {

  root  <- withr::local_tempdir()
  snaps <- file.path(root, "_snaps", "context")
  dir.create(snaps, recursive = TRUE)
  writeLines("reviewed", file.path(snaps, "changed.svg"))
  writeLines("reviewed", file.path(snaps, "kept.svg"))
  writeLines("reviewed", file.path(snaps, "deleted.md"))
  writeLines("failed comparison", file.path(snaps, "kept.new.svg"))
  before <- bayestools_snapshot_state(root)
  expect_setequal(basename(names(before)), c("changed.svg", "kept.svg", "deleted.md", "kept.new.svg"))

  # what a run can do: a new baseline, a changed one, and a cleanup deletion
  writeLines("unreviewed", file.path(snaps, "created.svg"))
  writeLines("rewritten", file.path(snaps, "changed.svg"))
  unlink(file.path(snaps, "deleted.md"))

  result <- bayestools_quarantine_snapshots(root, before)
  expect_setequal(basename(result$changed), c("created.svg", "changed.svg"))
  expect_identical(basename(result$deleted), "deleted.md")
  expect_identical(result$problems, character())
  expect_false(file.exists(file.path(snaps, "created.svg")))
  expect_identical(readLines(file.path(snaps, "created.new.svg")), "unreviewed")
  expect_identical(readLines(file.path(snaps, "changed.svg")), "reviewed")
  expect_identical(readLines(file.path(snaps, "changed.new.svg")), "rewritten")
  expect_identical(readLines(file.path(snaps, "deleted.md")), "reviewed")
  expect_identical(readLines(file.path(snaps, "kept.new.svg")), "failed comparison")

  # an unchanged run reports nothing
  result <- bayestools_quarantine_snapshots(root, bayestools_snapshot_state(root))
  expect_identical(result, list(changed = character(), deleted = character(), problems = character()))
})


test_that("the snapshot guard keeps testthat's fresh candidates and restores removed directories", {

  root  <- withr::local_tempdir()
  snaps <- file.path(root, "_snaps")
  dir.create(file.path(snaps, "sections"), recursive = TRUE)
  dir.create(file.path(snaps, "retired"))
  writeLines("old section", file.path(snaps, "sections", "file.md"))
  writeLines("stale candidate", file.path(snaps, "sections", "other.new.md"))
  writeLines("reviewed", file.path(snaps, "sections", "other.md"))
  writeLines("obsolete", file.path(snaps, "retired", "old.svg"))
  before <- bayestools_snapshot_state(root)

  # a changed and an added section: testthat writes the added section into the
  # baseline and its candidate holds the changed and the added section; a
  # stale candidate beside another changed file is replaced; cleanup removes a
  # retired directory
  writeLines(c("old section", "added section"), file.path(snaps, "sections", "file.md"))
  writeLines(c("changed section", "added section"), file.path(snaps, "sections", "file.new.md"))
  writeLines("rewritten", file.path(snaps, "sections", "other.md"))
  unlink(file.path(snaps, "retired"), recursive = TRUE)

  result <- bayestools_quarantine_snapshots(root, before)
  expect_identical(result$problems, character())
  expect_identical(readLines(file.path(snaps, "sections", "file.md")), "old section")
  expect_identical(readLines(file.path(snaps, "sections", "file.new.md")), c("changed section", "added section"))
  expect_identical(readLines(file.path(snaps, "sections", "other.md")), "reviewed")
  expect_identical(readLines(file.path(snaps, "sections", "other.new.md")), "rewritten")
  expect_identical(readLines(file.path(snaps, "retired", "old.svg")), "obsolete")
})


test_that("a failed candidate write still restores the reviewed baseline and is reported", {

  root  <- withr::local_tempdir()
  snaps <- file.path(root, "_snaps", "context")
  dir.create(snaps, recursive = TRUE)
  writeLines("reviewed", file.path(snaps, "plot.svg"))
  before <- bayestools_snapshot_state(root)

  writeLines("unreviewed", file.path(snaps, "plot.svg"))
  # a directory in place of the candidate makes writing it fail
  dir.create(file.path(snaps, "plot.new.svg"))

  result <- suppressWarnings(bayestools_quarantine_snapshots(root, before))
  expect_identical(readLines(file.path(snaps, "plot.svg")), "reviewed")
  expect_length(result$problems, 1L)
  expect_match(result$problems, "writing candidate", fixed = TRUE)
})


test_that("candidates of new snapshots are reviewed against an empty baseline", {

  root  <- withr::local_tempdir()
  snaps <- file.path(root, "_snaps", "context")
  dir.create(snaps, recursive = TRUE)
  writeLines("accept me", file.path(snaps, "accept.new.svg"))
  writeLines("reject me", file.path(snaps, "reject.new.md"))
  writeLines("skip me", file.path(snaps, "skip.new.svg"))
  writeLines("reviewed", file.path(snaps, "existing.svg"))
  writeLines("changed", file.path(snaps, "existing.new.svg"))

  expect_setequal(basename(bayestools_new_snapshot_candidates(root)$candidate),
                  c("accept.new.svg", "reject.new.md", "skip.new.svg"))

  # a stand-in for testthat's reviewer: accept moves a candidate over its
  # staged baseline, reject deletes it, skip leaves it
  review <- function(path) {
    staged <- file.path(path, "_snaps", "context")
    file.rename(file.path(staged, "accept.new.svg"), file.path(staged, "accept.svg"))
    file.remove(file.path(staged, "reject.new.md"))
  }
  result <- bayestools_review_new_snapshots(root, review = review)
  expect_identical(result$status[order(basename(result$candidate))], c("accepted", "rejected", "skipped"))
  expect_identical(readLines(file.path(snaps, "accept.svg")), "accept me")
  expect_false(file.exists(file.path(snaps, "accept.new.svg")))
  expect_false(file.exists(file.path(snaps, "reject.md")))
  expect_false(file.exists(file.path(snaps, "reject.new.md")))
  expect_true(file.exists(file.path(snaps, "skip.new.svg")))
  expect_false(file.exists(file.path(snaps, "skip.svg")))
  expect_identical(readLines(file.path(snaps, "existing.new.svg")), "changed")
})


test_that("missing reference baselines fail and cache a reviewable candidate", {

  print_dir <- withr::local_tempdir()
  table     <- data.frame(value = 1, row.names = "row")
  expected_candidate <- capture_output_lines(table, print = TRUE, width = 150)
  checks <- list(
    table = function() test_reference_table(table, "table.txt", print_dir = print_dir),
    stochastic = function() test_reference_table_stochastic(table, "stochastic.txt", print_dir = print_dir),
    text = function() test_reference_text("first\nsecond", "text.txt", print_dir = print_dir)
  )
  candidates <- list(table = expected_candidate, stochastic = expected_candidate, text = c("first", "second"))
  for (name in names(checks)) {
    reference_file <- file.path(print_dir, paste0(name, ".txt"))
    expect_failure(checks[[name]](), "not found")
    expect_false(file.exists(reference_file))
    expect_identical(readLines(reference_table_candidate_path(reference_file), warn = FALSE), candidates[[name]])
  }
  expect_identical(reference_table_candidates(print_dir)[["name"]], c("stochastic", "table", "text"))
})
