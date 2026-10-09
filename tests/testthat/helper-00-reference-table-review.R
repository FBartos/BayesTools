reference_table_candidate_path <- function(path) {

  sub("(\\.[^.]+)$", ".new\\1", path)
}


cache_reference_table_candidate <- function(actual_output, expected_output,
                                              reference_file,
                                              actual_comparison = actual_output,
                                              expected_comparison = expected_output) {

  candidate <- reference_table_candidate_path(reference_file)
  if (identical(actual_comparison, expected_comparison)) {
    unlink(candidate)
    return(NULL)
  }

  dir.create(dirname(candidate), recursive = TRUE, showWarnings = FALSE)
  writeLines(actual_output, candidate)
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}


reference_table_failure_info <- function(info_msg, candidate) {

  review_info <- paste0(
    "Changed reference table cached at ", candidate,
    " and will be offered for review after the test run."
  )
  if (is.null(info_msg) || !nzchar(info_msg)) {
    return(review_info)
  }

  paste(info_msg, review_info, sep = "\n")
}


reference_table_candidates <- function(root = file.path("tests", "results")) {

  empty <- data.frame(
    name      = character(),
    expected  = character(),
    candidate = character(),
    stringsAsFactors = FALSE
  )
  if (!dir.exists(root)) {
    return(empty)
  }

  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  candidates <- list.files(
    root,
    pattern    = "[.]new[.]txt$",
    recursive  = TRUE,
    full.names = TRUE
  )
  if (length(candidates) == 0L) {
    return(empty)
  }

  candidates <- normalizePath(candidates, winslash = "/", mustWork = TRUE)
  expected   <- sub("[.]new[.]txt$", ".txt", candidates)
  keep       <- file.exists(expected)
  candidates <- candidates[keep]
  expected   <- expected[keep]
  if (length(candidates) == 0L) {
    return(empty)
  }

  relative <- substring(candidates, nchar(root) + 2L)
  changes <- data.frame(
    name      = sub("[.]new[.]txt$", "", relative),
    expected  = expected,
    candidate = candidates,
    stringsAsFactors = FALSE
  )
  rownames(changes) <- NULL
  changes[order(changes[["name"]]), , drop = FALSE]
}


.copy_reference_table_file <- function(from, to, description) {

  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(from, to, overwrite = TRUE)
  if (!copied) {
    stop("Failed to write ", description, ": ", to, call. = FALSE)
  }

  invisible(to)
}


.apply_reference_table_review <- function(changes, staged_expected,
                                          staged_candidate) {

  changes[["status"]] <- "skipped"
  for (i in seq_len(nrow(changes))) {
    if (file.exists(staged_candidate[[i]])) {
      next
    }

    staged <- readLines(staged_expected[[i]], warn = FALSE, encoding = "UTF-8")
    old <- readLines(changes[["expected"]][[i]], warn = FALSE,
                     encoding = "UTF-8")
    new <- readLines(changes[["candidate"]][[i]], warn = FALSE,
                     encoding = "UTF-8")
    if (identical(staged, new)) {
      .copy_reference_table_file(
        changes[["candidate"]][[i]],
        changes[["expected"]][[i]],
        "accepted reference table"
      )
      changes[["status"]][[i]] <- "accepted"
      unlink(changes[["candidate"]][[i]])
    } else if (identical(staged, old)) {
      changes[["status"]][[i]] <- "rejected"
      unlink(changes[["candidate"]][[i]])
    } else {
      stop(
        "Reviewed reference table has an unexpected value: ",
        changes[["name"]][[i]], ".",
        call. = FALSE
      )
    }
  }

  changes
}


review_reference_tables <- function(root = file.path("tests", "results")) {

  changes <- reference_table_candidates(root)
  if (nrow(changes) == 0L) {
    return(invisible(changes))
  }

  message(
    "Changed reference tables:\n- ",
    paste(changes[["name"]], collapse = "\n- ")
  )
  if (!interactive()) {
    message("Cached candidates were kept for later interactive review.")
    return(invisible(changes))
  }

  required <- c("shiny", "diffviewer")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) {
    message(
      "Interactive reference-table review is unavailable; install: ",
      paste(required[!available], collapse = ", "), "."
    )
    return(invisible(changes))
  }

  review_root <- tempfile("bayestools-reference-table-review-")
  on.exit(unlink(review_root, recursive = TRUE), add = TRUE)
  staged_expected  <- character(nrow(changes))
  staged_candidate <- character(nrow(changes))
  for (i in seq_len(nrow(changes))) {
    relative <- paste0(changes[["name"]][[i]], ".txt")
    staged_expected[[i]] <- file.path(review_root, "_snaps", relative)
    staged_candidate[[i]] <- reference_table_candidate_path(
      staged_expected[[i]]
    )
    .copy_reference_table_file(
      changes[["expected"]][[i]],
      staged_expected[[i]],
      "review reference table"
    )
    .copy_reference_table_file(
      changes[["candidate"]][[i]],
      staged_candidate[[i]],
      "review candidate"
    )
  }

  testthat::snapshot_review(path = review_root)
  changes <- .apply_reference_table_review(
    changes,
    staged_expected,
    staged_candidate
  )
  message(
    "Reference-table review complete: ",
    paste0(
      changes[["name"]], " [", changes[["status"]], "]",
      collapse = ", "
    ),
    "."
  )

  invisible(changes)
}
