# Certification runner, dispatched by `tools/test-profile.R certification`:
#
#   certification              run the cases that are not current
#   certification --status     report the state of every case without running
#   certification --all        run every case
#   certification <case> ...   run the named cases
#
# A case is current when its marker holds its key
# (tests/testthat/helper-certification-cases.R). Each case runs in its own R
# process with the test controls of this one, writing its output to
# <BAYESTOOLS_TEST_FILES_DIR>/certification/<case>.log, within
# bayestools_certification_timeout seconds; that process writes the marker
# after a clean pass whose key did not change during the run.

# The package and the test helpers as testthat sources them (the helpers
# read files relative to the test directory while they are sourced); the
# certification functions of the helpers run from the project root.
.certification_helpers <- function() {

  pkgload::load_all(".", quiet = TRUE, helpers = FALSE, export_all = TRUE)
  old_testthat <- Sys.getenv("TESTTHAT", unset = NA_character_)
  on.exit(if (is.na(old_testthat)) Sys.unsetenv("TESTTHAT") else Sys.setenv(TESTTHAT = old_testthat), add = TRUE)
  Sys.setenv(TESTTHAT = "true")
  helpers <- new.env(parent = asNamespace("BayesTools"))
  invisible(testthat::source_test_helpers(file.path("tests", "testthat"), env = helpers))
  helpers
}

.certification_keys <- function(cases, helpers = .certification_helpers()) {
  stats::setNames(lapply(cases, function(case) {
    helpers$bayestools_certification_key(case, helpers, test_dir = file.path("tests", "testthat"))
  }), cases)
}

# One case in this process: the tests of its file under the certification
# profile, then the marker when they passed cleanly and the key is unchanged.
.run_certification_worker <- function(case, reporter) {

  Sys.setenv(BAYESTOOLS_TEST_PROFILE = "certification")
  helpers    <- .certification_helpers()
  key_before <- .certification_keys(case, helpers)[[1L]]
  unlink(helpers$bayestools_certification_marker_file(case))

  snapshot_root  <- file.path("tests", "testthat")
  snapshot_state <- bayestools_snapshot_state(snapshot_root)
  if (length(bayestools_snapshot_unreadable(snapshot_state)) > 0L) {
    stop("Cannot read the snapshot files, so a run could not restore them; no tests ran:\n- ",
         paste(bayestools_snapshot_unreadable(snapshot_state), collapse = "\n- "), call. = FALSE)
  }
  results <- tryCatch(
    as.data.frame(devtools::test(filter = paste0("^certification-", case, "$"),
                                 stop_on_failure = FALSE, reporter = reporter)),
    finally = {
      unreviewed <- bayestools_quarantine_snapshots(snapshot_root, snapshot_state)
    }
  )
  helpers   <- .certification_helpers()
  key_after <- .certification_keys(case, helpers)[[1L]]

  problems <- c(
    if (nrow(results) == 0L) "no tests ran",
    if (sum(results[["failed"]]) > 0L) paste(sum(results[["failed"]]), "failed expectation(s)"),
    if (any(results[["error"]])) "an error",
    if (any(results[["skipped"]])) "a skip",
    if (sum(results[["warning"]]) > 0L) paste(sum(results[["warning"]]), "warning(s)"),
    if (length(unlist(unreviewed)) > 0L) "unreviewed snapshot changes",
    if (!identical(key_before, key_after)) "its key changed during the run"
  )
  if (length(problems) > 0L) {
    stop("Certification case '", case, "' is not certified: ", paste(problems, collapse = ", "), ".", call. = FALSE)
  }
  helpers$bayestools_write_certification_marker(case, key_after)
  message("Certification case '", case, "' passed; marker written.")
  invisible(TRUE)
}

.run_certification_case <- function(case, script_path, helpers) {

  marker   <- helpers$bayestools_certification_marker_file(case)
  log_file <- file.path(dirname(marker), paste0(case, ".log"))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  # a rerun invalidates the case first, so a child that dies before its own
  # bookkeeping (a startup crash, the timeout) cannot leave the old marker
  unlink(marker)
  if (file.exists(marker)) {
    stop("Cannot remove the certification marker ", marker, ".", call. = FALSE)
  }
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  started <- Sys.time()
  # --vanilla: the child inherits this process's test controls, which a user
  # .Renviron read again would override
  status <- suppressWarnings(system2(
    rscript,
    c("--vanilla", shQuote(script_path), "certification", "--worker", case),
    stdout = log_file, stderr = log_file,
    timeout = helpers$bayestools_certification_timeout
  ))
  elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
  if (identical(as.integer(status), 0L)) {
    message("PASS ", case, " (", elapsed, " s)")
    return(TRUE)
  }
  message("FAIL ", case, " (", elapsed, " s, exit status ", status,
          if (identical(as.integer(status), 124L)) ", timed out" else "", "); log: ", log_file)
  if (file.exists(log_file)) {
    message(paste(utils::tail(readLines(log_file, warn = FALSE), 40L), collapse = "\n"))
  }
  FALSE
}

run_certification <- function(arguments, script_path, reporter) {

  options <- arguments[startsWith(arguments, "--")]
  named   <- arguments[!startsWith(arguments, "--")]
  unknown_options <- setdiff(options, c("--status", "--all", "--worker"))
  if (length(unknown_options) > 0L) {
    stop("Unknown certification option(s): ", paste(unknown_options, collapse = ", "), call. = FALSE)
  }
  helpers <- .certification_helpers()
  cases   <- names(helpers$bayestools_certification_cases())
  unknown <- setdiff(named, cases)
  if (length(unknown) > 0L) {
    stop("Unknown certification case(s): ", paste(unknown, collapse = ", "),
         ". Known cases: ", paste(cases, collapse = ", "), ".", call. = FALSE)
  }
  if ("--worker" %in% options) {
    if (length(named) != 1L) stop("A certification worker runs exactly one case.", call. = FALSE)
    return(.run_certification_worker(named, reporter))
  }

  keys   <- .certification_keys(cases, helpers)
  states <- lapply(cases, function(case) helpers$bayestools_certification_state(case, keys[[case]]))
  names(states) <- cases
  for (case in cases) {
    state <- states[[case]]
    message(sprintf(
      "%-28s %-8s %s", case, state$state,
      if (identical(state$state, "stale")) {
        paste0(length(state$changed), " changed: ", paste(utils::head(state$changed, 5L), collapse = ", "),
               if (length(state$changed) > 5L) ", ...")
      } else if (identical(state$state, "current")) {
        paste0("passed ", state$completed_at)
      } else {
        "never passed"
      }
    ))
  }
  if ("--status" %in% options) {
    return(invisible(states))
  }

  selected <- if (length(named) > 0L) {
    named
  } else if ("--all" %in% options) {
    cases
  } else {
    cases[vapply(states, function(state) !identical(state$state, "current"), logical(1))]
  }
  if (length(selected) == 0L) {
    message("Every certification case is current.")
    return(invisible(TRUE))
  }
  passed <- vapply(selected, .run_certification_case, logical(1), script_path = script_path, helpers = helpers)
  if (!all(passed)) {
    stop("Certification failed for: ", paste(selected[!passed], collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
