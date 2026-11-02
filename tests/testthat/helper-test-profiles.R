# The lanes of `all`; the certification lane (helper-certification-cases.R)
# runs only when it is named.
bayestools_known_test_profiles <- c("unit", "fixture", "visual", "visual-fixture", "fit", "cran", "meta")


bayestools_quiet_llm_reporter <- function(...) {

  reporter_class <- R6::R6Class(
    classname = "BayesToolsQuietLlmReporter",
    inherit   = testthat::LlmReporter,
    public    = list(
      add_result = function(context, test, result) {

        if (self$is_full()) {
          return(invisible())
        }
        if (inherits(result, "expectation_skip")) {
          self$n_skip <- self$n_skip + 1L
          return(invisible())
        }

        super$add_result(context, test, result)
      }
    )
  )

  reporter_class$new(...)
}

bayestools_test_profile_contexts <- list(
  unit = c(
    "model-averaging-contracts",
    "density-display-coherence",
    "backend-fingerprint",
    "bounded-plot-transformations",
    "distributions-mpoint",
    "distributions-point",
    "distributions-tools",
    "distributions-weightfunctions",
    "draws-metadata",
    "exported-schemas-and-aliases",
    "factor-interaction-coefficients",
    "formula-contribution-state",
    "JAGS-formula-selected-replay",
    "JAGS-formula-point-replay",
    "priors-bound-factor-rng",
    "formula-numerical-measures",
    "harrell-davis-quantile",
    "hypothesis-ast",
    "hypothesis-BF-inference",
    "hypothesis-BF",
    "hypothesis-BF-model",
    "hypothesis-ordered",
    "hypothesis-BF-parser-adversarial",
    "hypothesis-prior-region-grid",
    "interpret",
    "JAGS-bridge-compiler",
    "JAGS-bridge-marginal-random",
    "JAGS-bridgesampling-wrapper",
    "JAGS-convergence",
    "JAGS-deterministic-nodes",
    "JAGS-diagnostic-plot-data",
    "JAGS-diagnostics-controls",
    "JAGS-draw-geometry",
    "JAGS-draws-view",
    "JAGS-fit-contract",
    "JAGS-fit-settings",
    "JAGS-formula-coefficient-density",
    "JAGS-formula-coefficient-arithmetic",
    "JAGS-formula-contrasts",
    "JAGS-formula-default-priors",
    "JAGS-formula-fixed-design",
    "JAGS-formula-fixed-evaluation",
    "JAGS-formula-expressions",
    "JAGS-formula-random-design",
    "JAGS-formula-random-priors",
    "JAGS-formula-random-bridge",
    "JAGS-formula-random-evaluation",
    "JAGS-formula-random-scaling",
    "JAGS-formula-powers",
    "JAGS-formula-prediction-targets",
    "JAGS-formula-predictor-basis",
    "JAGS-indexed-parameters",
    "JAGS-lkj-cholesky",
    "JAGS-marginal-distributions",
    "JAGS-marglik-missing-monitors",
    "JAGS-parameter-catalog",
    "JAGS-parameter-coordinates",
    "JAGS-posterior-extraction",
    "JAGS-random-effects-compile",
    "JAGS-runtime-cache",
    "JAGS-runtime-setup",
    "JAGS-selection-inits",
    "JAGS-structured-rho-support",
    "JAGS-random-effect-scaling",
    "marginal-inference-conditioning",
    "model-averaging-edge-cases",
    "model-averaging-plots-edge-cases",
    "native-registration",
    "native-range",
    "parameter-labels",
    "prior-density-ordinate",
    "prior-ordered",
    "prior-ordered-direct-display",
    "priors-automatic-transformation",
    "priors-coverage",
    "priors-density-numeric",
    "priors-informed",
    "priors-linear-density",
    "priors-native-missing",
    "priors-nonlocal",
    "priors-numerical",
    "priors-plot-data",
    "priors-print",
    "priors-tools",
    "random-allocation-inclusion",
    "random-mean-parameterization",
    "random-parameterization",
    "random-effects-correlation-draws",
    "random-effects-diagonal-factor",
    "random-effects-factor-units",
    "random-effects-independent-scalability",
    "random-effects-marginal-covariance",
    "random-effects-marginal-update",
    "random-effects-memory",
    "random-effects-structured-local",
    "random-effects-summary-posterior",
    "robma-facing-exports",
    "selection-kernels",
    "simplex-finite-total",
    "summary-tables",
    "summary-tables-helpers",
    "tools-evaluation",
    "tools-input",
    "transformation-point-images",
    "weightfunction-plot-analytic",
    "weightfunction-redesign"
  ),
  fixture = c(
    "density-display-producers",
    "expression-priors-fixture",
    "formula-contribution-fixture",
    "fixture-integrity",
    "hypothesis-BF-fixture",
    "JAGS-deterministic-nodes-fixture",
    "JAGS-ensemble-tables",
    "JAGS-fit",
    "JAGS-formula-scale",
    "JAGS-formula-draws-fixture",
    "JAGS-formula",
    "JAGS-summary-tables",
    "model-averaging",
    "parameter-labels-fixture",
    "random-effects-summary-posterior-fixture",
    "selection-kernels",
    "weightfunction-redesign"
  ),
  visual = c(
    "JAGS-ensemble-plots",
    "model-averaging-plots",
    "priors",
    "priors-density",
    "priors-plot",
    "priors-print"
  ),
  `visual-fixture` = c(
    "JAGS-diagnostic-plots",
    "JAGS-ensemble-plots",
    "JAGS-marginal-distributions",
    "model-averaging-plots"
  ),
  fit = c(
    "00-model-fits",
    "fixture-integrity"
  ),
  cran = c(
    "cran-smoke"
  ),
  meta = c(
    "test-layout-policy",
    "fixture-catalog-static",
    "snapshot-path-portability",
    "interactive-test-runner",
    "reference-table-review",
    "stochastic-reference-policy"
  ),
  certification = c(
    "certification-count-marginal-likelihood",
    "certification-density-kernels",
    "certification-fit-assertions",
    "certification-hypothesis-bayes-factors",
    "certification-marginal-likelihood",
    "certification-model-averaging",
    "certification-predictions",
    "certification-savage-dickey-ordinates",
    "certification-scaled-priors",
    "certification-simplex-priors"
  )
)

bayestools_normalize_test_profiles <- function(profiles = NULL) {
  if (is.null(profiles) || length(profiles) == 0L) {
    profiles <- Sys.getenv("BAYESTOOLS_TEST_PROFILE", "unit")
  }

  profiles <- unlist(strsplit(as.character(profiles), "[,;[:space:]]+"))
  profiles <- tolower(profiles[nzchar(profiles)])

  if (length(profiles) == 0L) {
    profiles <- "all"
  }

  aliases <- c(
    fixtures = "fixture",
    plot = "visual",
    plots = "visual",
    snapshot = "visual",
    snapshots = "visual",
    `cached-visual` = "visual-fixture",
    cached_visual = "visual-fixture",
    heavy = "fit",
    `jags-visual` = "visual-fixture",
    jags_visual = "visual-fixture",
    slow = "fit",
    visual_fixture = "visual-fixture",
    `visual-fixtures` = "visual-fixture",
    visual_fixtures = "visual-fixture"
  )

  alias_matches <- profiles %in% names(aliases)
  profiles[alias_matches] <- unname(aliases[profiles[alias_matches]])

  if ("all" %in% profiles) {
    profiles <- bayestools_known_test_profiles
  }

  invalid <- setdiff(profiles, c(bayestools_known_test_profiles, "certification"))
  if (length(invalid) > 0L) {
    stop(
      "Unknown BAYESTOOLS_TEST_PROFILE value(s): ",
      paste(invalid, collapse = ", "),
      ". Expected one or more of: all, ",
      paste(c(bayestools_known_test_profiles, "certification"), collapse = ", "),
      call. = FALSE
    )
  }

  unique(profiles)
}

bayestools_test_profile_selected_contexts <- function(profiles = NULL) {
  unique(unlist(bayestools_test_profile_contexts[bayestools_normalize_test_profiles(profiles)]))
}

bayestools_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

bayestools_test_profile_filter <- function(profiles = NULL) {
  contexts <- bayestools_test_profile_selected_contexts(profiles)
  paste0("^(", paste(bayestools_escape_regex(contexts), collapse = "|"), ")$")
}

bayestools_test_profile <- function() {
  bayestools_normalize_test_profiles()
}

bayestools_test_profile_label <- function(profiles = bayestools_test_profile()) {
  paste(profiles, collapse = ",")
}

bayestools_test_profile_includes <- function(profiles) {
  any(bayestools_test_profile() %in% bayestools_normalize_test_profiles(profiles))
}

skip_if_not_test_profile <- function(profiles) {
  expected <- bayestools_normalize_test_profiles(profiles)

  if (!any(bayestools_test_profile() %in% expected)) {
    testthat::skip(sprintf(
      "Skipping %s-profile test file because BAYESTOOLS_TEST_PROFILE=%s.",
      paste(expected, collapse = "/"),
      bayestools_test_profile_label()
    ))
  }

  invisible(TRUE)
}

skip_if_not_visual_tests <- function() {
  skip_if_not_test_profile("visual")
  testthat::skip_if_not_installed("vdiffr")
  invisible(TRUE)
}

skip_if_not_visual_fixture_tests <- function() {
  skip_if_not_test_profile("visual-fixture")
  testthat::skip_if_not_installed("vdiffr")
  invisible(TRUE)
}

# testthat stores a snapshot without a baseline (a new SVG or snapshot file,
# or a new section of one) as if it were reviewed, and its cleanup deletes
# baselines it considers unused. The runner records the snapshot files before
# a run. Afterwards it first restores the reviewed baselines (recreating
# directories) and removes unreviewed new ones, then keeps the run's output for
# review: testthat's own '.new' candidate when it wrote one during the run (a
# changed section; testthat rewrites it on every failing run), otherwise the
# new or changed file as the candidate. A run starts only when every file was
# read (bayestools_snapshot_unreadable()); afterwards every file is attempted,
# also when another file cannot be read or written; failures are reported, and
# the run fails.
bayestools_snapshot_state <- function(root) {

  directory <- file.path(root, "_snaps")
  files <- if (dir.exists(directory)) {
    list.files(directory, recursive = TRUE, full.names = TRUE)
  } else {
    character()
  }

  stats::setNames(lapply(files, function(file) {
    tryCatch(
      list(bytes = readBin(file, "raw", file.size(file)), mtime = file.mtime(file)),
      error = function(error) list(bytes = NULL, mtime = file.mtime(file), error = conditionMessage(error))
    )
  }), files)
}

bayestools_snapshot_unreadable <- function(state) {
  unreadable <- names(state)[vapply(state, function(x) is.null(x$bytes), logical(1))]
  if (length(unreadable) == 0L) {
    return(character())
  }
  paste0(unreadable, ": ", vapply(state[unreadable], `[[`, character(1), "error"))
}

bayestools_snapshot_is_candidate <- function(files) {
  grepl("[.]new[.][^.]+$", files)
}

bayestools_snapshot_candidate_path <- function(files) {
  sub("([.][^.]+)$", ".new\\1", files)
}

bayestools_quarantine_snapshots <- function(root, before, after = bayestools_snapshot_state(root)) {

  baselines_before <- names(before)[!bayestools_snapshot_is_candidate(names(before))]
  baselines_after  <- names(after)[!bayestools_snapshot_is_candidate(names(after))]
  # an unreadable file counts as changed
  changed <- baselines_after[vapply(baselines_after, function(file) {
    is.null(after[[file]]$bytes) || !identical(after[[file]]$bytes, before[[file]]$bytes)
  }, logical(1))]
  deleted <- setdiff(baselines_before, baselines_after)

  problems <- bayestools_snapshot_unreadable(before)
  if (length(problems) > 0L) {
    problems <- paste0("reading before the run ", problems)
  }
  attempt  <- function(action, file, expr) {
    tryCatch(
      {
        force(expr)
        TRUE
      },
      error = function(error) {
        problems <<- c(problems, paste0(action, " ", file, ": ", conditionMessage(error)))
        FALSE
      }
    )
  }

  # the reviewed state first, so a failure below cannot leave unreviewed bytes
  for (file in c(changed, deleted)) {
    if (file %in% baselines_before) {
      attempt("restoring", file, {
        if (is.null(before[[file]]$bytes)) stop("its reviewed content could not be read before the run")
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        writeBin(before[[file]]$bytes, file)
      })
    } else {
      attempt("removing", file, if (!file.remove(file)) stop("the file could not be removed"))
    }
  }
  for (file in changed) {
    candidate <- bayestools_snapshot_candidate_path(file)
    written   <- after[[candidate]]
    # an unreadable candidate is judged by its modification time alone
    fresh     <- !is.null(written) && (
      is.null(before[[candidate]]) ||
        !identical(written$mtime, before[[candidate]]$mtime) ||
        (!is.null(written$bytes) && !identical(written$bytes, before[[candidate]]$bytes))
    )
    if (!fresh) {
      attempt("writing candidate", candidate, {
        if (is.null(after[[file]]$bytes)) stop("the run's output could not be read: ", after[[file]]$error)
        writeBin(after[[file]]$bytes, candidate)
      })
    }
  }

  list(changed = changed, deleted = deleted, problems = problems)
}

# testthat's reviewer lists only candidates next to an existing baseline, so
# candidates of new snapshots (no baseline) are reviewed against an empty
# placeholder baseline staged in a temporary tree (a blank SVG for figures,
# which the image views can display; never equal to the candidate, so that an
# accepted candidate is recognised by the staged baseline holding its bytes):
# accepting creates the baseline, rejecting removes the candidate, and
# skipping keeps it. A candidate that cannot be read, staged, or installed is
# kept and reported as failed.
bayestools_new_snapshot_candidates <- function(root) {

  state      <- names(bayestools_snapshot_state(root))
  candidates <- state[bayestools_snapshot_is_candidate(state)]
  baselines  <- sub("[.]new([.][^.]+)$", "\\1", candidates)
  missing    <- !file.exists(baselines)

  data.frame(candidate = candidates[missing], baseline = baselines[missing], stringsAsFactors = FALSE)
}

bayestools_review_new_snapshots <- function(root, review = testthat::snapshot_review, copy = file.copy) {

  changes <- bayestools_new_snapshot_candidates(root)
  if (nrow(changes) == 0L) {
    return(invisible(cbind(changes, status = character())))
  }

  review_root <- tempfile("bayestools-new-snapshot-review-")
  on.exit(unlink(review_root, recursive = TRUE), add = TRUE)
  snaps    <- normalizePath(file.path(root, "_snaps"), winslash = "/", mustWork = TRUE)
  relative <- substring(normalizePath(changes$baseline, winslash = "/", mustWork = FALSE), nchar(snaps) + 2L)
  staged_baseline  <- file.path(review_root, "_snaps", relative)
  staged_candidate <- bayestools_snapshot_candidate_path(staged_baseline)
  read_bytes <- function(file) {
    tryCatch(suppressWarnings(readBin(file, "raw", file.size(file))), error = function(error) NULL)
  }
  copied <- function(from, to) {
    isTRUE(tryCatch(suppressWarnings(copy(from, to, overwrite = TRUE)), error = function(error) FALSE))
  }

  changes$status  <- "skipped"
  candidate_bytes <- lapply(changes$candidate, read_bytes)
  placeholders    <- vector("list", nrow(changes))
  for (i in seq_len(nrow(changes))) {
    placeholder <- if (grepl("[.]svg$", changes$baseline[[i]], ignore.case = TRUE)) {
      charToRaw("<svg xmlns='http://www.w3.org/2000/svg' width='720' height='576' viewBox='0 0 720 576'></svg>\n")
    } else {
      raw()
    }
    if (identical(placeholder, candidate_bytes[[i]])) {
      placeholder <- c(placeholder, charToRaw("\n"))
    }
    placeholders[[i]] <- placeholder
    dir.create(dirname(staged_baseline[[i]]), recursive = TRUE, showWarnings = FALSE)
    writeBin(placeholder, staged_baseline[[i]])
    if (is.null(candidate_bytes[[i]]) || !copied(changes$candidate[[i]], staged_candidate[[i]])) {
      unlink(c(staged_baseline[[i]], staged_candidate[[i]]))
      changes$status[[i]] <- "failed to stage"
    }
  }

  review(path = review_root)

  for (i in which(changes$status == "skipped")) {
    if (file.exists(staged_candidate[[i]])) {
      next
    }
    # the decision is read from the staged baseline: the candidate's bytes
    # (accepted) or the placeholder (rejected); anything else, an unreadable
    # file included, keeps the candidate
    staged <- read_bytes(staged_baseline[[i]])
    if (!identical(staged, candidate_bytes[[i]]) && !identical(staged, placeholders[[i]])) {
      changes$status[[i]] <- "failed to read the review decision"
      next
    }
    accepted <- identical(staged, candidate_bytes[[i]])
    if (accepted) {
      dir.create(dirname(changes$baseline[[i]]), recursive = TRUE, showWarnings = FALSE)
      if (!copied(staged_baseline[[i]], changes$baseline[[i]])) {
        changes$status[[i]] <- "failed to install"
        next
      }
    }
    changes$status[[i]] <- if (accepted) "accepted" else "rejected"
    if (!isTRUE(tryCatch(file.remove(changes$candidate[[i]]), error = function(error) FALSE))) {
      changes$status[[i]] <- paste0(changes$status[[i]], "; failed to remove the candidate")
    }
  }

  failed <- grepl("failed", changes$status)
  if (any(failed)) {
    warning(
      "New-snapshot review could not complete for:\n- ",
      paste0(changes$candidate[failed], " (", changes$status[failed], ")", collapse = "\n- "),
      call. = FALSE
    )
  }

  invisible(changes)
}
