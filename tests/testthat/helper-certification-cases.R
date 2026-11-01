# Certification cases: deeper numerical verification outside the routine lanes
# and outside `all`. A case is the file test-certification-<case>.R, gated by
# skip_if_not_test_profile("certification"); it may fit models live.
# `tools/test-profile.R certification` runs each case in its own R process and,
# after a clean pass (no failure, error, warning, or skip), writes a marker
# holding the case's key (bayestools_certification_key()). A case whose key
# matches its marker is current and runs again only when requested by name or
# with --all. 'native' marks a case that runs native code (JAGS fits), whose
# key includes the native sources and the JAGS version.
bayestools_certification_cases <- function() {
  list(
    `density-kernels` = list(
      description = "Expensive reference sweeps of prior densities, ordinates, region grids, and quantile estimators moved out of the unit lane.",
      native      = FALSE
    ),
    `fit-assertions` = list(
      description = "Live-fit assertions of the post-fit interface, run when the code they reach changes.",
      native      = TRUE
    ),
    `hypothesis-bayes-factors` = list(
      description = "Point and region Bayes factors of hypothesis_BF() against a conjugate posterior.",
      native      = TRUE
    ),
    `marginal-likelihood` = list(
      description = "Bridge-sampling log marginal likelihoods of Gaussian models against exact integrals.",
      native      = TRUE
    ),
    `model-averaging` = list(
      description = "Model probabilities, inclusion Bayes factors, and averaged posteriors against exact enumeration.",
      native      = TRUE
    )
  )
}

# The limit of one case's process, in seconds.
bayestools_certification_timeout <- 3600

bayestools_certification_file <- function(case, test_dir = testthat::test_path()) {
  file.path(test_dir, paste0("test-certification-", case, ".R"))
}

bayestools_certification_marker_file <- function(case,
                                                 test_files_dir = Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")) {
  if (!nzchar(test_files_dir)) {
    stop("BAYESTOOLS_TEST_FILES_DIR must name the test cache directory.", call. = FALSE)
  }
  file.path(test_files_dir, "certification", paste0(case, ".txt"))
}

# The functions of the test helpers that 'references' (of the case file)
# reach: the closure under the names their code uses. Test code resolves a
# name to a helper before the package, so a helper that masks a package
# function stands for it. A helper is a function of 'helper_env' not enclosed
# by the package namespace: devtools::test() sources the helpers into the
# attached package environment, which also holds copies of the package's
# functions.
.bayestools_certification_helper_closure <- function(references, helper_env,
                                                     namespace = asNamespace("BayesTools")) {
  helper_names <- Filter(function(name) {
    value <- get(name, envir = helper_env, inherits = FALSE)
    is.function(value) && !identical(environment(value), namespace)
  }, ls(helper_env, all.names = TRUE))

  reached  <- character()
  used     <- references$used
  strings  <- references$strings
  pending  <- intersect(used, helper_names)
  while (length(pending) > 0L) {
    reached <- c(reached, pending)
    found   <- lapply(pending, function(name) {
      .test_cache_code_references(get(name, envir = helper_env, inherits = FALSE))
    })
    new_used <- unlist(lapply(found, `[[`, "used"))
    used     <- union(used, new_used)
    strings  <- union(strings, unlist(lapply(found, `[[`, "strings")))
    pending  <- setdiff(intersect(new_used, helper_names), reached)
  }

  list(helpers = sort(unique(reached), method = "radix"), used = used, strings = strings)
}

# The key of a certification case, as named hashes: the case file; the test
# helper functions and the package objects (.test_cache_reached_package_objects())
# that its code and the package's load hooks (which every run executes)
# reach; the package DESCRIPTION without its version
# (.test_cache_description_md5()); the R version; and, for a case that runs
# native code, the native sources and the JAGS version. 'helper_env' holds the
# test helpers as testthat sources them.
bayestools_certification_key <- function(case, helper_env,
                                         test_dir = testthat::test_path(),
                                         namespace = asNamespace("BayesTools")) {
  cases <- bayestools_certification_cases()
  if (!case %in% names(cases)) {
    stop("Unknown certification case: ", case, call. = FALSE)
  }
  case_file <- bayestools_certification_file(case, test_dir)
  if (!file.exists(case_file)) {
    stop("The certification case '", case, "' has no file ", basename(case_file), ".", call. = FALSE)
  }

  references <- .test_cache_code_references(parse(case_file, keep.source = FALSE, encoding = "UTF-8"))
  closure    <- .bayestools_certification_helper_closure(references, helper_env, namespace)
  roots      <- c(setdiff(closure$used, closure$helpers), ".onLoad", ".onAttach")
  reached    <- .test_cache_reached_package_objects(roots, closure$strings, namespace)

  helper_hashes <- vapply(closure$helpers, function(name) {
    .test_cache_object_hash(get(name, envir = helper_env, inherits = FALSE))
  }, character(1))
  package_objects <- lapply(reached, get, envir = namespace, inherits = FALSE)
  package_hashes  <- vapply(package_objects, .test_cache_object_hash, character(1))
  package_prefix  <- ifelse(vapply(package_objects, is.function, logical(1)), "package_fn_", "package_obj_")

  key <- c(
    case_file   = .test_cache_file_md5(case_file),
    description = .test_cache_description_md5(file.path(test_dir, "..", "..", "DESCRIPTION")),
    r_version   = R.version.string,
    stats::setNames(helper_hashes, paste0("helper_fn_", closure$helpers)),
    stats::setNames(package_hashes, paste0(package_prefix, reached))
  )
  if (isTRUE(cases[[case]]$native)) {
    native_files <- .test_cache_package_src_files(
      files    = c("Makevars.in", "Makevars.win", "Makevars.win.common", "Makevars.ucrt"),
      patterns = "\\.(c|cc|h)$"
    )
    key <- c(
      key,
      jags_version = as.character(rjags::jags.version()),
      vapply(native_files, .test_cache_file_md5, character(1))
    )
  }
  key
}

# The state of a case against its marker: "current" when the marker holds the
# key, "stale" with the names of the changed, added, and removed entries, or
# "missing" without a marker.
bayestools_certification_state <- function(case, key,
                                           marker_file = bayestools_certification_marker_file(case)) {
  marker <- .test_cache_marker_values(marker_file)
  if (length(marker) == 0L) {
    return(list(state = "missing", changed = character(), completed_at = NA_character_))
  }

  recorded <- marker[startsWith(names(marker), "key_")]
  names(recorded) <- substring(names(recorded), 5L)
  shared  <- intersect(names(key), names(recorded))
  changed <- c(
    shared[key[shared] != recorded[shared]],
    setdiff(names(key), names(recorded)),
    setdiff(names(recorded), names(key))
  )
  completed_at <- if ("completed_at" %in% names(marker)) unname(marker[["completed_at"]]) else NA_character_

  list(
    state        = if (length(changed) == 0L) "current" else "stale",
    changed      = sort(changed, method = "radix"),
    completed_at = completed_at
  )
}

bayestools_write_certification_marker <- function(case, key,
                                                  marker_file = bayestools_certification_marker_file(case)) {
  dir.create(dirname(marker_file), recursive = TRUE, showWarnings = FALSE)
  if (any(grepl("[\r\n]", key))) {
    stop("A certification key value spans lines.", call. = FALSE)
  }
  lines <- c(
    paste0("name: ", case),
    paste0("completed_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("key_", names(key), ": ", unname(key))
  )
  temporary <- paste0(marker_file, ".tmp")
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, marker_file)) {
    unlink(temporary)
    stop("Cannot write the certification marker ", marker_file, ".", call. = FALSE)
  }
  invisible(marker_file)
}
