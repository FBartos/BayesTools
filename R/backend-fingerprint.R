# Cached fingerprint of the implementation loaded in this R session.
.fit_backend_fingerprint_cache <- new.env(parent = emptyenv())


# Resolve the package tree that owns the loaded BayesTools namespace.
.fit_backend_root <- function() {

  root <- tryCatch(
    getNamespaceInfo(asNamespace("BayesTools"), "path"),
    error = function(e) ""
  )
  if (!is.character(root) || length(root) != 1L || is.na(root) ||
      !nzchar(root) || !file.exists(file.path(root, "DESCRIPTION"))) {
    root <- system.file(package = "BayesTools")
  }
  if (!is.character(root) || length(root) != 1L || is.na(root) ||
      !nzchar(root)) {
    return(NA_character_)
  }

  return(normalizePath(root, winslash = "/", mustWork = TRUE))
}


# Find the native library used by the loaded BayesTools namespace.
.fit_backend_dll <- function() {

  dlls <- getLoadedDLLs()
  if (!"BayesTools" %in% names(dlls)) {
    return(character())
  }
  path <- dlls[["BayesTools"]][["path"]]
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !file.exists(path)) {
    return(character())
  }

  return(normalizePath(path, winslash = "/", mustWork = TRUE))
}


# Enumerate the executable surface of a source or installed package tree.
.fit_backend_files <- function(package_root = .fit_backend_root(),
                               loaded_dll = .fit_backend_dll()) {

  if (!is.character(package_root) || length(package_root) != 1L ||
      is.na(package_root) || !nzchar(package_root)) {
    return(character())
  }
  package_root <- normalizePath(package_root, winslash = "/", mustWork = TRUE)

  r_root    <- file.path(package_root, "R")
  r_sources <- if (dir.exists(r_root)) {
    list.files(
      r_root,
      pattern    = "\\.[Rr]$",
      recursive  = TRUE,
      full.names = TRUE
    )
  } else {
    character()
  }
  if (length(r_sources) > 0L) {
    src_root  <- file.path(package_root, "src")
    src_files <- if (dir.exists(src_root)) {
      list.files(src_root, recursive = TRUE, full.names = TRUE)
    } else {
      character()
    }
    src_extensions <- c(
      "c", "cc", "cpp", "cxx", "f", "f77", "f90", "f95",
      "h", "hh", "hpp", "hxx", "inc"
    )
    src_files <- src_files[
      tolower(tools::file_ext(src_files)) %in% src_extensions |
        grepl("^Makevars", basename(src_files))
    ]
    root_files <- file.path(
      package_root,
      c(
        "DESCRIPTION", "NAMESPACE", "configure", "configure.win",
        "cleanup", "cleanup.win"
      )
    )
    files <- c(r_sources, src_files, root_files[file.exists(root_files)])
  } else {
    installed_roots <- file.path(package_root, c("R", "libs"))
    files <- unlist(lapply(installed_roots, function(path) {
      if (!dir.exists(path)) {
        return(character())
      }
      list.files(path, recursive = TRUE, full.names = TRUE)
    }), use.names = FALSE)
    root_files <- file.path(package_root, c("DESCRIPTION", "NAMESPACE"))
    files <- c(files, root_files[file.exists(root_files)])
  }
  files <- c(files, loaded_dll)
  files <- files[file.exists(files)]

  return(sort(unique(normalizePath(
    files,
    winslash = "/",
    mustWork = TRUE
  ))))
}


# Compute a stable fingerprint from an explicit package implementation.
.compute_fit_backend_fingerprint <- function(
    package_root = .fit_backend_root(),
    loaded_dll   = .fit_backend_dll()) {

  if (!is.character(package_root) || length(package_root) != 1L ||
      is.na(package_root) || !nzchar(package_root)) {
    return(NA_character_)
  }
  package_root <- normalizePath(package_root, winslash = "/", mustWork = TRUE)
  files        <- .fit_backend_files(package_root, loaded_dll)
  if (length(files) == 0L) {
    return(NA_character_)
  }
  hashes <- unname(tools::md5sum(files))
  labels <- ifelse(
    startsWith(files, paste0(package_root, "/")),
    substring(files, nchar(package_root) + 2L),
    paste0("loaded-dll/", basename(files))
  )
  input <- c(
    "BayesTools-fit-backend-schema:1",
    paste0(labels, ":", hashes)
  )
  normalized <- tempfile("BayesTools-fit-backend-", fileext = ".txt")
  on.exit(unlink(normalized), add = TRUE)
  writeLines(input, normalized, useBytes = TRUE)

  return(unname(tools::md5sum(normalized)))
}


# Clear memoized state after reloading BayesTools during development.
.clear_fit_backend_fingerprint_cache <- function() {

  rm(
    list  = ls(envir = .fit_backend_fingerprint_cache),
    envir = .fit_backend_fingerprint_cache
  )
  return(invisible(TRUE))
}


# Freeze the implementation fingerprint that belongs to the loaded namespace.
.freeze_fit_backend_fingerprint <- function(
    value = .compute_fit_backend_fingerprint()) {

  if (!exists(
    "value",
    envir    = .fit_backend_fingerprint_cache,
    inherits = FALSE
  )) {
    assign(
      "value",
      value,
      envir = .fit_backend_fingerprint_cache
    )
  }

  return(invisible(get(
    "value",
    envir    = .fit_backend_fingerprint_cache,
    inherits = FALSE
  )))
}


#' Fingerprint the loaded BayesTools fitting backend
#'
#' @description
#' Returns a deterministic fingerprint of the BayesTools implementation that
#' owns the currently loaded namespace. Downstream packages can store this
#' value with serialized fits and invalidate caches when the generic fitting
#' backend changes.
#'
#' @details
#' The fingerprint covers the package's executable R and native surface and an
#' explicit fit-backend schema version. It is captured when the namespace loads
#' and memoized for the R session so source files edited afterward cannot
#' silently restamp fits produced by the already loaded implementation.
#' Development workflows should reload BayesTools after changing its backend.
#'
#' The scope is intentionally conservative: any executable R/native or package
#' build artifact change can change the fingerprint. Installed lazy-load
#' databases and native libraries make the value build- and platform-local; it
#' is not a cross-platform cache key. External runtimes and dependencies such as
#' JAGS are outside this BayesTools implementation fingerprint.
#'
#' @return A length-one MD5 character string, or `NA_character_` when the
#' loaded implementation cannot be resolved.
#'
#' @export
fit_backend_fingerprint <- function() {

  .freeze_fit_backend_fingerprint()

  return(get(
    "value",
    envir    = .fit_backend_fingerprint_cache,
    inherits = FALSE
  ))
}
