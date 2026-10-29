vignette_cache <- function(name, objects, packages = character(), root = "../models") {

  stopifnot(
    is.character(name), length(name) == 1L,
    is.character(objects), length(objects) > 0L,
    is.character(packages),
    is.character(root), length(root) == 1L
  )

  object_names <- names(objects)
  file_names   <- unname(objects)

  if (is.null(object_names)) {
    object_names <- rep("", length(objects))
  }

  unnamed_objects <- !nzchar(object_names)
  object_names[unnamed_objects] <- sub(
    "\\.[Rr][Dd][Ss]$",
    "",
    basename(file_names[unnamed_objects])
  )

  missing_extension <- !grepl("\\.[Rr][Dd][Ss]$", file_names)
  file_names[missing_extension] <- paste0(file_names[missing_extension], ".RDS")

  paths <- file.path(root, name, file_names)
  names(paths) <- object_names

  structure(
    paths,
    class    = "vignette_cache",
    name     = name,
    packages = packages,
    root     = root
  )
}

vignette_cache_checking <- function() {

  ("CheckExEnv" %in% search()) ||
    any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
}

vignette_cache_missing <- function(cache) {

  names(cache)[!file.exists(unname(cache))]
}

vignette_cache_missing_packages <- function(cache) {

  packages <- attr(cache, "packages")
  if (is.null(packages)) {
    packages <- character()
  }

  has_package <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  packages[!has_package]
}

vignette_cache_refit <- function() {

  isTRUE(getOption("BayesTools.vignette_refit",
                   tolower(Sys.getenv("BAYESTOOLS_VIGNETTE_REFIT")) != "false"))
}

vignette_cache_eval <- function(cache) {

  !vignette_cache_checking() &&
    !length(vignette_cache_missing_packages(cache)) &&
    (!length(vignette_cache_missing(cache)) || vignette_cache_refit())
}

# A chunk's key covers its code without comments or formatting, the values of the
# workspace objects the code references, and the keys of the cached objects it uses.
vignette_cache_key <- function(cache, options, envir) {

  code  <- parse(text = options$code, keep.source = FALSE)
  input <- lapply(code, deparse)
  used  <- setdiff(unique(all.names(code)), options$cache_fits)
  used  <- used[vapply(used, exists, logical(1), envir = envir, inherits = FALSE)]

  for (name in sort(used)) {
    value <- if (name %in% names(cache)) {
      vignette_cache_read_key(cache[name])
    } else {
      deparse(get(name, envir = envir, inherits = FALSE))
    }
    input <- c(input, name, value)
  }

  key_file <- tempfile()
  on.exit(unlink(key_file))
  writeLines(as.character(unlist(input)), key_file)
  unname(tools::md5sum(key_file))
}

vignette_cache_key_paths <- function(cache) {

  sub("\\.[Rr][Dd][Ss]$", ".key", unname(cache))
}

vignette_cache_read_key <- function(cache) {

  vapply(vignette_cache_key_paths(cache), function(path) {
    if (file.exists(path)) readLines(path, n = 1L, warn = FALSE) else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

vignette_cache_write_key <- function(cache, key) {

  for (path in vignette_cache_key_paths(cache)) {
    writeLines(key, path)
  }
  invisible(cache)
}

vignette_cache_setup <- function(cache) {

  knitr::opts_hooks$set(cache_fits = function(options) {

    if (!isTRUE(options$eval) || !vignette_cache_eval(cache)) {
      options$eval <- FALSE
      return(options)
    }

    objects <- options$cache_fits
    stopifnot(is.character(objects), length(objects) > 0L, all(objects %in% names(cache)))
    chunk_cache <- cache[objects]
    options$cache_key <- vignette_cache_key(cache, options, envir = knitr::knit_global())

    if (!length(vignette_cache_missing(chunk_cache))) {
      stored_key <- vignette_cache_read_key(chunk_cache)
      changed    <- !is.na(stored_key) & stored_key != options$cache_key

      if (any(changed) && vignette_cache_refit()) {
        message("Refitting cached vignette objects with changed inputs: ",
                paste(objects, collapse = ", "))
        options$error <- FALSE
        return(options)
      }
      if (any(changed)) {
        warning("Cached vignette objects were created from different inputs: ",
                paste(objects[changed], collapse = ", "),
                ". Enable refitting or delete their RDS files.", call. = FALSE)
      }
      # Caches saved before keys existed are adopted as they are.
      if (anyNA(stored_key) && vignette_cache_refit()) {
        vignette_cache_write_key(chunk_cache[is.na(stored_key)], options$cache_key)
      }
      vignette_cache_load(chunk_cache, envir = knitr::knit_global())
      options$eval <- FALSE
    }
    options$error <- FALSE
    options
  })

  knitr::knit_hooks$set(cache_fits = function(before, options, envir) {

    if (!before && isTRUE(options$eval)) {
      vignette_cache_save(cache[options$cache_fits], envir = envir)
      vignette_cache_write_key(cache[options$cache_fits], options$cache_key)
    }
    invisible(NULL)
  })

  invisible(cache)
}

vignette_cache_load <- function(cache, envir = parent.frame()) {

  missing <- vignette_cache_missing(cache)
  if (length(missing)) {
    stop(
      "Missing cached vignette fits: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_along(cache)) {
    assign(
      names(cache)[i],
      readRDS(file = unname(cache[i])),
      envir = envir
    )
  }

  invisible(cache)
}

vignette_cache_save <- function(cache, envir = parent.frame(), compress = "xz") {

  for (cache_dir in unique(dirname(unname(cache)))) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  has_object <- vapply(
    names(cache),
    exists,
    logical(1),
    envir    = envir,
    inherits = FALSE
  )

  if (any(!has_object)) {
    stop(
      "Missing objects to cache: ",
      paste(names(cache)[!has_object], collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_along(cache)) {
    saveRDS(
      get(names(cache)[i], envir = envir, inherits = FALSE),
      file     = unname(cache[i]),
      compress = compress
    )
  }

  invisible(cache)
}
