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

vignette_cache_md5 <- function(lines) {

  file <- tempfile()
  on.exit(unlink(file))
  writeLines(as.character(lines), file)
  unname(tools::md5sum(file))
}

# A chunk's key has one "<input> <md5>" line for its code without comments or
# formatting, one for each workspace object the code references (its deparsed
# value), and one for each cached object it uses (that object's key).
vignette_cache_key <- function(cache, options, envir) {

  code <- parse(text = options$code, keep.source = FALSE)
  used <- setdiff(unique(all.names(code)), options$cache_fits)
  used <- sort(used[vapply(used, exists, logical(1), envir = envir, inherits = FALSE)])

  hashes <- vapply(used, function(name) {
    if (name %in% names(cache)) {
      vignette_cache_md5(vignette_cache_read_key(cache[name])[[1]])
    } else {
      vignette_cache_md5(deparse(get(name, envir = envir, inherits = FALSE)))
    }
  }, character(1))

  c(paste(".code", vignette_cache_md5(unlist(lapply(code, deparse)))),
    paste(used, hashes))
}

vignette_cache_key_paths <- function(cache) {

  sub("\\.[Rr][Dd][Ss]$", ".key", unname(cache))
}

# Returns one key per object; NULL when the key file is missing or predates
# the "<input> <md5>" format.
vignette_cache_read_key <- function(cache) {

  lapply(vignette_cache_key_paths(cache), function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    key <- readLines(path, warn = FALSE)
    if (!length(key) || !all(grepl(" ", key))) NULL else key
  })
}

vignette_cache_write_key <- function(cache, key) {

  for (path in vignette_cache_key_paths(cache)) {
    writeLines(key, path)
  }
  invisible(cache)
}

vignette_cache_changed_inputs <- function(stored_key, key) {

  inputs <- function(lines) sub(" .*$", "", lines)
  sort(unique(c(inputs(setdiff(key, stored_key)), inputs(setdiff(stored_key, key)))))
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
      unkeyed    <- vapply(stored_key, is.null, logical(1))
      changed    <- !unkeyed & !vapply(stored_key, identical, logical(1), options$cache_key)
      inputs     <- vignette_cache_changed_inputs(unlist(stored_key[changed]), options$cache_key)

      if (any(changed) && vignette_cache_refit()) {
        message("Refitting cached vignette objects ", paste(objects, collapse = ", "),
                "; changed inputs: ", paste(inputs, collapse = ", "))
        options$error <- FALSE
        return(options)
      }
      if (any(changed)) {
        warning("Cached vignette objects were created from different inputs: ",
                paste(objects[changed], collapse = ", "), " (changed: ",
                paste(inputs, collapse = ", "),
                "). Enable refitting or delete their RDS files.", call. = FALSE)
      }
      # Caches saved before keys existed are adopted as they are.
      if (any(unkeyed) && vignette_cache_refit()) {
        vignette_cache_write_key(chunk_cache[unkeyed], options$cache_key)
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
