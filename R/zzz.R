.BayesTools_private <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname){

  .BayesTools_private$module_location <- .BayesTools_module_location(libname, pkgname)
  .BayesTools_private$lib_name <- libname
  .BayesTools_private$native_routines_loaded <- FALSE
  .BayesTools_private$native_dll_path <- NULL
  .BayesTools_private$native_dll_loaded_manually <- FALSE

  if(requireNamespace("rjags", quietly = TRUE)){
    BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)
  }

  .BayesTools_load_native_routines(pkgname = pkgname, libname = libname, warn = FALSE)
}

.onUnload <- function(libpath){

  if(requireNamespace("rjags", quietly = TRUE)){
    tryCatch(
      {
        if("BayesTools" %in% rjags::list.modules()){
          rjags::unload.module("BayesTools")
        }
      },
      error = function(e) NULL
    )
  }

  tryCatch(
    library.dynam.unload("BayesTools", libpath),
    error = function(e) NULL
  )

  if(isTRUE(.BayesTools_private$native_dll_loaded_manually)){
    dll_path <- .BayesTools_private$native_dll_path
    if(!is.null(dll_path) && file.exists(dll_path)){
      tryCatch(dyn.unload(dll_path), error = function(e) NULL)
    }
  }

  .BayesTools_private$native_routines_loaded <- FALSE
  .BayesTools_private$native_dll_path <- NULL
  .BayesTools_private$native_dll_loaded_manually <- FALSE

  invisible(NULL)
}

#' Load the BayesTools JAGS module
#'
#' Loads the package-shipped JAGS extension module that provides compiled
#' distributions and helper functions used by generated BayesTools JAGS syntax.
#'
#' @param quiet logical scalar. Whether module-loading messages should be
#'   suppressed.
#' @param warn logical scalar. Whether a warning should be emitted when the
#'   module cannot be loaded.
#'
#' @return Invisibly returns \code{TRUE} when the module is loaded and
#'   \code{FALSE} otherwise.
#'
#' @export
BayesTools_load_JAGS_module <- function(quiet = TRUE, warn = TRUE){

  check_bool(quiet, "quiet", allow_NA = FALSE)
  check_bool(warn, "warn", allow_NA = FALSE)

  if(!requireNamespace("rjags", quietly = TRUE)){
    if(warn){
      warning(
        "The 'rjags' package is required to load the BayesTools JAGS module.",
        call. = FALSE
      )
    }
    return(invisible(FALSE))
  }

  loaded <- tryCatch("BayesTools" %in% rjags::list.modules(), error = function(e) FALSE)
  if(loaded){
    return(invisible(TRUE))
  }

  path <- .BayesTools_private$module_location
  if(is.null(path) || !dir.exists(path)){
    path <- .BayesTools_module_location(
      libname = .BayesTools_private$lib_name,
      pkgname = "BayesTools"
    )
    .BayesTools_private$module_location <- path
  }

  if(is.null(path) || !dir.exists(path)){
    if(warn){
      warning(
        "BayesTools JAGS module was not found in the package library. ",
        "Reinstall BayesTools after installing JAGS >= 4.3.0.",
        call. = FALSE
      )
    }
    return(invisible(FALSE))
  }

  load_error <- NULL
  tryCatch(
    rjags::load.module("BayesTools", path = path, quiet = quiet),
    error = function(e) load_error <<- conditionMessage(e)
  )

  loaded <- tryCatch("BayesTools" %in% rjags::list.modules(), error = function(e) FALSE)
  if(!loaded && warn){
    message <- paste0(
      "BayesTools JAGS module failed to load from '", path, "'. ",
      "Reinstall BayesTools after installing JAGS >= 4.3.0."
    )
    if(!is.null(load_error)){
      message <- paste0(message, " rjags error: ", load_error)
    }
    warning(message, call. = FALSE)
  }

  return(invisible(loaded))
}

.BayesTools_module_location <- function(libname, pkgname){

  arch <- if(.Platform$r_arch != "") .Platform$r_arch else ""

  if(!is.null(libname)){
    module_location <- file.path(libname, pkgname, "libs", arch)
    module_file <- file.path(module_location, paste0(pkgname, .Platform$dynlib.ext))
    if(file.exists(module_file)){
      return(normalizePath(module_location, winslash = "/", mustWork = TRUE))
    }
  }

  source_location <- .BayesTools_source_module_location(pkgname)
  if(!is.null(source_location)){
    return(source_location)
  }

  dll_path <- .BayesTools_loaded_dll_path(pkgname)
  if(is.null(dll_path)){
    return(NULL)
  }

  dirname(dll_path)
}

.BayesTools_source_module_location <- function(pkgname){

  source_path <- tryCatch(getNamespaceInfo(pkgname, "path"), error = function(e) NULL)
  if(is.null(source_path)){
    return(NULL)
  }

  arch <- if(.Platform$r_arch != "") .Platform$r_arch else ""
  paths <- unique(c(
    file.path(source_path, "src", arch),
    file.path(source_path, "src")
  ))

  for(path in paths){
    module_file <- file.path(path, paste0(pkgname, .Platform$dynlib.ext))
    if(file.exists(module_file)){
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }

  NULL
}

.BayesTools_loaded_dll_path <- function(pkgname){

  dlls <- getLoadedDLLs()
  if(!pkgname %in% names(dlls)){
    return(NULL)
  }

  dll_path <- dlls[[pkgname]][["path"]]
  if(is.null(dll_path) || !file.exists(dll_path)){
    return(NULL)
  }

  normalizePath(dll_path, winslash = "/", mustWork = TRUE)
}

.BayesTools_native_symbols <- function(){
  c(
    "BayesTools_lkj_cholesky_from_u",
    "BayesTools_lkj_corr_from_u",
    "BayesTools_lkj_log_prior_u",
    "BayesTools_lkj_alpha"
  )
}

.BayesTools_native_routines_loaded <- function(pkgname = "BayesTools"){

  if(identical(pkgname, "BayesTools") && isTRUE(.BayesTools_private$native_routines_loaded)){
    return(TRUE)
  }

  symbols <- .BayesTools_native_symbols()
  loaded <- all(vapply(symbols, is.loaded, logical(1), PACKAGE = pkgname))
  if(identical(pkgname, "BayesTools")){
    .BayesTools_private$native_routines_loaded <- loaded
  }

  loaded
}

.BayesTools_load_native_routines <- function(pkgname = "BayesTools",
                                             libname = .BayesTools_private$lib_name,
                                             warn = FALSE){

  if(isTRUE(.BayesTools_native_routines_loaded(pkgname = pkgname))){
    return(invisible(TRUE))
  }

  load_error <- NULL
  tryCatch(
    {
      library.dynam("BayesTools", pkgname, libname)
      .BayesTools_private$native_dll_path <- NULL
      .BayesTools_private$native_dll_loaded_manually <- FALSE
    },
    error = function(e) load_error <<- conditionMessage(e)
  )

  if(!isTRUE(.BayesTools_native_routines_loaded(pkgname = pkgname))){
    path <- .BayesTools_private$module_location
    if(is.null(path) || !dir.exists(path)){
      path <- .BayesTools_module_location(
        libname = .BayesTools_private$lib_name,
        pkgname = pkgname
      )
      .BayesTools_private$module_location <- path
    }
    if(!is.null(path) && dir.exists(path)){
      module_file <- file.path(path, paste0(pkgname, .Platform$dynlib.ext))
      if(file.exists(module_file)){
        tryCatch(
          {
            dyn.load(module_file)
            .BayesTools_private$native_dll_path <- normalizePath(
              module_file,
              winslash = "/",
              mustWork = TRUE
            )
            .BayesTools_private$native_dll_loaded_manually <- TRUE
          },
          error = function(e) load_error <<- conditionMessage(e)
        )
      }
    }
  }

  loaded <- .BayesTools_native_routines_loaded(pkgname = pkgname)
  if(!loaded){
    .BayesTools_private$native_dll_path <- NULL
    .BayesTools_private$native_dll_loaded_manually <- FALSE
  }
  if(!loaded && warn){
    message <- paste0(
      "BayesTools native routines failed to load from the package DLL. ",
      "Compiled LKJ helpers will be unavailable."
    )
    if(!is.null(load_error)){
      message <- paste0(message, " R loader error: ", load_error)
    }
    warning(message, call. = FALSE)
  }

  invisible(loaded)
}

.BayesTools_require_native_lkj <- function(){

  if(isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    return(invisible(TRUE))
  }

  .BayesTools_load_native_routines(pkgname = "BayesTools", warn = TRUE)
  if(!isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    stop("BayesTools native LKJ routines are not loaded.", call. = FALSE)
  }

  invisible(TRUE)
}
