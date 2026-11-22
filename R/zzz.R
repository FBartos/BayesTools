.BayesTools_private <- new.env(parent = emptyenv())

# Saved in the namespace during installation, before the binary is distributed.
.BayesTools_jags_build <- rjags::pkg.compile.flags("ALL")

.onLoad <- function(libname, pkgname){

  .BayesTools_private$module_location <- .BayesTools_module_location(libname, pkgname)
  .BayesTools_private$lib_name <- libname
  .BayesTools_private$native_routines_loaded <- FALSE
  .BayesTools_private$native_dll_path <- NULL
  .BayesTools_private$native_dll_loaded_manually <- FALSE
  loaded <- FALSE
  on.exit(if(!loaded) .onUnload(libname))

  # pkgload may source the R code without a compiled DLL during documentation.
  source_only <- is.null(.BayesTools_private$module_location) &&
    dir.exists(file.path(getNamespaceInfo(pkgname, "path"), "src"))
  if(!source_only){
    .BayesTools_load_module(quiet = TRUE)
    if(!.BayesTools_load_native_routines(pkgname = pkgname, libname = libname, warn = TRUE)){
      stop("The BayesTools native routines could not be loaded.", call. = FALSE)
    }
  }
  .freeze_fit_backend_fingerprint()
  loaded <- TRUE
}

.onUnload <- function(libpath){

  # Finalize unreferenced rjags models first: their samplers run code of this
  # DLL and must not outlive it.
  gc()

  if(requireNamespace("rjags", quietly = TRUE)){
    tryCatch(
      {
        if("BayesTools" %in% rjags::list.modules()){
          rjags::unload.module("BayesTools", quiet = TRUE)
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
  loaded <- tryCatch(
    .BayesTools_load_module(quiet),
    error = function(e){
      if(warn) warning(conditionMessage(e), call. = FALSE)
      FALSE
    }
  )
  invisible(loaded)
}

.BayesTools_load_module <- function(quiet){

  if("BayesTools" %in% rjags::list.modules()) return(TRUE)

  libname <- .BayesTools_private$lib_name
  path <- .BayesTools_module_location(libname, "BayesTools")
  if(is.null(path)){
    stop("The BayesTools JAGS module is missing. Reinstall BayesTools using JAGS 5 or later.",
      call. = FALSE)
  }
  .BayesTools_private$module_location <- path
  installed_dll <- file.path(libname, "BayesTools", "libs", .Platform$r_arch,
    paste0("BayesTools", .Platform$dynlib.ext))
  if(file.exists(installed_dll)){
    rjags::load.pkg.module("BayesTools", libname, .BayesTools_jags_build, quiet = quiet)
  }else{
    # pkgload builds in src/ instead of an installed package's libs/ directory.
    rjags::load.pkg.module("BayesTools", libname, .BayesTools_jags_build, check.only = TRUE)
    rjags::load.module("BayesTools", path = path, quiet = quiet)
  }
  TRUE
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
    "BayesTools_lkj_cholesky_from_u_ordered",
    "BayesTools_lkj_log_prior_u",
    "BayesTools_moment_d",
    "BayesTools_moment_p",
    "BayesTools_moment_q",
    "BayesTools_moment_r",
    "BayesTools_invmoment_d",
    "BayesTools_invmoment_p",
    "BayesTools_invmoment_q",
    "BayesTools_invmoment_r",
    "BayesTools_nonlocal_log_interval_mass",
    "BayesTools_nonlocal_truncated_quantile",
    "BayesTools_invgamma_d",
    "BayesTools_invgamma_p",
    "BayesTools_invgamma_q",
    "BayesTools_invgamma_r",
    "BayesTools_structured_cholesky",
    "BayesTools_draw_fingerprint",
    "BayesTools_native_range_environment",
    "BayesTools_decimal_read",
    "BayesTools_prepare_finite_readout"
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
  if(!loaded && warn){
    message <- paste0(
      "BayesTools native routines failed to load from the package DLL. ",
      "Compiled helper routines will be unavailable."
    )
    if(!is.null(load_error)){
      message <- paste0(message, " R loader error: ", load_error)
    }
    warning(message, call. = FALSE)
  }

  invisible(loaded)
}

# Stops with a clear error when the compiled routines of the package are not
# available; `label` names the routine family in the message.
.BayesTools_require_native <- function(label){

  if(isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    return(invisible(TRUE))
  }

  .BayesTools_load_native_routines(pkgname = "BayesTools", warn = TRUE)
  if(!isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    stop("BayesTools native ", label, " routines are not loaded.", call. = FALSE)
  }

  invisible(TRUE)
}

.BayesTools_require_native_lkj <- function() .BayesTools_require_native("LKJ")

.BayesTools_require_native_nonlocal <- function() .BayesTools_require_native("nonlocal-prior")

.BayesTools_require_native_invgamma <- function() .BayesTools_require_native("inverse-gamma")
