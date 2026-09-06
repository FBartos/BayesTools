.JAGS_require_packages <- function(required_packages, cl = NULL){

  if(length(required_packages) == 0)
    return(invisible(logical(0)))

  required_packages <- unique(required_packages)

  builds <- if(is.null(cl)){
    stats::setNames(lapply(required_packages, function(package){
      if(requireNamespace(package, quietly = TRUE)) TRUE else NULL
    }), required_packages)
  }else{
    .JAGS_package_builds(required_packages)
  }
  workers <- if(is.null(cl)){
    list()
  }else{
    # Send a self-contained probe: a worker may have an older BayesTools loaded.
    probe <- .JAGS_package_builds
    environment(probe) <- baseenv()
    parallel::clusterCall(cl, probe, required_packages)
  }
  package_loaded <- vapply(required_packages, function(package){
    !is.null(builds[[package]]) && all(vapply(workers, function(worker){
      !is.null(worker[[package]])
    }, logical(1)))
  }, logical(1))

  missing_packages <- names(package_loaded)[!package_loaded]
  if(length(missing_packages) > 0)
    stop(
      paste0(
        "Required packages are not available: '",
        paste0(missing_packages, collapse = "', '"),
        "'."
      ),
      call. = FALSE
    )

  mismatched <- required_packages[vapply(required_packages, function(package){
    any(vapply(workers, function(worker){
      !identical(worker[[package]], builds[[package]])
    }, logical(1)))
  }, logical(1))]
  if(length(mismatched) > 0L){
    stop(
      "Parallel JAGS fitting is unavailable with mismatched package versions, R code, or native builds: '",
      paste(mismatched, collapse = "', '"),
      "'. Install the parent-session builds into a library and set 'R_LIBS_USER' ",
      "to that library before starting R and its workers.",
      call. = FALSE
    )
  }

  invisible(package_loaded)
}

.JAGS_package_builds <- function(packages){

  out <- lapply(packages, function(package){
    if(!requireNamespace(package, quietly = TRUE)){
      return(NULL)
    }
    # JAGS modules may load a package DLL before its namespace loader does.
    dlls <- getLoadedDLLs()
    paths <- vapply(dlls, function(dll) dll[["path"]], character(1))
    # Development tools also register in-memory entries without a DLL file.
    paths <- paths[file.exists(paths)]
    package_path <- normalizePath(
      getNamespaceInfo(package, "path"), winslash = "/"
    )
    normalized_paths <- normalizePath(paths, winslash = "/")
    paths <- paths[names(paths) == package |
      startsWith(normalized_paths, paste0(package_path, "/libs/")) |
      startsWith(normalized_paths, paste0(package_path, "/src/"))]
    # Compare loaded definitions, not source files or lazy-load databases: the
    # latter differ between load_all() and an equivalent installed package.
    namespace <- asNamespace(package)
    functions <- Filter(function(value){
      is.function(value) && identical(environment(value), namespace)
    }, as.list(namespace, all.names = TRUE))
    functions <- functions[sort(names(functions), method = "radix")]
    definitions <- vapply(functions, function(value){
      paste(deparse(value, width.cutoff = 500L,
                    control = c("keepNA", "keepInteger", "niceNames")),
            collapse = "\n")
    }, character(1))
    code_file <- tempfile("JAGS-package-code-")
    on.exit(unlink(code_file), add = TRUE)
    saveRDS(definitions, code_file, version = 2L, compress = FALSE)
    list(
      version = as.character(utils::packageVersion(package)),
      r_code = unname(tools::md5sum(code_file)),
      dll = stats::setNames(unname(tools::md5sum(paths)), names(paths))
    )
  })
  names(out) <- packages
  out
}

.JAGS_load_modules <- function(jags_modules, cl = NULL, warn = TRUE){

  if(length(jags_modules) == 0){
    return(invisible(logical(0)))
  }

  jags_modules <- unique(jags_modules)

  if(is.null(cl)){
    loaded <- vapply(
      jags_modules,
      function(module){
        if(module == "BayesTools"){
          return(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = warn)))
        }
        if(!requireNamespace("rjags", quietly = TRUE)){
          if(warn){
            warning(
              "The 'rjags' package is required to load JAGS modules.",
              call. = FALSE
            )
          }
          return(FALSE)
        }
        if(module %in% rjags::list.modules()){
          return(TRUE)
        }
        load_error <- NULL
        tryCatch(
          rjags::load.module(module, quiet = TRUE),
          error = function(e) load_error <<- conditionMessage(e)
        )
        loaded <- module %in% rjags::list.modules()
        if(!loaded && warn){
          message <- paste0("JAGS module '", module, "' failed to load.")
          if(!is.null(load_error)){
            message <- paste0(message, " rjags error: ", load_error)
          }
          warning(message, call. = FALSE)
        }
        loaded
      },
      logical(1)
    )
  }else{
    loaded <- parallel::clusterCall(
      cl,
      function(jags_modules){
        vapply(
          jags_modules,
          function(module){
            if(module == "BayesTools"){
              return(isTRUE(BayesTools::BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))
            }
            if(!requireNamespace("rjags", quietly = TRUE)){
              return(FALSE)
            }
            if(module %in% rjags::list.modules()){
              return(TRUE)
            }
            tryCatch(rjags::load.module(module, quiet = TRUE), error = function(e) NULL)
            module %in% rjags::list.modules()
          },
          logical(1)
        )
      },
      jags_modules
    )
    loaded <- Reduce("&", loaded)
  }

  failed <- names(loaded)[!loaded]
  if(length(failed) > 0){
    stop(
      paste0(
        "Required JAGS modules failed to load: '",
        paste0(failed, collapse = "', '"),
        "'."
      ),
      call. = FALSE
    )
  }

  invisible(loaded)
}

.JAGS_prior_uses_BayesTools_module <- function(prior){

  if(!is.prior(prior)){
    return(FALSE)
  }

  if(isTRUE(prior[["distribution"]] %in% c("invgamma", "moment", "invmoment"))){
    return(TRUE)
  }

  if(is.prior.ordered(prior)){
    return(.JAGS_prior_uses_BayesTools_module(prior$total))
  }

  if(is.prior.spike_and_slab(prior)){
    return(
      .JAGS_prior_uses_BayesTools_module(.get_spike_and_slab_variable(prior)) ||
        .JAGS_prior_uses_BayesTools_module(.get_spike_and_slab_inclusion(prior))
    )
  }

  if(is.prior.mixture(prior)){
    return(any(vapply(as.list(prior), .JAGS_prior_uses_BayesTools_module, logical(1))))
  }

  if(is.prior.weightfunction(prior) && identical(prior$weights$type, "independent")){
    return(.JAGS_prior_uses_BayesTools_module(prior$weights$prior))
  }

  if(is_prior_phacking(prior)){
    return(.JAGS_prior_uses_BayesTools_module(prior$alpha))
  }

  if(is_prior_bias(prior)){
    uses_selection <- !is.null(prior$selection) && .JAGS_prior_uses_BayesTools_module(prior$selection)
    uses_phacking  <- !is.null(prior$phacking)  && .JAGS_prior_uses_BayesTools_module(prior$phacking)
    return(uses_selection || uses_phacking)
  }

  FALSE
}
.JAGS_prior_list_uses_BayesTools_module <- function(prior_list){

  length(prior_list) > 0L && any(vapply(prior_list, .JAGS_prior_uses_BayesTools_module, logical(1)))
}
