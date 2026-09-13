.JAGS_validate_runtime_setup <- function(runtime_setup){

  if(!is.null(runtime_setup) && !is.function(runtime_setup)){
    stop("'runtime_setup' must be NULL or a function accepting one context argument.", call. = FALSE)
  }
  invisible(runtime_setup)
}

.JAGS_validate_runtime_cache <- function(runtime_cache){

  if(!is.null(runtime_cache) && !is.function(runtime_cache)){
    stop("'runtime_cache' must be NULL or a function accepting 'context' and 'state' arguments.", call. = FALSE)
  }
  invisible(runtime_cache)
}

# Cached values are optional computational state. Keep payloads out of callback
# environments and send each old process shard to at most one current process.
.JAGS_run_runtime_cache <- function(runtime_cache, phase, chains, cl = NULL,
                                    state = NULL){

  if(is.null(runtime_cache)) return(NULL)
  if(is.null(state)) state <- list()
  if(phase == "restore" && !is.list(state)){
    warning("Runtime cache restore is unavailable because the retained state is not a list.", call. = FALSE)
    state <- list()
  }
  processes <- if(is.null(cl)) 1L else length(cl)
  assigned <- rep(list(list()), processes)
  if(phase == "restore"){
    for(index in seq_along(state)){
      if(is.null(state[[index]])) next
      process <- (index - 1L) %% processes + 1L
      assigned[[process]] <- c(assigned[[process]], list(state[[index]]))
    }
  }
  tasks <- lapply(seq_len(processes), function(process){
    list(context = .JAGS_runtime_context(chains, processes, !is.null(cl), process, phase),
      state = if(phase == "restore") assigned[[process]] else NULL)
  })
  worker <- function(task, callback){

    messages <- character()
    value <- tryCatch(withCallingHandlers(
      callback(task$context, state = task$state),
      warning = function(condition){
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }), error = function(condition){
        messages <<- c(messages, conditionMessage(condition))
        NULL
      })
    list(value = value, messages = messages)
  }
  environment(worker) <- baseenv()
  results <- tryCatch({
    if(is.null(cl)) list(worker(tasks[[1L]], runtime_cache)) else
      parallel::clusterApply(cl, tasks, worker, callback = runtime_cache)
  }, error = function(condition){
    warning("Runtime cache ", phase, " could not be completed: ", conditionMessage(condition),
      call. = FALSE)
    NULL
  })
  if(is.null(results)) return(NULL)
  for(process in seq_along(results)){
    location <- if(is.null(cl)) "the local process" else paste("worker", process)
    for(message in results[[process]]$messages){
      warning("Runtime cache ", phase, " in ", location, ": ", message, call. = FALSE)
    }
  }
  if(phase != "capture") return(NULL)
  captured <- lapply(results, `[[`, "value")
  if(all(vapply(captured, is.null, logical(1L)))) NULL else captured
}

.JAGS_runtime_context <- function(chains, processes = 1L, parallel = FALSE,
                                  process_id = 1L, phase = "start"){

  # runjags assigns chains round-robin to its simulations, with at most one
  # simulation per worker. The first workers receive any remaining chains.
  process_chains <- if(!parallel){
    chains
  }else if(process_id == 0L){
    0L
  }else{
    chains %/% processes + as.integer(process_id <= chains %% processes)
  }
  list(
    phase          = phase,
    role           = if(!parallel) "local" else if(process_id == 0L) "coordinator" else "worker",
    chains         = as.integer(chains),
    processes      = as.integer(processes),
    process_id     = as.integer(process_id),
    process_chains = as.integer(process_chains),
    parallel       = parallel
  )
}

.JAGS_run_runtime_setup <- function(runtime_setup, chains, cl = NULL){

  if(is.null(runtime_setup)){
    return(invisible(NULL))
  }
  context <- .JAGS_runtime_context(chains,
    processes = if(is.null(cl)) 1L else length(cl),
    parallel = !is.null(cl), process_id = if(is.null(cl)) 1L else 0L)
  runtime_setup(context)
  if(!is.null(cl)){
    contexts <- lapply(seq_along(cl), function(process_id){
      .JAGS_runtime_context(chains, length(cl), TRUE, process_id)
    })
    run_setup <- function(context, setup){

      setup(context)
      NULL
    }
    environment(run_setup) <- baseenv()
    parallel::clusterApply(cl, contexts, run_setup, setup = runtime_setup)
  }
  invisible(NULL)
}

.JAGS_finish_runtime_setup <- function(runtime_setup, chains, cl = NULL,
                                        operation = "Parallel JAGS"){

  # Release workers before restoring resources to the calling process.
  if(!is.null(cl)){
    cleanup_error <- tryCatch({
      parallel::stopCluster(cl)
      NULL
    }, error = identity)
    if(!is.null(cleanup_error)){
      warning(
        operation, " worker cleanup failed: ", conditionMessage(cleanup_error),
        ". The runtime finish callback was not run.", call. = FALSE
      )
      return(invisible(NULL))
    }
  }
  if(!is.null(runtime_setup)){
    runtime_setup(.JAGS_runtime_context(chains,
      processes = if(is.null(cl)) 1L else length(cl),
      parallel = !is.null(cl), process_id = if(is.null(cl)) 1L else 0L,
      phase = "finish"))
  }
  invisible(NULL)
}

.JAGS_require_packages <- function(required_packages, cl = NULL,
                                   operation = "Parallel JAGS fitting"){

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
      operation, " is unavailable with mismatched package versions, R code, or native builds: '",
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
    values <- as.list(namespace, all.names = TRUE)
    functions <- Filter(function(value){
      is.function(value) && identical(environment(value), namespace)
    }, values)
    functions <- functions[sort(names(functions), method = "radix")]
    definitions <- vapply(functions, function(value){
      paste(deparse(value, width.cutoff = 500L,
                    control = c("keepNA", "keepInteger", "niceNames")),
            collapse = "\n")
    }, character(1))
    # Numerical rule vectors and immutable metadata are executable settings
    # too. Reject environments/pointers/closures/language anywhere, including
    # attributes, so mutable caches and installation-specific DLL records do
    # not enter this source-versus-installed comparison.
    immutable <- function(value){

      if(!typeof(value) %in% c("NULL", "logical", "integer", "double", "complex", "character", "raw", "list")) return(FALSE)
      if(is.list(value) && !all(vapply(value, immutable, logical(1L)))) return(FALSE)
      value_attributes <- attributes(value)
      is.null(value_attributes) || all(vapply(value_attributes, immutable, logical(1L)))
    }
    constants <- Filter(immutable, values)
    constants <- constants[sort(names(constants), method = "radix")]
    code_file <- tempfile("JAGS-package-code-")
    on.exit(unlink(code_file), add = TRUE)
    saveRDS(list(functions = definitions, constants = constants), code_file,
      version = 2L, compress = FALSE)
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
