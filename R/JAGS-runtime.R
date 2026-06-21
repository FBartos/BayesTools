.JAGS_require_packages <- function(required_packages, cl = NULL){

  if(length(required_packages) == 0)
    return(invisible(logical(0)))

  required_packages <- unique(required_packages)

  if(is.null(cl)){
    package_loaded <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  }else{
    package_loaded <- vapply(required_packages, function(package){
      all(unlist(parallel::clusterCall(
        cl,
        function(package) requireNamespace(package, quietly = TRUE),
        package
      ), use.names = FALSE))
    }, logical(1))
  }

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

  invisible(package_loaded)
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
.JAGS_prior_uses_nonlocal <- function(prior){

  .JAGS_prior_uses_BayesTools_module(prior)
}
.JAGS_prior_list_uses_BayesTools_module <- function(prior_list){

  length(prior_list) > 0L && any(vapply(prior_list, .JAGS_prior_uses_BayesTools_module, logical(1)))
}
.JAGS_prior_list_uses_nonlocal <- function(prior_list){

  .JAGS_prior_list_uses_BayesTools_module(prior_list)
}

