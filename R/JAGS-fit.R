#' @title Fits a 'JAGS' model
#'
#' @description A wrapper around
#' \link[runjags]{run.jags}  that simplifies fitting 'JAGS' models
#' with usage with pre-specified model part of the 'JAGS' syntax, data and list
#' of prior distributions.
#' @param model_syntax jags syntax for the model part
#' @param data list containing data to fit the model (not including data for the formulas)
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names) of parameters not specified within the
#' \code{formula_list}
#' @param formula_list named list of formulas to be added to the model
#' (names correspond to the parameter name created by each of the formula)
#' @param formula_data_list named list of data frames containing data for each formula
#' (names of the lists correspond to the parameter name created by each of the formula)
#' @param formula_prior_list named list of named lists of prior distributions
#' (names of the lists correspond to the parameter name created by each of the formula and
#' the names of the prior distribution correspond to the parameter names) of parameters specified
#' within the \code{formula}
#' @param formula_scale_list named list of named lists for standardizing continuous predictors
#' (names of the lists correspond to the parameter name created by each of the formula).
#' Each entry should be a named list where continuous predictors with \code{TRUE} values will
#' be standardized. Defaults to \code{NULL} (no standardization).
#' @param chains number of chains to be run, defaults to \code{4}
#' @param adapt number of samples used for adapting the MCMC chains, defaults to \code{500}
#' @param burnin number of burnin iterations of the MCMC chains, defaults to \code{1000}
#' @param sample number of sampling iterations of the MCMC chains, defaults to \code{4000}
#' @param thin thinning interval for the MCMC samples, defaults to \code{1}
#' @param autofit whether the models should be refitted until convergence criteria
#' specified in \code{autofit_control}. Defaults to \code{FALSE}.
#' @param autofit_control a list of arguments controlling the autofit function.
#' Possible options are:
#' \describe{
#'   \item{max_Rhat}{maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}.}
#'   \item{min_ESS}{minimum effective sample size. Defaults to \code{500}.}
#'   \item{max_error}{maximum MCMC error. Defaults to \code{0.01}.}
#'   \item{max_SD_error}{maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.}
#'   \item{max_time}{list specifying the time \code{time} and \code{units}
#'   after which the automatic fitting function is stopped. The units arguments
#'   need to correspond to \code{units} passed to \link[base]{difftime} function.}
#'   \item{max_extend}{number of times after which the automatic fitting function is stopped.}
#'   \item{sample_extend}{number of samples between each convergence check. Defaults to
#'   \code{1000}.}
#'   \item{restarts}{number of times new initial values should be generated in case the model
#'   fails to initialize. Defaults to \code{10}.}
#'   \item{check_indicators}{whether model indicator variables should be included
#'   in convergence checks. Defaults to \code{FALSE}.}
#' }
#' @param parallel whether the chains should be run in parallel \code{FALSE}
#' @param cores number of cores used for multithreading if \code{parallel = TRUE},
#'  defaults to \code{chains}
#' @param silent whether the function should proceed silently, defaults to \code{TRUE}
#' @param seed seed for random number generation
#' @param add_parameters vector of additional parameter names that should be used
#' monitored but were not specified in the \code{prior_list}
#' @param required_packages character vector specifying list of packages containing
#' JAGS models required for sampling (in case that the function is run in parallel or in
#' detached R session). Defaults to \code{NULL}.
#' @param fit a 'BayesTools_fit' object (created by \code{JAGS_fit()} function) to be
#' extended
#' @param ... additional hidden arguments
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#' }
#'
#' @return \code{JAGS_fit} returns an object of class 'runjags' and 'BayesTools_fit'.
#'
#' @seealso [JAGS_check_convergence()]
#'
#' @export JAGS_fit
#' @export JAGS_extend
#' @name JAGS_fit
NULL

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

#' @rdname JAGS_fit
JAGS_fit <- function(model_syntax, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL, formula_scale_list = NULL,
                     chains = 4, adapt = 500, burnin = 1000, sample = 4000, thin = 1,
                     autofit = FALSE, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10, check_indicators = FALSE),
                     parallel = FALSE, cores = chains, silent = TRUE, seed = NULL,
                     add_parameters = NULL, required_packages = NULL, ...){

  .check_runjags()
  dots <- list(...)

  ### check input
  .check_JAGS_syntax(model_syntax)
  JAGS_check_and_list_fit_settings(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed)
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE)
  check_char(required_packages, "required_packages", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_list(formula_list, "formula_list", allow_NULL = TRUE)
  check_list(formula_data_list, "formula_data_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))
  check_list(formula_prior_list, "formula_prior_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))
  check_list(formula_scale_list, "formula_scale_list", allow_NULL = TRUE)

  ### add formulas
  if(!is.null(formula_list)){

    # obtain settings for each formula
    formula_output <- list()
    for(parameter in names(formula_list)){
      formula_output[[parameter]] <- JAGS_formula(
        formula        = formula_list[[parameter]],
        parameter      = parameter,
        data           = formula_data_list[[parameter]],
        prior_list     = formula_prior_list[[parameter]],
        formula_scale  = if(!is.null(formula_scale_list)) formula_scale_list[[parameter]] else NULL)
    }

    # merge with the rest of the input
    prior_list     <- c(do.call(c, unname(lapply(formula_output, function(output) output[["prior_list"]]))), prior_list)
    data           <- c(do.call(c, unname(lapply(formula_output, function(output) output[["data"]]))),       data)
    formula_syntax <- paste0(lapply(formula_output, function(output) output[["formula_syntax"]]), collapse = "")
    
    # collect formula_scale information
    formula_scale_info <- lapply(formula_output, function(output) output[["formula_scale"]])
    formula_scale_info <- formula_scale_info[!sapply(formula_scale_info, is.null)]
    if(length(formula_scale_info) == 0) formula_scale_info <- NULL

    # add the formula syntax to the model syntax
    opening_bracket <- regexpr("{", model_syntax, fixed = TRUE)[1]
    syntax_start    <- substr(model_syntax, 1, opening_bracket)
    syntax_end      <- substr(model_syntax, opening_bracket + 1, nchar(model_syntax))
    model_syntax    <- paste0(syntax_start, "\n", formula_syntax, "\n", syntax_end)
  }else{
    formula_scale_info <- NULL
  }


  ### create the model call
  model_call <- list(
    model     = JAGS_add_priors(syntax = model_syntax, prior_list = prior_list),
    data      = data,
    inits     = JAGS_get_inits(prior_list, chains = chains, seed = seed),
    monitor   = c(JAGS_to_monitor(prior_list), add_parameters),
    n.chains  = chains,
    adapt     = adapt,
    burnin    = burnin,
    sample    = sample,
    thin      = thin,
    summarise = FALSE
  )

  # parallel vs. not
  if(parallel){
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    .JAGS_require_packages(required_packages, cl)
    model_call <- c(
      model_call,
      method = "rjparallel",
      cl     = list(cl)
    )
  }else{
    .JAGS_require_packages(required_packages)
    model_call <- c(
      model_call,
      method = "rjags"
    )
  }


  if(!is.null(seed)){
    set.seed(seed)
  }

  # set silent mode
  if(silent){
    user_silent.jags    <- runjags::runjags.getOption("silent.jags")
    user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
    on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  }

  start_time <- Sys.time()
  # special fitting procedure for JASP
  # singlcore interrupted fits allowing for bar progression
  if(isTRUE(dots[["is_JASP"]])){

    model_call_adapt  <- model_call
    model_call_adapt[["sample"]] <- 1 # at least one burnin & adapt need to be specified

    # adapt & burnin
    .JASP_progress_bar_start(n = 1, label = paste0(if(!is.null(dots[["is_JASP_prefix"]])) paste0(dots[["is_JASP_prefix"]], ": "), "Adapting and burnin the model"))
    fit <- tryCatch(do.call(runjags::run.jags, model_call_adapt), error = function(e) e)
    .JASP_progress_bar_tick()

    # sample
    .JASP_progress_bar_start(n = 5, label = paste0(if(!is.null(dots[["is_JASP_prefix"]])) paste0(dots[["is_JASP_prefix"]], ": "), "Sampling the model"))
    for(i in 1:5){
      if(!inherits(fit, "error")){
        fit <- tryCatch(runjags::extend.jags(fit, burnin = 0, sample = floor((model_call[["sample"]])/5)), error = function(e)e)
        .JASP_progress_bar_tick()
      }
    }

  }else{
    if(is.null(autofit_control[["restarts"]])){
      fit <- tryCatch(do.call(runjags::run.jags, model_call), error = function(e) e)
    }else{
      for(i in 1:autofit_control[["restarts"]]){
        fit <- tryCatch(do.call(runjags::run.jags, model_call), error = function(e) e)
        if(!inherits(fit, "error")){
          break
        }else{
          # restart with different inits
          model_call$inits <- JAGS_get_inits(prior_list, chains = chains, seed = if(!is.null(seed)) seed + i)
        }
      }
    }
  }



  if(inherits(fit, "error") & !silent)
    warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)

  if(autofit && !inherits(fit, "error")){

    converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)
    itteration <- 1

    if(!converged && isTRUE(dots[["is_JASP"]]))
      .JASP_progress_bar_start(n = if (!is.null(autofit_control[["max_extend"]])) autofit_control[["max_extend"]] else 10, label = paste0(if(!is.null(dots[["is_JASP_prefix"]])) paste0(dots[["is_JASP_prefix"]], ": "), "Extending the model (autofit)"))

    while(!converged){

      if(!is.null(autofit_control[["max_time"]]) && difftime(Sys.time(), start_time, units = autofit_control[["max_time"]][["unit"]]) > autofit_control[["max_time"]][["time"]]){
        if(!silent){
          attr(fit, "warning") <- "The automatic model fitting was terminated due to the 'max_time' constraint."
          warning(attr(fit, "warning"), immediate. = TRUE)
        }
        break
      }
      if(!is.null(autofit_control[["max_extend"]]) && itteration > autofit_control[["max_extend"]]){
        if(!silent){
          attr(fit, "warning") <- "The automatic model fitting was terminated due to the 'max_extend' constraint."
          warning(attr(fit, "warning"), immediate. = TRUE)
        }
        break
      }

      fit <- tryCatch(runjags::extend.jags(fit, sample = autofit_control[["sample_extend"]]), error = function(e)e)

      if(inherits(fit, "error")){
        if(!silent)
          warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)
        break
      }

      fit <- runjags::add.summary(fit)

      converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)
      itteration <- itteration + 1

      if(isTRUE(dots[["is_JASP"]]))
        .JASP_progress_bar_tick()
    }
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "required_packages") <- required_packages
  if(!is.null(formula_scale_info)){
    # Keep formula_scale as a nested list keyed by parameter name
    # Each element contains the scaling info for that parameter's predictors
    attr(fit, "formula_scale") <- formula_scale_info
  }

  class(fit) <- c(class(fit), "BayesTools_fit")

  return(fit)
}

#' @rdname JAGS_fit
JAGS_extend <- function(fit, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10, check_indicators = FALSE),
                        parallel = FALSE, cores = NULL, silent = TRUE, seed = NULL){

  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a 'BayesTools_fit'")

  # extract fitting information
  prior_list        <- attr(fit, "prior_list")
  model_syntax      <- attr(fit, "model_syntax")
  required_packages <- attr(fit, "required_packages")
  formula_scale     <- attr(fit, "formula_scale")
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)

  # parallel vs. not
  if(parallel){
    if(is.null(cores)){
      cores <- length(fit[["mcmc"]])
    }
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    .JAGS_require_packages(required_packages, cl)
    refit_call <- list(
      runjags.object = fit,
      sample         = autofit_control[["sample_extend"]],
      method         = "rjparallel",
      cl             = cl,
      summarise      = FALSE
    )
  }else{
    .JAGS_require_packages(required_packages)
    refit_call <- list(
      runjags.object = fit,
      sample         = autofit_control[["sample_extend"]],
      method         = "rjags",
      summarise      = FALSE
    )
  }


  if(!is.null(seed)){
    set.seed(seed)
  }

  # set silent mode
  if(silent){
    user_silent.jags    <- runjags::runjags.getOption("silent.jags")
    user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
    on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  }

  start_time <- Sys.time()
  itteration <- 0
  converged  <- FALSE

  while(!converged){

    if(!is.null(autofit_control[["max_time"]]) && difftime(Sys.time(), start_time, units = autofit_control[["max_time"]][["unit"]]) > autofit_control[["max_time"]][["time"]]){
      if(!silent){
        attr(fit, "warning") <- "The automatic model fitting was terminated due to the 'max_time' constraint."
        warning(attr(fit, "warning"), immediate. = TRUE)
      }
      break
    }
    if(!is.null(autofit_control[["max_extend"]]) && itteration >= autofit_control[["max_extend"]]){
      if(!silent){
        attr(fit, "warning") <- "The automatic model fitting was terminated due to the 'max_extend' constraint."
        warning(attr(fit, "warning"), immediate. = TRUE)
      }
      break
    }

    fit <- tryCatch(do.call(runjags::extend.jags, refit_call), error = function(e)e)

    if(inherits(fit, "error")){
      if(!silent)
        warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)

      break
    }

    converged <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)

    # update the refit call
    if(!converged){
      itteration <- itteration + 1
      refit_call$runjags.object <- fit
    }
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "required_packages") <- required_packages
  if(!is.null(formula_scale)){
    attr(fit, "formula_scale") <- formula_scale
  }

  class(fit) <- c(class(fit), "BayesTools_fit")

  return(fit)
}


#' @title Assess convergence of a runjags model
#'
#' @description Checks whether the supplied \link[runjags]{runjags-package} model
#' satisfied convergence criteria.
#' @param fit a runjags model
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param max_Rhat maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}.
#' @param min_ESS minimum effective sample size. Defaults to \code{500}.
#' @param max_error maximum MCMC error. Defaults to \code{0.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of additional parameter names that should be used
#' (only allows removing last, fixed, omega element if omega is tracked manually).
#' @param fail_fast whether the function should stop after the first failed convergence check.
#' @param check_indicators whether model indicator variables should be included
#' in convergence checks. Defaults to \code{FALSE}.
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#' JAGS_check_convergence(fit, priors_list)
#' }
#' @return \code{JAGS_check_convergence} returns a boolean
#' indicating whether the model converged or not, with an
#' attribute 'errors' carrying the failed convergence checks (if any).
#'
#' @seealso [JAGS_fit()]
#' @export
JAGS_check_convergence <- function(fit, prior_list, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, add_parameters = NULL, fail_fast = FALSE, check_indicators = FALSE){

  # check input
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  check_list(prior_list, "prior_list", allow_NULL = TRUE)
  if(!is.null(prior_list) && any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_real(max_Rhat,     "max_Rhat",     lower = 1, allow_NULL = TRUE)
  check_real(min_ESS,      "min_ESS",      lower = 0, allow_NULL = TRUE)
  check_real(max_error,    "max_error",    lower = 0, allow_NULL = TRUE)
  check_real(max_SD_error, "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE)
  check_bool(check_indicators, "check_indicators")

  # extract samples and parameter information
  mcmc_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  mcmc_samples      <- do.call(rbind, mcmc_samples_list)
  
  # build remove_parameters list: point priors, spike priors, indicators, inclusions
  remove_params <- c(
    # point priors
    names(prior_list)[sapply(prior_list, is.prior.point)],
    # mixture with single point prior
    names(prior_list)[sapply(prior_list, function(p) {
      is.prior.mixture(p) && length(p) == 1 && is.prior.point(p[[1]])
    })],
    # add_parameters that should be excluded
    add_parameters
  )
  
  # use helper to remove auxiliary parameters
  cleaned <- .remove_auxiliary_parameters(mcmc_samples, prior_list, remove_params)
  mcmc_samples <- cleaned$model_samples
  
  # remove auxiliary inclusion probabilities and, by default, model indicators
  indicator_cols <- grepl("_indicator", colnames(mcmc_samples))
  inclusion_cols <- grepl("_inclusion", colnames(mcmc_samples))
  mcmc_samples <- mcmc_samples[, !(inclusion_cols | (!check_indicators & indicator_cols)), drop = FALSE]
  
  if(ncol(mcmc_samples) == 0){
    return(TRUE)
  }
  
  # convert back to mcmc.list for convergence checks
  n_chains <- length(mcmc_samples_list)
  samples_per_chain <- nrow(mcmc_samples) / n_chains
  mcmc_samples_list_cleaned <- lapply(1:n_chains, function(i) {
    start_idx <- (i - 1) * samples_per_chain + 1
    end_idx <- i * samples_per_chain
    coda::as.mcmc(mcmc_samples[start_idx:end_idx, , drop = FALSE])
  })
  mcmc_samples <- coda::as.mcmc.list(mcmc_samples_list_cleaned)

  ### check the convergence
  fails <- NULL

  # assess R-hat
  if(!is.null(max_Rhat)){
    if(length(fit$mcmc) == 1){
      warning("Only one chain was run. R-hat cannot be computed.", immediate. = TRUE)
    }else{
      temp_Rhat <- coda::gelman.diag(mcmc_samples, multivariate = FALSE, autoburnin = FALSE)$psrf
      temp_Rhat[is.na(temp_Rhat)] <- 1
      temp_Rhat <- max(temp_Rhat)
      if(temp_Rhat > max_Rhat){
        fails <- c(fails, paste0("R-hat ", round(temp_Rhat, 3), " is larger than the set target (", max_Rhat, ")."))
        if(fail_fast){
          return(FALSE)
        }
      }
    }
  }

  if(!is.null(min_ESS)){
    temp_ESS <- coda::effectiveSize(mcmc_samples)
    temp_ESS[is.nan(temp_ESS) | temp_ESS == 0] <- Inf
    temp_ESS <- min(temp_ESS)
    if(temp_ESS < min_ESS){
      fails <- c(fails, paste0("ESS ", round(temp_ESS), " is lower than the set target (", min_ESS, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  # compute the MCMC error and & SD error
  if(!(is.null(max_error) && is.null(max_SD_error))){
    temp_summary <- summary(mcmc_samples, quantiles = NULL)$statistics
    if(is.null(dim(temp_summary))){
      temp_summary <- t(temp_summary)
    }
  }


  if(!is.null(max_error)){
    temp_error    <- temp_summary[,"Time-series SE"]
    temp_error[is.na(temp_error)] <- 0
    temp_error    <- max(temp_error)
    if(temp_error > max_error){
      fails <- c(fails, paste0("MCMC error ", round(temp_error, 5), " is larger than the set target (", max_error, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  if(!is.null(max_SD_error)){
    temp_error_SD <- temp_summary[,"Time-series SE"] / temp_summary[,"SD"]
    temp_error_SD[is.na(temp_error_SD)] <- 0
    temp_error_SD <- max(temp_error_SD)
    if(temp_error_SD > max_SD_error){
      fails <- c(fails, paste0("MCMC SD error ", round(temp_error_SD, 3), " is larger than the set target (", max_SD_error, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  converged <- length(fails) == 0
  attr(converged, "errors") <- fails
  return(converged)
}


#' @title Add 'JAGS' prior
#'
#' @description Adds priors to a 'JAGS' syntax.
#'
#' @param syntax JAGS model syntax
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#'
#' @return \code{JAGS_add_priors} returns a JAGS syntax.
#'
#' @export
JAGS_add_priors           <- function(syntax, prior_list){

  # return the original syntax in case that no prior was specified
  if(length(prior_list) == 0){
    return(syntax)
  }

  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_JAGS_syntax(syntax)

  # identify parts of the syntax
  opening_bracket <- regexpr("{", syntax, fixed = TRUE)[1]
  syntax_start    <- substr(syntax, 1, opening_bracket)
  syntax_end      <- substr(syntax, opening_bracket + 1, nchar(syntax))

  # create the priors relevant syntax
  syntax_priors <- .JAGS_add_priors.fun(prior_list)

  # merge everything back together
  syntax <- paste0(syntax_start, "\n", syntax_priors, "\n", syntax_end)

  return(syntax)
}

.JAGS_add_priors.fun       <- function(prior_list){

  syntax_priors <- ""
  syntax_attributes <- NULL

  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.factor(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(syntax_priors)
}
.JAGS_prior.simple         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # parse expressions to test
  if(.is_prior_expression(prior)){
    prior <- .prior_expression_to_character(prior)
  }

  # distribution
  syntax <- switch(
    prior[["distribution"]],
    "point"     = paste0(parameter_name," = ",prior$parameter[["location"]]),
    "normal"    = paste0(parameter_name," ~ dnorm(",prior$parameter[["mean"]],",", .JAGS_parameter_to_precision(prior$parameter[["sd"]]),")"),
    "lognormal" = paste0(parameter_name," ~ dlnorm(",prior$parameter[["meanlog"]],",", .JAGS_parameter_to_precision(prior$parameter[["sdlog"]]),")"),
    "t"         = paste0(parameter_name," ~ dt(",prior$parameter[["location"]],",", .JAGS_parameter_to_precision(prior$parameter[["scale"]]),",", prior$parameter[["df"]],")"),
    "gamma"     = paste0(parameter_name," ~ dgamma(",prior$parameter[["shape"]],",",prior$parameter[["rate"]],")"),
    "invgamma"  = paste0("inv_",parameter_name," ~ dgamma(",prior$parameter[["shape"]],",",prior$parameter[["scale"]],")"),
    "exp"       = paste0(parameter_name," ~ dexp(",prior$parameter[["rate"]],")"),
    "beta"      = paste0(parameter_name," ~ dbeta(",prior$parameter[["alpha"]],",",prior$parameter[["beta"]],")"),
    "bernoulli" = paste0(parameter_name," ~ dbern(",prior$parameter[["probability"]],")"),
    "uniform"   = paste0(parameter_name," ~ dunif(",prior$parameter[["a"]],",",prior$parameter[["b"]],")")
  )

  # add truncation
  if(!.is_prior_default_range(prior)){
    # the truncation for invgamma needs to be done in reverse since we sample from gamma
    if(prior[["distribution"]] == "invgamma"){
      syntax <- paste0(syntax, "T(",
                       ifelse(is.infinite(prior$truncation[["upper"]]^-1),"",prior$truncation[["upper"]]^-1),
                       ",",
                       ifelse(is.infinite(prior$truncation[["lower"]]^-1),"",prior$truncation[["lower"]]^-1),
                       ")")
    }else{
      syntax <- paste0(syntax, "T(",
                       ifelse(is.infinite(prior$truncation[["lower"]]),"",prior$truncation[["lower"]]),
                       ",",
                       ifelse(is.infinite(prior$truncation[["upper"]]),"",prior$truncation[["upper"]]),
                       ")")
    }
  }

  # finish the line
  syntax <- paste0(syntax, "\n")

  # transform the parameter in case of inverse-gamma
  if(prior[["distribution"]] == "invgamma"){
    syntax <- paste0(syntax, "  ", parameter_name," = pow(inv_",parameter_name,", -1)\n")
  }

  return(syntax)
}
.JAGS_prior.vector         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  if(prior[["distribution"]] != "mpoint")
    .check_vector_truncation_unsupported(prior$truncation)

  # parse expressions to test
  if(.is_prior_expression(prior)){
    prior <- .prior_expression_to_character(prior)
  }

  if(prior[["distribution"]] %in% c("mnormal", "mt")){
    # create the location/means vector the sigma matrix

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]]
    )
    par2 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["sd"]],
      "mt"      = prior$parameter[["scale"]]
    )

    par2 <- .JAGS_parameter_to_precision(par2)

    # TODO: beautify this code by specific JAGS distributions?
    if(prior[["distribution"]] == "mt"){
      # using the chisq * covariance parametrization since the mt fails with 1 df
      # (using a common df parameter as in Rouder et al. 2012)
      syntax <- paste0("prior_par1_", parameter_name, " = rep(0,", prior$parameter[["K"]], ")\n")
      syntax <- paste0(syntax, "prior_par_s_", parameter_name, " ~ dgamma(", prior$parameter[["df"]]/2, ", ", prior$parameter[["df"]]/2,")\n")
      syntax <- paste0(
        syntax,
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  prior_par2_", parameter_name, "[i,i] <- ", par2, "\n",
        "  for(j in 1:(i-1)){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "  for (j in (i+1):", prior$parameters[["K"]], "){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "}\n",
        "prior_par_z_", parameter_name, " ~ dmnorm(prior_par1_", parameter_name, ",prior_par2_", parameter_name, ")\n",
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  ", parameter_name, "[i] <- prior_par_z_", parameter_name, "[i]/sqrt(prior_par_s_", parameter_name, ") + ", par1, " \n",
        "}\n")
    }else if(prior[["distribution"]] == "mnormal"){
      syntax <- paste0("prior_par1_", parameter_name, " = rep(", par1, ",", prior$parameter[["K"]], ")\n")
      syntax <- paste0(
        syntax,
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  prior_par2_", parameter_name, "[i,i] <- ", par2, "\n",
        "  for(j in 1:(i-1)){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "  for (j in (i+1):", prior$parameters[["K"]], "){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "}\n")
      syntax <- paste0(syntax, parameter_name," ~ dmnorm(prior_par1_", parameter_name, ",prior_par2_", parameter_name, ")\n")
    }

  }else if(prior[["distribution"]] == "mpoint"){

    syntax <- paste0(
      "for(i in 1:", prior$parameters[["K"]], "){\n",
      "  ", parameter_name, "[i] = ", prior$parameter[["location"]], " \n",
      "}\n")

  }


  return(syntax)
}
.JAGS_prior.factor         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(.get_prior_factor_levels(prior), "levels", lower = 1)

  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    syntax <- paste0(
      "for(i in 1:", .get_prior_factor_levels(prior), "){\n",
      "  ", .JAGS_prior.simple(prior, paste0(parameter_name, "[i]")),
      "}\n")

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    syntax <- .JAGS_prior.vector(prior, parameter_name)

  }

  return(syntax)
}
.JAGS_prior.PP             <- function(prior){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    syntax <- .JAGS_prior.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    syntax <- .JAGS_prior.simple(prior, "PEESE")
  }

  return(syntax)
}
.JAGS_prior.weightfunction <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}
.JAGS_prior.phacking      <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}
.JAGS_prior.bias          <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}

.JAGS_selection_backend_syntax <- function(spec){

  code <- c(spec$prior_code, spec$transform_code)
  code <- code[nzchar(code)]
  if(length(code) == 0L){
    return("")
  }
  code <- sub("[\r\n]+$", "", code)

  paste0(paste0(code, collapse = "\n"), "\n")
}

.JAGS_weightfunction_component_syntax <- function(prior, component_id = NULL, global_cuts = NULL, force_one_sided = FALSE){

  J <- .weightfunction_n_bins(prior)
  syntax <- character()

  expansion     <- .weightfunction_mapping_expansion(prior, force_one_sided)
  all_cuts      <- if(is.null(global_cuts)) expansion$cuts else global_cuts
  needs_mapping <- !identical(all_cuts, .weightfunction_local_cuts(prior)) ||
    !identical(expansion$index, seq_len(J))

  omega_local <- if(is.null(component_id) && !needs_mapping) "omega" else if(is.null(component_id)) "omega_local" else paste0("omega_local_component_", component_id)
  omega_target <- if(is.null(component_id)) "omega" else paste0("omega_component_", component_id)

  if(prior$weights$type == "cumulative"){
    eta_name <- if(is.null(component_id)) "eta" else paste0("eta_component_", component_id)
    std_eta_name <- if(is.null(component_id)) "std_eta" else paste0("std_eta_component_", component_id)

    for(i in seq_len(J)){
      syntax <- paste0(syntax, eta_name, "[", i, "] ~ dgamma(", prior$weights$alpha[i], ", 1)\n")
    }
    syntax <- paste0(syntax,
      "for(j in 1:", J, "){\n",
      "  ", std_eta_name, "[j] <- ", eta_name, "[j] / sum(", eta_name, ")\n",
      "}\n",
      omega_local, "[1] <- 1\n"
    )
    if(J > 1L){
      syntax <- paste0(syntax,
        "for(j in 2:", J, "){\n",
        "  ", omega_local, "[j] <- sum(", std_eta_name, "[j:", J, "])\n",
        "}\n"
      )
    }

  }else if(prior$weights$type == "fixed"){
    for(i in seq_len(J)){
      syntax <- paste0(syntax, omega_local, "[", i, "] <- ", prior$weights$omega[i], "\n")
    }

  }else if(prior$weights$type == "independent"){
    syntax <- paste0(syntax, omega_local, "[1] <- 1\n")
    if(J > 1L){
      for(i in 2:J){
        if(prior$weights$scale == "omega"){
          syntax <- paste0(syntax, .JAGS_prior.simple(prior$weights$prior, paste0(omega_local, "[", i, "]")))
        }else if(prior$weights$scale == "log_omega"){
          log_omega_name <- if(is.null(component_id)) "log_omega" else paste0("log_omega_component_", component_id)
          syntax <- paste0(
            syntax,
            .JAGS_prior.simple(prior$weights$prior, paste0(log_omega_name, "[", i, "]")),
            omega_local, "[", i, "] <- exp(", log_omega_name, "[", i, "])\n"
          )
        }
      }
    }
  }

  if(!is.null(component_id) || needs_mapping){
    global_bin_indices <- .weightfunction_global_bin_indices(all_cuts, expansion)
    for(i in seq_len(length(all_cuts) - 1L)){
      ind <- global_bin_indices[i]
      syntax <- paste0(syntax, omega_target, "[", i, "] <- ", omega_local, "[", expansion$index[ind], "]\n")
    }
  }

  syntax
}
.JAGS_weightfunction_none_component_syntax <- function(component_id, n_bins){

  syntax <- character()
  omega_target <- if(is.null(component_id)) "omega" else paste0("omega_component_", component_id)
  for(i in seq_len(n_bins)){
    syntax <- paste0(syntax, omega_target, "[", i, "] <- 1\n")
  }
  syntax
}
.JAGS_prior.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable_list  <- list(.get_spike_and_slab_variable(prior))
  prior_inclusion_list <- list(.get_spike_and_slab_inclusion(prior))
  names(prior_variable_list)  <- paste0(parameter_name, "_variable")
  names(prior_inclusion_list) <- paste0(parameter_name, "_inclusion")

  syntax <- paste0(
    .JAGS_add_priors.fun(prior_variable_list),
    .JAGS_add_priors.fun(prior_inclusion_list),
    parameter_name, "_indicator ~ dbern(",   paste0(parameter_name, "_inclusion"), ")\n",
    parameter_name, " = ",  parameter_name, "_variable * ", parameter_name, "_indicator\n"
  )

  return(syntax)
}
.JAGS_prior.mixture        <- function(prior_list, parameter_name){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_phacking       <- sapply(prior_list, is_prior_phacking)
    is_bias           <- sapply(prior_list, is_prior_bias)
    is_none           <- sapply(prior_list, is.prior.none)
    branch_info       <- lapply(prior_list, .selection_branch_info)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    prior_weights <- attr(prior_list, "prior_weights")
    if(any(has_selection) || any(has_phacking)){
      spec <- selection_backend_spec(prior_list)
      syntax <- .JAGS_selection_backend_syntax(spec)
    }else{
      syntax <- paste0(" bias_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))\n")
    }

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      named_prior_PET <- prior_list[[which(is_PET)]]
      class(named_prior_PET) <- class(named_prior_PET)[!class(named_prior_PET) %in% "prior.PET"]
      named_prior_PET <- list("PET_1" = named_prior_PET)

      syntax <- paste0(
        syntax,
        .JAGS_add_priors.fun(named_prior_PET),
        " PET <- PET_1 * equals(bias_indicator, ", which(is_PET), ")\n"
      )
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      named_prior_PEESE <- prior_list[[which(is_PEESE)]]
      class(named_prior_PEESE) <- class(named_prior_PEESE)[!class(named_prior_PEESE) %in% "prior.PEESE"]
      named_prior_PEESE <- list("PEESE_1" = named_prior_PEESE)

      syntax <- paste0(
        syntax,
        .JAGS_add_priors.fun(named_prior_PEESE),
        " PEESE <- PEESE_1 * equals(bias_indicator, ", which(is_PEESE), ")\n"
      )
    }

  }else{

    prior_weights    <- attr(prior_list, "prior_weights")
    prior_components <- as.list(prior_list)
    class(prior_components) <- "list"
    names(prior_components) <- paste0(parameter_name, "_component_", seq_along(prior_components))

    syntax <- paste0(
      " ", parameter_name, "_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))\n",
      sapply(.JAGS_add_priors.fun(prior_components), paste, collapse = "\n"),
      " ", parameter_name, " = ",  paste0(names(prior_components), " * ", paste0("(", parameter_name, "_indicator == ", seq_along(prior_components), ")"), collapse = " + "), "\n"
    )
  }

  return(syntax)
}


.add_JAGS_vector   <- function(name, vector){

  if(!is.vector(vector))
    stop("vector must be a vector")
  check_char(name, "name")

  syntax <- paste0(" ", name, " = c(", paste0(vector, collapse = ", "), ")\n")

  return(syntax)
}
.add_JAGS_matrix   <- function(name, matrix){

  if(!is.matrix(matrix))
    stop("matrix must be a matrix")
  check_char(name, "name")

  syntax <- ""

  # this unfortunatelly cannot be defined on row/column basis
  # I tried simplifying this before but only possible initialization is elementwise
  for(i in 1:nrow(matrix)){
   syntax <- paste0(
     syntax, " ",
     paste0(name,"[", i, ",", seq_len(ncol(matrix)), "] = ", matrix[i,], collapse = "; "), "\n"
   )
  }

  return(syntax)
}
.check_JAGS_syntax <- function(syntax){

  check_char(syntax, "syntax", allow_NULL = TRUE)
  if(is.null(syntax)){
    syntax <- "model{}"
  }
  if(!grepl("model", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("{", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("}", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
}
.JAGS_parameter_to_precision <- function(parameter){

  if(is.character(parameter)){
    parameter <- paste0("1/pow(", parameter, ", 2)")
  }else{
    parameter <- 1/parameter^2
  }

  return(parameter)
}

#' @title Create initial values for 'JAGS' model
#'
#' @description Creates initial values for priors in
#' a 'JAGS' model.
#'
#' @param chains number of chains
#' @param seed seed for random number generation
#'
#' @inheritParams JAGS_add_priors
#'
#' @return \code{JAGS_add_priors} returns a list of JAGS
#' initial values.
#'
#' @export
JAGS_get_inits            <- function(prior_list, chains, seed){

  # return empty list in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  check_int(chains, "chains", lower = 1)
  check_real(seed, "seed", allow_NULL = TRUE)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # select seed at random if none was specified
  if(is.null(seed)){
    seed <- sample(666666, 1)
  }
  set.seed(seed)


  # create the starting values
  inits <- vector("list", chains)
  for(j in 1:chains){

    temp_inits <- .JAGS_get_inits.fun(prior_list)

    temp_inits[[".RNG.seed"]] <- seed + j
    temp_inits[[".RNG.name"]] <- if(chains > 4) "lecuyer::RngStream" else "base::Super-Duper"

    inits[[j]] <- temp_inits
  }

  return(inits)
}

.JAGS_get_inits.fun        <- function(prior_list){

  temp_inits <- list()

  for(i in seq_along(prior_list)){

    if(is.prior.point(prior_list[[i]])){

      next

    }else if(is.prior.weightfunction(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.factor(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(temp_inits)
}
.JAGS_init.simple          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(prior[["distribution"]] == "point"){

    return()

  }else{
    init <- list()

    if(prior[["distribution"]] == "invgamma"){

      sampling_prior <- prior(
        "distribution" = "gamma",
        "parameters"   = list("shape" = prior$parameters[["shape"]], "rate" = prior$parameters[["scale"]]),
        "truncation"   = list("lower" = prior$truncation[["upper"]]^-1, "upper" = prior$truncation[["lower"]]^-1))
      init[[paste0("inv_", parameter_name)]] <- rng(sampling_prior, 1)

    }else{

      init[[parameter_name]] <- rng(prior, 1)

    }
  }

  return(init)
}
.JAGS_init.vector          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(prior[["distribution"]] == "point"){

    return()

  }else{

    init <- list()


    if(prior[["distribution"]] == "mt"){
      init[[paste0("prior_par_s_", parameter_name)]] <- rng(prior("gamma", list(shape = prior$parameters[["df"]]/2, rate = prior$parameters[["df"]]/2)), 1)
      init[[paste0("prior_par_z_", parameter_name)]] <- rng(prior("mnormal", list(mean = 0, sd = prior$parameters[["scale"]], K = prior$parameters[["K"]])), 1)[1,]
    }else{
      init[[parameter_name]] <- rng(prior, 1)[1,]
    }

  }

  return(init)
}
.JAGS_init.factor          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(.get_prior_factor_levels(prior), "levels", lower = 1)

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    init <- list()
    init[[parameter_name]] <- rng(prior, .get_prior_factor_levels(prior))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    # remove the orthonormal/meandif class, otherwise samples from the transformed distributions are generated
    class(prior) <- class(prior)[!class(prior) %in% c("prior.orthonormal", "prior.meandif")]

    init <- .JAGS_init.vector(prior, parameter_name)

  }

  return(init)
}
.JAGS_init.PP              <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    init <- .JAGS_init.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    init <- .JAGS_init.simple(prior, "PEESE")
  }

  return(init)
}
.JAGS_init.weightfunction  <- function(prior, component_id = NULL){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  if(is.null(component_id)){
    return(selection_backend_spec(prior)$init)
  }

  init <- list()
  if(prior$weights$type == "fixed"){
    return()
  }else if(prior$weights$type == "cumulative"){
    eta_name <- paste0("eta_component_", component_id)
    init[[eta_name]] <- stats::rgamma(length(prior$weights[["alpha"]]), shape = prior$weights[["alpha"]], rate = 1)
  }

  return(init)
}
.JAGS_init.phacking       <- function(prior, component_id = NULL){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  if(is.null(component_id)){
    return(selection_backend_spec(prior)$init)
  }

  alpha_name <- paste0("alpha_component_", component_id)
  init <- .JAGS_init.simple(prior$alpha, alpha_name)

  return(init)
}
.JAGS_init.bias           <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  return(selection_backend_spec(prior)$init)
}
.JAGS_init.spike_and_slab  <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")

  prior_variable        <- list(.get_spike_and_slab_variable(prior))
  names(prior_variable) <- paste0(parameter_name, "_variable")
  init <- .JAGS_get_inits.fun(prior_variable)

  if(!is.prior.point(.get_spike_and_slab_inclusion(prior))){
    init[[paste0(parameter_name, "_inclusion")]] <- rng(.get_spike_and_slab_inclusion(prior), 1)
  }


  return(init)
}
.JAGS_init.mixture         <- function(prior_list, parameter_name){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_phacking       <- sapply(prior_list, is_prior_phacking)
    is_bias           <- sapply(prior_list, is_prior_bias)
    is_none           <- sapply(prior_list, is.prior.none)
    branch_info       <- lapply(prior_list, .selection_branch_info)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

    init <- list()

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      named_prior_PET <- prior_list[[which(is_PET)]]
      class(named_prior_PET) <- class(named_prior_PET)[!class(named_prior_PET) %in% "prior.PET"]
      named_prior_PET <- list("PET_1" = named_prior_PET)

      init <- c(init, .JAGS_get_inits.fun(named_prior_PET))
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      named_prior_PEESE <- prior_list[[which(is_PEESE)]]
      class(named_prior_PEESE) <- class(named_prior_PEESE)[!class(named_prior_PEESE) %in% "prior.PEESE"]
      named_prior_PEESE <- list("PEESE_1" = named_prior_PEESE)

      init <- c(init, .JAGS_get_inits.fun(named_prior_PEESE))
    }
    if(any(has_selection) || any(has_phacking)){
      init <- c(init, selection_backend_spec(prior_list)$init)
    }else{
      init[["bias_indicator"]] <- rng(prior_list, 1, sample_components = TRUE)
    }

  }else{

    prior_components <- as.list(prior_list)
    class(prior_components) <- "list"
    names(prior_components) <- paste0(parameter_name, "_component_", seq_along(prior_components))

    init <- .JAGS_get_inits.fun(prior_components)
    init[[paste0(parameter_name, "_indicator")]] <- rng(prior_list, 1, sample_components = TRUE)

  }

  return(init)
}


#' @title Create list of monitored parameters for 'JAGS' model
#'
#' @description Creates a vector of parameter names to be
#' monitored in a 'JAGS' model.
#'
#' @inheritParams JAGS_add_priors
#'
#' @return \code{JAGS_to_monitor} returns a character vector of
#' parameter names.
#'
#' @export
JAGS_to_monitor             <- function(prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    return("")
  }

  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # add the monitored parameters
  monitor <- character()
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.factor(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(monitor)
}


.JAGS_monitor.simple         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!(is.prior.simple(prior) | is.prior.vector(prior) | is.prior.factor(prior)))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(prior[["distribution"]] == "invgamma"){
    monitor <- c(parameter_name, paste0("inv_", parameter_name))
  }else{
    monitor <- parameter_name
  }

  return(monitor)
}
.JAGS_monitor.vector         <- function(prior, parameter_name){

  monitor <- .JAGS_monitor.simple(prior, parameter_name)

  return(monitor)
}
.JAGS_monitor.factor         <- function(prior, parameter_name){

  monitor <- .JAGS_monitor.simple(prior, parameter_name)

  return(monitor)
}
.JAGS_monitor.PP             <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    monitor <- .JAGS_monitor.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    monitor <- .JAGS_monitor.simple(prior, "PEESE")
  }

  return(monitor)
}
.JAGS_monitor.weightfunction <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  return(selection_backend_spec(prior)$monitor)
}
.JAGS_monitor_private.weightfunction <- function(prior){

  if(prior$weights$type == "cumulative"){
    return("eta")
  }
  if(prior$weights$type == "independent" && prior$weights$scale == "log_omega"){
    return("log_omega")
  }

  character()
}
.JAGS_monitor.phacking      <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)$monitor
}
.JAGS_monitor.bias          <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)$monitor
}
.JAGS_monitor.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable  <- list(.get_spike_and_slab_variable(prior))
  prior_inclusion <- list(.get_spike_and_slab_inclusion(prior))
  names(prior_variable)  <- paste0(parameter_name, "_variable")
  names(prior_inclusion) <- paste0(parameter_name, "_inclusion")

  monitor <- c(
    paste0(parameter_name, "_indicator"),
    JAGS_to_monitor(prior_inclusion),
    parameter_name,
    JAGS_to_monitor(prior_variable)
  )

  return(monitor)
}
.JAGS_monitor.mixture        <- function(prior_list, parameter_name){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_phacking       <- sapply(prior_list, is_prior_phacking)
    is_bias           <- sapply(prior_list, is_prior_bias)
    is_none           <- sapply(prior_list, is.prior.none)
    branch_info       <- lapply(prior_list, .selection_branch_info)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    monitor <- if(any(has_selection) || any(has_phacking)){
      selection_backend_spec(prior_list)$monitor
    }else{
      "bias_indicator"
    }

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PET")
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PEESE")
    }
  }else{

    monitor <- c(paste0(parameter_name, "_indicator"), parameter_name)
  }

  return(unique(monitor))
}

#' @title Check and list 'JAGS' fitting settings
#'
#' @description Checks and lists settings for the
#' [JAGS_fit] function.
#'
#' @param check_mins named list of minimal values for which
#' should some input be checked. Defaults to:
#' \describe{
#'   \item{chains}{\code{1}}
#'   \item{adapt}{\code{50}}
#'   \item{burnin}{\code{50}}
#'   \item{sample}{\code{100}}
#'   \item{thin}{\code{1}}
#' }
#' @param skip_sample_extend whether \code{sample_extend}
#' is allowed to be NULL and skipped in the check
#'
#' @inheritParams JAGS_fit
#' @inheritParams check_input
#'
#' @return \code{JAGS_check_and_list_fit_settings} invisibly returns a
#' list of checked fit settings. \code{JAGS_check_and_list_autofit_settings}
#' invisibly returns a list of checked autofit settings.
#' parameter names.
#'
#' @export JAGS_check_and_list_fit_settings
#' @export JAGS_check_and_list_autofit_settings
#' @name JAGS_check_and_list
NULL

#' @rdname JAGS_check_and_list
JAGS_check_and_list_fit_settings     <- function(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed, check_mins = list(chains = 1, adapt = 50, burnin = 50, sample = 100, thin = 1), call = ""){

  check_int(chains, "chains", lower = check_mins[["chains"]], call = call)
  check_int(adapt,  "adapt",  lower = check_mins[["adapt"]],  call = call)
  check_int(burnin, "burnin", lower = check_mins[["burnin"]], call = call)
  check_int(sample, "sample", lower = check_mins[["sample"]], call = call)
  check_int(thin,   "thin",   lower = check_mins[["thin"]],   call = call)
  check_bool(parallel, "parallel",                call = call)
  check_int(cores,     "cores", lower = 1,        call = call)
  check_bool(autofit,  "autofit",                 call = call)
  check_bool(silent,   "silent",                  call = call)
  check_int(seed,      "seed", allow_NULL = TRUE, call = call)

  return(invisible(list(
    chains   = chains,
    adapt    = adapt,
    burnin   = burnin,
    sample   = sample,
    thin     = thin,
    autofit  = autofit,
    parallel = parallel,
    cores    = cores,
    silent   = silent,
    seed     = seed
  )))
}

#' @rdname JAGS_check_and_list
JAGS_check_and_list_autofit_settings <- function(autofit_control, skip_sample_extend = FALSE, call = ""){

  check_list(autofit_control, "autofit_control", check_names = c("max_Rhat", "min_ESS", "max_error", "max_SD_error",  "max_time", "sample_extend", "restarts", "max_extend", "check_indicators"), call = call)
  if(is.null(autofit_control[["check_indicators"]])){
    autofit_control[["check_indicators"]] <- FALSE
  }
  check_real(autofit_control[["max_Rhat"]],     "max_Rhat",     lower = 1, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["min_ESS"]],      "min_ESS",      lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_error"]],    "max_error",    lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_SD_error"]], "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE, call = call)
  check_bool(autofit_control[["check_indicators"]], "check_indicators", call = call)
  check_list(autofit_control[["max_time"]],     "max_time", check_names = c("time", "unit"), check_length = 2, allow_NULL = TRUE, call = call)
  if(!is.null(autofit_control[["max_time"]])){
    if(is.null(names(autofit_control[["max_time"]]))){
      names(autofit_control[["max_time"]]) <- c("time", "unit")
    }
    check_real(autofit_control[["max_time"]][["time"]], "max_time:time", lower = 0, call = call)
    check_char(autofit_control[["max_time"]][["unit"]], "max_time:unit", allow_values = c("secs", "mins", "hours", "days", "weeks"), call = call)
  }
  check_int(autofit_control[["sample_extend"]], "sample_extend", lower = 1, allow_NULL = skip_sample_extend, call = call)
  check_int(autofit_control[["restarts"]], "restarts", lower = 1, allow_NULL = TRUE, call = call)
  check_int(autofit_control[["max_extend"]], "max_extend", lower = 1, allow_NULL = TRUE, call = call)

  return(invisible(autofit_control))
}
