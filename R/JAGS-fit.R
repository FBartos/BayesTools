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
#'   \item{max_error}{maximum MCMC error. Defaults to \code{1.01}.}
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

#' @rdname JAGS_fit
JAGS_fit <- function(model_syntax, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL,
                     chains = 4, adapt = 500, burnin = 1000, sample = 4000, thin = 1,
                     autofit = FALSE, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10),
                     parallel = FALSE, cores = chains, silent = TRUE, seed = NULL,
                     add_parameters = NULL, required_packages = NULL, ...){

  .check_runjags()
  dots <- list(...)

  ### check input
  .check_JAGS_syntax(model_syntax)
  JAGS_check_and_list_fit_settings(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed)
  JAGS_check_and_list_autofit_settings(autofit_control)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE)
  check_char(required_packages, "required_packages", check_length = 0, allow_NULL = TRUE)
  check_list(formula_list, "formula_list", allow_NULL = TRUE)
  check_list(formula_data_list, "formula_data_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))
  check_list(formula_prior_list, "formula_prior_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))

  ### add formulas
  if(!is.null(formula_list)){

    # obtain settings for each formula
    formula_output <- list()
    for(parameter in names(formula_list)){
      formula_output[[parameter]] <- JAGS_formula(
        formula    = formula_list[[parameter]],
        parameter  = parameter,
        data       = formula_data_list[[parameter]],
        prior_list = formula_prior_list[[parameter]])
    }

    # merge with the rest of the input
    prior_list     <- c(do.call(c, unname(lapply(formula_output, function(output) output[["prior_list"]]))), prior_list)
    data           <- c(do.call(c, unname(lapply(formula_output, function(output) output[["data"]]))),       data)
    formula_syntax <- paste0(lapply(formula_output, function(output) output[["formula_syntax"]]), collapse = "")

    # add the formula syntax to the model syntax
    opening_bracket <- regexpr("{", model_syntax, fixed = TRUE)[1]
    syntax_start    <- substr(model_syntax, 1, opening_bracket)
    syntax_end      <- substr(model_syntax, opening_bracket + 1, nchar(model_syntax))
    model_syntax    <- paste0(syntax_start, "\n", formula_syntax, "\n", syntax_end)
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
    on.exit(try(parallel::stopCluster(cl)))
    for(i in seq_along(required_packages)){
      parallel::clusterCall(cl, function(x) requireNamespace(required_packages[i]))
    }
    model_call <- c(
      model_call,
      method = "rjparallel",
      cl     = list(cl)
    )
  }else{
    for(i in seq_along(required_packages)){
      requireNamespace(required_packages[i])
    }
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
    on.exit(runjags::runjags.options(silent.jags = runjags::runjags.getOption("silent.jags"), silent.runjags = runjags::runjags.getOption("silent.runjags")))
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

    converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], fail_fast = TRUE)
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

      converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], fail_fast = TRUE)
      itteration <- itteration + 1

      if(isTRUE(dots[["is_JASP"]]))
        .JASP_progress_bar_tick()
    }
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "required_packages") <- required_packages

  class(fit) <- c(class(fit), "BayesTools_fit")

  return(fit)
}

#' @rdname JAGS_fit
JAGS_extend <- function(fit, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10),
                        parallel = FALSE, cores = NULL, silent = TRUE, seed = NULL){

  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a 'BayesTools_fit'")

  # extract fitting information
  prior_list        <- attr(fit, "prior_list")
  model_syntax      <- attr(fit, "model_syntax")
  required_packages <- attr(fit, "required_packages")
  JAGS_check_and_list_autofit_settings(autofit_control)

  # parallel vs. not
  if(parallel){
    if(is.null(cores)){
      cores <- length(fit[["mcmc"]])
    }
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(try(parallel::stopCluster(cl)))
    for(i in seq_along(required_packages)){
      parallel::clusterCall(cl, function(x) requireNamespace(required_packages[i]))
    }
    refit_call <- list(
      runjags.object = fit,
      sample         = autofit_control[["sample_extend"]],
      method         = "rjparallel",
      cl             = cl,
      summarise      = FALSE
    )
  }else{
    for(i in seq_along(required_packages)){
      requireNamespace(required_packages[i])
    }
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
    on.exit(runjags::runjags.options(silent.jags = runjags::runjags.getOption("silent.jags"), silent.runjags = runjags::runjags.getOption("silent.runjags")))
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  }

  start_time <- Sys.time()
  itteration <- 0
  converged  <- FALSE

  while(!converged & itteration < autofit_control[["restarts"]]){

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

    fit <- tryCatch(do.call(runjags::extend.jags, refit_call), error = function(e)e)

    if(inherits(fit, "error")){
      if(!silent)
        warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)

      break
    }

    converged <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], fail_fast = TRUE)

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
#' @param max_error maximum MCMC error. Defaults to \code{1.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of additional parameter names that should be used
#' (only allows removing last, fixed, omega element if omega is tracked manually).
#' @param fail_fast whether the function should stop after the first failed convergence check.
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
JAGS_check_convergence <- function(fit, prior_list, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, add_parameters = NULL, fail_fast = FALSE){

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

  # extract samples and parameter information
  mcmc_samples    <- coda::as.mcmc.list(fit)
  parameter_names <- colnames(mcmc_samples[[1]])
  parameters_keep <- rep(TRUE, length(parameter_names))

  # remove auxiliary and support parameters from the summary
  for(i in seq_along(prior_list)){
    if(is.prior.weightfunction(prior_list[[i]])){
      if(prior_list[[i]][["distribution"]] %in% c("one.sided", "two.sided")){
        parameters_keep[grepl("eta", parameter_names)] <- FALSE
      }
      parameter_names[max(grep("omega", parameter_names))] <- FALSE
    }else if(is.prior.mixture(prior_list[[i]]) && any(sapply(prior_list[[i]], is.prior.weightfunction))){
      parameters_keep[max(grep("omega", parameter_names))] <- FALSE
    }else if(is.prior.point(prior_list[[i]])){
      parameters_keep[parameter_names == names(prior_list)[i]] <- FALSE
    }else if(is.prior.simple(prior_list[[i]]) && prior_list[[i]][["distribution"]] == "invgamma"){
      parameters_keep[parameter_names == paste0("inv_",names(prior_list)[i])] <- FALSE
    }else if(is.prior.mixture(prior_list[[i]]) && length(prior_list[[i]]) == 1 && is.prior.point(prior_list[[i]][[1]])){
      parameters_keep[parameter_names == names(prior_list)[i]] <- FALSE
    }
  }

  # remove indicators/inclusions
  parameters_keep[grepl("_indicator", parameter_names)] <- FALSE
  parameters_keep[grepl("_inclusion", parameter_names)] <- FALSE

  if(all(!parameters_keep)){
    return(TRUE)
  }

  # remove parameters that are not monitored
  for(i in seq_along(mcmc_samples)){
    mcmc_samples[[i]] <- mcmc_samples[[i]][,parameters_keep,drop=FALSE]
  }

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

  # create an empty attribute holder if any data need to be passed with the syntax
  syntax_attributes <- NULL

  # identify parts of the syntax
  opening_bracket <- regexpr("{", syntax, fixed = TRUE)[1]
  syntax_start    <- substr(syntax, 1, opening_bracket)
  syntax_end      <- substr(syntax, opening_bracket + 1, nchar(syntax))

  # create the priors relevant syntax
  syntax_priors <- .JAGS_add_priors.fun(prior_list)
  if(!is.null(attr(syntax_priors, "auxiliary_data"))){
    syntax_attributes <- attr(syntax_priors, "auxiliary_data")
  }

  # merge everything back together
  syntax <- paste0(syntax_start, "\n", syntax_priors, "\n", syntax_end)
  if(!is.null(attr(syntax_priors, "auxiliary_data"))){
    attr(syntax, "auxiliary_data") <- syntax_attributes
  }

  return(syntax)
}

.JAGS_add_priors.fun       <- function(prior_list){

  syntax_priors <- ""
  syntax_attributes <- NULL

  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.weightfunction(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      syntax_mixture <- .JAGS_prior.mixture(prior_list[[i]], names(prior_list)[i])
      syntax_priors  <- paste(syntax_priors, syntax_mixture)
      if(!is.null(attr(syntax_mixture, "auxiliary_data"))){
        syntax_attributes <- attr(syntax_mixture, "auxiliary_data")
      }

    }else if(is.prior.factor(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  if(length(syntax_attributes) > 0){
    attr(syntax_priors, "auxiliary_data") <- syntax_attributes
  }

  return(syntax_priors)
}
.JAGS_prior.simple         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # distribution
  syntax <- switch(
    prior[["distribution"]],
    "point"     = paste0(parameter_name," = ",prior$parameter[["location"]]),
    "normal"    = paste0(parameter_name," ~ dnorm(",prior$parameter[["mean"]],",",1/prior$parameter[["sd"]]^2,")"),
    "lognormal" = paste0(parameter_name," ~ dlnorm(",prior$parameter[["meanlog"]],",",1/prior$parameter[["sdlog"]]^2,")"),
    "t"         = paste0(parameter_name," ~ dt(",prior$parameter[["location"]],",",1/prior$parameter[["scale"]]^2,",", prior$parameter[["df"]],")"),
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

  .check_prior(prior)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")


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

    # TODO: beautify this code by specific JAGS distributions?
    if(prior[["distribution"]] == "mt"){
      # using the chisq * covariance parametrization since the mt fails with 1 df
      # (using a common df parameter as in Rouder et al. 2012)
      syntax <- paste0("prior_par1_", parameter_name, " = rep(0,", prior$parameter[["K"]], ")\n")
      syntax <- paste0(syntax, "prior_par_s_", parameter_name, " ~ dgamma(", prior$parameter[["df"]]/2, ", ", prior$parameter[["df"]]/2,")\n")
      syntax <- paste0(
        syntax,
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  prior_par2_", parameter_name, "[i,i] <- ", 1/par2^2, "\n",
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
        "  prior_par2_", parameter_name, "[i,i] <- ", 1/par2^2, "\n",
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

  .check_prior(prior)
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

  .check_prior(prior)
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

  # creating cummulative dirichlet distribution using gammas (in order to bypass bugs in bridgesampling)
  if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){
    syntax <- character()
    for(i in 1:length(prior$parameters[["alpha"]])){
      syntax <- paste0(syntax, "eta[",i,"] ~ dgamma(",prior$parameters[["alpha"]][i],", 1)\n")
    }
    syntax <- paste0(syntax,
                     "for(j in 1:",length(prior$parameters[["alpha"]]),"){\n",
                     "  std_eta[j]  = eta[j] / sum(eta)\n",
                     "  omega[j]    = sum(std_eta[1:j])\n",
                     "}\n")
  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){
    syntax <- character()
    for(i in 1:length(prior$parameters[["alpha1"]])){
      syntax <- paste0(syntax, "eta1[",i,"] ~ dgamma(",prior$parameters[["alpha1"]][i],", 1)\n")
    }
    for(i in 1:length(prior$parameters[["alpha2"]])){
      syntax <- paste0(syntax, "eta2[",i,"] ~ dgamma(",prior$parameters[["alpha2"]][i],", 1)\n")
    }
    syntax <- paste0(syntax,
                     "for(j1 in 1:",length(prior$parameters[["alpha1"]]),"){\n",
                     "  std_eta1[j1]      = eta1[j1] / sum(eta1)\n",
                     "  omega[",length(prior$parameters[["alpha2"]])," - 1 + j1] = sum(std_eta1[1:j1])\n",
                     "}\n",
                     "for(j2 in 1:",length(prior$parameters[["alpha2"]]),"){\n",
                     "  std_eta2[j2]  = (eta2[j2] / sum(eta2)) * (1 - std_eta1[1])\n",
                     "}\n",
                     "for(j2 in 2:",length(prior$parameters[["alpha2"]]),"){\n",
                      "  omega[j2-1] = sum(std_eta2[j2:",length(prior$parameters[["alpha2"]]),"]) + std_eta1[1]\n",
                     "}\n")
  }else if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){
    syntax <- character()
    for(i in 1:length(prior$parameters[["omega"]])){
      syntax <- paste0(syntax, "omega[",i,"] = ",prior$parameters[["omega"]][i],"\n")
    }
  }

  return(syntax)
}
.JAGS_prior.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable_list  <- prior["variable"]
  prior_inclusion_list <- prior["inclusion"]
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

  .check_prior_list(prior_list)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_none           <- sapply(prior_list, is.prior.none)

    prior_weights <- attr(prior_list, "prior_weights")
    syntax <- paste0(" bias_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))\n")

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      named_prior_PET <- prior_list[[which(is_PET)]]
      class(named_prior_PET) <- class(named_prior_PET)[!class(named_prior_PET) %in% "prior.PET"]
      named_prior_PET <- list("PET_1" = named_prior_PET)

      syntax <- paste0(
        syntax,
        .JAGS_add_priors.fun(named_prior_PET),
        " PET = PET_1 * (bias_indicator == ", which(is_PET), ")\n"
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
        " PEESE = PEESE_1 * (bias_indicator == ", which(is_PEESE), ")\n"
      )
    }
    if(any(is_weightfunction)){
      # we cannot simulate weights from the mixture distribution directly because
      # JAGS does not allow complex support for the cumulative simplex parameter
      # (we could make it on the non-cumulative simplex but it does not give more advantage)

      # create a vector of the alpha parameters
      alpha <- lapply(prior_list[is_weightfunction], function(x){
        if(grepl("fixed", x[["distribution"]])){
          return(x$parameters[["omega"]])
        }else{
          return(x$parameters[["alpha"]])
        }
      })

      # dispatch the prior distribution on weight parameters
      syntax <- paste0(
        syntax,
        " for(i in 1:", max(lengths(alpha)), "){\n",
        "   eta[i] ~ dgamma(eta_shape[i, bias_indicator], 1)\n",
        " }\n"
      )

      # transform etas into weights (eta2omega JAGS function is in the RoBMA package)
      syntax <- paste0(syntax, " omega = eta2omega(eta, omega_index[,bias_indicator], eta_index[,bias_indicator], eta_index_max[bias_indicator])\n")

      # add the necessary auxiliary data: omega_index, eta_index, eta_index_max, eta_shape
      # create the weightfunction mapping for weights
      omega_index_weighfunction <- weightfunctions_mapping(prior_list[is_weightfunction], one_sided = TRUE)
      omega_index_weighfunction <- lapply(omega_index_weighfunction, rev)
      omega_index_weighfunction <- do.call(rbind, omega_index_weighfunction)
      omega_index <- matrix(0, ncol = ncol(omega_index_weighfunction), length(prior_list))
      omega_index[is_weightfunction,] <- omega_index_weighfunction

      # in case of fixed weight functions, the omega_index direly corresponds to the fixed weights
      for(i in seq_along(prior_list)){
        if(is.prior.weightfunction(prior_list[[i]])){
          if(grepl("fixed", prior_list[[i]]$distribution)){
            omega_index[i,] <- prior_list[[i]]$parameters[["omega"]][omega_index[i,]]
          }
        }
      }

      # create the eta to omega mapping
      # eta_index_max helps dispatching within the eta2omega function
      # 0  = non-weightfunction, all weights are set to 0
      # >1 = indicates how many alpha parameters are needed to construct the weightfunction based on the eta_index
      # -1 = indicates fixed weightfunction, omega index already encoded all weights
      eta_index     <- matrix(0, nrow = length(prior_list), ncol = max(lengths(alpha)))
      eta_index_max <- rep(0, length(prior_list))
      for(i in seq_along(prior_list)){
        if(is.prior.weightfunction(prior_list[[i]])){
          if(grepl("fixed", prior_list[[i]]$distribution)){
            eta_index_max[i] <- -1
            eta_index[i,]    <- -1
          }else{
            temp_index       <- unique(omega_index[i,])
            eta_index[i,1:length(temp_index)] <- sort(temp_index)
            eta_index_max[i] <- length(temp_index)
          }
        }else{
          eta_index_max[i] <- 0
        }
      }

      # create priors for eta (set alpha to 1 for non-weightfunctions to keep the sampling in the expected area)
      eta_shape <- matrix(1, nrow = length(prior_list), ncol = max(lengths(alpha)))
      for(i in seq_along(prior_list)){
        if(is.prior.weightfunction(prior_list[[i]])){
          if(!grepl("fixed", prior_list[[i]]$distribution)){
            temp_shape <- prior_list[[i]]$parameters[["alpha"]]
            eta_shape[i,1:length(temp_shape)] <- temp_shape
          }
        }
      }

      # paste the matricies directly into JAGS code (simplifies data handling)
      syntax <- paste0(
        syntax,
        .add_JAGS_matrix("omega_index",   t(omega_index)),
        .add_JAGS_matrix("eta_index",     t(eta_index)),
        .add_JAGS_vector("eta_index_max", eta_index_max),
        .add_JAGS_matrix("eta_shape",     t(eta_shape))
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

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

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

  .check_prior(prior)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

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

  .check_prior(prior)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(.get_prior_factor_levels(prior), "levels", lower = 1)

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
.JAGS_init.weightfunction  <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  init <- list()
  if(prior[["distribution"]] %in% c("one.sided.fixed", "two.sided.fixed")){

    return()

  }else if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){

    init[["eta"]] <- stats::rgamma(length(prior$parameters[["alpha"]]), shape = prior$parameters[["alpha"]], rate = 1)

  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){

    init[["eta1"]] <- stats::rgamma(length(prior$parameters[["alpha1"]]), shape = prior$parameters[["alpha1"]], rate = 1)
    init[["eta2"]] <- stats::rgamma(length(prior$parameters[["alpha2"]]), shape = prior$parameters[["alpha2"]], rate = 1)

  }

  return(init)
}
.JAGS_init.spike_and_slab  <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")

  prior_variable        <- prior["variable"]
  names(prior_variable) <- paste0(parameter_name, "_variable")
  init <- .JAGS_get_inits.fun(prior_variable)

  if(!is.prior.point(prior[["inclusion"]])){
    init[[paste0(parameter_name, "_inclusion")]] <- rng(prior[["inclusion"]], 1)
  }


  return(init)
}
.JAGS_init.mixture         <- function(prior_list, parameter_name){

  .check_prior_list(prior_list)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_none           <- sapply(prior_list, is.prior.none)

    init <- list()

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_none)))
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
    if(any(is_weightfunction)){

      # find prior with the most alpha parameters and simulate initial values from it
      alpha <- sapply(prior_list[is_weightfunction], function(x){
        if(grepl("fixed", x[["distribution"]])){
          return(length(x$parameters[["omega"]]))
        }else{
          return(length(x$parameters[["alpha"]]))
        }
      })

      init <- c(init, .JAGS_get_inits.fun(prior_list[is_weightfunction][which.max(alpha)]))
    }

    init[["bias_indicator"]] <- rng(prior_list, 1, sample_components = TRUE)

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

  .check_prior(prior)
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

  monitor <- "omega"
  if(all(names(prior[["parameters"]]) %in% c("alpha", "steps"))){
    monitor <- c(monitor, "eta")
  }else if(all(names(prior[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){
    monitor <- c(monitor, "eta1", "eta2")
  }

  return(monitor)
}
.JAGS_monitor.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable  <- prior["variable"]
  prior_inclusion <- prior["inclusion"]
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

  .check_prior_list(prior_list)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_none           <- sapply(prior_list, is.prior.none)

    monitor <- "bias_indicator"

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PET")
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PEESE")
    }
    if(any(is_weightfunction)){

      monitor <- c(monitor, "omega")
    }

  }else{

    monitor <- c(paste0(parameter_name, "_indicator"), parameter_name)
  }

  return(monitor)
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

  check_list(autofit_control, "autofit_control", check_names = c("max_Rhat", "min_ESS", "max_error", "max_SD_error",  "max_time", "sample_extend", "restarts", "max_extend"), call = call)
  check_real(autofit_control[["max_Rhat"]],     "max_Rhat",     lower = 1, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["min_ESS"]],      "min_ESS",      lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_error"]],    "max_error",    lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_SD_error"]], "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE, call = call)
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
