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
#' @param formula_random_prior_list optional named list of `prior_random()`
#' objects for random effects in `formula_list`. Required for any formula that
#' contains random effects.
#' @param formula_random_effects_compile_list optional named list of
#' `random_effects_compile()` objects controlling which formula random-effect
#' blocks are sampled and which are compiled as marginalized structural blocks.
#' Defaults to \code{NULL}, which preserves the current all-sampled behavior.
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
#' @param jags_modules character vector specifying JAGS modules required by the
#' generated model syntax. Defaults to \code{NULL}.
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
JAGS_fit <- function(model_syntax, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL, formula_scale_list = NULL, formula_random_prior_list = NULL, formula_random_effects_compile_list = NULL,
                     chains = 4, adapt = 500, burnin = 1000, sample = 4000, thin = 1,
                     autofit = FALSE, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10, check_indicators = FALSE),
                     parallel = FALSE, cores = chains, silent = TRUE, seed = NULL,
                     add_parameters = NULL, required_packages = NULL, jags_modules = NULL, ...){

  .check_runjags()
  dots <- list(...)

  ### check input
  model_syntax <- .check_JAGS_syntax(model_syntax)
  JAGS_check_and_list_fit_settings(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed)
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(required_packages, "required_packages", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(jags_modules, "jags_modules", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_list(formula_list, "formula_list", allow_NULL = TRUE)
  check_list(formula_data_list, "formula_data_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))
  check_list(formula_prior_list, "formula_prior_list", check_names = names(formula_list), allow_other = FALSE, all_objects = TRUE, allow_NULL = is.null(formula_list))
  check_list(formula_random_prior_list, "formula_random_prior_list", check_names = names(formula_list), allow_other = FALSE, all_objects = FALSE, allow_NULL = TRUE)
  check_list(formula_random_effects_compile_list, "formula_random_effects_compile_list", check_names = names(formula_list), allow_other = FALSE, all_objects = FALSE, allow_NULL = TRUE)
  check_list(formula_scale_list, "formula_scale_list", check_names = names(formula_list), allow_other = FALSE, all_objects = FALSE, allow_NULL = TRUE)
  if(!is.null(formula_random_prior_list)){
    for(parameter in names(formula_random_prior_list)){
      .bt_check_prior_random(formula_random_prior_list[[parameter]])
    }
  }
  if(!is.null(formula_random_effects_compile_list)){
    for(parameter in names(formula_random_effects_compile_list)){
      .bt_check_random_effects_compile(formula_random_effects_compile_list[[parameter]])
    }
  }
  if(!is.null(formula_list)){
    for(parameter in names(formula_list)){
      if((is.language(formula_list[[parameter]]) ||
          inherits(formula_list[[parameter]], "BayesTools_random_effects")) &&
         .has_random_effects(formula_list[[parameter]]) &&
         (is.null(formula_random_prior_list) || is.null(formula_random_prior_list[[parameter]]))){
        stop(
          "JAGS_fit() requires 'formula_random_prior_list' with a prior_random() object for formula random effects in parameter '",
          parameter,
          "'.",
          call. = FALSE
        )
      }
    }
  }

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
        formula_scale  = if(!is.null(formula_scale_list)) formula_scale_list[[parameter]] else NULL,
        prior_random   = if(!is.null(formula_random_prior_list)) formula_random_prior_list[[parameter]] else NULL,
        random_effects_compile = if(!is.null(formula_random_effects_compile_list)) formula_random_effects_compile_list[[parameter]] else NULL)
    }

    # merge with the rest of the input
    prior_list     <- c(do.call(c, unname(lapply(formula_output, function(output) output[["prior_list"]]))), prior_list)
    data           <- c(do.call(c, unname(lapply(formula_output, function(output) output[["data"]]))),       data)
    formula_syntax <- paste0(lapply(formula_output, function(output) output[["formula_syntax"]]), collapse = "")

    # collect formula_scale information
    formula_scale_info <- lapply(formula_output, function(output) output[["formula_scale"]])
    formula_scale_info <- formula_scale_info[!sapply(formula_scale_info, is.null)]
    if(length(formula_scale_info) == 0) formula_scale_info <- NULL
    formula_design_info <- lapply(formula_output, function(output) output[["formula_design"]])
    formula_add_parameters <- unique(unlist(lapply(formula_output, function(output) output[["add_parameters"]])), use.names = FALSE)
    formula_jags_modules <- unique(unlist(lapply(formula_output, function(output) output[["jags_modules"]])), use.names = FALSE)
    formula_required_packages <- unique(unlist(lapply(formula_output, function(output) output[["required_packages"]])), use.names = FALSE)

    # add the formula syntax to the model syntax
    opening_bracket <- regexpr("{", model_syntax, fixed = TRUE)[1]
    syntax_start    <- substr(model_syntax, 1, opening_bracket)
    syntax_end      <- substr(model_syntax, opening_bracket + 1, nchar(model_syntax))
    model_syntax    <- paste0(syntax_start, "\n", formula_syntax, "\n", syntax_end)
  }else{
    formula_scale_info <- NULL
    formula_design_info <- NULL
    formula_add_parameters <- character()
    formula_jags_modules <- character()
    formula_required_packages <- character()
  }

  add_parameters <- unique(c(add_parameters, formula_add_parameters))
  jags_modules <- unique(c(jags_modules, formula_jags_modules))
  required_packages <- unique(c(required_packages, formula_required_packages))

  if(.JAGS_prior_list_uses_BayesTools_module(prior_list)){
    jags_modules <- unique(c(jags_modules, "BayesTools"))
    required_packages <- unique(c(required_packages, "BayesTools"))
  }
  prior_list <- .complete_factor_metadata_prior_list(prior_list)
  .bt_validate_jags_add_parameters(add_parameters, prior_list)

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
    .JAGS_load_modules(jags_modules, cl, warn = !silent)
    model_call <- c(
      model_call,
      method = "rjparallel",
      cl     = list(cl)
    )
  }else{
    .JAGS_require_packages(required_packages)
    .JAGS_load_modules(jags_modules, warn = !silent)
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

    converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], add_parameters = add_parameters, check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)
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

      converged  <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], add_parameters = add_parameters, check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)
      itteration <- itteration + 1

      if(isTRUE(dots[["is_JASP"]]))
        .JASP_progress_bar_tick()
    }
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "add_parameters") <- add_parameters
  attr(fit, "required_packages") <- required_packages
  attr(fit, "jags_modules") <- jags_modules
  if(!is.null(formula_scale_info)){
    # Keep formula_scale as a nested list keyed by parameter name
    # Each element contains the scaling info for that parameter's predictors
    attr(fit, "formula_scale") <- formula_scale_info
  }
  if(!is.null(formula_design_info)){
    attr(fit, "formula_design") <- formula_design_info
  }

  class(fit) <- c(class(fit), "BayesTools_fit")

  return(fit)
}

.bt_validate_jags_add_parameters <- function(add_parameters, prior_list){

  if(length(add_parameters) == 0L){
    return(invisible(TRUE))
  }
  if(any(!nzchar(add_parameters))){
    stop(
      "The 'add_parameters' argument cannot contain empty parameter names.",
      call. = FALSE
    )
  }
  if(length(prior_list) == 0L){
    return(invisible(TRUE))
  }

  prior_parameters <- unique(c(names(prior_list), JAGS_to_monitor(prior_list)))
  parameter_base <- function(x){
    sub("\\[.*$", "", x)
  }
  overlap <- add_parameters[
    parameter_base(add_parameters) %in% parameter_base(prior_parameters)
  ]
  if(length(overlap) > 0L){
    stop(
      "The 'add_parameters' argument must not include parameters already ",
      "monitored through 'prior_list': ",
      paste(unique(overlap), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
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
  jags_modules      <- attr(fit, "jags_modules")
  add_parameters    <- attr(fit, "add_parameters")
  formula_scale     <- attr(fit, "formula_scale")
  formula_design    <- attr(fit, "formula_design")
  prior_list        <- .complete_factor_metadata_prior_list(prior_list)
  if(is.null(add_parameters)){
    add_parameters <- character()
  }
  .bt_validate_jags_add_parameters(add_parameters, prior_list)
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)

  # parallel vs. not
  if(parallel){
    if(is.null(cores)){
      cores <- length(fit[["mcmc"]])
    }
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    .JAGS_require_packages(required_packages, cl)
    .JAGS_load_modules(jags_modules, cl, warn = !silent)
    refit_call <- list(
      runjags.object = fit,
      sample         = autofit_control[["sample_extend"]],
      method         = "rjparallel",
      cl             = cl,
      summarise      = FALSE
    )
  }else{
    .JAGS_require_packages(required_packages)
    .JAGS_load_modules(jags_modules, warn = !silent)
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

    converged <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]], add_parameters = add_parameters, check_indicators = autofit_control[["check_indicators"]], fail_fast = TRUE)

    # update the refit call
    if(!converged){
      itteration <- itteration + 1
      refit_call$runjags.object <- fit
    }
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "add_parameters") <- add_parameters
  attr(fit, "required_packages") <- required_packages
  attr(fit, "jags_modules") <- jags_modules
  if(!is.null(formula_scale)){
    attr(fit, "formula_scale") <- formula_scale
  }
  if(!is.null(formula_design)){
    attr(fit, "formula_design") <- formula_design
  }

  class(fit) <- c(class(fit), "BayesTools_fit")

  return(fit)
}
