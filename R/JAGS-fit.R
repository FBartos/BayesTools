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
#'   \item{sample_extend}{number of samples between each convergence check. Defaults to
#'   \code{1000}.}
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
#' @return \code{JAGS_fit} returns an object of class 'runjags'.
#'
#' @seealso [JAGS_check_convergence()]
#' @export
JAGS_fit <- function(model_syntax, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL,
                     chains = 4, adapt = 500, burnin = 1000, sample = 4000, thin = 1,
                     autofit = FALSE, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000),
                     parallel = FALSE, cores = chains, silent = TRUE, seed = NULL,
                     add_parameters = NULL, required_packages = NULL){

  .check_runjags()

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
    user_silent.jags    <- runjags::runjags.getOption("silent.jags")
    user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  }

  start_time <- Sys.time()
  fit <- tryCatch(do.call(runjags::run.jags, model_call), error = function(e)e)

  if(inherits(fit, "error") & !silent)
    warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)

  if(autofit & !inherits(fit, "error")){

    converged <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]])

    while(!converged){

      if(!is.null(autofit_control[["max_time"]]) && difftime(Sys.time(), start_time, units = autofit_control[["max_time"]][["unit"]]) > autofit_control[["max_time"]][["time"]]){
        if(!silent){
          attr(fit, "warning") <- "The automatic model fitting was terminated due to the 'max_time' constraint."
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

      converged <- JAGS_check_convergence(fit, prior_list, autofit_control[["max_Rhat"]], autofit_control[["min_ESS"]], autofit_control[["max_error"]], autofit_control[["max_SD_error"]])
    }
  }

  # return user settings
  if(silent){
    runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags)
  }

  if(parallel){
    parallel::stopCluster(cl)
  }

  # add information to the fitted object
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax

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
JAGS_check_convergence <- function(fit, prior_list, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05){

  # check input
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_real(max_Rhat,     "max_Rhat",     lower = 1, allow_NULL = TRUE)
  check_real(min_ESS,      "min_ESS",      lower = 0, allow_NULL = TRUE)
  check_real(max_error,    "max_error",    lower = 0, allow_NULL = TRUE)
  check_real(max_SD_error, "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE)

  fails         <- NULL
  invisible(utils::capture.output(temp_summary <- suppressWarnings(summary(fit, silent.jags = TRUE))))

  # remove auxiliary and support parameters from the summary
  for(i in seq_along(prior_list)){
    if(is.prior.weightfunction(prior_list[[i]])){
      if(prior_list[[i]][["distribution"]] %in% c("one.sided", "two.sided")){
        temp_summary <- temp_summary[!grepl("eta", rownames(temp_summary)),,drop=FALSE]
      }
      temp_summary <- temp_summary[-max(grep("omega", rownames(temp_summary))),,drop=FALSE]
    }else if(is.prior.point(prior_list[[i]])){
      temp_summary <- temp_summary[rownames(temp_summary) != names(prior_list)[i],,drop=FALSE]
    }else if(is.prior.simple(prior_list[[i]]) && prior_list[[i]][["distribution"]] == "invgamma"){
      temp_summary <- temp_summary[rownames(temp_summary) != paste0("inv_",names(prior_list)[i]),,drop=FALSE]
    }
  }

  # check the convergence
  if(!is.null(max_Rhat)){
    temp_Rhat <- max(ifelse(is.na(temp_summary[, "psrf"]), 1, temp_summary[, "psrf"]))
    if(temp_Rhat > max_Rhat){
      fails <- c(fails, paste0("R-hat ", round(temp_Rhat, 3), " is larger than the set target (", max_Rhat, ")."))
    }
  }

  if(!is.null(min_ESS)){
    temp_ESS <- min(ifelse(is.na(temp_summary[, "SSeff"]), Inf, temp_summary[, "SSeff"]))
    if(temp_ESS < min_ESS){
      fails <- c(fails, paste0("ESS ", round(temp_ESS), " is lower than the set target (", min_ESS, ")."))
    }
  }

  if(!is.null(max_error)){
    temp_error    <- max(ifelse(is.na(temp_summary[, "MCerr"]), 0, temp_summary[, "MCerr"]))
    if(temp_error > max_error){
      fails <- c(fails, paste0("MCMC error ", round(temp_error, 5), " is larger than the set target (", max_error, ")."))
    }
  }

  if(!is.null(max_SD_error)){
    temp_error_SD <- max(ifelse(is.na(temp_summary[, "MC%ofSD"]),   0, temp_summary[, "MC%ofSD"]))
    if(temp_error_SD/100 > max_SD_error){
      fails <- c(fails, paste0("MCMC SD error ", round(temp_error_SD/100, 4), " is larger than the set target (", max_SD_error, ")."))
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

  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.weightfunction(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

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

  if(is.prior.PET(prior[["variable"]]) | is.prior.PEESE(prior[["variable"]]) | is.prior.weightfunction(prior[["variable"]]))
    stop("Spike and slab functionality is not implemented for publication bias prior distributions.")
  if(is.prior.spike_and_slab(prior[["variable"]]))
     stop("Spike and slab prior distribution cannot be nested inside of a spike and slab prior distribution.")


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
    temp_inits[[".RNG.name"]] <- "base::Super-Duper"

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
    parameter_name,
    JAGS_to_monitor(prior_variable),
    JAGS_to_monitor(prior_inclusion),
    paste0(parameter_name, "_indicator")
  )

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

  check_list(autofit_control, "autofit_control", check_names = c("max_Rhat", "min_ESS", "max_error", "max_SD_error",  "max_time", "sample_extend"), call = call)
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

  return(invisible(autofit_control))
}
