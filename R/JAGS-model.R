#' @title Add JAGS prior
#'
#' @description Adds priors to JAGS syntax or
#' generates syntax for a single prior.
#'
#' @param syntax JAGS model syntax
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param prior prior distribution
#' @param parameter_name name of a parameter
#'
#' @export JAGS_add_priors
#' @export JAGS_prior.simple
#' @export JAGS_prior.PP
#' @export JAGS_prior.weightfunction
#' @name JAGS_add_priors
NULL

#' @rdname JAGS_add_priors
JAGS_add_priors           <- function(syntax, prior_list){

  # return the original syntax in case that no prior was specified
  if(length(prior_list) == 0){
    return(syntax)
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  if(!is.character(syntax) | length(syntax) != 1)
    stop("'syntax' must be a character vector of length 1.")
  if(!grepl("model", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("{", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("}", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")


  # identify parts of the syntax
  opening_bracket <- regexpr("{", syntax, fixed = TRUE)[1]
  syntax_start    <- substr(syntax, 1, opening_bracket)
  syntax_end      <- substr(syntax, opening_bracket + 1, nchar(syntax))

  # create the priors relevant syntax
  syntax_priors <- ""
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, JAGS_prior.weightfunction(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, JAGS_prior.PP(prior_list[[i]]))

    }else if(is.prior.simple(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, JAGS_prior.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  # merge everything back together
  syntax <- paste0(syntax_start, "\n", syntax_priors, "\n", syntax_end)

  return(syntax)
}
#' @rdname JAGS_add_priors
JAGS_prior.simple         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  if(!is.character(parameter_name) | length(parameter_name) != 1)
    stop("'parameter_name' must be a character vector of length 1.")

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
#' @rdname JAGS_add_priors
JAGS_prior.PP             <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    syntax <- JAGS_prior.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    syntax <- JAGS_prior.simple(prior, "PEESE")
  }

  return(syntax)
}
#' @rdname JAGS_add_priors
JAGS_prior.weightfunction <- function(prior){

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


#' @title Create initial values for JAGS model
#'
#' @description Creates initial values for priors.
#'
#' @param chains number of chains
#' @param seed seed for random number generation
#'
#' @inheritParams JAGS_add_priors
#' @export JAGS_get_inits
#' @export JAGS_init.simple
#' @export JAGS_init.PP
#' @export JAGS_init.weightfunction
#' @name JAGS_get_inits
NULL

#' @rdname JAGS_get_inits
JAGS_get_inits            <- function(prior_list, chains, seed){

  # return empty list in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  .check_n(chains, "chains")
  .check_x(seed, name = "seed")
  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
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

    temp_inits <- list()

    for(i in seq_along(prior_list)){

      if(is.prior.weightfunction(prior_list[[i]])){

        temp_inits <- c(temp_inits, JAGS_init.weightfunction(prior_list[[i]]))

      }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

        temp_inits <- c(temp_inits, JAGS_init.PP(prior_list[[i]]))

      }else if(is.prior.simple(prior_list[[i]])){

        temp_inits <- c(temp_inits, JAGS_init.simple(prior_list[[i]], names(prior_list)[i]))

      }
    }

    temp_inits[[".RNG.seed"]] <- seed + j
    temp_inits[[".RNG.name"]] <- "base::Super-Duper"

    inits[[j]] <- temp_inits
  }

  return(inits)
}
#' @rdname JAGS_get_inits
JAGS_init.simple          <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")

  if(prior[["distribution"]] == "point"){

    return()

  }else{
    init <- list()

    if(prior[["distribution"]] == "invgamma"){

      sampling_prior <- prior(
        "distribution" = "gamma",
        "parameters"   = list("shape" = prior$parameters[["shape"]], "rate" = prior$parameters[["scale"]]),
        "truncation"   = list("lower" = prior$truncation[["upper"]]^-1, "upper" = prior$truncation[["lower"]]^-1))
      init[[paste0("inv_", parameter_name)]] <- rng(1, sampling_prior)

    }else{

      init[[parameter_name]] <- rng(1, prior)

    }
  }

  return(init)
}
#' @rdname JAGS_get_inits
JAGS_init.PP              <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    init <- JAGS_init.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    init <- JAGS_init.simple(prior, "PEESE")
  }

  return(init)
}
#' @rdname JAGS_get_inits
JAGS_init.weightfunction  <- function(prior){

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


#' @title Creates list of monitored parameters for JAGS model
#'
#' @description Creates a vector of parameter names to be
#' monitored.
#'
#'
#' @inheritParams JAGS_add_priors
#' @export JAGS_to_monitor
#' @export JAGS_monitor.simple
#' @export JAGS_monitor.PP
#' @export JAGS_monitor.weightfunction
#' @name JAGS_to_monitor
NULL

#' @rdname JAGS_to_monitor
JAGS_to_monitor             <- function(prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    return("")
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # add the monitored parameters
  monitor <- character()
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      monitor <- c(monitor, JAGS_monitor.weightfunction(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      monitor <- c(monitor, JAGS_monitor.PP(prior_list[[i]]))

    }else if(is.prior.simple(prior_list[[i]])){

      monitor <- c(monitor, JAGS_monitor.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(monitor)
}
#' @rdname JAGS_to_monitor
JAGS_monitor.simple         <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  if(!is.character(parameter_name) | length(parameter_name) != 1)
    stop("'parameter_name' must be a character vector of length 1.")

  if(prior[["distribution"]] == "invgamma"){
    monitor <- c(parameter_name, paste0("inv_", parameter_name))
  }else{
    monitor <- parameter_name
  }

  return(monitor)
}
#' @rdname JAGS_to_monitor
JAGS_monitor.PP             <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    monitor <- JAGS_monitor.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    monitor <- JAGS_monitor.simple(prior, "PEESE")
  }

  return(monitor)
}
#' @rdname JAGS_to_monitor
JAGS_monitor.weightfunction <- function(prior){

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
