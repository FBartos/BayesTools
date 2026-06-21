#' @title Creates a prior distribution
#'
#' @description \code{prior} creates a prior distribution.
#' The prior can be visualized by the \code{plot} function.
#'
#' @param distribution name of the prior distribution. The
#' possible options are
#' \describe{
#'   \item{\code{"point"}}{for a point density characterized by a
#'   \code{location} parameter.}
#'   \item{\code{"normal"}}{for a normal distribution characterized
#'   by a \code{mean} and \code{sd} parameters.}
#'   \item{\code{"lognormal"}}{for a lognormal distribution characterized
#'   by a \code{meanlog} and \code{sdlog} parameters.}
#'   \item{\code{"cauchy"}}{for a Cauchy distribution characterized
#'   by a \code{location} and \code{scale} parameters. Internally
#'   converted into a generalized t-distribution with \code{df = 1}.}
#'   \item{\code{"t"}}{for a generalized t-distribution characterized
#'   by a \code{location}, \code{scale}, and \code{df} parameters.}
#'   \item{\code{"gamma"}}{for a gamma distribution characterized
#'   by either \code{shape} and \code{rate}, or \code{shape} and
#'   \code{scale} parameters. The latter is internally converted to
#'   the \code{shape} and \code{rate} parametrization}
#'   \item{\code{"invgamma"}}{for an inverse-gamma distribution
#'   characterized by a \code{shape} and \code{scale} parameters. The
#'   JAGS part uses a 1/gamma distribution with a shape and rate
#'   parameter.}
#'   \item{\code{"beta"}}{for a beta distribution
#'   characterized by an \code{alpha} and \code{beta} parameters.}
#'   \item{\code{"exp"}}{for an exponential distribution
#'   characterized by either \code{rate} or \code{scale}
#'   parameter. The later is internally converted to
#'   \code{rate}.}
#'   \item{\code{"uniform"}}{for a uniform distribution defined on a
#'   range from \code{a} to \code{b}}
#'   \item{\code{"moment"}}{for a Johnson-Rossell moment prior characterized
#'   by exactly one of \code{mode} or \code{tau}, an optional \code{order}
#'   that defaults to 1, and an optional \code{location} that defaults to 0.}
#'   \item{\code{"invmoment"}}{for a Johnson-Rossell inverse-moment prior
#'   characterized by exactly one of \code{mode} or \code{tau}, a \code{df},
#'   an optional \code{order} that defaults to 1, and an optional
#'   \code{location} that defaults to 0.}
#'   \item{\code{"dirichlet"}}{for a Dirichlet distribution over a simplex,
#'   characterized by a positive concentration vector \code{alpha}. Density
#'   and plot methods show one beta marginal per simplex coordinate.}
#' }
#' @param parameters list of appropriate parameters for a given
#' \code{distribution}.
#' @param truncation list with two elements, \code{lower} and
#' \code{upper}, that define the lower and upper truncation of the
#' distribution. Defaults to \code{list(lower = -Inf, upper = Inf)}.
#' The truncation is automatically set to the bounds of the support.
#' @param prior_weights prior odds associated with a given distribution.
#' The value is passed into the model fitting function, which creates models
#' corresponding to all combinations of prior distributions for each of
#' the model parameters and sets the model priors odds to the product
#' of its prior distributions.
#'
#' @details Moment and inverse-moment priors are symmetric nonlocal priors
#' with zero density at \code{location}. The \code{mode} parameter is the
#' positive distance from \code{location} to each symmetric mode; supplied
#' negative values are converted to their absolute value. Positional parameters
#' are interpreted as \code{mode} for \code{"moment"} and \code{mode, df}
#' for \code{"invmoment"}. Use named parameters for \code{tau}, \code{order},
#' \code{location}, or the inverse-moment \code{nu} alias. Positional input
#' always treats the first parameter as \code{mode}. For moment priors,
#' \code{tau = mode^2 / (2 * order)}. For inverse-moment priors,
#' \code{tau = mode^2 * ((df + 1) / (2 * order))^(1 / order)}.
#' Aliases \code{"pmom"}, \code{"pimom"}, and \code{"inversemoment"} are
#' accepted, and inverse-moment \code{df} may also be supplied as \code{nu}.
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # create a half-normal standard normal prior distribution
#' p2 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1),
#' truncation = list(lower = 0, upper = Inf))
#'
#' # the prior distribution can be visualized using the plot function
#' # (see ?plot.prior for all options)
#' plot(p1)
#'
#' # create nonlocal priors centered on zero
#' p3 <- prior("moment", list(mode = 0.5))
#' p4 <- prior("invmoment", list(mode = 0.5, df = 3))
#'
#' # tau can be specified directly when named
#' p5 <- prior("pmom", list(tau = 0.125))
#' p6 <- prior("pimom", list(tau = 0.5, nu = 3))
#'
#' @return \code{prior} and \code{prior_none} return an object of class 'prior'.
#' A named list containing the distribution name, parameters, and prior weights.
#'
#' @name prior
#' @export prior
#' @export prior_none
#' @seealso [plot.prior()], \link[stats]{Normal}, \link[stats]{Lognormal}, \link[stats]{Cauchy},
#' \link[stats]{Beta}, \link[stats]{Exponential},
#' \link[extraDistr]{LocationScaleT}, \link[extraDistr]{InvGamma}.

#' @rdname prior
prior <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_weights = 1){

  # general input check (detailed checks are performed withing the constructors)
  check_char(distribution, "distribution")
  check_list(parameters, "parameters")
  #sapply(seq_along(parameters), function(i)check_real(parameters[[i]], names(parameters[[i]]), check_length = 0))
  check_list(truncation, "truncation")
  .check_prior_weight(prior_weights)

  # clean the input name
  distribution <- .prior_clean_input_name(distribution)

  if(distribution %in% c("norm", "normal")){
    distribution <- "normal"
  }else if(distribution %in% c("lnorm", "lognormal")){
    distribution <- "lognormal"
  }else if(distribution %in% c("t", "student", "studentt")){
    distribution <- "t"
  }else if(distribution %in% c("cauchy")){
    distribution <- "cauchy"
  }else if(distribution %in% c("invgamma", "inversegamma")){
    distribution <- "invgamma"
  }else if(distribution %in% c("gamma")){
    distribution <- "gamma"
  }else if(distribution %in% c("beta")){
    distribution <- "beta"
  }else if(distribution %in% c("bernoulli", "bernouli", "bern")){
    distribution <- "bernoulli"
  }else if(distribution %in% c("exp", "exponential")){
    distribution <- "exp"
  }else if(distribution %in% c("uniform", "unif")){
    distribution <- "uniform"
  }else if(distribution %in% c("moment", "pmom")){
    distribution <- "moment"
  }else if(distribution %in% c("invmoment", "inversemoment", "pimom")){
    distribution <- "invmoment"
  }else if(distribution %in% c("point", "spike")){
    distribution <- "point"
  }else if(distribution %in% c("multivariatenorm", "multivariatenormal", "mnorm", "mnormal")){
    distribution <- "mnormal"
  }else if(distribution %in% c("multivariatet", "multivariatestudent", "mt", "mstudent")){
    distribution <- "mt"
  }else if(distribution %in% c("multivariatecauchy", "mcauchy")){
    distribution <- "mcauchy"
  }else if(distribution %in% c("mpoint", "mspike")){
    distribution <- "mpoint"
  }else if(distribution %in% c("dirichlet", "simplex")){
    distribution <- "dirichlet"
  }else{
    stop(paste0("The specified distribution name '", distribution,"' is not known. Please, see '?prior' for more information about supported prior distributions."))
  }

  # check the passed settings
  output <- do.call(paste0(".prior_", distribution), list(parameters = parameters, truncation = truncation))

  # add the prior odds
  output$prior_weights <- prior_weights

  return(output)
}

#' @rdname prior
prior_none <- function(prior_weights = 1){

  .check_prior_weight(prior_weights)

  out <- list()
  out$distribution <- "none"
  out$prior_weights   <- prior_weights
  class(out)       <- c("prior", "prior.none")

  return(out)
}


#' @title Creates a prior distribution for PET or PEESE models
#'
#' @description \code{prior} creates a prior distribution for fitting a PET or
#' PEESE style models in RoBMA. The prior distribution can be visualized
#' by the \code{plot} function.
#'
#' @examples
#' # create a half-Cauchy prior distribution
#' # (PET and PEESE specific functions automatically set lower truncation at 0)
#' p1 <- prior_PET(distribution = "Cauchy", parameters = list(location = 0, scale = 1))
#'
#' plot(p1)
#'
#' @return \code{prior_PET} and \code{prior_PEESE} return an object of class 'prior'.
#'
#' @inheritParams prior
#' @export prior_PET
#' @export prior_PEESE
#' @seealso [plot.prior()], [prior()]
#' @name prior_PP
NULL

#' @rdname prior_PP
prior_PET   <- function(distribution, parameters, truncation = list(lower = 0, upper = Inf), prior_weights = 1){

  output <- prior(distribution, parameters, truncation, prior_weights)

  class(output) <- c(class(output), "prior.PET")

  return(output)
}
#' @rdname prior_PP
prior_PEESE <- function(distribution, parameters, truncation = list(lower = 0, upper = Inf), prior_weights = 1){

  output <- prior(distribution, parameters, truncation, prior_weights)

  class(output) <- c(class(output), "prior.PEESE")

  return(output)
}

#' @title Creates a prior distribution for factors
#'
#' @description \code{prior_factor} creates a prior distribution for fitting
#' models with factor predictors. (Note that results across different operating
#' systems might vary due to differences in JAGS numerical precision.)
#'
#' @param contrast type of contrast for the prior distribution. The possible options are
#' \describe{
#'   \item{\code{"meandif"}}{for contrast centered around the grand mean
#'   with equal marginal distributions, making the prior distribution exchangeable
#'   across factor levels. In contrast to \code{"orthonormal"}, the marginal distributions
#'   are identical regardless of the number of factor levels and the specified prior
#'   distribution corresponds to the difference from grand mean for each factor level.
#'   Only supports \code{distribution = "mnormal"} and \code{distribution = "mt"}
#'   which generates the corresponding multivariate normal/t distributions.}
#'   \item{\code{"orthonormal"}}{for contrast centered around the grand mean
#'   with equal marginal distributions, making the prior distribution exchangeable
#'   across factor levels. Only supports \code{distribution = "mnormal"} and
#'   \code{distribution = "mt"} which generates the corresponding multivariate normal/t
#'   distributions.}
#'   \item{\code{"treatment"}}{for contrasts using the first level as a comparison
#'   group and setting equal prior distribution on differences between the individual
#'   factor levels and the comparison level.}
#'   \item{\code{"independent"}}{for contrasts specifying dependent prior distribution
#'   for each factor level (note that this leads to an overparameterized model if the
#'   intercept is included).}
#' }
#'
#'
#' @examples
#' # create an orthonormal prior distribution
#' p1 <- prior_factor(distribution = "mnormal", contrast = "orthonormal",
#'                    parameters = list(mean = 0, sd = 1))
#'
#' @return return an object of class 'prior'.
#'
#' @inheritParams prior
#' @export  prior_factor
#' @seealso [prior()]
prior_factor <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_weights = 1, contrast = "meandif"){

  # general input check (detailed checks are performed withing the constructors)
  check_char(contrast, "contrast", allow_values = c("meandif", "orthonormal", "treatment", "dummy", "independent"))

  # check its compatibility with the contrasts
  if(contrast %in% c("meandif", "orthonormal")){

    # add the (yet unspecified) dimensions parameter
    if(is.null(names(parameters))){
      parameters <- c(parameters, NA)
    }else{
      parameters[["K"]] <- NA
    }

    # change spike/point into mpoint dispatch
    if(distribution %in% c("spike", "mspike", "point", "mpoint"))
      distribution <- "mpoint"

    if(!distribution %in% c("multivariatenorm", "multivariatenormal", "mnorm", "mnormal",
                            "multivariatet", "multivariatestudent", "mt", "mstudent",
                            "multivariatecauchy", "mcauchy",
                            "mpoint"))
      stop(paste0("'", contrast,"' contrasts require multivariate prior disribution."))

    # generate the prior object
    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    if(!is.prior.vector(output))
      stop(paste0("'", contrast,"' contrasts require vector prior distribution."))
    if(output[["distribution"]] != "mpoint" && !all(sapply(output[["truncation"]], is.infinite)))
      stop(paste0("'", contrast,"' contrasts do not support truncation."))

    class(output) <- c(class(output), "prior.factor", paste0("prior.", contrast))

  }else if(contrast %in% c("treatment", "dummy")){

    # generate the prior object
    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    if(!is.prior.simple(output))
      stop("'treatment' contrasts require univariate prior distribution.")

    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    class(output) <- c(class(output), "prior.factor", "prior.treatment")

  }else if(contrast  == "independent"){

    # generate the prior object
    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    if(!is.prior.simple(output))
      stop("'independent' contrasts require univariate prior distribution.")

    output <- prior(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights)

    class(output) <- c(class(output), "prior.factor", "prior.independent")
  }

  return(output)
}


#' @title Creates a spike and slab prior distribution
#'
#' @description \code{prior_spike_and_slab} creates a spike and slab prior
#' distribution corresponding to the specification in
#' \insertCite{kuo1998variable;textual}{BayesTools} (see
#' \insertCite{ohara2009review;textual}{BayesTools} for further details). I.e.,
#' a prior distribution is multiplied by an independent indicator with values
#' either zero or one.
#'
#' @param prior_parameter a prior distribution for the parameter
#' @param prior_inclusion a prior distribution for the inclusion probability. The
#' inclusion probability must be bounded within 0 and 1 range. Defaults to
#' \code{prior("spike", parameters = list(location = 0.5))} which corresponds to 1/2
#' prior probability of including the slab prior distribution (but other prior
#' distributions, like beta etc can be also specified).
#'
#'
#' @examples
#' # create a spike and slab prior distribution
#' p1 <- prior_spike_and_slab(
#'    prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
#'    prior_inclusion = prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))
#' )
#'
#' @return return an object of class 'prior'.
#'
#' @inheritParams prior
#' @seealso [prior()]
#' @export
prior_spike_and_slab <- function(prior_parameter,
                                 prior_inclusion = prior(distribution = "spike", parameters = list(location = 0.5)),
                                 prior_weights = 1){
  if(!is.prior(prior_parameter))
    stop("'prior_parameter' must be a prior distribution")
  if(!is.prior(prior_inclusion))
    stop("'prior_inclusion' must be a prior distribution")
  .check_spike_and_slab_inclusion_prior(prior_inclusion)
  .check_prior_weight(prior_weights)
  if(is.prior.point(prior_inclusion) && (prior_inclusion$parameters[["location"]] < 0 | prior_inclusion$parameters[["location"]] > 1))
    stop("The probability parameter of 'prior_inclusion' must be within 0 and 1.")
  if(!is.prior.point(prior_inclusion) && (prior_inclusion$truncation[["lower"]] < 0 | prior_inclusion$truncation[["upper"]] > 1))
    stop("The range of the probability parameter (set via the 'truncation' argument) of 'prior_inclusion' must be within 0 and 1.")

  # Create the spike component (point at 0)
  if(is.prior.factor(prior_parameter)){
    # For factor priors, create a factor spike
    priors_type <- .get_prior_factor_list_type(list(prior_parameter))
    contrast_type <- gsub("prior.", "", priors_type[["class"]], fixed = TRUE)
    
    spike_component <- prior_factor(
      distribution = "point",
      parameters   = list(location = 0),
      contrast     = contrast_type
    )
  } else {
    # For simple priors, create a simple spike
    spike_component <- prior(
      distribution = "point", 
      parameters   = list(location = 0)
    )
  }
  
  # Create the mixture using the mixture backend
  mixture_output <- prior_mixture(
    prior_list = list(prior_parameter, spike_component),
    components = c("alternative", "null")
  )
  
  # Store inclusion prior as attribute so it can be retrieved by helper functions
  attr(mixture_output, "inclusion_prior") <- prior_inclusion
  attr(mixture_output, "model_prior_weights") <- prior_weights
  
  # Add spike_and_slab classes for specialized behavior while keeping mixture functionality
  if(is.prior.factor(prior_parameter)){
    # obtain and store the contrast type
    priors_type <- .get_prior_factor_list_type(list(prior_parameter))
    
    attr(prior_parameter, "K") <- priors_type[["K"]]
    class(mixture_output) <- c("prior", "prior.spike_and_slab", "prior.factor_spike_and_slab", 
                              class(mixture_output)[-1], priors_type[["class"]])
  }else if(is.prior.simple(prior_parameter)){
    class(mixture_output) <- c("prior", "prior.spike_and_slab", "prior.simple_spike_and_slab", 
                              class(mixture_output)[-1])
  }else{
    stop("The 'prior_parameter' must be either a simple or factor prior distribution.")
  }

  return(mixture_output)
}

.check_spike_and_slab_inclusion_prior <- function(prior_inclusion){

  scalar_probability_prior <- is.prior.simple(prior_inclusion) &&
    !is.prior.vector(prior_inclusion) &&
    !is.prior.factor(prior_inclusion) &&
    !is.prior.simplex(prior_inclusion) &&
    !is.prior.weightfunction(prior_inclusion) &&
    !is.prior.mixture(prior_inclusion) &&
    !is.prior.spike_and_slab(prior_inclusion) &&
    !is.prior.PET(prior_inclusion) &&
    !is.prior.PEESE(prior_inclusion) &&
    !is_prior_phacking(prior_inclusion) &&
    !is_prior_bias(prior_inclusion)

  if(!scalar_probability_prior){
    stop("'prior_inclusion' must be a scalar probability prior.", call. = FALSE)
  }

  invisible(TRUE)
}

# Helper functions to extract variable and inclusion from spike_and_slab mixture structure
.get_spike_and_slab_variable <- function(spike_and_slab_prior) {
  if (!is.prior.spike_and_slab(spike_and_slab_prior)) {
    stop("This function only works with spike_and_slab priors")
  }
 
  # Find the alternative component (this is the variable/slab part)
  components    <- attr(spike_and_slab_prior, "components") 
  alternative_idx <- which(components == "alternative")
  
  return(spike_and_slab_prior[[alternative_idx]])
}

.get_spike_and_slab_inclusion <- function(spike_and_slab_prior) {
  if (!is.prior.spike_and_slab(spike_and_slab_prior)) {
    stop("This function only works with spike_and_slab priors")
  }
  
  # For backward compatibility, use stored inclusion if available
  if (!is.null(spike_and_slab_prior[["inclusion"]])) {
    return(spike_and_slab_prior[["inclusion"]])
  }
  
  # Get inclusion prior from attribute
  inclusion_prior <- attr(spike_and_slab_prior, "inclusion_prior")
  return(inclusion_prior)
}

# Setter functions to allow modifying variable and inclusion components
.set_spike_and_slab_variable_attr <- function(spike_and_slab_prior, attr_name, value) {
  if (!is.prior.spike_and_slab(spike_and_slab_prior)) {
    stop("This function only works with spike_and_slab priors")
  }
  
  # Find the alternative component (this is the variable/slab part)
  components <- attr(spike_and_slab_prior, "components") 
  alternative_idx <- which(components == "alternative")

  # Set attribute on the variable component
  attr(spike_and_slab_prior[[alternative_idx]], attr_name) <- value
  
  return(spike_and_slab_prior)
}




#' @title Creates a mixture of prior distributions
#' @description \code{prior_mixture} creates a mixture of prior distributions.
#' This is a more generic version of the \code{prior_spike_and_slab} function.
#'
#' @param prior_list a list of prior distributions to be mixed.
#' @param is_null a logical vector indicating which of the prior distributions
#' should be considered as a null distribution. Defaults to \code{rep(FALSE, length(prior_list))}.
#' @param components a character vector indicating which of the prior distributions
#' belong to the same mixture component (this is an alternative specification to the \code{is_null} argument).
#' Defaults to \code{NULL} (i.e., \code{is_null} is used.
#'
#' @seealso [prior()]
#' @export
prior_mixture <- function(prior_list, is_null = rep(FALSE, length(prior_list)), components = NULL){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  check_bool(is_null, "is_null",       check_length = length(prior_list), allow_NULL = TRUE)
  check_char(components, "components", check_length = length(prior_list), allow_NULL = TRUE)
  if(is.null(is_null) && is.null(components))
    stop("Either 'is_null' or 'components' must be specified.")

  if(is.null(components)){
    components <- ifelse(is_null, "null", "alternative")
  }

  for(i in seq_along(prior_list)){
    attr(prior_list[[i]], "component") <- components[i]
  }


  # distinguish normal, factor, and publication bias mixture priors
  if(any(sapply(prior_list, is.prior.factor))){

    # test that the prior is either a factor prior or a spike prior
    if(!all(sapply(prior_list, is.prior.factor) | sapply(prior_list, is.prior.point) | sapply(prior_list, is.prior.none)))
      stop("Factor prior mixture requires that all priors are either factor priors or spike prior distributions")

    # obtain and store the contrast type
    priors_type <- .get_prior_factor_list_type(prior_list)

    # change prior none/spikes into factor prior spikes
    for(i in seq_along(prior_list)){
      if(is.prior.point(prior_list[[i]])){
        # Save mixture metadata before recreating the factor-compatible prior.
        component_attr <- attr(prior_list[[i]], "component")
        prior_weight   <- .prior_model_weight(prior_list[[i]])
        prior_list[[i]] <- prior_factor(
          distribution = "point",
          parameters   = list(location = prior_list[[i]][["parameters"]][["location"]]),
          prior_weights = prior_weight,
          contrast     = gsub("prior.", "", priors_type[["class"]], fixed = TRUE)
        )
        attr(prior_list[[i]], "component") <- component_attr
      }else if(is.prior.none(prior_list[[i]])){
        # Save mixture metadata before recreating the factor-compatible prior.
        component_attr <- attr(prior_list[[i]], "component")
        prior_weight   <- .prior_model_weight(prior_list[[i]])
        prior_list[[i]] <- prior_factor(
          distribution = "point",
          parameters   = list(location = 0),
          prior_weights = prior_weight,
          contrast     = gsub("prior.", "", priors_type[["class"]], fixed = TRUE)
        )
        attr(prior_list[[i]], "component") <- component_attr
      }
    }

    attr(prior_list, "K")  <- priors_type[["K"]]
    class(prior_list)      <- c("prior", "prior.factor_mixture", "prior.mixture", priors_type[["class"]])


  }else if(any(sapply(prior_list, is.prior.PET)) || any(sapply(prior_list, is.prior.PEESE)) ||
           any(sapply(prior_list, is.prior.weightfunction)) || any(sapply(prior_list, is_prior_phacking)) ||
           any(sapply(prior_list, is_prior_bias))){

    # test that the prior is either a PET, PEESE, weightfunction, p-hacking,
    # composed bias, or none prior
    if(!all(sapply(prior_list, is.prior.PET) | sapply(prior_list, is.prior.PEESE) |
            sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, is_prior_phacking) |
            sapply(prior_list, is_prior_bias) | sapply(prior_list, is.prior.none)))
      stop("Publication-bias prior mixtures require PET, PEESE, weightfunction, p-hacking, composed bias, or none prior distributions.")

    class(prior_list) <- c("prior", "prior.bias_mixture", "prior.mixture")


  }else if(any(sapply(prior_list, is.prior.simple))){

    # test that all priors are simple priors
    if(!all(sapply(prior_list, is.prior.simple) | sapply(prior_list, is.prior.none)))
      stop("Simple prior mixture requires that all priors are simple prior distributions")

    # change none into prior spikes
    for(i in seq_along(prior_list)){
      if(is.prior.none(prior_list[[i]])){
        prior_weight <- .prior_model_weight(prior_list[[i]])
        prior_list[[i]] <- prior(
          distribution = "point",
          parameters   = list(location = 0),
          prior_weights = prior_weight
        )
      }
    }

    class(prior_list) <- c("prior", "prior.simple_mixture", "prior.mixture")

  }else{
    stop("The prior mixture must contain either factors, publication bias components, or simple prior distributions.")
  }

  attr(prior_list, "components")    <- components
  attr(prior_list, "prior_weights") <- sapply(prior_list, .prior_model_weight)

  return(prior_list)
}

