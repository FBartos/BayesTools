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
#'   \code{scale} parameters. The later is internally converted to
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
  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

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

  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

  out <- list()
  out$distribution <- "none"
  out$prior_weights   <- prior_weights
  class(out)       <- c("prior", "prior.none")

  return(out)
}


#' @title Creates a prior distribution for a weight function
#'
#' @description \code{prior_weightfunction} creates a prior distribution for
#' fitting a RoBMA selection model. The \code{side} and \code{steps} arguments
#' define the p-value bins, and the \code{weights} argument defines the prior on
#' the publication weights in those bins.
#'
#' @param side side geometry. Either \code{"one-sided"} or \code{"two-sided"}.
#' @param steps increasing p-value cut points between 0 and 1.
#' @param weights a weight-prior object created by \code{wf_cumulative()},
#' \code{wf_fixed()}, or \code{wf_independent()}.
#' @param reference reference bin. Currently only \code{"most_significant"} is
#' supported and fixes the most significant bin to \code{omega = 1}.
#' @param prior_weights prior odds associated with a given distribution.
#'
#' @examples
#' p1 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05, .10),
#'   weights = wf_cumulative(alpha = c(1, 1, 1))
#' )
#'
#' p2 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05),
#'   weights = wf_independent(prior("beta", list(1, 1)))
#' )
#'
#' @return \code{prior_weightfunction} returns an object of class 'prior'.
#'
#' @export prior_weightfunction
#' @seealso [plot.prior()]
prior_weightfunction <- function(side = "one-sided", steps = c(.025, .05),
                                 weights = wf_cumulative(),
                                 reference = "most_significant",
                                 prior_weights = 1){

  check_char(side, "side")
  check_real(steps, "steps", check_length = 0)
  check_char(reference, "reference", allow_values = "most_significant")
  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

  side <- .weightfunction_normalize_side(side)
  steps <- .weightfunction_validate_steps(steps)
  bins <- .weightfunction_bins(steps)
  weights <- .weightfunction_validate_weights(weights, n_bins = nrow(bins), reference = reference)
  truncation <- .weightfunction_weights_truncation(weights)

  output <- list(
    distribution  = "weightfunction",
    side          = side,
    steps         = steps,
    bins          = bins,
    reference     = reference,
    weights       = weights,
    parameters    = list(steps = steps),
    truncation    = truncation,
    prior_weights = prior_weights
  )

  class(output) <- c("prior", "prior.weightfunction")

  return(output)
}

#' @rdname prior_weightfunction
#' @param alpha positive cumulative-Dirichlet concentration parameters. If
#' omitted, a flat Dirichlet prior is used with one concentration parameter per
#' bin.
#' @export
wf_cumulative <- function(alpha = NULL){

  if(!is.null(alpha)){
    check_real(alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = 0)
  }

  out <- list(type = "cumulative", alpha = alpha)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.cumulative")
  return(out)
}

#' @rdname prior_weightfunction
#' @param omega fixed publication weights, one per bin.
#' @export
wf_fixed <- function(omega){

  check_real(omega, "omega", lower = 0, upper = 1, check_length = 0)

  out <- list(type = "fixed", omega = omega)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.fixed")
  return(out)
}

#' @rdname prior_weightfunction
#' @param prior prior distribution for each non-reference weight.
#' @param scale latent scale for independent weights. \code{"omega"} places the
#' prior directly on the non-negative publication weight. \code{"log_omega"}
#' places the prior on \code{log(omega)} and transforms with
#' \code{omega = exp(log_omega)}, allowing weights above one whenever the log
#' prior assigns mass above zero.
#' @export
wf_independent <- function(prior, scale = "omega"){

  .check_prior(prior)
  if(!is.prior.simple(prior) || is.prior.point(prior) || is.prior.discrete(prior)){
    stop("'prior' must be a continuous simple prior distribution.", call. = FALSE)
  }

  check_char(scale, "scale", allow_values = c("omega", "log_omega", "log"))
  if(scale == "log"){
    scale <- "log_omega"
  }

  if(scale == "omega"){
    if(prior$truncation[["lower"]] < 0){
      stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
    }
  }

  out <- list(type = "independent", scale = scale, prior = prior)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.independent")
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
        # Save the component attribute before recreating
        component_attr <- attr(prior_list[[i]], "component")
        prior_list[[i]] <- prior_factor(
          distribution = "point",
          parameters   = list(location = prior_list[[i]][["parameters"]][["location"]]),
          contrast     = gsub("prior.", "", priors_type[["class"]], fixed = TRUE)
        )
        # Restore the component attribute
        attr(prior_list[[i]], "component") <- component_attr
      }else if(is.prior.none(prior_list[[i]])){
        # Save the component attribute before recreating
        component_attr <- attr(prior_list[[i]], "component")
        prior_list[[i]] <- prior_factor(
          distribution = "point",
          parameters   = list(location = 0),
          contrast     = gsub("prior.", "", priors_type[["class"]], fixed = TRUE)
        )
        # Restore the component attribute
        attr(prior_list[[i]], "component") <- component_attr
      }
    }

    attr(prior_list, "K")  <- priors_type[["K"]]
    class(prior_list)      <- c("prior", "prior.factor_mixture", "prior.mixture", priors_type[["class"]])


  }else if(any(sapply(prior_list, is.prior.PET)) || any(sapply(prior_list, is.prior.PEESE)) || any(sapply(prior_list, is.prior.weightfunction))){

    # test that the prior is either a PET, PEESE, or weightfunction prior
    if(!all(sapply(prior_list, is.prior.PET) | sapply(prior_list, is.prior.PEESE) | sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, is.prior.none)))
      stop("PET/PEESE/weightfunction prior mixture requires that all priors are either PET, PEESE, or weightfunction prior distributions")

    class(prior_list) <- c("prior", "prior.bias_mixture", "prior.mixture")


  }else if(any(sapply(prior_list, is.prior.simple))){

    # test that all priors are simple priors
    if(!all(sapply(prior_list, is.prior.simple) | sapply(prior_list, is.prior.none)))
      stop("Simple prior mixture requires that all priors are simple prior distributions")

    # change none into prior spikes
    for(i in seq_along(prior_list)){
      if(is.prior.none(prior_list[[i]])){
        prior_list[[i]] <- prior(
          distribution = "point",
          parameters   = list(location = 0)
        )
      }
    }

    class(prior_list) <- c("prior", "prior.simple_mixture", "prior.mixture")

  }else{
    stop("The prior mixture must contain either factors, publication bias components, or simple prior distributions.")
  }

  attr(prior_list, "components")    <- components
  attr(prior_list, "prior_weights") <- sapply(prior_list, function(p) p[["prior_weights"]])

  return(prior_list)
}

#### functions for constructing prior distributions ####
.prior_normal    <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("mean", "sd"), "normal")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$mean, "mean")
  .check_parameter(parameters$sd,   "sd")
  .check_parameter_positive(parameters$sd, "sd")

  # add the values to the output
  output$distribution <- "normal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_lognormal <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("meanlog", "sdlog"), "lognormal")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$meanlog, "meanlog")
  .check_parameter(parameters$sdlog,   "sdlog")
  .check_parameter_positive(parameters$sdlog, "sdlog")

  # add the values to the output
  output$distribution <- "lognormal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_cauchy    <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale"), "Cauchy")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter_positive(parameters$scale, "scale")

  # deal with as with a t-distribution
  parameters$df <- 1

  output$distribution <- "t"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_t         <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "df"), "student-t")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter(parameters$df,       "df")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_positive(parameters$df,    "df")

  # add the values to the output
  output$distribution <- "t"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_gamma     <- function(parameters, truncation){

  output <- list()

  # deal with possible scale parametrization
  if(!is.null(names(parameters))){
    if("scale" %in% names(parameters)){
      parameters <- .check_and_name_parameters(parameters, c("shape", "scale"), "gamma")
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
  }

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("shape", "rate"), "gamma")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$shape, "shape")
  .check_parameter(parameters$rate,  "rate")
  .check_parameter_positive(parameters$shape, "shape")
  .check_parameter_positive(parameters$rate,  "rate")

  # add the values to the output
  output$distribution <- "gamma"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_invgamma  <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("shape", "scale"), "invgamma")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$shape, "shape")
  .check_parameter(parameters$scale, "scale")
  .check_parameter_positive(parameters$shape, "shape")
  .check_parameter_positive(parameters$scale, "scale")

  # add the values to the output
  output$distribution <- "invgamma"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_exp       <- function(parameters, truncation){

  output <- list()

  # deal with possible scale parametrization
  if(!is.null(names(parameters))){
    if("scale" %in% names(parameters)){
      parameters <- .check_and_name_parameters(parameters, c("scale"), "exp")
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
  }

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("rate"), "exp")
  truncation <- .check_and_set_truncation(truncation, lower = 0)

  # check individual parameters
  .check_parameter(parameters$rate, "rate")
  .check_parameter_positive(parameters$rate, "rate")

  # add the values to the output
  output$distribution <- "exp"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_beta      <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("alpha", "beta"), "beta")
  truncation <- .check_and_set_truncation(truncation, lower = 0, upper = 1)

  # check individual parameters
  .check_parameter(parameters$alpha, "alpha")
  .check_parameter(parameters$beta,  "beta")
  .check_parameter_positive(parameters$alpha, "alpha")
  .check_parameter_positive(parameters$beta,  "beta")

  # add the values to the output
  output$distribution <- "beta"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_bernoulli <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, "probability", "bernoulli")
  truncation <- .check_and_set_truncation(truncation, lower = 0, upper = 1)

  # check individual parameters
  .check_parameter(parameters$probability, "probability")
  .check_parameter_range(parameters$probability, "probability", lower = 0, upper = 1, include_bounds = TRUE)

  # add the values to the output
  output$distribution <- "bernoulli"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple", "prior.discrete")

  return(output)
}
.prior_uniform   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("a", "b"), "uniform")

  # check individual parameters
  .check_parameter(parameters$a, "a")
  .check_parameter(parameters$b, "b")

  if(parameters$a >= parameters$b)
    stop("Parameter 'a' must be lower than the parameter 'b'.")

  # add the values to the output
  output$distribution <- "uniform"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$a, upper = parameters$b)

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_point     <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location"), "point")

  # check individual parameters
  .check_parameter(parameters$location, "location")

  # add the values to the output
  output$distribution <- "point"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$location, upper = parameters$location)

  class(output) <- c("prior", "prior.simple", "prior.point")

  return(output)
}
.prior_mnormal   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("mean", "sd", "K"), "multivariate normal")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$mean, "mean")
  .check_parameter(parameters$sd,   "sd")
  .check_parameter_positive(parameters$sd, "sd")
  .check_parameter_dimensions(parameters$K, "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor

  # add the values to the output
  output$distribution <- "mnormal"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}
.prior_mcauchy   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "K"), "multivariate Cauchy")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_dimensions(parameters$K,   "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor

  # deal with as with a t-distribution
  parameters$df <- 1

  output$distribution <- "mt"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}
.prior_mt        <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "scale", "df", "K"), "multivariate student-t")
  truncation <- .check_and_set_truncation(truncation)

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter(parameters$scale,    "scale")
  .check_parameter(parameters$df,       "df")
  .check_parameter_positive(parameters$scale, "scale")
  .check_parameter_positive(parameters$df,    "df")
  .check_parameter_dimensions(parameters$K,   "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor


  # add the values to the output
  output$distribution <- "mt"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector")

  return(output)
}
.prior_mpoint    <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location", "K"), "multivariate point")

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_parameter_dimensions(parameters$K, "K", allow_NA = TRUE)   # allow undetermined dimensions if called by prior_factor

  # add the values to the output
  output$distribution <- "mpoint"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$location, upper = parameters$location)

  class(output) <- c("prior", "prior.vector", "prior.point")

  return(output)
}
.weightfunction_normalize_side <- function(side){

  side_clean <- .prior_clean_input_name(side)
  if(side_clean %in% c("onesided", "one")){
    return("one-sided")
  }
  if(side_clean %in% c("twosided", "two")){
    return("two-sided")
  }

  stop("'side' must be either 'one-sided' or 'two-sided'.", call. = FALSE)
}
.weightfunction_validate_steps <- function(steps){

  check_real(steps, "steps", check_length = 0)

  if(length(steps) == 0){
    stop("'steps' must contain at least one p-value cut point.", call. = FALSE)
  }
  if(any(steps >= 1) || any(steps <= 0)){
    stop("'steps' must be higher than 0 and lower than 1.", call. = FALSE)
  }
  if(anyDuplicated(steps)){
    stop("'steps' must not contain duplicate cut points.", call. = FALSE)
  }
  if(!all(steps == cummax(steps))){
    stop("'steps' must be monotonically increasing.", call. = FALSE)
  }

  steps
}
.weightfunction_bins <- function(steps){

  cuts <- c(0, steps, 1)
  data.frame(
    lower     = cuts[-length(cuts)],
    upper     = cuts[-1],
    reference = c(TRUE, rep(FALSE, length(cuts) - 2L))
  )
}
.weightfunction_validate_weights <- function(weights, n_bins, reference){

  if(!inherits(weights, "weightfunction_weights")){
    stop("'weights' must be created by wf_cumulative(), wf_fixed(), or wf_independent().", call. = FALSE)
  }

  if(weights$type == "cumulative"){
    if(is.null(weights$alpha)){
      weights$alpha <- rep(1, n_bins)
    }
    check_real(weights$alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = n_bins)

  }else if(weights$type == "fixed"){
    check_real(weights$omega, "omega", lower = 0, upper = 1, check_length = n_bins)
    if(reference == "most_significant" && !isTRUE(all.equal(weights$omega[1], 1))){
      stop("The reference-bin fixed weight must be exactly 1.", call. = FALSE)
    }

  }else if(weights$type == "independent"){
    .check_prior(weights$prior)
    if(!is.prior.simple(weights$prior) || is.prior.point(weights$prior) || is.prior.discrete(weights$prior)){
      stop("'weights$prior' must be a continuous simple prior distribution.", call. = FALSE)
    }
    if(weights$scale == "omega"){
      if(weights$prior$truncation[["lower"]] < 0){
        stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
      }
    }else if(weights$scale != "log_omega"){
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

  }else{
    stop("Unsupported weightfunction weight prior type.", call. = FALSE)
  }

  weights
}
.weightfunction_weights_truncation <- function(weights){

  if(weights$type %in% c("fixed", "cumulative")){
    return(list(lower = 0, upper = 1))
  }

  if(weights$type == "independent"){
    if(weights$scale == "omega"){
      upper <- weights$prior$truncation[["upper"]]
    }else if(weights$scale == "log_omega"){
      upper <- weights$prior$truncation[["upper"]]
      upper <- if(is.infinite(upper)) Inf else exp(upper)
    }else{
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

    return(list(lower = 0, upper = max(1, upper)))
  }

  stop("Unsupported weightfunction weight prior type.", call. = FALSE)
}
.weightfunction_n_bins <- function(prior){

  if(!is.prior.weightfunction(prior)){
    stop("'prior' must be a weightfunction prior.", call. = FALSE)
  }

  nrow(prior$bins)
}
.weightfunction_is_fixed <- function(prior){
  prior$weights$type == "fixed"
}
.weightfunction_local_cuts <- function(prior){

  c(prior$bins$lower[1], prior$bins$upper)
}
.weightfunction_reference_index <- function(prior){
  which(prior$bins$reference)
}
.weightfunction_alpha_marginal <- function(alpha, index){

  if(index <= 1L){
    return(list(type = "point", location = 1))
  }

  list(
    type  = "beta",
    alpha = sum(alpha[index:length(alpha)]),
    beta  = sum(alpha[seq_len(index - 1L)])
  )
}
.weightfunction_rng <- function(prior, n){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    out <- matrix(rep(prior$weights$omega, each = n), nrow = n)

  }else if(prior$weights$type == "cumulative"){
    theta <- extraDistr::rdirichlet(n, alpha = prior$weights$alpha)
    out   <- t(apply(theta[,J:1, drop = FALSE], 1, cumsum))[,J:1, drop = FALSE]

  }else if(prior$weights$type == "independent"){
    out <- matrix(1, nrow = n, ncol = J)
    if(J > 1L){
      draws <- matrix(rng(prior$weights$prior, n * (J - 1L)), nrow = n, ncol = J - 1L)
      if(prior$weights$scale == "log_omega"){
        draws <- exp(draws)
      }
      out[,2:J] <- draws
    }
  }

  colnames(out) <- paste0("omega[", seq_len(J), "]")
  out
}
.weightfunction_marginal_components <- function(prior){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    return(lapply(prior$weights$omega, function(x){
      list(type = "point", location = x)
    }))
  }

  if(prior$weights$type == "cumulative"){
    return(lapply(seq_len(J), function(j){
      .weightfunction_alpha_marginal(prior$weights$alpha, j)
    }))
  }

  if(prior$weights$type == "independent"){
    return(lapply(seq_len(J), function(j){
      if(j == 1L){
        list(type = "point", location = 1)
      }else{
        list(type = "prior", prior = prior$weights$prior, scale = prior$weights$scale)
      }
    }))
  }
}
.prior_weightfunction_component_range <- function(component, quantiles = .005){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      if(component$scale == "omega"){
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) mquant(component$prior, quantiles) else lower
        upper <- if(is.infinite(upper)) mquant(component$prior, 1 - quantiles) else upper
        c(lower, upper)
      }else{
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) 0 else exp(lower)
        upper <- if(is.infinite(upper)) exp(mquant(component$prior, 1 - quantiles)) else exp(upper)
        c(lower, upper)
      }
    }
  )
}
.weightfunction_range <- function(prior, quantiles = .005){

  ranges <- do.call(rbind, lapply(
    .weightfunction_marginal_components(prior),
    .prior_weightfunction_component_range,
    quantiles = quantiles
  ))

  x_range <- range(as.vector(ranges), finite = TRUE)
  if(x_range[1] == x_range[2]){
    x_range <- range(c(0, 1, x_range), finite = TRUE)
  }

  x_range
}
.prior_weightfunction_component_cdf <- function(component, q){

  switch(
    component$type,
    "point" = ppoint(q, location = component$location),
    "beta"  = stats::pbeta(q, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mcdf(component$prior, q)
      }else{
        p <- numeric(length(q))
        p[q <= 0] <- 0
        inside <- q > 0
        p[inside] <- mcdf(component$prior, log(q[inside]))
        p
      }
    }
  )
}
.prior_weightfunction_component_ccdf <- function(component, q){
  if(component$type == "point"){
    return(ppoint(q, location = component$location, lower.tail = FALSE))
  }
  1 - .prior_weightfunction_component_cdf(component, q)
}
.prior_weightfunction_component_lpdf <- function(component, x){

  switch(
    component$type,
    "point" = dpoint(x, location = component$location, log = TRUE),
    "beta"  = stats::dbeta(x, shape1 = component$alpha, shape2 = component$beta, log = TRUE),
    "prior" = {
      if(component$scale == "omega"){
        mlpdf(component$prior, x)
      }else{
        out <- rep(-Inf, length(x))
        inside <- x > 0
        out[inside] <- mlpdf(component$prior, log(x[inside])) - log(x[inside])
        out
      }
    }
  )
}
.prior_weightfunction_component_quant <- function(component, p){

  switch(
    component$type,
    "point" = qpoint(p, location = component$location),
    "beta"  = stats::qbeta(p, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mquant(component$prior, p)
      }else{
        exp(mquant(component$prior, p))
      }
    }
  )
}
.prior_weightfunction_component_mean <- function(component){

  switch(
    component$type,
    "point" = component$location,
    "beta"  = component$alpha / (component$alpha + component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mean(component$prior)
      }else{
        stats::integrate(
          f     = function(x, prior) {
            y <- exp(x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
      }
    }
  )
}
.prior_weightfunction_component_var <- function(component){

  switch(
    component$type,
    "point" = 0,
    "beta"  = (component$alpha * component$beta) /
      ((component$alpha + component$beta)^2 * (component$alpha + component$beta + 1)),
    "prior" = {
      if(component$scale == "omega"){
        var(component$prior)
      }else{
        m1 <- .prior_weightfunction_component_mean(component)
        m2 <- stats::integrate(
          f     = function(x, prior) {
            y <- exp(2 * x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
        m2 - m1^2
      }
    }
  )
}

#### elementary prior related functions ####
#' @title Elementary prior related functions
#'
#' @description Density (pdf / lpdf), distribution
#' function (cdf / ccdf), quantile function (quant),
#' random generation (rng), mean, standard deviation (sd),
#' and marginal variants of the functions (mpdf, mlpf, mcdf,
#' mccdf, mquant) for prior distributions.
#'
#' @param x prior distribution
#' @param y vector of observations
#' @param q vector or matrix of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # generate a random sample from the prior
#' rng(p1, 10)
#'
#' # compute cumulative density function
#' cdf(p1, 0)
#'
#' # obtain quantile
#' quant(p1, .5)
#'
#' # compute probability density
#' pdf(p1, c(0, 1, 2))
#'
#' @return \code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
#' the (marginal) density and the log of (marginal) density,
#' \code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
#' (marginal) distribution and the complement of (marginal) distribution function,
#' \code{quant} (\code{mquant}) give the (marginal) quantile function,
#' and \code{rng} generates random deviates for an object of class 'prior'.
#'
#' @exportS3Method rng prior
#' @exportS3Method cdf prior
#' @exportS3Method ccdf prior
#' @exportS3Method quant prior
#' @exportS3Method lpdf prior
#' @exportS3Method pdf prior
#' @exportS3Method mcdf prior
#' @exportS3Method mccdf prior
#' @exportS3Method mquant prior
#' @exportS3Method mlpdf prior
#' @exportS3Method mpdf prior
#' @name prior_functions
NULL


#### joint distribution functions ####
#' @rdname prior_functions
rng.prior   <- function(x, n, ...){

  prior <- x

  .check_n(n)
  .check_prior(prior)

  dots  <- list(...)
  if(!is.null(dots[["transform_factor_samples"]])){
    check_bool(dots[["transform_factor_samples"]], "transform_factor_samples")
    transform_factor_samples <- dots[["transform_factor_samples"]]
  }else{
    transform_factor_samples <- TRUE
  }
  if(!is.null(dots[["sample_components"]])){
    check_bool(dots[["sample_components"]], "sample_components")
    sample_components <- dots[["sample_components"]]
  }else{
    sample_components <- FALSE
  }

  if(is.prior.spike_and_slab(prior)){

    inclusion <- stats::rbinom(n, size = 1, prob = rng(.get_spike_and_slab_inclusion(prior), n))

    if(sample_components)
      return(inclusion)

    x         <- rng(.get_spike_and_slab_variable(prior), n) * inclusion
    attr(x, "inclusion") <- inclusion

  }else if(is.prior.mixture(prior)){

    component_probabilities <- attr(prior, "prior_weights")
    components              <- sample(seq_along(component_probabilities), size = n, replace = TRUE, prob = component_probabilities)

    if(sample_components)
      return(components)

    if(inherits(prior, "prior.factor_mixture")){

      prior_type <- .get_prior_factor_list_type(prior)

      if(transform_factor_samples){
        x <- matrix(NA, nrow = n, ncol = prior_type[["K"]] + 1)
      }else{
        x <- matrix(NA, nrow = n, ncol = prior_type[["K"]])
      }

      for(component in unique(components)){
        x[component == components,] <- rng(prior[[component]], sum(component == components), transform_factor_samples = transform_factor_samples)
      }

    }else if(inherits(prior, "prior.simple_mixture")){

      x <- rep(NA, n)
      for(component in unique(components)){
        x[component == components] <- rng(prior[[component]], sum(component == components))
      }

    }else{
      stop("unsupported prior mixture type")
    }


    attr(x, "components") <- components

  }else if(is.prior.simple(prior)){

    x <- .prior_simple_rng(prior, n)

  }else if(transform_factor_samples && (is.prior.orthonormal(prior) | is.prior.meandif(prior))){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    par1 <- rep(0, prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
      "mpoint"     = rmpoint(n, location = par1)
    )


    if(is.prior.orthonormal(prior)){
      x <- x %*% t(contr.orthonormal(1:(prior$parameters[["K"]] + 1)))
    }else if(is.prior.meandif(prior)){
      x <- x %*% t(contr.meandif(1:(prior$parameters[["K"]] + 1)))
    }

  }else if(is.prior.vector(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    # TODO: generalize this to priors with covariances
    par1 <- rep(par1, length = prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    x <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::rmvnorm(n, mean = par1, sigma = par2),
      "mt"         = mvtnorm::rmvt(n, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted"),
      "mpoint"     = rmpoint(n, location = par1)
    )

  }else if(is.prior.weightfunction(prior)){

    x <- .weightfunction_rng(prior, n)

  }

  return(x)
}
#' @rdname prior_functions
cdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No cdfs are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No cdfs are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_cdf(prior, q)

  }else if(is.prior.vector(prior)){

    stop("No cdfs are implemented for vector priors.")

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal cdfs are implemented for prior weightfunctions.")

  }

  return(p)
}
#' @rdname prior_functions
ccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No ccdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No ccdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_ccdf(prior, q)

  }else if(is.prior.vector(prior)){

    stop("No cdfs are implemented for vector priors.")

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal ccdf functions are implemented for prior weightfunctions.")

  }

  return(p)
}
#' @rdname prior_functions
lpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No lpdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No lpdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    log_lik <- .prior_simple_lpdf(prior, x)

  }else if(is.prior.vector(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    # TODO: generalize this to priors with covariances
    par1 <- rep(par1, length = prior$parameter[["K"]])

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      par2 <- diag(switch(
        prior[["distribution"]],
        "mnormal" = prior$parameter[["sd"]]^2,
        "mt"      = prior$parameter[["scale"]]^2
      ), ncol = prior$parameter[["K"]], nrow = prior$parameter[["K"]])
    }

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal"    = mvtnorm::dmvnorm(x, mean = par1, sigma = par2, log = TRUE),
      "mt"         = mvtnorm::dmvt(x, delta = par1, sigma = par2, df = prior$parameter[["df"]], type = "shifted", log = TRUE),
      "mpoint"     = dmpoint(x, location = par1, log = TRUE)
    )

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal lpdf are implemented for prior weightfunctions.")

  }

  return(log_lik)
}
#' @rdname prior_functions
pdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.simple(prior)){
    lik <- .prior_simple_pdf(prior, x)
  }else{
    log_lik <- lpdf.prior(prior, x)
    lik     <- exp(log_lik)
  }

  return(lik)
}
#' @rdname prior_functions
quant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No quantile functions are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No quantile functions are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    q <- .prior_simple_quant(prior, p)

  }else if(is.prior.weightfunction(prior)){

    stop("Only marginal quantile functions are implemented for prior weightfunctions.")

  }

  return(q)
}

.prior_simple_base_d <- function(prior, x, log = FALSE){

  switch(
    prior[["distribution"]],
    "normal"    = stats::dnorm(x, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], log = log),
    "lognormal" = stats::dlnorm(x, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], log = log),
    "t"         = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], log = log),
    "gamma"     = stats::dgamma(x, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], log = log),
    "invgamma"  = extraDistr::dinvgamma(x, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], log = log),
    "beta"      = stats::dbeta(x, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], log = log),
    "bernoulli" = stats::dbinom(x, size = 1, prob = prior$parameters[["probability"]], log = log),
    "exp"       = stats::dexp(x, rate = prior$parameters[["rate"]], log = log),
    "uniform"   = stats::dunif(x, min = prior$parameters[["a"]], max = prior$parameters[["b"]], log = log),
    "point"     = dpoint(x, location = prior$parameters[["location"]], log = log)
  )
}

.prior_simple_base_p <- function(prior, q, lower.tail = TRUE){

  switch(
    prior[["distribution"]],
    "normal"    = stats::pnorm(q, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = lower.tail, log.p = FALSE),
    "lognormal" = stats::plnorm(q, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = lower.tail, log.p = FALSE),
    "t"         = extraDistr::plst(q, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "gamma"     = stats::pgamma(q, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "invgamma"  = extraDistr::pinvgamma(q, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = lower.tail, log.p = FALSE),
    "beta"      = stats::pbeta(q, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = lower.tail, log.p = FALSE),
    "bernoulli" = stats::pbinom(q, size = 1, prob = prior$parameters[["probability"]], lower.tail = lower.tail, log.p = FALSE),
    "exp"       = stats::pexp(q, rate = prior$parameters[["rate"]], lower.tail = lower.tail, log.p = FALSE),
    "uniform"   = stats::punif(q, min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = lower.tail, log.p = FALSE),
    "point"     = ppoint(q, location = prior$parameters[["location"]], lower.tail = lower.tail, log.p = FALSE)
  )
}

.prior_simple_base_q <- function(prior, p){

  switch(
    prior[["distribution"]],
    "normal"    = stats::qnorm(p, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
    "lognormal" = stats::qlnorm(p, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
    "t"         = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
    "gamma"     = stats::qgamma(p, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
    "invgamma"  = extraDistr::qinvgamma(p, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
    "beta"      = stats::qbeta(p, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
    "bernoulli" = stats::qbinom(p, size = 1, prob = prior$parameters[["probability"]], lower.tail = TRUE, log.p = FALSE),
    "exp"       = stats::qexp(p, rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
    "uniform"   = stats::qunif(p, min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
    "point"     = qpoint(p, location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
  )
}

.prior_simple_base_r <- function(prior, n){

  switch(
    prior[["distribution"]],
    "normal"    = stats::rnorm(n, mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]]),
    "lognormal" = stats::rlnorm(n, meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]]),
    "t"         = extraDistr::rlst(n, df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]]),
    "gamma"     = stats::rgamma(n, shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]]),
    "invgamma"  = extraDistr::rinvgamma(n, alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]]),
    "beta"      = stats::rbeta(n, shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]]),
    "bernoulli" = stats::rbinom(n, size = 1, prob = prior$parameters[["probability"]]),
    "exp"       = stats::rexp(n, rate = prior$parameters[["rate"]]),
    "uniform"   = stats::runif(n, min = prior$parameters[["a"]], max = prior$parameters[["b"]]),
    "point"     = rpoint(n, location = prior$parameters[["location"]])
  )
}

.prior_simple_cdf <- function(prior, q){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_p(prior, q, lower.tail = TRUE))
  }

  p        <- numeric(length(q))
  q_lower  <- q < prior$truncation[["lower"]]
  q_higher <- q > prior$truncation[["upper"]]
  q_inside <- !q_lower & !q_higher

  p[q_lower]  <- 0
  p[q_higher] <- 1

  if(any(q_inside)){
    p[q_inside] <- .prior_simple_base_p(prior, q[q_inside], lower.tail = TRUE)

    if(prior[["distribution"]] != "point"){
      C1          <- .prior_C1(prior)
      p[q_inside] <- (p[q_inside] - C1) / (.prior_C2(prior) - C1)
    }
  }

  return(p)
}

.prior_simple_ccdf <- function(prior, q){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_p(prior, q, lower.tail = FALSE))
  }

  p        <- numeric(length(q))
  q_lower  <- q < prior$truncation[["lower"]]
  q_higher <- q > prior$truncation[["upper"]]
  q_inside <- !q_lower & !q_higher

  p[q_lower]  <- 1
  p[q_higher] <- 0

  if(any(q_inside)){
    p[q_inside] <- .prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE)

    if(prior[["distribution"]] != "point"){
      C2          <- .prior_C2(prior)
      p[q_inside] <- (p[q_inside] - (1 - C2)) / (C2 - .prior_C1(prior))
    }
  }

  return(p)
}

.prior_simple_lpdf <- function(prior, x){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_d(prior, x, log = TRUE))
  }

  log_lik <- .prior_simple_base_d(prior, x, log = TRUE)
  log_lik[x < prior$truncation[["lower"]] | x > prior$truncation[["upper"]]] <- -Inf

  if(prior[["distribution"]] != "point"){
    log_lik <- log_lik - log(.prior_C(prior))
  }

  return(log_lik)
}

.prior_simple_pdf <- function(prior, x){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_d(prior, x, log = FALSE))
  }

  lik <- .prior_simple_base_d(prior, x, log = FALSE)
  lik[x < prior$truncation[["lower"]] | x > prior$truncation[["upper"]]] <- 0

  if(prior[["distribution"]] != "point"){
    lik <- lik / .prior_C(prior)
  }

  return(lik)
}

.prior_simple_quant <- function(prior, p){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_q(prior, p))
  }

  if(is.prior.discrete(prior)){
    return(.prior_simple_quant_optim(prior, p))
  }

  C1 <- .prior_C1(prior)
  .prior_simple_base_q(prior, C1 + p * (.prior_C2(prior) - C1))
}

.prior_simple_quant_optim <- function(prior, p){

  if(!is.infinite(prior$truncation[["lower"]]) & !is.infinite(prior$truncation[["upper"]])){
    start_value <- prior$truncation[["lower"]] + (prior$truncation[["upper"]]  - prior$truncation[["lower"]]) / 2
  }else if(!is.infinite(prior$truncation[["upper"]])){
    start_value <- prior$truncation[["upper"]] - 1
  }else if(!is.infinite(prior$truncation[["lower"]])){
    start_value <- prior$truncation[["lower"]] + 1
  }else{
    start_value <- 0
  }

  sapply(p, function(p_i){
    stats::optim(
      par     = start_value,
      fn      = function(x, prior, p_i)(.prior_simple_cdf(prior, x) - p_i)^2,
      lower   = prior$truncation[["lower"]],
      upper   = prior$truncation[["upper"]],
      prior   = prior,
      p_i     = p_i,
      method  = "L-BFGS-B",
      control = list(
        factr = 1e3
      )
    )$par
  })
}

.prior_simple_rng <- function(prior, n){

  if(.is_prior_default_range(prior)){
    return(.prior_simple_base_r(prior, n))
  }

  if(is.prior.discrete(prior)){
    return(.prior_simple_rng_rejection(prior, n))
  }

  C1 <- .prior_C1(prior)
  .prior_simple_base_q(prior, stats::runif(n, min = C1, max = .prior_C2(prior)))
}

.prior_simple_rng_rejection <- function(prior, n){

  x  <- NULL
  nn <- round(n * 1 / .prior_C(prior) * 1.10)

  while(length(x) < n){
    temp_x <- .prior_simple_base_r(prior, nn)
    x      <- c(x, temp_x[temp_x >= prior$truncation[["lower"]] & temp_x <= prior$truncation[["upper"]]])
  }

  x[1:n]
}

.prior_C1 <- function(prior){

  if(is.prior.simple(prior)){

    C1 <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(prior$truncation[["lower"]], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(prior$truncation[["lower"]], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(prior$truncation[["lower"]], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(prior$truncation[["lower"]], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(prior$truncation[["lower"]], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["lower"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "bernoulli" = stats::pbinom(prior$truncation[["lower"]] - 1, size = 1, prob = prior$parameters[["probability"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["lower"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["lower"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "point"     = ppoint(prior$truncation[["lower"]], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

  }

  return(C1)
}
.prior_C2 <- function(prior){

  if(is.prior.simple(prior)){

    C2 <- switch(
      prior[["distribution"]],
      "normal"    = stats::pnorm(prior$truncation[["upper"]], mean = prior$parameters[["mean"]], sd = prior$parameters[["sd"]], lower.tail = TRUE, log.p = FALSE),
      "lognormal" = stats::plnorm(prior$truncation[["upper"]], meanlog = prior$parameters[["meanlog"]], sdlog = prior$parameters[["sdlog"]], lower.tail = TRUE, log.p = FALSE),
      "t"         = extraDistr::plst(prior$truncation[["upper"]], df = prior$parameters[["df"]], mu = prior$parameters[["location"]], sigma = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "gamma"     = stats::pgamma(prior$truncation[["upper"]], shape = prior$parameters[["shape"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "invgamma"  = extraDistr::pinvgamma(prior$truncation[["upper"]], alpha = prior$parameters[["shape"]], beta = prior$parameters[["scale"]], lower.tail = TRUE, log.p = FALSE),
      "beta"      = stats::pbeta(prior$truncation[["upper"]], shape1 = prior$parameters[["alpha"]], shape2 = prior$parameters[["beta"]], lower.tail = TRUE, log.p = FALSE),
      "bernoulli" = stats::pbinom(prior$truncation[["upper"]] + 1, size = 1, prob = prior$parameters[["probability"]], lower.tail = TRUE, log.p = FALSE),
      "exp"       = stats::pexp(prior$truncation[["upper"]], rate = prior$parameters[["rate"]], lower.tail = TRUE, log.p = FALSE),
      "uniform"   = stats::punif(prior$truncation[["upper"]], min = prior$parameters[["a"]], max = prior$parameters[["b"]], lower.tail = TRUE, log.p = FALSE),
      "point"     = ppoint(prior$truncation[["upper"]], location = prior$parameters[["location"]], lower.tail = TRUE, log.p = FALSE)
    )

  }

  return(C2)
}
.prior_C  <- function(prior){

  if(is.prior.simple(prior)){

    C <- .prior_C2(prior) - .prior_C1(prior)

  }

  return(C)
}

# tools
.is_prior_default_range   <- function(prior){
  default_range <- switch(
    prior[["distribution"]],
    "normal"    = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "lognormal" = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "t"         = is.infinite(prior$truncation[["lower"]])          & is.infinite(prior$truncation[["upper"]]),
    "gamma"     = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "invgamma"  = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "beta"      = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & isTRUE(all.equal(prior$truncation[["upper"]], 1)),
    "bernoulli" = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & isTRUE(all.equal(prior$truncation[["upper"]], 1)),
    "exp"       = isTRUE(all.equal(prior$truncation[["lower"]], 0)) & is.infinite(prior$truncation[["upper"]]),
    "uniform"   = TRUE,
    "point"     = TRUE,
    "mpoint"          = TRUE,
    "weightfunction"  = TRUE,
    "none"            = TRUE
  )
}

#### marginal distribution functions ####
# weightfunctions require just marginals for plotting and etc
# (also, the joint are not implemented in general, since they differ between JAGS/Stan implementations)
#' @rdname prior_functions
mcdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mcdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mcdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_cdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    p <- do.call(cbind, lapply(components, .prior_weightfunction_component_cdf, q = q))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2),
      "mpoint"  = ppoint(q, location = 0)
    )

  }

  return(p)
}
#' @rdname prior_functions
mccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(q)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mccdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mccdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_ccdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    p <- do.call(cbind, lapply(components, .prior_weightfunction_component_ccdf, q = q))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2, lower.tail = FALSE),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2, lower.tail = FALSE),
      "mpoint"  = ppoint(q, location = 0, lower.tail = FALSE),
    )

  }

  return(p)
}
#' @rdname prior_functions
mlpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mlpdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mlpdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    log_lik <- .prior_simple_lpdf(prior, x)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    log_lik <- do.call(cbind, lapply(components, .prior_weightfunction_component_lpdf, x = x))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal" = stats::dnorm(x, mean = 0, sd = par2, log = TRUE),
      "mt"      = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = 0, sigma = par2, log = TRUE),
      "mpoint"  = dpoint(x, location = 0, log = TRUE),
    )

  }

  return(log_lik)
}
#' @rdname prior_functions
mpdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(x)
  .check_prior(prior)

  if(is.prior.simple(prior)){
    lik <- .prior_simple_pdf(prior, x)
  }else{
    log_lik <- mlpdf.prior(prior, x)
    lik     <- exp(log_lik)
  }

  return(lik)
}
#' @rdname prior_functions
mquant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No quantile functions are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No quantile functions are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    q <- .prior_simple_quant(prior, p)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    q <- do.call(cbind, lapply(components, .prior_weightfunction_component_quant, p = p))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    q <- switch(
      prior[["distribution"]],
      "mnormal" = stats::qnorm(p, mean = 0, sd = par2),
      "mt"      = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = 0, sigma = par2),
      "mpoint"  = qpoint(p, location = 0)
    )

  }

  return(q)
}

### create generic method
#' @title Creates generics for common statistical functions
#'
#' @description Density (pdf / lpdf), distribution
#' function (cdf / ccdf), quantile function (quant),
#' random generation (rng), mean, standard deviation (sd),
#' and marginal variants of the functions (mpdf, mlpf, mcdf,
#' mccdf, mquant).
#'
#' @param x main argument
#' @param ... unused arguments
#'
#' @return \code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
#' the (marginal) density and the log of (marginal) density,
#' \code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
#' (marginal) distribution and the complement of (marginal) distribution function,
#' \code{quant} (\code{mquant}) give the (marginal) quantile function,
#' and \code{rng} generates random deviates for an object of class 'prior'.
#'
#' The \code{pdf} function proceeds to PDF graphics device if \code{x} is a character.
#'
#' @export rng
#' @export cdf
#' @export ccdf
#' @export quant
#' @export lpdf
#' @export pdf
#' @export mcdf
#' @export mccdf
#' @export mquant
#' @export mlpdf
#' @export mpdf
#' @importFrom grDevices pdf
#' @name prior_functions_methods
NULL

#' @rdname prior_functions_methods
rng    <- function(x, ...){
  UseMethod("rng")
}
#' @rdname prior_functions_methods
cdf    <- function(x, ...){
  UseMethod("cdf")
}
#' @rdname prior_functions_methods
ccdf   <- function(x, ...){
  UseMethod("ccdf")
}
#' @rdname prior_functions_methods
quant  <- function(x, ...){
  UseMethod("quant")
}
#' @rdname prior_functions_methods
lpdf   <- function(x, ...){
  UseMethod("lpdf")
}
#' @rdname prior_functions_methods
pdf    <- function(x, ...){
  UseMethod("pdf")
}
#' @rdname prior_functions_methods
mcdf   <- function(x, ...){
  UseMethod("mcdf")
}
#' @rdname prior_functions_methods
mccdf  <- function(x, ...){
  UseMethod("mccdf")
}
#' @rdname prior_functions_methods
mquant <- function(x, ...){
  UseMethod("mquant")
}
#' @rdname prior_functions_methods
mlpdf  <- function(x, ...){
  UseMethod("mlpdf")
}
#' @rdname prior_functions_methods
mpdf   <- function(x, ...){
  UseMethod("mpdf")
}

#' @export
pdf.default  <- function(x, ...){
  grDevices::pdf(x, ...)
}


#### additional prior related function ####
#' @title Prior mean
#'
#' @description Computes mean of a prior
#' distribution. (In case of orthonormal prior distributions
#' for factors, the mean of for the deviations from intercept
#' is returned.)
#'
#' @param x a prior
#' @param ... unused
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute mean of the prior distribution
#' mean(p1)
#'
#' @return a mean of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @rdname mean.prior
#' @export
mean.prior   <- function(x, ...){

  .check_prior(x, "x")

  if(is.prior.spike_and_slab(x)){

    m <- mean(.get_spike_and_slab_variable(x)) * mean(.get_spike_and_slab_inclusion(x))

  }else if(is.prior.mixture(x)){

    stop("No mean is implemented for prior mixtures.")

  }else if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      m <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["mean"]],
        "lognormal" = exp(x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2/2),
        "t"         = ifelse(x$parameters[["df"]] > 1, x$parameters[["location"]], NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]],
        "invgamma"  = ifelse(x$parameters[["shape"]] > 1, x$parameters[["scale"]]/(x$parameters[["shape"]] - 1), NaN),
        "beta"      = x$parameters[["alpha"]] / (x$parameters[["alpha"]] + x$parameters[["beta"]]),
        "bernoulli" = x$parameters[["probability"]],
        "exp"       = 1 / x$parameters[["rate"]],
        "uniform"   = (x$parameters[["a"]] + x$parameters[["b"]]) / 2,
        "point"     = x$parameters[["location"]]
      )

    }else{

      # check for undefined values
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 1){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 1){
          return(NaN)
        }
      }

      m <- stats::integrate(
        f       = function(x, prior) x * pdf(prior, x),
        lower   = x$truncation[["lower"]],
        upper   = x$truncation[["upper"]],
        prior   = x
      )$value

    }


  }else if(is.prior.weightfunction(x)){

    m <- .mean.weightfunction(x)

  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){

    par1 <- switch(
      x[["distribution"]],
      "mnormal" = x$parameter[["mean"]],
      "mt"      = x$parameter[["location"]],
      "mpoint"  = x$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(x[["distribution"]] == "mt" && x$parameters[["df"]] <= 1){
      m <- NaN
    }else{
      # orthonormal prior distributions must be centered
      m <- 0
    }

  }

  return(m)
}


### create generic methods from sd and var
#' @title Creates generic for sd function
#'
#' @param x main argument
#' @param ... additional arguments
#'
#' @return \code{sd} returns a standard deviation
#' of the supplied object (if it is either a numeric vector
#' or an object of class 'prior').
#'
#' @seealso \link[stats]{sd}
#' @export
sd  <- function(x, ...){
  UseMethod("sd")
}

#' @title Creates generic for var function
#'
#' @param x main argument
#' @param ... additional arguments
#'
#' @return \code{var} returns a variance
#' of the supplied object (if it is either a numeric vector
#' or an object of class 'prior').
#'
#' @seealso \link[stats]{cor}
#' @export
var <- function(x, ...){
  UseMethod("var")
}

#' @export
sd.default  <- function(x, ...){
  stats::sd(x, ...)
}
#' @export
var.default <- function(x, ...){
  stats::var(x, ...)
}


#' @title Prior var
#'
#' @description Computes variance
#' of a prior distribution.
#'
#' @param x a prior
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute variance of the prior distribution
#' var(p1)
#'
#' @return a variance of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @importFrom stats var
#' @rdname var.prior
#' @export
var.prior   <- function(x, ...){

  .check_prior(x, "x")

  if(is.prior.spike_and_slab(x)){

    inclusion_prior <- .get_spike_and_slab_inclusion(x)
    variable_prior <- .get_spike_and_slab_variable(x)
    
    # Handle different inclusion prior types
    if(inclusion_prior[["distribution"]] == "beta") {
      # the inclusion is beta -> indicators are betabinom
      var_inclusion <- with(inclusion_prior[["parameters"]], (alpha * beta * (alpha + beta + 1) ) / ( (alpha + beta)^2 * (alpha + beta + 1)  ) )
    } else if(inclusion_prior[["distribution"]] == "point") {
      # the inclusion is a spike (point) -> no variance
      var_inclusion <- 0
    } else {
      # for other distributions, compute variance directly
      var_inclusion <- var(inclusion_prior)
    }
    
    var <-
      (var(variable_prior)  + mean(variable_prior)^2) *
      (var_inclusion         + mean(inclusion_prior)^2) -
      (mean(variable_prior)^2 * mean(inclusion_prior)^2)

  }else if(is.prior.mixture(x)){

    stop("No var is implemented for prior mixtures.")

  }else if(is.prior.simple(x)){

    if(.is_prior_default_range(x)){

      var <- switch(
        x[["distribution"]],
        "normal"    = x$parameters[["sd"]]^2,
        "lognormal" = (exp(x$parameters[["sdlog"]]^2) - 1) * exp(2 * x$parameters[["meanlog"]] + x$parameters[["sdlog"]]^2),
        "t"         = ifelse(x$parameters[["df"]] > 2, x$parameters[["scale"]]^2 * x$parameters[["df"]] / (x$parameters[["df"]] - 2), NaN),
        "gamma"     = x$parameters[["shape"]] / x$parameters[["rate"]]^2,
        "invgamma"  = ifelse(x$parameters[["shape"]] > 2, x$parameters[["scale"]]^2 / (x$parameters[["shape"]] - 1)^2 * (x$parameters[["shape"]] - 2), NaN),
        "beta"      = (x$parameters[["alpha"]] * x$parameters[["beta"]]) / ((x$parameters[["alpha"]] + x$parameters[["beta"]])^2 * (x$parameters[["alpha"]] + x$parameters[["beta"]] + 1)),
        "bernoulli" = (x$parameters[["probability"]] * (1 - x$parameters[["probability"]]) ),
        "exp"       = 1 / x$parameters[["rate"]]^2,
        "uniform"   = (x$parameters[["b"]] - x$parameters[["a"]])^2 / 12,
        "point"     = 0
      )

    }else{

      # check for undefined values
      if(x[["distribution"]] == "t"){
        if(x$parameters[["df"]] <= 2){
          return(NaN)
        }
      }
      if(x[["distribution"]] == "invgamma"){
        if(x$parameters[["shape"]] <= 2){
          return(NaN)
        }
      }

      E2 <- stats::integrate(
        f       = function(x, prior) x^2 * pdf(prior, x),
        lower   = x$truncation[["lower"]],
        upper   = x$truncation[["upper"]],
        prior   = x
      )$value

      var <- E2 - mean(x)^2
    }


  }else if(is.prior.weightfunction(x)){

    var <- .var.weightfunction(x)

  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){

    par1 <- switch(
      x[["distribution"]],
      "mnormal" = x$parameter[["mean"]],
      "mt"      = x$parameter[["location"]],
      "mpoint"  = x$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(x[["distribution"]] == "mpoint"){
      var <- 0
    }else if(x[["distribution"]] == "mt" && x$parameters[["df"]] <= 2){
      var <- NaN
    }else{


      if(is.na(x$parameters[["K"]]) && !is.null(attr(x, "levels"))){
        x$parameters[["K"]] <- .get_prior_factor_levels(x)
      }else if(is.na(x$parameters[["K"]])){
        x$parameters[["K"]] <- 1
        warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
      }


      if(is.prior.orthonormal(x)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(x$parameters[["K"]] + 1))[1,] * switch(
          x[["distribution"]],
          "mnormal" = x$parameters[["sd"]],
          "mt"      = x$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(x)){
        par2 <- sqrt(sum( (contr.meandif(1:(x$parameters[["K"]] + 1))[1,] * switch(
          x[["distribution"]],
          "mnormal" = x$parameters[["sd"]],
          "mt"      = x$parameters[["scale"]]) )^2 ))
      }

      # use the univariate functions
      var <- switch(
        x[["distribution"]],
        "mnormal" = var.prior(prior("normal", parameters = list(mean = 0, sd = par2))),
        "mt"      = var.prior(prior("t",      parameters = list(location = 0, scale = par2, df = x$parameters[["df"]]))))
    }

  }

  return(var)
}

#' @title Prior sd
#'
#' @description Computes standard deviation
#' of a prior distribution.
#'
#' @param x a prior
#' @param ... unused arguments
#'
#' @examples
#' # create a standard normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # compute sd of the prior distribution
#' sd(p1)
#'
#' @return a standard deviation of an object of class 'prior'.
#'
#' @seealso [prior()]
#' @importFrom stats sd
#' @rdname sd.prior
#' @export
sd.prior     <- function(x, ...){

  .check_prior(x, "x")

  sd <- sqrt(var(x))

  return(sd)
}

.mean.weightfunction <- function(prior){

  components <- .weightfunction_marginal_components(prior)
  m <- vapply(components, .prior_weightfunction_component_mean, numeric(1))

  return(m)
}
.var.weightfunction  <- function(prior){

  components <- .weightfunction_marginal_components(prior)
  m <- vapply(components, .prior_weightfunction_component_var, numeric(1))

  return(m)
}
