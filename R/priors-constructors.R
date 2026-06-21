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
  bernoulli_support <- c(0, 1)
  bernoulli_keep    <- bernoulli_support >= truncation[["lower"]] & bernoulli_support <= truncation[["upper"]]
  if(!any(bernoulli_keep))
    stop("Bernoulli truncation must contain at least one support point: 0 or 1.", call. = FALSE)
  if(!is.expression(parameters$probability)){
    bernoulli_prob <- c(1 - parameters$probability, parameters$probability)
    if(sum(bernoulli_prob[bernoulli_keep]) <= 0)
      stop("Bernoulli truncation must retain positive probability mass.", call. = FALSE)
  }

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

  truncation <- .check_and_set_truncation(truncation, lower = parameters$a, upper = parameters$b)

  # add the values to the output
  output$distribution <- "uniform"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_point     <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("location"), "point")

  # check individual parameters
  .check_parameter(parameters$location, "location")
  .check_point_truncation(truncation, parameters$location, "point")

  # add the values to the output
  output$distribution <- "point"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$location, upper = parameters$location)

  class(output) <- c("prior", "prior.simple", "prior.point")

  return(output)
}
.prior_moment    <- function(parameters, truncation){

  output <- list()

  parameters <- .nonlocal_parameters_moment(parameters)
  truncation <- .check_and_set_truncation(truncation)

  output$distribution <- "moment"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_invmoment <- function(parameters, truncation){

  output <- list()

  parameters <- .nonlocal_parameters_invmoment(parameters)
  truncation <- .check_and_set_truncation(truncation)

  output$distribution <- "invmoment"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.simple")

  return(output)
}
.prior_mnormal   <- function(parameters, truncation){

  output <- list()

  # check overall settings
  parameters <- .check_and_name_parameters(parameters, c("mean", "sd", "K"), "multivariate normal")
  truncation <- .check_and_set_truncation(truncation)
  .check_vector_truncation_unsupported(truncation)

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
  .check_vector_truncation_unsupported(truncation)

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
  .check_vector_truncation_unsupported(truncation)

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
  .check_point_truncation(truncation, parameters$location, "multivariate point")

  # add the values to the output
  output$distribution <- "mpoint"
  output$parameters   <- parameters
  output$truncation   <- list(lower = parameters$location, upper = parameters$location)

  class(output) <- c("prior", "prior.vector", "prior.point")

  return(output)
}
.prior_dirichlet <- function(parameters, truncation){

  output <- list()

  if(!is.null(names(parameters))){
    names(parameters)[names(parameters) == "concentration"] <- "alpha"
  }
  parameters <- .check_and_name_parameters(parameters, "alpha", "Dirichlet")
  truncation <- .check_and_set_truncation(truncation)
  .check_vector_truncation_unsupported(truncation)

  .check_parameter(parameters$alpha, "alpha", length = 0)
  .check_parameter_positive(parameters$alpha, "alpha")
  if(is.numeric(parameters$alpha) && any(!is.finite(parameters$alpha))){
    stop("The 'alpha' concentration parameters must be finite.", call. = FALSE)
  }
  if(length(parameters$alpha) < 2L){
    stop("The Dirichlet 'alpha' concentration vector must contain at least two values.", call. = FALSE)
  }

  parameters$K <- length(parameters$alpha)

  output$distribution <- "dirichlet"
  output$parameters   <- parameters
  output$truncation   <- truncation

  class(output) <- c("prior", "prior.vector", "prior.simplex")

  return(output)
}
.check_point_truncation <- function(truncation, location, distribution){

  if(is.expression(location)){
    return(invisible(NULL))
  }

  if(length(truncation) > 2)
    stop("More than two truncation points were supplied.", call. = FALSE)

  if(!is.null(names(truncation))){
    if(!all(names(truncation) %in% c("lower", "upper")))
      stop("Truncation points must be named 'lower' and 'upper'.", call. = FALSE)
    if(length(truncation) == 1){
      if(names(truncation) == "lower"){
        truncation$upper <- Inf
      }else if(names(truncation) == "upper"){
        truncation$lower <- -Inf
      }
    }
  }else{
    if(length(truncation) == 2){
      names(truncation) <- c("lower", "upper")
    }else if(length(truncation) == 1){
      names(truncation) <- "lower"
      truncation$upper <- Inf
    }else{
      truncation <- list(lower = -Inf, upper = Inf)
    }
  }

  truncation <- truncation[c("lower", "upper")]

  if(!is.numeric(truncation[["lower"]]) || length(truncation[["lower"]]) != 1 ||
     !is.numeric(truncation[["upper"]]) || length(truncation[["upper"]]) != 1){
    stop("Truncation points must be numeric values.", call. = FALSE)
  }

  if(truncation[["lower"]] > truncation[["upper"]]){
    stop("The lower truncation point must be lower or equal to the upper truncation point.", call. = FALSE)
  }

  if(truncation[["lower"]] > location || truncation[["upper"]] < location){
    stop(paste0("The ", distribution, " prior truncation must contain the point location."), call. = FALSE)
  }

  invisible(NULL)
}
.check_vector_truncation_unsupported <- function(truncation){

  if(!is.infinite(truncation[["lower"]]) || !is.infinite(truncation[["upper"]])){
    stop("Vector priors do not support truncation.", call. = FALSE)
  }
}
