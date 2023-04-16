# functions for checking input into prior distributions related functions
.check_and_name_parameters <- function(parameters, names, distribution){

  if(length(parameters) != length(names))
    stop(paste0(distribution, " prior distribution requires ", length(names), " parameters."), call. = FALSE)
  if(!is.null(names(parameters))){
    if(!all(names(parameters) %in% names))
      stop(paste0("Parameters ", paste(paste0("'", names(parameters)[!names(parameters) %in% c(names, "")], "'"), sep = ", ", collapse = ""), " are not supported for a ", distribution," distribution."), call. = FALSE)
  }else{
    names(parameters) <- names
  }

  # order
  parameters <- parameters[names]

  return(parameters)
}
.check_and_set_truncation  <- function(truncation, lower = -Inf, upper = Inf){

  if(length(truncation) > 2)
    stop("More than two truncation points were supplied.")

  if(!is.null(names(truncation))){
    if(!all(names(truncation) %in% c("lower", "upper")))
      stop("Truncation points must be named 'lower' and 'upper'.")
    if(length(truncation) == 1){
      if(names(truncation) == "lower"){
        truncation$upper  <- upper
      }else if(names(truncation) == "upper"){
        truncation$lower  <- lower
      }
    }
  }else{
    if(length(truncation) == 2){
      names(truncation) <- c("lower", "upper")
    }else if(length(truncation) == 1){
      names(truncation) <- "lower"
      truncation$upper  <- upper
    }else{
      truncation <- list(
        lower = lower,
        upper = upper
      )
    }
  }


  if(lower != -Inf){
    # change the default value to a distribution specific value or throw an error of misspecified by user
    if(truncation$lower == -Inf){
      truncation$lower <- lower
    }else if(truncation$lower < lower)
      stop(paste0("Lower truncation point must be larger or equal to ", lower, "."))
  }

  if(upper != Inf){
    # change the default value to a distribution specific value or throw an error of misspecified by user
    if(truncation$upper == Inf){
      truncation$upper <- upper
    }else if(truncation$upper > upper)
      stop(paste0("Upper truncation point must be smaller or equal to ", upper, "."))
  }

  if(truncation$lower >= truncation$upper)
    stop("The lower truncation point must be lower than the upper truncation points.")

  # order
  truncation <- truncation[c("lower", "upper")]

  return(truncation)
}
.check_parameter           <- function(parameter, name, length = 1){
  if(length == -1){
    if(!is.numeric(parameter) | !is.vector(parameter) | length(parameter) > 1)
      stop(paste0("The '", name, "' must be a numeric vector of length at least 2."), call. = FALSE)
  }else if(length == 0){
    if(!is.numeric(parameter) | !is.vector(parameter))
      stop(paste0("The '", name, "' must be a numeric vector."), call. = FALSE)
  }else{
    if(!is.numeric(parameter) | !is.vector(parameter) | length(parameter) != length)
      stop(paste0("The '", name, "' must be a numeric vector of length ", length, "."), call. = FALSE)
  }
}
.check_parameter_positive  <- function(parameter, name, include_zero = FALSE){
  if(include_zero){
    if(any(parameter < 0))
      stop(paste0("The '", name, "' must be non-negative."), call. = FALSE)
  }else{
    if(any(parameter <= 0))
      stop(paste0("The '", name, "' must be positive."), call. = FALSE)
  }
}
.check_parameter_negative  <- function(parameter, name, include_zero = FALSE){
  if(include_zero){
    if(any(parameter > 0))
      stop(paste0("The '", name, "' must be non-positive."), call. = FALSE)
  }else{
    if(any(parameter >= 0))
      stop(paste0("The '", name, "' must be negative."), call. = FALSE)
  }
}
.check_parameter_dimensions<- function(parameter, name, allow_NA = FALSE){
  if(!allow_NA && is.na(parameter)){
    stop(paste0("The '", name, "' must be defined."), call. = FALSE)
  }else if(!allow_NA){
    if(length(parameter) != 1 || !is.numeric(parameter))
      stop(paste0("The '", name, "' must be a numeric vector of length 1."), call. = FALSE)
    if(!.is.wholenumber(parameter))
      stop(paste0(call, "The '", name ,"' must be an integer"), call. = FALSE)
  }
}
.check_parameter_weigthfunction <- function(steps, alpha = NULL, alpha1 = NULL, alpha2 = NULL, omega = NULL){

  if(!is.null(alpha)){

    .check_parameter(steps, "steps", length = 0)
    .check_parameter(alpha, "alpha", length = 0)

    if(any(steps >= 1) | any(steps <= 0))
      stop("Parameter 'steps' must be higher than 0 and lower than 1.", call. = FALSE)
    if(!(all(steps == cummax(steps))))
      stop("Parameters 'steps' must be monotonically increasing.", call. = FALSE)
    if(length(steps) != length(alpha) - 1)
      stop("The parameter alpha needs to have one more argument then there are steps.", call. = FALSE)

    .check_parameter_positive(alpha, "alpha")

  }else if(!is.null(alpha1) & !is.null(alpha2)){

    .check_parameter(steps,  "steps",  length = 0)
    .check_parameter(alpha1, "alpha1", length = 0)
    .check_parameter(alpha2, "alpha2", length = 0)

    if(any(steps >= 1) | any(steps <= 0))
      stop("Parameter 'steps' must be higher than 0 and lower than 1.", call. = FALSE)
    if(!(all(steps == cummax(steps))))
      stop("Parameters 'steps' must be monotonically increasing.", call. = FALSE)
    if(sum(steps <= .5) != length(alpha1) - 1)
      stop("The parameter alpha1 needs to have one more argument then there are steps <= .5.", call. = FALSE)
    if(sum(steps > .5)  != length(alpha2) - 1)
      stop("The parameter alpha2 needs to have one more argument then there are steps > .5.", call. = FALSE)

    .check_parameter_positive(alpha1, "alpha1")
    .check_parameter_positive(alpha2, "alpha2")

  }else if(!is.null(omega)){

    .check_parameter(steps, "steps", length = 0)
    .check_parameter(omega, "omega", length = 0)

    if(any(steps >= 1) | any(steps <= 0))
      stop("Parameter 'steps' must be higher than 0 and lower than 1.", call. = FALSE)
    if(!(all(steps == cummax(steps))))
      stop("Parameters 'steps' must be monotonically increasing.", call. = FALSE)
    if(length(steps) != length(omega) - 1)
      stop("The parameter omega needs to have one more argument then there are steps.", call. = FALSE)
    if(any(omega < 0) | any(omega > 1))
      stop("Parameter 'omega' must be within 0 to 1 range.", call. = FALSE)
    if(!any(sapply(omega, function(omega_i)isTRUE(all.equal(omega_i, 1)))))
      stop("At least one 'omega' parameter must be equal to 1.", call. = FALSE)

  }
}
.check_parameter_range          <- function(parameter, name, lower, upper, include_bounds = FALSE){
  if(include_bounds){
    if(any(parameter < lower) | any(parameter > upper))
      stop(paste0("The '", name, "' must be higher than ", lower, " and lower than ", upper, "."), call. = FALSE)
  }else{
    if(any(parameter <= lower) | any(parameter >= upper))
      stop(paste0("The '", name, "' must be higher or equal to than ", lower, " and lower or equal to than ", upper, "."), call. = FALSE)
  }
}

.get_prior_factor_levels      <- function(prior){
  if(is.prior.independent(prior)){
    return(attr(prior, "levels"))
  }else if(is.prior.dummy(prior)){
    return(attr(prior, "levels") - 1)
  }else if(is.prior.orthonormal(prior)){
    return(attr(prior, "levels") - 1)
  }else if(is.prior.meandif(prior)){
    return(attr(prior, "levels") - 1)
  }else if(is.prior.point(prior)){
    # allow to deal with spike priors assuming there is an intercept & the prior is not independent
    return(attr(prior, "levels") - 1)
  }
}
.get_prior_factor_level_names <- function(prior){
  if(is.null(attr(prior, "level_names"))){
    return(1:.get_prior_factor_levels(prior))
  }else{
    return(attr(prior, "level_names"))
  }
}
.is_prior_interaction         <- function(prior){
  if(is.null(attr(prior, "interaction"))){
    return(FALSE)
  }else{
    return(attr(prior, "interaction"))
  }
}

#' @title Reports whether x is a a prior object
#'
#' @description Reports whether x is a a prior object. Note that
#' point priors inherit the prior.simple property
#'
#' @param x an object of test
#'
#' @examples
#' # create some prior distributions
#' p0 <- prior(distribution = "point",  parameters = list(location = 0))
#' p1 <- prior_PET(distribution = "normal", parameters = list(mean = 0, sd = 1))
#'
#' is.prior(p0)
#' is.prior.simple(p0)
#' is.prior.point(p0)
#' is.prior.PET(p0)
#'
#' is.prior(p1)
#' is.prior.simple(p1)
#' is.prior.point(p1)
#' is.prior.PET(p1)
#'
#' @return returns a boolean indicating whether the test object
#' is a prior (of specific type).
#'
#' @export is.prior
#' @export is.prior.simple
#' @export is.prior.discrete
#' @export is.prior.vector
#' @export is.prior.point
#' @export is.prior.none
#' @export is.prior.PET
#' @export is.prior.PEESE
#' @export is.prior.weightfunction
#' @export is.prior.factor
#' @export is.prior.orthonormal
#' @export is.prior.meandif
#' @export is.prior.dummy
#' @export is.prior.independent
#' @export is.prior.spike_and_slab
#' @name is.prior
NULL

#' @rdname is.prior
is.prior                 <- function(x){
  inherits(x, "prior")
}
#' @rdname is.prior
is.prior.point           <- function(x){
  inherits(x, "prior.point")
}
#' @rdname is.prior
is.prior.none            <- function(x){
  inherits(x, "prior.none")
}
#' @rdname is.prior
is.prior.simple          <- function(x){
  inherits(x, "prior.simple")
}
#' @rdname is.prior
is.prior.discrete        <- function(x){
  inherits(x, "prior.discrete")
}
#' @rdname is.prior
is.prior.vector          <- function(x){
  inherits(x, "prior.vector")
}
#' @rdname is.prior
is.prior.PET             <- function(x){
  inherits(x, "prior.PET")
}
#' @rdname is.prior
is.prior.PEESE           <- function(x){
  inherits(x, "prior.PEESE")
}
#' @rdname is.prior
is.prior.weightfunction  <- function(x){
  inherits(x, "prior.weightfunction")
}
#' @rdname is.prior
is.prior.factor          <- function(x){
  inherits(x, "prior.factor")
}
#' @rdname is.prior
is.prior.orthonormal     <- function(x){
  inherits(x, "prior.orthonormal")
}
#' @rdname is.prior
is.prior.dummy           <- function(x){
  inherits(x, "prior.dummy")
}
#' @rdname is.prior
is.prior.independent     <- function(x){
  inherits(x, "prior.independent")
}
#' @rdname is.prior
is.prior.spike_and_slab  <- function(x){
  inherits(x, "prior.spike_and_slab")
}
#' @rdname is.prior
is.prior.meandif         <- function(x){
  inherits(x, "prior.meandif")
}

.check_prior <- function(prior, name = "prior"){
  if(!is.prior(prior))
    stop(paste0("The '", name, "' argument must be a valid prior object."))
}


# clean user input
.prior_clean_input_name <- function(name){

  name <- tolower(name)
  name <- gsub(".", "", name, fixed = TRUE)
  name <- gsub(",", "", name, fixed = TRUE)
  name <- gsub(";", "", name, fixed = TRUE)
  name <- gsub("/", "", name, fixed = TRUE)
  name <- gsub("-", "", name, fixed = TRUE)
  name <- gsub("_", "", name, fixed = TRUE)
  name <- gsub(" ", "", name, fixed = TRUE)

  return(name)
}
