#' @title Check input
#'
#' @description A set of convenience functions for checking
#' objects/arguments to a function passed by a user.
#'
#' @param x object to be checked
#' @param name name of the object that will be print in the error
#' message.
#' @param check_length length of the object to be checked. Defaults to
#' \code{1}. Set to \code{0} in order to not check object length.
#' @param allow_NULL whether the object can be \code{NULL}.
#' If so, no checks are executed.
#' @param check_names names of entries allowed in a list. Defaults to
#' \code{NULL} (do not check).
#' @param all_objects whether all entries in \code{check_names} must be
#' present. Defaults to \code{FALSE}.
#' @param allow_other whether additional entries then the specified in
#' \code{check_names} might be present
#' @param lower lower bound of allowed values.
#' Defaults to \code{-Inf} (do not check).
#' @param upper upper bound of allowed values.
#' Defaults to \code{Inf} (do not check).
#' @param allow_bound whether the values at the boundary are allowed.
#' Defaults to \code{TRUE}.
#'
#' @name check_input
#' @export check_bool
#' @export check_char
#' @export check_int
#' @export check_real
#' @export check_list

#' @rdname check_input
check_bool   <- function(x, name, check_length = 1, allow_NULL = FALSE){

  if(is.null(x)){
    if(allow_NULL){
      return()
    }else{
      stop(paste0("The '", name, "' argument cannot be NULL."))
    }
  }

  if(!is.logical(x) | !is.vector(x))
    stop(paste0("The '", name, "' argument must be a logical vector."))

  if(check_length != 0  && length(x) != check_length)
    stop(paste0("The '", name, "' argument must have length '", check_length, "'."))

  return()
}

#' @rdname check_input
check_char   <- function(x, name, check_length = 1, allow_NULL = FALSE){

  if(is.null(x)){
    if(allow_NULL){
      return()
    }else{
      stop(paste0("The '", name, "' argument cannot be NULL."))
    }
  }

  if(!is.character(x) | !is.vector(x))
    stop(paste0("The '", name, "' argument must be a character vector."))

  if(check_length != 0 && length(x) != check_length)
    stop(paste0("The '", name, "' argument must have length '", check_length, "'."))

  return()
}

#' @rdname check_input
check_real   <- function(x, name, lower = -Inf, upper = Inf, allow_bound = TRUE, check_length = 1, allow_NULL = FALSE){

  if(is.null(x)){
    if(allow_NULL){
      return()
    }else{
      stop(paste0("The '", name, "' argument cannot be NULL."))
    }
  }

  if(!is.numeric(x) | !is.vector(x))
    stop(paste0("The '", name, "' argument must be a numeric vector."))

  if(!is.infinite(lower)){
    if(!allow_bound){
      if(any(x <= lower))
        stop(paste0("All '", name ,"' must be higher than ", lower,"."))
    }else{
      if(any(x < lower))
        stop(paste0("All '", name ,"' must be equal or higher than ", lower,"."))
    }
  }

  if(!is.infinite(upper)){
    if(!allow_bound){
      if(any(x >= upper))
        stop(paste0("All '", name ,"' must be lower than ", upper,"."))
    }else{
      if(any(x > upper))
        stop(paste0("All '", name ,"' must be equal or lower than ", upper,"."))
    }
  }

  if(check_length != 0 && length(x) != check_length)
    stop(paste0("The '", name, "' argument must have length '", check_length, "'."))

  return()
}

#' @rdname check_input
check_int    <- function(x, name, lower = -Inf, upper = Inf, allow_bound = TRUE, check_length = 1, allow_NULL = FALSE){

  if(is.null(x)){
    if(allow_NULL){
      return()
    }else{
      stop(paste0("The '", name, "' argument cannot be NULL."))
    }
  }

  check_real(x, name, lower, upper, allow_bound, check_length, allow_NULL)

  if(!all(.is.wholenumber(x)))
    stop(paste0("The '", name ,"' argument must be an integer vector."))

  return()
}

#' @rdname check_input
check_list   <- function(x, name, check_length = 0, check_names = NULL, all_objects = FALSE, allow_other = FALSE, allow_NULL = FALSE){

  if(is.null(x)){
    if(allow_NULL){
      return()
    }else{
      stop(paste0("The '", name, "' argument cannot be NULL."))
    }
  }

  if(!is.list(x))
    stop(paste0("The '", name, "' argument must be a list."))

  if(check_length != 0 && length(x) != check_length)
    stop(paste0("The '", name, "' argument must have length '", check_length, "'."))

  if(!is.null(check_names)){

    if(all_objects && any(!check_names %in% names(x)))
      stop(paste0("The '", paste0(check_names[!check_names %in% names(x)], collapse = "', '") ,"' objects are missing in the '", name, "' argument."))

    if(!allow_other && any(!names(x) %in% check_names))
      stop(paste0("The '", paste0(names(x)[!names(x) %in% check_names], collapse = "', '") ,"' objects are not recognized by the '", name, "' argument."))
  }

  return()
}

# helper functions
.is.wholenumber  <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}
