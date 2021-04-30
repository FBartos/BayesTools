# functions for checking input into custom distributions related functions
.check_log        <- function(log){
  if(length(log) != 1 & !is.logical(log))
    stop("'log' must be a logical vector of length one.")
}
.check_log.p      <- function(log.p){
  if(length(log.p) != 1 & !is.logical(log.p))
    stop("'log.p' must be a logical vector of length one.")
}
.check_lower.tail <- function(lower.tail){
  if(length(lower.tail) != 1 & !is.logical(lower.tail))
    stop("'lower.tail' must be a logical vector of length one.")
}
.check_x <- function(x, lower = -Inf, upper = Inf, name = "x"){
  if(!is.numeric(x) | !is.vector(x))
    stop(paste0("'", name ,"' must be a numeric vector."))
  if(!is.infinite(lower))if(any(x < lower))
    stop(paste0("all '", name ,"' must be higher than ", lower))
  if(!is.infinite(upper))if(any(x > upper))
    stop(paste0("all '", name ,"' must be lower than ", upper))
}
.check_n <- function(n, name = "n"){
  if(!is.numeric(n) | length(n) != 1)
    stop(paste0("'", name ,"' must be a numeric vector of length one."))
  if(n <= 0)
    stop(paste0("'", name ,"' must be larger than 0"))
  if(!.is.wholenumber(n))
    stop(paste0("'", name ,"' must be an integer"))
}
.check_q <- function(q, lower = -Inf, upper = Inf, name = "n"){
  if(!is.numeric(q) | !is.vector(q))
    stop(paste0("'", name ,"' must be a numeric vector."))
  if(!is.infinite(lower))if(any(q < lower))
    stop(paste0("all '", name ,"' must be higher than ", lower))
  if(!is.infinite(upper))if(any(q > upper))
    stop(paste0("all '", name ,"' must be lower than ", upper))
}
.check_p <- function(p, log.p, name = "n"){
  if(!is.numeric(p) | !is.vector(p))
    stop(paste0("'", name ,"' must be a numeric vector."))
  if(!log.p){
    if(any(p < 0) | any(p > 1))
      stop(paste0("all '", name ,"' must be between 0 and 1."))
  }else{
    if(any(p > 0))
      stop(paste0("all log('", name ,"') must be lower or equal to 0."))
  }
}
.check_logical    <- function(x, name = "x"){
  if(!is.logical(x) | length(x) != 1)
    stop(paste0("The '", x, "' argument must be a logical vector of length 1."))
}
.is.wholenumber   <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

