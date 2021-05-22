# functions for checking input into custom distributions related functions
.check_log        <- function(log){
  check_bool(log, "log")
}
.check_log.p      <- function(log.p){
  check_bool(log.p, "log.p")
}
.check_lower.tail <- function(lower.tail){
  check_bool(lower.tail, "lower.tail")
}
.check_x <- function(x, lower = -Inf, upper = Inf){
  check_real(x, "x", lower = lower, upper = upper, check_length = 0)
}
.check_n <- function(n){
  check_int(n, "n", lower = 1)
}
.check_q <- function(q, lower = -Inf, upper = Inf){
  check_real(q, "q", lower = lower, upper = upper, check_length = 0)
}
.check_p <- function(p, log.p){
  if(!log.p){
    check_real(p, "p", lower = 0, upper = 1, check_length = 0)
  }else{
    check_real(p, "log(p)", upper = 0, check_length = 0)
  }
}
