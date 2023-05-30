#' @title Prints a prior object
#'
#' @param x a prior
#' @param short_name whether prior distribution names should be
#' shorted
#' @param parameter_names whether parameter names should be
#' printed
#' @param plot to return \link[base]{bquote} formatted
#'   prior name for plotting.
#' @param digits_estimates number of decimals to be displayed
#'   for printed parameters.
#' @param silent to silently return the print message.
#' @param ... additional arguments
#'
#' @examples
#' # create some prior distributions
#' p0 <- prior(distribution = "point",  parameters = list(location = 0))
#' p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
#'
#' # print them
#' p0
#' p1
#'
#' # use short names
#' print(p1, short_name = TRUE)
#'
#' # print parameter names
#' print(p1, parameter_names = TRUE)
#'
#' # generate bquote plotting syntax
#' plot(0, main = print(p1, plot = TRUE))
#'
#' @return \code{print.prior} invisibly returns the print statement.
#'
#' @seealso [prior()]
#' @rdname print.prior
#' @export
print.prior <- function(x, short_name = FALSE, parameter_names = FALSE, plot = FALSE,
                        digits_estimates = 2, silent = FALSE, ...){

  .check_prior(x, "x")
  check_bool(short_name, "short_name")
  check_bool(parameter_names, "parameter_names")
  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(plot, "plot")
  check_bool(silent, "silent")

  if(plot){
    silent <- TRUE
  }


  if(is.prior.none(x)){
    output <- .print.prior.none(x, short_name, parameter_names, plot, digits_estimates, silent)
  }else if(is.prior.simple(x) || is.prior.vector(x)){
    output <- .print.prior.simple(x, short_name, parameter_names, plot, digits_estimates, silent)
  }else if(is.prior.weightfunction(x)){
    output <- .print.prior.weightfunction(x, short_name, parameter_names, plot, digits_estimates, silent)
  }else if(is.prior.spike_and_slab(x)){
    output <- .print.prior.spike_and_slab(x, short_name, parameter_names, plot, digits_estimates, silent)
  }


  if(!silent){
    cat(output)
  }
  return(invisible(output))
}

.print.prior.simple         <- function(x, short_name, parameter_names, plot, digits_estimates, silent){

  # check whether the range was truncated (before the object is modified)
  if(is.prior.vector(x)){
    needs_truncation <- FALSE
  }else{
    needs_truncation <- !.is_prior_default_range(x)
  }

  # deal with exceptions - Cauchy is passed as T
  if(x[["distribution"]] == "t" && x[["parameters"]][["df"]] == 1){
      x[["distribution"]] <- "Cauchy"
      x[["parameters"]]   <- list(
        location = x[["parameters"]][["location"]],
        scale    = x[["parameters"]][["scale"]])
  }else if(x[["distribution"]] == "mt" && x[["parameters"]][["df"]] == 1){
    x[["distribution"]] <- "mCauchy"
    x[["parameters"]]   <- list(
      location = x[["parameters"]][["location"]],
      scale    = x[["parameters"]][["scale"]])
  }

  ### prepare prior name
  if(short_name){
    out_name <- switch(
      x[["distribution"]],
      "normal"       = "N",
      "lognormal"    = "Ln",
      "t"            = "T",
      "Cauchy"       = "C",
      "gamma"        = "G",
      "invgamma"     = "Ig",
      "point"        = "S",
      "beta"         = "B",
      "bernoulli"    = "Br",
      "exp"          = "E",
      "uniform"      = "U",
      "mnormal"      = "mN",
      "mt"           = "mT",
      "mCauchy"      = "mC",
      "mpoint"       = "mS"
    )
  }else{
    out_name <- switch(
      x[["distribution"]],
      "normal"       = "Normal",
      "lognormal"    = "Lognormal",
      "t"            = "Student-t",
      "Cauchy"       = "Cauchy",
      "gamma"        = "Gamma",
      "invgamma"     = "InvGamma",
      "point"        = "Spike",
      "beta"         = "Beta",
      "bernoulli"    = "Bernoulli",
      "exp"          = "Exponential",
      "uniform"      = "Uniform",
      "mnormal"      = "mNormal",
      "mt"           = "mStudent-t",
      "mCauchy"      = "mCauchy",
      "mpoint"       = "mSpike"
    )
  }

  # add prefix
  if(is.prior.PET(x)){
    out_prefix <- "PET ~ "
  }else if(is.prior.PEESE(x)){
    out_prefix <- "PEESE ~ "
  }else if(is.prior.treatment(x)){
    out_prefix <- "treatment contrast: "
  }else if(is.prior.orthonormal(x)){
    out_prefix <- "orthonormal contrast: "
  }else if(is.prior.meandif(x)){
    out_prefix <- "mean difference contrast: "
  }else if(is.prior.independent(x)){
    out_prefix <- "independent contrast: "
  }else{
    out_prefix <- NULL
  }

  # remove the dimensions parameter from multivariate prior distributions
  if(is.prior.vector(x)){
    x[["parameters"]] <- x[["parameters"]][names(x[["parameters"]]) != "K"]
  }

  ### prepare prior parameters
  # round the parameters and truncation for printing
  for(i in seq_along(x[["parameters"]])){
    x[["parameters"]][[i]] <- round(x[["parameters"]][[i]], digits_estimates)
  }
  for(i in seq_along(x[["truncation"]])){
    x[["truncation"]][[i]] <- round(x[["truncation"]][[i]], digits_estimates)
  }

  if(parameter_names){
    out_parameters <- paste(sapply(seq_along(x[["parameters"]]), function(i){
      paste0(names(x[["parameters"]])[i], " = ", x[["parameters"]][[i]])
    }), collapse = ", ")
  }else{
    out_parameters <- paste(sapply(seq_along(x[["parameters"]]), function(i){
      x[["parameters"]][[i]]
    }), collapse = ", ")
  }


  ### prepare prior truncation
  if(needs_truncation){
    out_truncation <- paste(x[["truncation"]], collapse = ", ")
  }else{
    out_truncation <- NULL
  }


  ### combine the results together
  if(!plot){
    output <- paste0(
      if(!is.null(out_prefix)) out_prefix,
      out_name,
      if(!is.null(out_parameters)) paste0("(", out_parameters, ")"),
      if(!is.null(out_truncation)) paste0("[", out_truncation, "]"))
  }else{
    output <- bquote(
      .(if(!is.null(out_prefix)){ out_prefix }else{ bquote() } )*
        italic(.(out_name))*
        .(if(!is.null(out_parameters)){ paste0("(", out_parameters, ")") }else{ bquote() } )*
        .(if(!is.null(out_truncation)){
          if(is.infinite(x[["truncation"]][["lower"]])){
            bquote(""["["*-infinity*", "*.(x[["truncation"]][["upper"]])*"]"])
          }else if(is.infinite(x[["truncation"]][["upper"]])){
            bquote(""["["*.(x[["truncation"]][["lower"]])*", "*infinity*"]"])
          }else{
            bquote(""[.(paste0("[", out_truncation, "]"))])
          }
        }else{ bquote() })
    )
  }

  return(output)
}
.print.prior.weightfunction <- function(x, short_name, parameter_names, plot, digits_estimates, silent){

  ### prepare prior name
  # add prefix
  if(plot){
    out_prefix <- bquote(omega)
  }else{
    out_prefix <- "omega"
  }

  # type of steps
  if(short_name){
    steps_name <- switch(
      x[["distribution"]],
      "two.sided"    = "2s: ",
      "one.sided"    = "1s: ",
      "two.sided.fixed" = "2s: ",
      "one.sided.fixed" = "1s: "
    )
  }else{
    steps_name <- switch(
      x[["distribution"]],
      "two.sided"    = "two-sided: ",
      "one.sided"    = "one-sided: ",
      "two.sided.fixed" = "two-sided: ",
      "one.sided.fixed" = "one-sided: "
    )
  }

  # add steps
  out_steps  <- paste(trimws(x$parameters[["steps"]], which = "left", whitespace = "0"), collapse = ", ")

  # distribution
  if(all(names(x[["parameters"]]) %in% c("alpha", "steps"))){

    out_parameters <- paste(round(x$parameters[["alpha"]], digits_estimates), collapse = ", ")
    if(parameter_names){
      out_parameters <- paste0("alpha = ", out_parameters)
    }
    if(short_name){
      out_distribution <- paste0("CumD")
    }else{
      out_distribution <- paste0("CumDirichlet")
    }

    if(!plot){
      output <- paste0(out_prefix, "[", steps_name, out_steps, "]", " ~ ", out_distribution, "(", out_parameters, ")")
    }else{
      output <- bquote(italic(.(out_prefix))[.(steps_name)*.(out_steps)]~"~"~italic(.(out_distribution))*(.(out_parameters)))
    }

  }else if(all(names(x[["parameters"]]) %in% c("alpha1", "alpha2", "steps"))){

    out_parameters1 <- paste(round(x$parameters[["alpha1"]],      digits_estimates), collapse = ", ")
    out_parameters2 <- paste(round(rev(x$parameters[["alpha2"]]), digits_estimates), collapse = ", ")
    if(parameter_names){
      out_parameters1 <- paste0("alpha1 = ", out_parameters1)
      out_parameters2 <- paste0("alpha2 = ", out_parameters2)
    }
    if(short_name){
      out_distribution1 <- paste0("CumD")
      out_distribution2 <- paste0("rCumD")
    }else{
      out_distribution1 <- paste0("CumDirichlet")
      out_distribution2 <- paste0("revCumDirichlet")
    }

    if(!plot){
      output <- paste0(out_prefix, "[", steps_name, out_steps, "]", " ~ ", out_distribution1, "(", out_parameters1, "), ", out_distribution2, "(", out_parameters2, ")")
    }else{
      output <- bquote(italic(.(out_prefix))[.(steps_name)*.(out_steps)]~"~"~italic(.(out_distribution1))*(.(out_parameters1))~","~~italic(.(out_distribution2))*(.(out_parameters2)))
    }

  }else if(all(names(x[["distribution"]]) %in% c("one.sided.fixed", "two.sided.fixed"))){

    out_parameters <- paste0(round(x$parameters[["omega"]], digits_estimates), collapse = ", ")

    if(!plot){
      output <- paste0(out_prefix, "[", steps_name, out_steps, "]", " = ", "(", out_parameters, ")")
    }else{
      output <- bquote(italic(.(out_prefix))[.(steps_name)*.(out_steps)]~"="~(.(out_parameters)))
    }
  }

  return(output)
}
.print.prior.none           <- function(x, short_name, parameter_names, plot, digits_estimates, silent){

  ### prepare prior name
  out_name <- "None"


  ### combine the results together
  if(!plot){
    output <- out_name
  }else{
    output <- bquote(italic(.(out_name)))
  }

  return(output)
}
.print.prior.spike_and_slab <- function(x, short_name, parameter_names, plot, digits_estimates, silent){

  variable  <- print(x[["variable"]],  short_name, parameter_names, plot, digits_estimates, silent = TRUE)
  inclusion <- print(x[["inclusion"]], short_name, parameter_names, plot, digits_estimates, silent = TRUE)

  ### combine the results together
  if(!plot){
    output <- paste0(variable, " * ", inclusion)
  }else{
    output <- bquote(.(variable) ~ "*" ~ .(inclusion))
  }

  return(output)
}
