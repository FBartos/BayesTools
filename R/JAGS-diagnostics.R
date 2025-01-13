#' @title Plot diagnostics of a 'JAGS' model
#'
#' @description Creates density plots, trace plots, and autocorrelation plots
#' for a given parameter of a JAGS model.
#'
#' @param fit a JAGS model fitted via [JAGS_fit()]
#' @param parameter parameter to be plotted
#' @param type what type of model diagnostic should be plotted. The available
#' options are \code{"density"}, \code{"trace"}, and \code{"autocorrelation"}
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting.
#' @param xlim x plotting range
#' @param ylim y plotting range
#' @param lags number of lags to be shown for the autocorrelation plot.
#' Defaults to \code{30}.
#' @param ... additional arguments
#' @inheritParams density.prior
#' @inheritParams plot.prior
#' @inheritParams print.prior
#' @inheritParams BayesTools_model_tables
#'
#' @seealso [JAGS_fit()] [JAGS_check_convergence()]
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @name JAGS_diagnostics
#' @aliases JAGS_diagnostics_density JAGS_diagnostics_autocorrelation JAGS_diagnostics_trace
#' @export JAGS_diagnostics
#' @export JAGS_diagnostics_density
#' @export JAGS_diagnostics_autocorrelation
#' @export JAGS_diagnostics_trace

#' @rdname JAGS_diagnostics
JAGS_diagnostics                 <- function(fit, parameter, type, plot_type = "base",
                                             xlim = NULL, ylim = NULL, lags = 30, n_points = 1000,
                                             transformations = NULL, transform_factors = FALSE, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){


  # check fits
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit")

  check_char(type, "type", allow_values = c("density", "trace", "autocorrelation"))
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(parameter, "parameter", allow_values = c(names(prior_list), if(any(names(prior_list) == "bias")) c("PET", "PEESE", "omega")))
  if(!is.null(transformations))
    check_char(names(transformations), "names(transformations)", allow_values = parameter)

  # do not produce diagnostics for a spike prior
  if(is.prior.point(prior_list[[parameter]])){
    message("No diagnostic plots are produced for a spike prior distribution")
    return(NULL)
  }

  # depreciate
  transform_factors <- .depreciate.transform_orthonormal(transform_orthonormal, transform_factors)

  # prepare the plot data
  plot_data <- .diagnostics_plot_data(fit = fit, parameter = parameter, prior_list = prior_list, transformations = transformations, transform_factors = transform_factors)
  plot_data <- switch(
    type,
    "density"         = .diagnostics_plot_data_density(plot_data, n_points, xlim),
    "trace"           = .diagnostics_plot_data_trace(plot_data, n_points, ylim),
    "autocorrelation" = .diagnostics_plot_data_autocorrelation(plot_data, n_points, lags)
  )


  # prepare nice parameter names
  if(!is.null(attr(prior_list[[parameter]], "parameter"))){
    parameter_name <- format_parameter_names(attr(plot_data, "parameter_name"), attr(prior_list[[parameter]], "parameter"), formula_prefix = formula_prefix)
  }else{
    parameter_name <- attr(plot_data, "parameter_name")
  }

  # get default plot settings
  dots      <- list(...)
  main      <- print(prior_list[[parameter]], plot = TRUE, short_name = short_name, parameter_names = parameter_names)
  xlab      <- switch(
    type,
    "density"         = parameter_name,
    "trace"           = "Iteration",
    "autocorrelation" = "Lag"
  )
  ylab      <- switch(
    type,
    "density"         = "Density",
    "trace"           = parameter_name,
    "autocorrelation" = paste0("Autocorrelation(", parameter_name, ")")
  )


  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab

  if(is.null(dots[["col"]]))  dots$col  <-  "black"
  if(is.null(dots[["lty"]]))  dots$lty  <-  1
  if(is.null(dots[["lwd"]]))  dots$lwd  <-  1

  if(length(dots[["col"]]) == 1)  dots$col  <- rep(dots$col, attr(plot_data, "chains"))
  if(length(dots[["lty"]]) == 1)  dots$lty  <- rep(dots$lty, attr(plot_data, "chains"))
  if(length(dots[["lwd"]]) == 1)  dots$lwd  <- rep(dots$lwd, attr(plot_data, "chains"))

  plots <- list()

  for(i in seq_along(plot_data)){

    temp_dots <- dots

    temp_dots[["xlab"]] <- if(length(temp_dots[["xlab"]]) > 1) temp_dots[["xlab"]][i] else temp_dots[["xlab"]]
    temp_dots[["ylab"]] <- if(length(temp_dots[["ylab"]]) > 1) temp_dots[["ylab"]][i] else temp_dots[["ylab"]]

    if(is.null(temp_dots[["xlim"]])) temp_dots$xlim <-  attr(plot_data[[i]], "x_range")
    if(is.null(temp_dots[["ylim"]])) temp_dots$ylim <-  attr(plot_data[[i]], "y_range")


    if(plot_type == "ggplot"){
      plots[[i]] <- .ggplot.prior_empty("simple", temp_dots)
    }else{
      .plot.prior_empty("simple", temp_dots)
    }

    for(j in seq_along(plot_data[[i]])){

      temp_args <- list()
      temp_args[["plot_data"]] <- plot_data[[i]][[j]]
      temp_args[["col"]]       <- dots[["col"]][j]
      temp_args[["lwd"]]       <- dots[["lwd"]][j]
      temp_args[["lty"]]       <- dots[["lty"]][j]

      if(plot_type == "ggplot"){
        plots[[i]] <- plots[[i]] + do.call(switch(
          type,
          "density"         = .geom_diagnostics.density,
          "trace"           = .geom_diagnostics.trace,
          "autocorrelation" = .geom_diagnostics.autocorrelation
        ), temp_args)
      }else{
        do.call(switch(
          type,
          "density"         = .lines_diagnostics.density,
          "trace"           = .lines_diagnostics.trace,
          "autocorrelation" = .lines_diagnostics.autocorrelation
        ), temp_args)
      }

    }
  }


  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    if(length(plots) == 1){
      plots <- plots[[1]]
    }
    return(plots)
  }


}

#' @rdname JAGS_diagnostics
JAGS_diagnostics_density         <- function(fit, parameter, plot_type = "base",
                                             xlim = NULL, n_points = 1000,
                                             transformations = NULL, transform_factors = FALSE, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  JAGS_diagnostics(fit = fit, parameter = parameter, plot_type = plot_type, type = "density",
                   xlim = xlim, n_points = n_points,
                   transformations = transformations, transform_factors = transform_factors, transform_orthonormal = transform_orthonormal,
                   short_name = short_name, parameter_names = parameter_names, formula_prefix = formula_prefix, ...)
}

#' @rdname JAGS_diagnostics
JAGS_diagnostics_trace           <- function(fit, parameter, plot_type = "base",
                                             ylim = NULL,
                                             transformations = NULL, transform_factors = FALSE, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  JAGS_diagnostics(fit = fit, parameter = parameter, plot_type = plot_type, type = "trace",
                   ylim = ylim,
                   transformations = transformations, transform_factors = transform_factors, transform_orthonormal = transform_orthonormal,
                   short_name = short_name, parameter_names = parameter_names, formula_prefix = formula_prefix, ...)
}

#' @rdname JAGS_diagnostics
JAGS_diagnostics_autocorrelation <- function(fit, parameter, plot_type = "base",
                                             lags = 30,
                                             transformations = NULL, transform_factors = FALSE, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  JAGS_diagnostics(fit = fit, parameter = parameter, plot_type = plot_type, type = "autocorrelation",
                   lags = lags,
                   transformations = transformations, transform_factors = transform_factors, transform_orthonormal = transform_orthonormal,
                   short_name = short_name, parameter_names = parameter_names, formula_prefix = formula_prefix, ...)
}


.diagnostics_plot_data                 <- function(fit, parameter, prior_list, transformations, transform_factors){

  check_list(transformations, "transformations", allow_NULL = TRUE)
  if(!is.null(transformations) && any(!sapply(transformations, function(trans)is.function(trans[["fun"]]))))
    stop("'transformations' must be list of functions in the 'fun' element.")

  model_samples <- coda::as.mcmc.list(fit)
  samples_chain <- lapply(seq_along(model_samples), function(i) {
    return(rep(i, nrow(model_samples[[i]])))
  })
  samples_iter  <- lapply(seq_along(model_samples), function(i) {
    return(1:nrow(model_samples[[i]]))
  })
  model_samples <- do.call(rbind, model_samples)


  # extract the relevant parameters
  if(is.prior.spike_and_slab(prior_list[[parameter]])){

    if(sum(model_samples[,colnames(model_samples) == paste0(parameter, "_indicator")] != 0) < 10)
      stop("The parameter with a spike and slab prior did not result in enough samples under the slab for producing a diagnostic figure.")

    # change the samples between conditional
    model_samples[
      model_samples[,colnames(model_samples) == paste0(parameter, "_indicator")] == 0,
      colnames(model_samples) == parameter] <- NA

    # modify the parameter list
    prior_list[[parameter]] <- prior_list[[parameter]]$variable
  }

  if(is.prior.factor(prior_list[[parameter]])){
    if(.get_prior_factor_levels(prior_list[[parameter]]) > 1){
      model_samples <- model_samples[,paste0(parameter, "[", 1:.get_prior_factor_levels(prior_list[[parameter]]), "]"),drop = FALSE]
    }else{
      model_samples <- model_samples[,parameter,drop = FALSE]
    }
  }else if(is.prior.weightfunction(prior_list[[parameter]])){
    model_samples <- model_samples[,paste0("omega", "[", (length(weightfunctions_mapping(list(prior_list[[parameter]]), cuts_only = TRUE)) - 2):1, "]"),drop = FALSE]
  }else if(parameter == "omega" && !is.null(prior_list[["bias"]])){
    model_samples <- model_samples[,paste0("omega", "[", 2:(length(weightfunctions_mapping(prior_list[["bias"]][sapply(prior_list[["bias"]], is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)) - 1), "]"),drop = FALSE]
  }else if(is.prior.vector(prior_list[[parameter]])){
    if(prior_list[[parameter]]$parameters[["K"]] > 1){
      model_samples <- model_samples[,paste0(parameter, "[", 1:prior_list[[parameter]]$parameters[["K"]], "]"),drop = FALSE]
    }else{
      model_samples <- model_samples[,parameter,drop = FALSE]
    }
  }else{
    if(is.prior.PET(prior_list[[parameter]])){
      parameter <- "PET"
    }else if(is.prior.PEESE(prior_list[[parameter]])){
      parameter <- "PEESE"
    }
    model_samples <- model_samples[,parameter,drop = FALSE]
  }

  # deal with bias mixture dispatching
  if(parameter %in% c("PET", "PEESE", "omega") && !is.null(prior_list[["bias"]])){
    prior_list      <- prior_list[["bias"]]
    parameter_names <- parameter
  }else{
    prior_list      <- prior_list[parameter]
    parameter_names <- parameter
  }


  # mostly adapted from runjags_estimates_table
  # apply transformations
  if(!is.null(transformations)){
    for(par in names(transformations)){
      model_samples[,par] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[,par]), transformations[[par]][["arg"]]))
    }
  }

  # transform meandif and orthonormal factors to differences from runjags_estimates_table
  if(transform_factors & any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))){
    for(par in names(prior_list)[sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x))]){

      if(.get_prior_factor_levels(prior_list[[par]]) == 1){
        par_names <- par
      }else{
        par_names <- paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")
      }

      original_samples <- model_samples[,par_names,drop = FALSE]

      if(is.prior.orthonormal(prior_list[[par]])){
        model_samples <- original_samples %*% t(contr.orthonormal(1:(.get_prior_factor_levels(prior_list[[par]])+1)))
      }else if(is.prior.meandif(prior_list[[par]])){
        model_samples <- original_samples %*% t(contr.meandif(1:(.get_prior_factor_levels(prior_list[[par]])+1)))
      }


      if(attr(prior_list[[par]], "interaction")){
        if(length(.get_prior_factor_level_names(prior_list[[par]])) == 1){
          parameter_names <- paste0(par, " [dif: ", .get_prior_factor_level_names(prior_list[[par]])[[1]],"]")
        }else{
          stop("orthonormal/meandif de-transformation for interaction of multiple factors is not implemented.")
        }
      }else{
        parameter_names <- paste0(par, " [dif: ", .get_prior_factor_level_names(prior_list[[par]]),"]")
      }
    }
  }else if(any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))){
    for(par in names(prior_list)[sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x))]){
      if(.get_prior_factor_levels(prior_list[[par]]) == 1){
        parameter_names <- par
      }else{
        parameter_names <- paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")
      }
    }
  }

  # rename treatment factor levels
  if(any(sapply(prior_list, is.prior.treatment))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.treatment)]){
      if(!.is_prior_interaction(prior_list[[par]])){
        if(.get_prior_factor_levels(prior_list[[par]]) == 1){
          parameter_names <- par
        }else{
          parameter_names <- paste0(par,"[",.get_prior_factor_level_names(prior_list[[par]])[-1], "]")
        }
      }else if(length(attr(prior_list[[par]], "levels")) == 1){
        parameter_names <- paste0(par,"[",.get_prior_factor_level_names(prior_list[[par]])[[1]][-1], "]")
      }
    }
  }

  # rename independent factor levels
  if(any(sapply(prior_list, is.prior.independent))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.independent)]){
      if(!.is_prior_interaction(prior_list[[par]])){
        if(.get_prior_factor_levels(prior_list[[par]]) == 1){
          parameter_names <- par
        }else{
          parameter_names <- paste0(par,"[",.get_prior_factor_level_names(prior_list[[par]]), "]")
        }
      }else if(length(attr(prior_list[[par]], "levels")) == 1){
        parameter_names <- paste0(par,"[",.get_prior_factor_level_names(prior_list[[par]]), "]")
      }
    }
  }

  # rename weightfunctions factor levels
  if(any(sapply(prior_list, is.prior.weightfunction)) && !is.prior.mixture(prior_list)){
    for(par in names(prior_list)[sapply(prior_list, is.prior.weightfunction)]){
      omega_cuts      <- weightfunctions_mapping(prior_list[par], cuts_only = TRUE)
      parameter_names <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
      parameter_names <- parameter_names[-1]
    }
  }

  if(is.prior.mixture(prior_list) && parameter == "omega"){
    omega_cuts      <- weightfunctions_mapping(prior_list[sapply(prior_list, is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)
    parameter_names <- sapply(2:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
  }

  # attach the relevant attributes
  colnames(model_samples)          <- parameter_names
  attr(model_samples, "chain")     <- do.call(c, samples_chain)
  attr(model_samples, "iter")      <- do.call(c, samples_iter)
  attr(model_samples, "parameter") <- parameter
  attr(model_samples, "prior")     <- if(is.prior.mixture(prior_list)) prior_list else prior_list[[parameter]]

  return(model_samples)
}
.diagnostics_plot_data_density         <- function(plot_data, n_points, xlim){

  chain <- attr(plot_data, "chain")
  prior <- attr(plot_data, "prior")

  if(is.prior.weightfunction(prior)){
    prior_lower <- 0
    prior_upper <- 1
  }else{
    prior_lower <- prior$truncation[["lower"]]
    prior_upper <- prior$truncation[["upper"]]
  }

  out   <- list()

  for(i in 1:ncol(plot_data)){

    if(is.null(xlim)){
      x_range <- range(plot_data[,i], na.rm = TRUE)
    }else{
      x_range <- xlim
    }

    for(j in seq_along(unique(chain))){

      temp_args    <- list(x = plot_data[chain == j,i], n = n_points, from = x_range[1], to = x_range[2], na.rm = TRUE)
      temp_density <- do.call(stats::density, temp_args)

      x_den    <- temp_density$x
      y_den    <- temp_density$y

      # check for truncation
      if(isTRUE(all.equal(prior_lower, x_den[1])) | prior_lower >= x_den[1]){
        y_den <- c(0, y_den)
        x_den <- c(x_den[1], x_den)
      }
      if(isTRUE(all.equal(prior_upper, x_den[length(x_den)])) | prior_upper <= x_den[length(x_den)]){
        y_den <- c(y_den, 0)
        x_den <- c(x_den, x_den[length(x_den)])
      }

      temp_density    <- list(
        call    = call("density"),
        bw      = NULL,
        n       = n_points,
        x       = x_den,
        y       = y_den
      )

      class(temp_density) <- c("density")
      attr(temp_density, "x_range")        <- range(x_den)
      attr(temp_density, "y_range")        <- c(0, max(y_den))
      attr(temp_density, "chain")          <- j
      attr(temp_density, "parameter")      <- attr(plot_data, "parameter")
      attr(temp_density, "parameter_name") <- colnames(plot_data)[i]

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_density
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- x_range
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- c(0, max(sapply(out[[i]], function(x) attr(x, "y_range"))))
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter")
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter")
  attr(out, "parameter_name") <- colnames(plot_data)

  return(out)
}
.diagnostics_plot_data_trace           <- function(plot_data, n_points, ylim){

  chain <- attr(plot_data, "chain")
  iter  <- attr(plot_data, "iter")
  prior <- attr(plot_data, "prior")

  out   <- list()

  for(i in 1:ncol(plot_data)){

    if(is.null(ylim)){
      y_range <- range(plot_data[,i], na.rm = TRUE)
    }else{
      y_range <- ylim
    }
    x_range <- range(iter)

    for(j in seq_along(unique(chain))){

      temp_x       <- iter[chain == j]
      temp_y       <- plot_data[chain == j,i]

      temp_chain <- list(
        x       = temp_x,
        y       = temp_y
      )

      class(temp_chain) <- c("BayesTools_chain")
      attr(temp_chain, "x_range")        <- x_range
      attr(temp_chain, "y_range")        <- y_range
      attr(temp_chain, "chain")          <- j
      attr(temp_chain, "parameter")      <- attr(plot_data, "parameter")
      attr(temp_chain, "parameter_name") <- colnames(plot_data)[i]

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_chain
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- x_range
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- y_range
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter")
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter")
  attr(out, "parameter_name") <- colnames(plot_data)

  return(out)
}
.diagnostics_plot_data_autocorrelation <- function(plot_data, n_points, lags){

  chain <- attr(plot_data, "chain")
  iter  <- attr(plot_data, "iter")
  prior <- attr(plot_data, "prior")

  out   <- list()

  for(i in 1:ncol(plot_data)){

    x_range <-

    for(j in seq_along(unique(chain))){

      temp_x  <- 0:lags
      temp_y  <- stats::acf(plot_data[chain == j,i], lag.max = lags, plot = FALSE, na.action = stats::na.pass)$acf[, , 1L]


      temp_autocor <- list(
        x       = temp_x,
        y       = temp_y
      )

      class(temp_autocor) <- c("BayesTools_autocorrelation")
      attr(temp_autocor, "x_range")        <- c(0, lags)
      attr(temp_autocor, "y_range")        <- range(c(0, temp_y))
      attr(temp_autocor, "chain")          <- j
      attr(temp_autocor, "parameter")      <- attr(plot_data, "parameter")
      attr(temp_autocor, "parameter_name") <- colnames(plot_data)[i]

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_autocor
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- c(0, lags)
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- c(0, max(sapply(out[[i]], function(x) attr(x, "y_range"))))
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter")
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter")
  attr(out, "parameter_name") <- colnames(plot_data)

  return(out)
}


.lines_diagnostics.density         <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines_diagnostics.trace           <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines_diagnostics.autocorrelation <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]])) dots[["col"]] else .plot.prior_settings()[["col"]]

  graphics::rect(
    xleft   = plot_data$x + 0.075 - 0.5,
    ybottom = 0,
    xright  = plot_data$x + 0.925 - 0.5,
    ytop    = plot_data$y,
    col     = col)

  return(invisible())
}

.geom_diagnostics.density         <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  geom <- ggplot2::geom_line(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["y"]]),
    linewidth = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_diagnostics.trace           <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  geom <- ggplot2::geom_line(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["y"]]),
    linewidth = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_diagnostics.autocorrelation <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]])) dots[["col"]] else .plot.prior_settings()[["col"]]

  geom <- ggplot2::geom_bar(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes(
      x      = .data[["x"]],
      weight = .data[["y"]]),
    color = col, fill = col)

  return(geom)
}
