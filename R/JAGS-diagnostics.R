#' @title Plot diagnostics of a 'JAGS' model
#'
#' @description Creates density plots, trace plots, and autocorrelation plots
#' for a given parameter of a JAGS model.
#'
#' @param fit a JAGS model
#' @param parameter parameter to be plotted
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting.
#' @param xlim x plotting range
#' @param par_name a type of parameter for which the prior is
#' specified. Only relevant if the prior corresponds to a mu
#' parameter that needs to be transformed.
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when priors for the omega
#' parameter are plotted and \code{individual = TRUE}.
#' @param lags number of lags to be shown for
#' \code{type = "autocorrelation"}. Defaults to \code{30}.
#' @param ... additional arguments
#' @inheritParams density.prior
#'
#'
#' @return \code{JAGS_fit} returns an object of class 'runjags'.
#'
#' @seealso [JAGS_check_convergence()]
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @rdname JAGS_diagnostics
#' @aliases JAGS_diagnostics_density JAGS_diagnostics_autocorrelation JAGS_diagnostics_trace
#' @export JAGS_diagnostics_density
#' @export JAGS_diagnostics_autocorrelation
#' @export JAGS_diagnostics_trace

JAGS_diagnostics_density         <- function(fit, parameter, plot_type = "base",
                                             xlim = NULL, n_points = 1000,
                                             transformations = NULL, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  # check fits
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit")

  check_char(plot_type, "plot_type")
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(parameter, "parameter", allow_values = names(prior_list))
  if(!is.null(transformations))
    check_char(names(transformations), "names(transformations)", allow_values = parameter)

  # do not produce diagnostics for a spike prior
  if(is.prior.point(prior_list[[parameter]])){
    message("No diagnostic plots are produced for a spike prior distribution")
    return(NULL)
  }

  # prepare the plot data
  plot_data <- .diagnostics_plot_data(fit = fit, parameter = parameter, prior_list = prior_list, transformations = transformations, transform_orthonormal = transform_orthonormal)
  plot_data <- .diagnostics_plot_data_density(plot_data, n_points, xlim)

  # prepare nice parameter names
  if(formula_prefix && !is.null(attr(prior_list[[parameter]], "parameter"))){
    parameter_name <- format_parameter_names(attr(plot_data, "parameter_name"), attr(prior_list[[parameter]], "parameter"), formula_prefix = formula_prefix)
  }else{
    parameter_name <- attr(plot_data, "parameter_name")
  }

  # get default plot settings
  dots      <- list(...)
  main      <- print(prior_list[[parameter]], plot = TRUE, short_name = short_name, parameter_names = parameter_names)
  xlab      <- parameter_name
  ylab      <- "Density"

  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab

  if(is.null(dots[["col"]]))  dots$col  <-  "black"
  if(is.null(dots[["lty"]]))  dots$lty  <-  1
  if(is.null(dots[["lwd"]]))  dots$lwd  <-  1

  if(length(dots[["xlab"]]) == 1) dots$xlab <- rep(dots[["xlab"]], attr(plot_data, "chains"))
  if(length(dots[["col"]]) == 1)  dots$col  <- rep(dots$col, attr(plot_data, "chains"))
  if(length(dots[["lty"]]) == 1)  dots$lty  <- rep(dots$lty, attr(plot_data, "chains"))
  if(length(dots[["lwd"]]) == 1)  dots$lwd  <- rep(dots$lwd, attr(plot_data, "chains"))

  plots <- list()

  for(i in seq_along(plot_data)){

    temp_dots <- dots

    temp_dots[["xlab"]] <- temp_dots[["xlab"]][i]
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
        plots[[i]] <- plots[[i]] + do.call(.geom_prior.simple, temp_args)
      }else{
        do.call(.lines.prior.simple, temp_args)
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

JAGS_diagnostics_trace           <- function(fit, parameter, plot_type = "base",
                                             ylim = NULL,
                                             transformations = NULL, transform_orthonormal = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  # check fits
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit")

  check_char(plot_type, "plot_type")
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(parameter, "parameter", allow_values = names(prior_list))
  if(!is.null(transformations))
    check_char(names(transformations), "names(transformations)", allow_values = parameter)

  # do not produce diagnostics for a spike prior
  if(is.prior.point(prior_list[[parameter]])){
    message("No diagnostic plots are produced for a spike prior distribution")
    return(NULL)
  }

  # prepare the plot data
  plot_data <- .diagnostics_plot_data(fit = fit, parameter = parameter, prior_list = prior_list, transformations = transformations, transform_orthonormal = transform_orthonormal)
  plot_data <- .diagnostics_plot_data_trace(plot_data, n_points, ylim)

  # prepare nice parameter names
  if(formula_prefix && !is.null(attr(prior_list[[parameter]], "parameter"))){
    parameter_name <- format_parameter_names(attr(plot_data, "parameter_name"), attr(prior_list[[parameter]], "parameter"), formula_prefix = formula_prefix)
  }else{
    parameter_name <- attr(plot_data, "parameter_name")
  }

  # get default plot settings
  dots      <- list(...)
  main      <- print(prior_list[[parameter]], plot = TRUE, short_name = short_name, parameter_names = parameter_names)
  xlab      <- parameter_name
  ylab      <- "Density"

  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab

  if(is.null(dots[["col"]]))  dots$col  <-  "black"
  if(is.null(dots[["lty"]]))  dots$lty  <-  1
  if(is.null(dots[["lwd"]]))  dots$lwd  <-  1

  if(length(dots[["xlab"]]) == 1) dots$xlab <- rep(dots[["xlab"]], attr(plot_data, "chains"))
  if(length(dots[["col"]]) == 1)  dots$col  <- rep(dots$col, attr(plot_data, "chains"))
  if(length(dots[["lty"]]) == 1)  dots$lty  <- rep(dots$lty, attr(plot_data, "chains"))
  if(length(dots[["lwd"]]) == 1)  dots$lwd  <- rep(dots$lwd, attr(plot_data, "chains"))

  plots <- list()

  for(i in seq_along(plot_data)){

    temp_dots <- dots

    temp_dots[["xlab"]] <- temp_dots[["xlab"]][i]
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
        plots[[i]] <- plots[[i]] + do.call(.geom_diagnostics.chain, temp_args)
      }else{
        do.call(.lines_diagnostics.chain, temp_args)
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


JAGS_diagnostics_autocorrelation <- function(x) x



.diagnostics_plot_data         <- function(fit, parameter, prior_list, transformations, transform_orthonormal){

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
  if(is.prior.factor(prior_list[[parameter]])){
    if(attr(prior_list[[parameter]], "levels") > 2){
      model_samples <- model_samples[,paste0(parameter, "[", 1:(attr(prior_list[[parameter]], "levels")-1), "]"),drop = FALSE]
    }else{
      model_samples <- model_samples[,parameter,drop = FALSE]
    }
  }else if(is.prior.weightfunction(prior_list[[parameter]])){
    model_samples <- model_samples[,paste0("omega", "[", (length(weightfunctions_mapping(list(prior_list[[parameter]]), cuts_only = TRUE)) - 2):1, "]"),drop = FALSE]
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

  prior_list      <- prior_list[parameter]
  parameter_names <- parameter

  # mostly adapted from runjags_estimates_table
  # apply transformations
  if(!is.null(transformations)){
    for(par in names(transformations)){
      model_samples[,par] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[,par]), transformations[[par]][["arg"]]))
    }
  }

  # transform orthonormal factors to differences from runjags_estimates_table
  if(transform_orthonormal & any(sapply(prior_list, is.prior.orthonormal))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.orthonormal)]){

      if((attr(prior_list[[par]], "levels") - 1) == 1){
        par_names <- par
      }else{
        par_names <- paste0(par, "[", 1:(attr(prior_list[[par]], "levels") - 1), "]")
      }

      orthonormal_samples <- model_samples[,par_names,drop = FALSE]
      model_samples       <- orthonormal_samples %*% t(contr.orthonormal(1:attr(prior_list[[par]], "levels")))

      if(attr(prior_list[[par]], "interaction")){
        if(length(attr(prior_list[[par]], "level_names")) == 1){
          parameter_names <- paste0(par, " [dif: ", attr(prior_list[[par]], "level_names")[[1]],"]")
        }else{
          stop("orthonormal de-transformation for interaction of multiple factors is not implemented.")
        }
      }else{
        parameter_names <- paste0(par, " [dif: ", attr(prior_list[[par]], "level_names"),"]")
      }
    }
  }else if(any(sapply(prior_list, is.prior.orthonormal))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.orthonormal)]){
      if((attr(prior_list[[par]], "levels") - 1) == 1){
        parameter_names <- par
      }else{
        parameter_names <- paste0(par, "[", 1:(attr(prior_list[[par]], "levels") - 1), "]")
      }
    }
  }

  # rename treatment factor levels
  if(any(sapply(prior_list, is.prior.dummy))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.dummy)]){
      if(!attr(prior_list[[par]], "interaction")){
        if(attr(prior_list[[par]], "levels") == 2){
          parameter_names <- paste0(par,"[",attr(prior_list[[par]], "level_names")[-1], "]")
        }else{
          parameter_names <- paste0(par,"[",attr(prior_list[[par]], "level_names")[-1], "]")
        }
      }else if(length(attr(prior_list[[par]], "levels")) == 1){
        parameter_names <- paste0(par,"[",attr(prior_list[[par]], "level_names")[[1]][-1], "]")
      }
    }
  }

  # rename weightfunctions factor levels
  if(any(sapply(prior_list, is.prior.weightfunction))){
    for(par in names(prior_list)[sapply(prior_list, is.prior.weightfunction)]){
      omega_cuts      <- weightfunctions_mapping(prior_list[par], cuts_only = TRUE)
      parameter_names <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
      parameter_names <- parameter_names[-1]
    }
  }

  # attach the relevant attributes
  colnames(model_samples)          <- parameter_names
  attr(model_samples, "chain")     <- do.call(c, samples_chain)
  attr(model_samples, "iter")      <- do.call(c, samples_iter)
  attr(model_samples, "parameter") <- parameter
  attr(model_samples, "prior")     <- prior_list[[parameter]]

  return(model_samples)
}
.diagnostics_plot_data_density <- function(plot_data, n_points, xlim){

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
      x_range <- range(plot_data[,i])
    }else{
      x_range <- xlim
    }

    for(j in seq_along(unique(chain))){

      temp_args    <- list(x = plot_data[chain == j,i], n = n_points, from = x_range[1], to = x_range[2])
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
.diagnostics_plot_data_trace   <- function(plot_data, n_points, ylim){

  chain <- attr(plot_data, "chain")
  iter  <- attr(plot_data, "iter")
  prior <- attr(plot_data, "prior")

  out   <- list()

  for(i in 1:ncol(plot_data)){

    if(is.null(ylim)){
      y_range <- range(plot_data[,i])
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

.lines_diagnostics.chain <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.geom_diagnostics.chain  <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  geom <- ggplot2::geom_line(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes_string(
      x = "x",
      y = "y"),
    size = lwd, linetype = lty, color = col)

  return(geom)
}


#
#
# diagnostics <- function(fit, parameter, type, plot_type = "base", show_models = NULL,
#                         , title = is.null(show_models) | length(show_models) > 1, ...){
#
#   # check settings
#   if(class(fit) != "RoBMA")
#     stop("Diagnostics are available only for RoBMA models.")
#   if(fit$add_info$save == "min")
#     stop("Diagnostics cannot be produced because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' while fitting the model (see ?RoBMA for more details).")
#   BayesTools::check_char(parameter, "parameter")
#   BayesTools::check_char(type, "type")
#   BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
#
#   # deal with bad type names
#   if(substr(type, 1, 1) == "c"){
#     type <- "chains"
#   }else if(substr(type, 1, 1) == "t"){
#     type <- "chains" # for trace
#   }else if(substr(type, 1, 1) == "d"){
#     type <- "densities"
#   }else if(substr(type, 1, 1) == "a"){
#     type <- "autocorrelation"
#   }else{
#     stop("Unsupported diagnostics type. See '?diagnostics' for more details.")
#   }
#
#   # deal with bad parameter names for PET-PEESE, weightfunction
#   if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
#     parameter <- "omega"
#   }else if(parameter %in% c("mu", "tau", "PET", "PEESE")){
#     parameters <- parameter
#   }else{
#     stop("The passed parameter is not supported for plotting. See '?plot.RoBMA' for more details.")
#   }
#
#   # omit the first figure for publication bias weights (it is constants for all interesting weightfunctions)
#   if(parameter == "omega"){
#     show_figures <-  -1
#   }else{
#     show_figures <-  NULL
#   }
#
#
#   # do the plotting
#   out <- list()
#
#   models_ind <- 1:length(fit$models)
#   if(!is.null(show_models)){
#     models_ind <- models_ind[show_models]
#   }
#
#   # a message with info about multiple plots
#   if(plot_type == "base" & (length(models_ind) > 1 | parameter == "omega"))
#     message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")
#
#
#   for(m in models_ind){
#
#     temp_out  <- NULL
#     # add ability to transform the estimates with 'par_transform'
#     temp_data <- .diagnostics_plot_data(fit, m, parameter, NULL)
#
#     # deal with no parameter in model
#     if(is.null(temp_data)){
#       if(length(models_ind) == 1){
#         message("Selected model does not containt the parameter of interest.")
#         return(invisible(NULL))
#       }else{
#         out[m] <- temp_out
#         next
#       }
#     }
#
#     # make the plots
#     par_ind <- 1:length(temp_data)
#     if(!is.null(show_figures)){
#       par_ind <- par_ind[show_figures]
#     }
#
#     for(i in par_ind){
#       if(type == "chains"){
#         temp_out <- c(temp_out, list(.diagnostics_plot_trace(temp_data[[i]], plot_type, if(title) m else NULL, ...)))
#       }else if(type == "densities"){
#         temp_out <- c(temp_out, list(.diagnostics_plot_density(temp_data[[i]], plot_type, if(title) m else NULL, parameter, ...)))
#       }else if(type == "autocorrelation"){
#         temp_out <- c(temp_out, list(.diagnostics_plot_ac(temp_data[[i]], plot_type, if(title) m else NULL, lags, ...)))
#       }
#     }
#
#     if(length(temp_out) == 1){
#       temp_out <- temp_out[[1]]
#     }
#
#     if(length(models_ind) == 1){
#       out    <- temp_out
#     }else{
#       out[m] <- list(temp_out)
#     }
#   }
#
#   # return the plots
#   if(plot_type == "base"){
#     return(invisible(NULL))
#   }else if(plot_type == "ggplot"){
#     return(out)
#   }
# }
#
#
# .diagnostics_plot_trace   <- function(plot_data, plot_type, title, ...){
#
#   if(plot_type == "base"){
#
#     # save plotting settings
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(mar = oldpar[["mar"]]))
#
#     # set up margins
#     if(length(list(...)) == 0){
#       graphics::par(mar = c(4, 4, if(!is.null(title)) 3 else 1, 1))
#     }else{
#       graphics::par(list(...))
#     }
#
#     graphics::plot(NA, type = "n", xlim = range(plot_data$samp$iteration), ylim = range(plot_data$samp$value),
#                    xlab = "", ylab = "", bty = "n", las = 1)
#     for(i in as.numeric(unique(plot_data$samp$chain))){
#       graphics::lines(plot_data$samp$iteration[plot_data$samp$chain == i], plot_data$samp$value[plot_data$samp$chain == i],
#                       col = .diagnostics_color(plot_data$nchains)[i])
#     }
#     if(!is.null(title)){
#       graphics::mtext(paste0("Model ",title), side = 3, line = 1, cex = 1.25)
#     }
#     graphics::mtext(plot_data$parameter, side = 2, line = 2.5, cex = 1.25)
#
#     graph <- NULL
#
#   }else if(plot_type == "ggplot"){
#
#     graph <- ggplot2::ggplot(plot_data$samp, ggplot2::aes_string(x = "iteration", y = "value", color = "chain")) +
#       ggplot2::geom_path() +
#       ggplot2::scale_color_manual(values = .diagnostics_color(plot_data$nchains))
#     temp_x_range <- range(plot_data$samp$iteration)
#     temp_y_range <- range(plot_data$samp$value)
#     graph <- graph + ggplot2::scale_x_continuous(
#       name   = "Iterations",
#       limits = temp_x_range,
#       breaks = pretty(temp_x_range, n = 3),
#       labels = pretty(temp_x_range, n = 3)
#     ) +
#       ggplot2::scale_y_continuous(
#         name   = plot_data$parameter,
#         limits = temp_y_range,
#         breaks = pretty(temp_y_range),
#         labels = pretty(temp_y_range)
#       )
#     if(!is.null(title)){
#       graph <- graph + ggplot2::ggtitle(paste0("Model ",title))
#     }
#   }
#
#   return(graph)
# }
# .diagnostics_plot_density <- function(plot_data, plot_type, title, par, ...){
#
#   if(plot_type == "base"){
#
#     # save plotting settings
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(mar = oldpar[["mar"]]))
#
#     # set up margins
#     if(length(list(...)) == 0){
#       graphics::par(mar = c(4, 4, if(!is.null(title)) 3 else 1, 1))
#     }else{
#       graphics::par(list(...))
#     }
#     with_trunc <- list()
#     if(!is.infinite(plot_data$lower)){
#       with_trunc$from <- plot_data$lower
#     }
#     if(!is.infinite(plot_data$upper)){
#       with_trunc$to   <- plot_data$upper
#     }
#
#
#     temp_den <- vector(mode = "list", length = length(unique(plot_data$samp$chain)))
#     for(i in as.numeric(unique(plot_data$samp$chain))){
#       # deal with first weights if requested
#       if(all(plot_data$samp$value[plot_data$samp$chain == i] == 1) & par == "omega"){
#         temp_den[[i]] <- NULL
#       }else{
#         temp_den[[i]] <- do.call(stats::density, c(list(x = plot_data$samp$value[plot_data$samp$chain == i]), with_trunc))
#       }
#     }
#
#     graphics::plot(
#       NA, type = "n",
#       xlim = if(all(sapply(temp_den, is.null))) c(0, 1) else range(sapply(1:length(temp_den), function(i)temp_den[[i]]$x)),
#       ylim = if(all(sapply(temp_den, is.null))) c(0, 1) else c(0, max(sapply(1:length(temp_den), function(i)temp_den[[i]]$y))),
#       xlab = "", ylab = "", bty = "n", las = 1)
#     for(i in 1:length(temp_den)){
#       if(is.null(temp_den[[i]]) & par == "omega"){
#         graphics::arrows(x0 = 1, y0 = 0, y1 = 1, lwd = 2, lty = 1, col = .diagnostics_color(plot_data$nchains)[i])
#       }else{
#         graphics::lines(temp_den[[i]], col = .diagnostics_color(plot_data$nchains)[i])
#         graphics::polygon(
#           x = c(if(!is.infinite(plot_data$lower)) plot_data$lower, temp_den[[i]]$x, if(!is.infinite(plot_data$upper)) plot_data$upper),
#           y = c(if(!is.infinite(plot_data$lower)) 0,               temp_den[[i]]$y, if(!is.infinite(plot_data$upper)) 0),
#           border = .diagnostics_color(plot_data$nchains)[i],
#           col    = scales::alpha(.diagnostics_color(plot_data$nchains)[i], alpha = .5))
#       }
#     }
#     if(!is.null(title)){
#       graphics::mtext(paste0("Model ",title), side = 3, line = 1, cex = 1.25)
#     }
#     graphics::mtext(if(all(sapply(temp_den, is.null))) "Probability" else "Density", side = 2, line = 2.5, cex = 1.25)
#     graphics::mtext(plot_data$parameter,                                             side = 1, line = 2.5, cex = 1.25)
#
#     graph <- NULL
#
#   }else if(plot_type == "ggplot"){
#
#     graph <-  ggplot2::ggplot(plot_data$samp, ggplot2::aes_string(x = "value")) +
#       ggplot2::geom_density(mapping = ggplot2::aes_string(fill = "chain"), color = "black", alpha = .5) +
#       ggplot2::scale_fill_manual(values = .diagnostics_color(plot_data$nchains))
#     temp_y_max   <- max(ggplot2::ggplot_build(graph)$data[[1]]$density)
#     temp_x_range <- if(par == "omega") c(0, 1) else range(plot_data$samp$value)
#     graph <- graph +  ggplot2::scale_y_continuous(
#       name   = "Density",
#       limits = range(pretty(c(0, temp_y_max))),
#       breaks = pretty(c(0, temp_y_max)),
#       labels = pretty(c(0, temp_y_max))
#     ) +
#       ggplot2::scale_x_continuous(
#         name   = plot_data$parameter,
#         limits = range(pretty(temp_x_range)),
#         breaks = pretty(temp_x_range),
#         labels = pretty(temp_x_range)
#       )
#     if(!is.null(title)){
#       graph <- graph + ggplot2::ggtitle(paste0("Model ",title))
#     }
#
#   }
#
#   return(graph)
# }
# .diagnostics_plot_ac      <- function(plot_data, plot_type, title, lags = 30, ...){
#
#   ac_dat   <- .diagnostics_ac_data(dat = plot_data$samp, lags = lags)
#
#   if(plot_type == "base"){
#
#     # save plotting settings
#     oldpar <- graphics::par(no.readonly = TRUE)
#     on.exit(graphics::par(mar = oldpar[["mar"]]))
#
#     # set up margins
#     if(length(list(...)) == 0){
#       graphics::par(mar = c(4,4,3,1))
#     }else{
#       graphics::par(list(...))
#     }
#
#     temp_dat <- as.numeric(by(ac_dat$ac, ac_dat$lag, mean))
#     temp_dat[is.nan(temp_dat)] <- 1
#     graphics::barplot(temp_dat, names.arg = unique(ac_dat$lag), col = "#B2001D", las = 1)
#     graphics::mtext("Lag",                  side = 1, line = 2.5, cex = 1.25)
#     graphics::mtext("Avg. autocorrelation", side = 2, line = 2.5, cex = 1.25)
#     if(!is.null(title)){
#       graphics::mtext(bquote(paste("Model"," ", .(title),":"," ", .(eval(plot_data$parameter)))), side = 3, line = 1, cex = 1.25)
#     }else{
#       graphics::mtext(plot_data$parameter, side = 3, line = 1, cex = 1.25)
#     }
#
#     graph <- NULL
#
#   }else if(plot_type == "ggplot"){
#     graph     <- ggplot2::ggplot(ac_dat, ggplot2::aes_string(x = "lag", y = "ac")) +
#       ggplot2::geom_bar(size = .5, color = "black", fill = "#B2001D", position = "dodge", stat = "summary", fun = "mean") +
#       ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
#       ggplot2::labs(x = "Lag", y = "Avg. autocorrelation")
#     if(!is.null(title)){
#       graph <- graph + ggplot2::ggtitle(bquote(paste("Model"," ", .(title),":"," ", .(eval(plot_data$parameter)))))
#     }else{
#       graph <- graph + ggplot2::ggtitle(plot_data$parameter)
#     }
#   }
#
#   return(graph)
# }
# .diagnostics_ac_data      <- function(dat, lags){
#   ch <- dat[, grep("chain", colnames(dat))]
#   nc <- length(unique(ch))
#
#   ac_list <- tapply(dat$value, INDEX = ch, FUN = function(x)stats::acf(x, lag.max = lags, plot = FALSE)$acf[, , 1L], simplify = FALSE)
#
#   nl <- lags + 1
#   ch <- factor(rep(1:nc, each = nl), labels = paste0("chain:", 1:nc))
#   ll <- rep(seq(0, lags), nc)
#
#   return(data.frame(chains = ch, ac = do.call(c, ac_list), lag = ll))
# }
# .diagnostics_color        <- function(n){
#   return(rep_len(c("#E66101", "#998EC3", "#542788", "#F1A340", "#D8DAEB", "#FEE0B6"), n))
# }
# .diagnostics_plot_data    <- function(fit, model, par, par_transform){
#
#   if(length(fit$models[[model]]$fit) == 0){
#
#     return(NULL)
#
#   }else{
#
#     # do not plot spike priors
#     if(is.prior.point(fit$models[[model]]$priors[[par]]))
#       return(NULL)
#
#     samples <- coda::as.array.mcmc.list(fit$models[[model]]$fit$mcmc, drop = FALSE)
#
#     if(!any(grepl(par, dimnames(samples)$var)))
#       return(NULL)
#
#     # create parameter names and get parameter indexes
#     if(par %in% c("mu", "tau", "PET", "PEESE")){
#       ind       <- c(1:length(dimnames(samples)$var))[par == dimnames(samples)$var]
#       par_names <- .plot.RoBMA_par_names(par, fit, fit$add_info$prior_scale)
#     }else{
#       ind <- c(1:length(dimnames(samples)$var))[grepl(par, dimnames(samples)$var)]
#       ind <- rev(ind)
#       summary_info <- summary(fit, "individual")
#       summary_info <- summary_info[["models"]][[model]][["estimates"]]
#       omega_names  <- rownames(summary_info)[grepl(par, rownames(summary_info))]
#       par_names    <- vector("list", length = length(omega_names))
#       for(i in 1:length(par_names)){
#         par_names[[i]] <- bquote(~omega[~.(substr(omega_names[i],6,nchar(omega_names[i])))])
#       }
#     }
#
#     plot_data <- list()
#     for(i in 1:length(ind)){
#       plot_data[[dimnames(samples)$var[ind[i]]]] <- list(
#         samp = data.frame(
#           value     = as.vector(samples[,ind[i],]),
#           parameter = dimnames(samples)$var[ind[i]],
#           chain     = as.factor(c(unlist(sapply(1:dim(samples)[3], function(x)rep(x,dim(samples)[1]))))),
#           iteration = rep(1:dim(samples)[1], dim(samples)[3])
#         ),
#         nchains   = dim(samples)[3],
#         nparams   = 1,
#         warmup    = 0,
#         parameter = par_names[[i]],
#         lower     = if(par == "omega") 0 else fit$models[[model]]$priors[[par]]$truncation[["lower"]],
#         upper     = if(par == "omega") 1 else fit$models[[model]]$priors[[par]]$truncation[["upper"]]
#       )
#     }
#
#     # TODO: implement later
#     # transform the values if requested
#     # if(par_transform){
#     #   if(par %in% c("mu", "theta") & fit$add_info$effect_size %in% c("r", "OR")){
#     #     plot_data[[1]]$samp$value <- .transform(plot_data[[1]]$samp$value, fit$add_info$effect_size, fit$add_info$transformation)
#     #   }
#     # }
#
#   }
#
#   return(plot_data)
# }
